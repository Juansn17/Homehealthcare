"""Long-Term Fundamental Agent — Claude Sonnet.

Horizon: 2–8 weeks (10–40 trading days).
Focus:   earnings acceleration, sector rotation leaders, macro tailwinds,
         technical breakouts from multi-week consolidation.
LLM:     Claude Sonnet — depth of reasoning matters more than speed.
"""
from __future__ import annotations

import logging
from datetime import datetime, timezone

import pandas as pd
import yfinance as yf
from pydantic import BaseModel, Field, field_validator

from trading.config.strategy_config import LONG_TERM
from trading.llm import claude_client
from trading.tools.news_macro import (
    get_earnings_data,
    get_valuation,
    get_news_headlines,
    get_macro_indicators,
    get_sector_performance,
)

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = f"""You are a portfolio manager specializing in 2-8 week positional trades in US equities.

Your edge is identifying stocks entering multi-week trends driven by:
1. Earnings acceleration — EPS growth rate improving quarter-over-quarter
2. Sector rotation — sector ETF showing relative strength vs SPY over 3+ months
3. Macro tailwinds — current rate/inflation/growth environment favoring this sector
4. Technical structure — weekly chart showing breakout from consolidation (not extended)

Decision framework (ALL must hold for BUY):
- Fundamental: forward P/E < sector average OR strong EPS growth (>20% YoY)
- Sector: sector ETF 3-month return > SPY 3-month return (relative strength)
- Macro: macro environment not a headwind (e.g., avoid rate-sensitive names if rates rising)
- Technical: price above SMA200 AND SMA50 trending upward (not near resistance)
- Risk/reward: stop at ATR×{LONG_TERM.stop_atr_multiplier:.1f} below current price, target ≥{LONG_TERM.reward_risk_ratio:.1f}:1 R:R

Output HOLD if:
- Earnings report within next 5 trading days (binary event risk)
- Strong technical but poor fundamentals (or vice versa — both required)
- Sector is rotating OUT (underperforming SPY on 3-month basis)

Strength >0.8 only when fundamentals AND sector AND macro AND technicals all align."""


class LTSignalOutput(BaseModel):
    ticker: str
    signal: str = Field(pattern="^(BUY|SELL|HOLD)$")
    strength: float = Field(ge=0.0, le=1.0)
    fundamental_score: str        # "strong" | "neutral" | "weak"
    sector_momentum: str          # "outperforming" | "inline" | "underperforming"
    macro_stance: str             # "tailwind" | "neutral" | "headwind"
    earnings_risk: bool           # True = earnings within 5 days → avoid
    pe_ratio: float | None = None
    reasons: list[str]

    @field_validator("reasons")
    @classmethod
    def limit_reasons(cls, v):
        return v[:5]


def _weekly_trend(ticker: str) -> dict:
    """Fetch weekly OHLCV and return simple trend context."""
    try:
        df = yf.download(ticker, period="6mo", interval="1wk", progress=False, auto_adjust=True)
        if df.empty or len(df) < 4:
            return {}
        close = df["Close"] if "Close" in df.columns else df["close"]
        close = close.dropna()
        ret_4w = float(close.iloc[-1] / close.iloc[-4] - 1) if len(close) >= 4 else None
        ret_12w = float(close.iloc[-1] / close.iloc[-12] - 1) if len(close) >= 12 else None
        return {
            "weekly_return_4w": ret_4w,
            "weekly_return_12w": ret_12w,
            "weekly_bars": len(close),
        }
    except Exception as e:
        logger.debug("Weekly trend fetch failed for %s: %s", ticker, e)
        return {}


def analyze_ticker(snapshot: dict, macro: dict, sector_perf: dict) -> dict:
    """Deep fundamental analysis for one ticker."""
    ticker = snapshot["ticker"]
    daily_ind = snapshot.get("indicators", {})
    valuation = get_valuation(ticker)
    earnings = get_earnings_data(ticker)
    news = get_news_headlines(ticker, max_items=5)
    weekly = _weekly_trend(ticker)

    # Sector context
    from trading.config.universe import SECTOR_MAP
    sector_name = SECTOR_MAP.get(ticker, valuation.get("sector", "Unknown"))
    sector_etf_map = {
        "Technology": "XLK", "Financials": "XLF", "Energy": "XLE",
        "Healthcare": "XLV", "Consumer Discretionary": "XLY",
        "Communication": "XLC", "Industrials": "XLI",
    }
    etf = sector_etf_map.get(sector_name)
    sector_data = sector_perf.get(etf, {}) if etf else {}
    spy_data = sector_perf.get("SPY", {})
    spy_3m = spy_data.get("return_3m")
    sector_3m = sector_data.get("return_3m")
    relative_str = (
        f"{sector_3m:.1%} vs SPY {spy_3m:.1%}"
        if sector_3m is not None and spy_3m is not None
        else "N/A"
    )

    news_text = "\n".join(f"- {n['title']} ({n['publisher']})" for n in news) or "No recent news"

    def _fmt(v):
        return f"{v:.3f}" if isinstance(v, float) else str(v)

    user_prompt = f"""Long-term positional analysis for {ticker} (target hold: {LONG_TERM.horizon_days_min}-{LONG_TERM.horizon_days_max} trading days):

FUNDAMENTALS:
- Sector: {sector_name}
- P/E (trailing): {valuation.get('pe_ratio', 'N/A')}
- Forward P/E: {valuation.get('forward_pe', 'N/A')}
- PEG: {valuation.get('peg_ratio', 'N/A')}
- EPS growth: {valuation.get('earnings_growth', 'N/A')}
- Revenue growth: {valuation.get('revenue_growth', 'N/A')}
- ROE: {valuation.get('return_on_equity', 'N/A')}
- Debt/Equity: {valuation.get('debt_to_equity', 'N/A')}
- Next earnings: {earnings.get('next_earnings_date', 'N/A')}  ← flag if within 5 trading days

SECTOR MOMENTUM ({sector_name} / {etf or 'N/A'}):
- 1-month return: {_fmt(sector_data.get('return_1m'))}
- 3-month return: {relative_str}

MACRO ENVIRONMENT:
- Fed Funds Rate: {_fmt(macro.get('fed_funds_rate'))}%
- 10Y-2Y Spread: {_fmt(macro.get('t10y2y_spread'))}%  (negative = inverted = caution)
- CPI YoY: {_fmt(macro.get('cpi_yoy'))}
- Unemployment: {_fmt(macro.get('unemployment_rate'))}%

TECHNICAL STRUCTURE (daily + weekly):
- Price: ${snapshot['price']:.2f}
- Above SMA50: {daily_ind.get('above_sma50', 'N/A')}
- Above SMA200: {daily_ind.get('above_sma200', 'N/A')}
- Golden cross: {daily_ind.get('golden_cross', 'N/A')}
- Daily ATR(14): {_fmt(daily_ind.get('atr_14'))} → stop at {LONG_TERM.stop_atr_multiplier}×ATR = ${(daily_ind.get('atr_14') or 0)*LONG_TERM.stop_atr_multiplier:.2f} below
- 4-week price change: {_fmt(weekly.get('weekly_return_4w'))}
- 12-week price change: {_fmt(weekly.get('weekly_return_12w'))}

NEWS:
{news_text}

Return JSON: ticker, signal, strength, fundamental_score, sector_momentum,
macro_stance, earnings_risk (bool), pe_ratio, reasons."""

    result = claude_client.call_structured(
        system_prompt=SYSTEM_PROMPT,
        user_prompt=user_prompt,
        output_schema=LTSignalOutput,
        max_tokens=600,
    )

    return {
        "ticker": result.ticker,
        "signal": result.signal,
        "strength": result.strength,
        "reasons": result.reasons,
        "source": "lt_agent",
        "strategy": LONG_TERM.name,
        "fundamental_score": result.fundamental_score,
        "sector_momentum": result.sector_momentum,
        "macro_stance": result.macro_stance,
        "earnings_risk": result.earnings_risk,
        "pe_ratio": result.pe_ratio,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def run(snapshots: list[dict]) -> list[dict]:
    """Run long-term agent. Returns strategy-tagged signals."""
    macro = get_macro_indicators()
    sector_perf = get_sector_performance()
    signals = []

    for snapshot in snapshots:
        ticker = snapshot.get("ticker", "")
        try:
            sig = analyze_ticker(snapshot, macro, sector_perf)
            signals.append(sig)
        except Exception as e:
            logger.error("LT agent failed for %s: %s", ticker, e)
            signals.append({
                "ticker": ticker,
                "signal": "HOLD",
                "strength": 0.3,
                "reasons": [f"Analysis failed: {str(e)[:100]}"],
                "source": "lt_agent_error",
                "strategy": LONG_TERM.name,
                "fundamental_score": "neutral",
                "sector_momentum": "inline",
                "macro_stance": "neutral",
                "earnings_risk": False,
                "pe_ratio": None,
                "timestamp": datetime.now(timezone.utc).isoformat(),
            })

    return signals
