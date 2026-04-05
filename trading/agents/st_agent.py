"""Short-Term Momentum Agent — Gemma 4 via Ollama.

Horizon: 1–5 trading days.
Focus:   intraday momentum breakouts, RSI extremes on 15m charts,
         MACD crossovers with volume confirmation, tight setups only.
LLM:     Gemma 4 local (speed matters — runs every 30 min intraday).
Fallback: rule-based Python signal if Ollama is unavailable.
"""
from __future__ import annotations

import logging
from datetime import datetime, timezone

import pandas as pd
import yfinance as yf
from pydantic import BaseModel, Field, field_validator

from trading.config.strategy_config import SHORT_TERM
from trading.llm import ollama_client
from trading.tools.indicators import rule_based_signal, compute_all

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = f"""You are a short-term momentum trader analyzing US equities for 1-5 day swing trades.

You receive BOTH daily indicators (trend context) and 15-minute intraday indicators (entry timing).

Decision framework:
1. Daily trend filter (required): Price must be above SMA50 for BUY, below for SELL.
   Without trend alignment, output HOLD regardless of intraday signal.
2. Intraday momentum (primary signal):
   - 15m RSI crossing above 50 from below → bullish momentum entry
   - 15m MACD histogram turning positive after being negative → entry confirmation
   - Price reclaiming 15m VWAP after pullback → intraday strength
   - 15m BB squeeze breaking out (bandwidth contracting then expanding) → breakout
3. Volume confirmation: current 15m volume must be > 1.5x the 15m average
4. Risk/reward: only BUY if setup allows stop at ATR×{SHORT_TERM.stop_atr_multiplier:.1f} from entry
   with target at {SHORT_TERM.reward_risk_ratio:.1f}:1 reward/risk.

Output HOLD if: no clear 15m momentum, trend against direction, or conflicting signals.
Conviction (strength) should be HIGH (>0.7) only when all 3 conditions align:
trend + intraday momentum + volume confirmation."""


class STSignalOutput(BaseModel):
    ticker: str
    signal: str = Field(pattern="^(BUY|SELL|HOLD)$")
    strength: float = Field(ge=0.0, le=1.0)
    trend_aligned: bool
    intraday_momentum: str         # brief description of the intraday setup
    volume_confirmed: bool
    reasons: list[str]

    @field_validator("reasons")
    @classmethod
    def limit_reasons(cls, v):
        return v[:4]


def _fetch_intraday_indicators(ticker: str) -> dict:
    """Fetch 15-minute bars and compute indicators. Returns dict or {} on failure."""
    try:
        df = yf.download(ticker, period="5d", interval="15m", progress=False, auto_adjust=True)
        if df.empty or len(df) < 20:
            return {}
        df.columns = [c.lower() for c in df.columns]
        indicators = compute_all(df)
        # VWAP approximation: (H+L+C)/3 rolling cumulative
        typical = (df["high"] + df["low"] + df["close"]) / 3
        vwap = (typical * df["volume"]).cumsum() / df["volume"].cumsum()
        indicators["vwap"] = float(vwap.iloc[-1])
        indicators["price_vs_vwap"] = float(df["close"].iloc[-1] - vwap.iloc[-1])
        # 15m volume ratio
        avg_vol = df["volume"].tail(20).mean()
        indicators["intraday_volume_ratio"] = float(df["volume"].iloc[-1] / avg_vol) if avg_vol > 0 else 1.0
        return indicators
    except Exception as e:
        logger.debug("Intraday fetch failed for %s: %s", ticker, e)
        return {}


def analyze_ticker(snapshot: dict) -> dict:
    """Analyze one ticker for short-term momentum entry."""
    ticker = snapshot["ticker"]
    daily_ind = snapshot.get("indicators", {})
    intraday_ind = _fetch_intraday_indicators(ticker)

    def _fmt(v):
        return f"{v:.3f}" if isinstance(v, float) else str(v)

    user_prompt = f"""Short-term momentum analysis for {ticker}:

DAILY CONTEXT (trend filter):
- Price: ${snapshot['price']:.2f}
- Above SMA50: {daily_ind.get('above_sma50', 'N/A')}  |  Above SMA200: {daily_ind.get('above_sma200', 'N/A')}
- Daily RSI(14): {_fmt(daily_ind.get('rsi_14'))}
- Daily MACD histogram: {_fmt(daily_ind.get('macd_hist'))}
- Daily ATR(14): {_fmt(daily_ind.get('atr_14'))} (stop distance at {SHORT_TERM.stop_atr_multiplier}×ATR = ${(daily_ind.get('atr_14') or 0)*SHORT_TERM.stop_atr_multiplier:.2f})

INTRADAY 15-MINUTE SIGNALS (entry timing):
- 15m RSI(14): {_fmt(intraday_ind.get('rsi_14'))}
- 15m MACD hist: {_fmt(intraday_ind.get('macd_hist'))}
- 15m BB%: {_fmt(intraday_ind.get('bb_pct'))}
- Price vs VWAP: {_fmt(intraday_ind.get('price_vs_vwap'))} (positive = above VWAP)
- 15m Volume ratio: {_fmt(intraday_ind.get('intraday_volume_ratio'))}x

Return JSON: ticker, signal, strength, trend_aligned (bool), intraday_momentum (string),
volume_confirmed (bool), reasons (list)."""

    def _fallback(snap: dict) -> dict:
        base = rule_based_signal(snap)
        above_sma50 = snap.get("indicators", {}).get("above_sma50", False)
        return {
            **base,
            "trend_aligned": bool(above_sma50) if base["signal"] == "BUY" else not bool(above_sma50),
            "intraday_momentum": "rule-based fallback",
            "volume_confirmed": snap.get("volume_ratio", 1.0) > 1.5,
        }

    result = ollama_client.call_structured(
        system_prompt=SYSTEM_PROMPT,
        user_prompt=user_prompt,
        output_schema=STSignalOutput,
        max_tokens=300,
        fallback_fn=_fallback,
        fallback_input=snapshot,
    )

    return {
        "ticker": result.ticker,
        "signal": result.signal,
        "strength": result.strength,
        "reasons": result.reasons,
        "source": "st_agent",
        "strategy": SHORT_TERM.name,
        "intraday_momentum": result.intraday_momentum,
        "trend_aligned": result.trend_aligned,
        "volume_confirmed": result.volume_confirmed,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def run(snapshots: list[dict]) -> list[dict]:
    """Run short-term agent on all snapshots. Returns strategy-tagged signals."""
    signals = []
    for snapshot in snapshots:
        ticker = snapshot.get("ticker", "")
        try:
            sig = analyze_ticker(snapshot)
            signals.append(sig)
        except Exception as e:
            logger.error("ST agent failed for %s: %s", ticker, e)
            fallback = rule_based_signal(snapshot)
            fallback.update({
                "strategy": SHORT_TERM.name,
                "source": "st_agent_fallback",
                "timestamp": datetime.now(timezone.utc).isoformat(),
            })
            signals.append(fallback)
    return signals
