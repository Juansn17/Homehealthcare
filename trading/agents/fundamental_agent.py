"""Fundamental / Macro Analysis Agent — Claude Sonnet.

Analyzes earnings momentum, valuation, macro environment, and sector
context for tickers that cleared the technical pre-filter.
"""
from __future__ import annotations

import logging
from datetime import datetime, timezone

from pydantic import BaseModel, Field, field_validator

from trading.llm import claude_client
from trading.tools.news_macro import (
    get_earnings_data,
    get_valuation,
    get_news_headlines,
    get_macro_indicators,
    get_sector_performance,
)

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are a fundamental equity analyst specializing in US large-cap stocks.
You analyze earnings momentum, valuation, macroeconomic context, and recent news to determine
if a stock is fundamentally attractive for a swing trade (3-10 trading day horizon).

Consider:
1. Valuation: Is the stock cheap/fair/expensive vs peers and historical averages?
2. Earnings: Recent EPS surprises? Upcoming earnings risk? Revenue trend?
3. Macro: Is the current macro environment (rates, inflation, growth) a tailwind or headwind?
4. Sector: Is the sector currently in favor (positive momentum) or out of favor?
5. News: Any near-term catalysts or red flags?

Output a BUY/SELL/HOLD signal with conviction strength (0.0-1.0) and concise reasons."""


class FundamentalSignalOutput(BaseModel):
    ticker: str
    signal: str = Field(pattern="^(BUY|SELL|HOLD)$")
    strength: float = Field(ge=0.0, le=1.0)
    pe_ratio: float | None = None
    earnings_surprise: float | None = None
    macro_context: str = ""
    reasons: list[str]

    @field_validator("reasons")
    @classmethod
    def limit_reasons(cls, v):
        return v[:5]


def analyze_ticker(snapshot: dict, macro: dict, sector_perf: dict) -> dict:
    """Run fundamental analysis for one ticker."""
    ticker = snapshot["ticker"]
    valuation = get_valuation(ticker)
    earnings = get_earnings_data(ticker)
    news = get_news_headlines(ticker, max_items=5)

    # Build sector context
    from trading.config.universe import SECTOR_MAP
    sector_name = SECTOR_MAP.get(ticker, "Unknown")
    sector_etf_map = {v: k for k, v in {
        "XLK": "Technology", "XLF": "Financials", "XLE": "Energy",
        "XLV": "Healthcare", "XLY": "Consumer Discretionary",
        "XLC": "Communication", "XLI": "Industrials",
    }.items()}
    etf = sector_etf_map.get(sector_name)
    sector_data = sector_perf.get(etf, {}) if etf else {}

    news_text = "\n".join(f"- {n['title']} ({n['publisher']})" for n in news) or "No recent news"

    user_prompt = f"""Fundamental analysis for {ticker}:

VALUATION:
- P/E (trailing): {valuation.get('pe_ratio', 'N/A')}
- Forward P/E: {valuation.get('forward_pe', 'N/A')}
- PEG Ratio: {valuation.get('peg_ratio', 'N/A')}
- Price/Book: {valuation.get('price_to_book', 'N/A')}
- Profit Margin: {valuation.get('profit_margin', 'N/A')}
- Debt/Equity: {valuation.get('debt_to_equity', 'N/A')}
- Sector: {valuation.get('sector', sector_name)}

EARNINGS & GROWTH:
- EPS growth: {valuation.get('earnings_growth', 'N/A')}
- Revenue growth: {valuation.get('revenue_growth', 'N/A')}
- Forward EPS: {earnings.get('forward_eps', 'N/A')}
- Next earnings date: {earnings.get('next_earnings_date', 'N/A')}

SECTOR MOMENTUM ({sector_name}):
- 1-month return: {sector_data.get('return_1m', 'N/A')}
- 3-month return: {sector_data.get('return_3m', 'N/A')}

MACRO ENVIRONMENT:
- Fed Funds Rate: {macro.get('fed_funds_rate', 'N/A')}%
- 10Y-2Y Yield Spread: {macro.get('t10y2y_spread', 'N/A')}%
- CPI YoY: {macro.get('cpi_yoy', 'N/A')}
- Unemployment: {macro.get('unemployment_rate', 'N/A')}%

RECENT NEWS:
{news_text}

Return JSON: ticker, signal (BUY/SELL/HOLD), strength (0.0-1.0), pe_ratio, earnings_surprise, macro_context (one sentence), reasons (list)."""

    result = claude_client.call_structured(
        system_prompt=SYSTEM_PROMPT,
        user_prompt=user_prompt,
        output_schema=FundamentalSignalOutput,
        max_tokens=512,
    )

    return {
        "ticker": result.ticker,
        "signal": result.signal,
        "strength": result.strength,
        "pe_ratio": result.pe_ratio,
        "earnings_surprise": result.earnings_surprise,
        "macro_context": result.macro_context,
        "reasons": result.reasons,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def run(snapshots: list[dict]) -> list[dict]:
    """Run fundamental analysis for a filtered list of snapshots."""
    macro = get_macro_indicators()
    sector_perf = get_sector_performance()
    signals = []

    for snapshot in snapshots:
        ticker = snapshot.get("ticker", "")
        try:
            sig = analyze_ticker(snapshot, macro, sector_perf)
            signals.append(sig)
        except Exception as e:
            logger.error("Fundamental agent failed for %s: %s", ticker, e)
            # Produce a neutral HOLD on failure
            signals.append({
                "ticker": ticker,
                "signal": "HOLD",
                "strength": 0.3,
                "pe_ratio": None,
                "earnings_surprise": None,
                "macro_context": "Analysis unavailable",
                "reasons": [f"Analysis failed: {str(e)[:100]}"],
                "timestamp": datetime.now(timezone.utc).isoformat(),
            })

    return signals
