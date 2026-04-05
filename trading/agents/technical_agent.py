"""Technical Analysis Agent — Gemma 4 via Ollama with rule-based fallback.

Screens the full stock universe using RSI, MACD, Bollinger Bands, and
moving averages. Runs concurrently with the fundamental agent in LangGraph.
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timezone

from pydantic import BaseModel, Field, field_validator

from trading.llm import ollama_client
from trading.tools.indicators import rule_based_signal

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are an expert technical analyst specializing in US equities.
You will analyze market data and technical indicators for a stock and produce a trading signal.

Your signal must be:
- "BUY": price likely to increase in the next 3-10 trading days
- "SELL": consider reducing/exiting if held; short opportunity
- "HOLD": no clear directional edge

Base your analysis on:
1. RSI(14): oversold <30 (bullish), overbought >70 (bearish)
2. MACD histogram: positive trending (bullish), negative trending (bearish)
3. Bollinger Band %B: <0.1 near lower band (potential reversal), >0.9 near upper (overbought)
4. Price vs SMA50/SMA200: above both is bullish trend confirmation
5. Golden cross (SMA50 > SMA200): bullish; death cross: bearish
6. Volume ratio >1.5: confirms breakout/breakdown

Be concise and decisive. Strength is 0.0-1.0 (higher = more conviction)."""


class TechnicalSignalOutput(BaseModel):
    ticker: str
    signal: str = Field(pattern="^(BUY|SELL|HOLD)$")
    strength: float = Field(ge=0.0, le=1.0)
    reasons: list[str]

    @field_validator("reasons")
    @classmethod
    def limit_reasons(cls, v):
        return v[:5]  # max 5 reasons


def analyze_ticker(snapshot: dict) -> dict:
    """Analyze one ticker. Returns a TechnicalSignal dict."""
    ticker = snapshot["ticker"]
    ind = snapshot.get("indicators", {})

    user_prompt = f"""Analyze {ticker}:

Price: ${snapshot['price']:.2f}
Volume ratio (vs 20d avg): {snapshot.get('volume_ratio', 1.0):.2f}x

Technical Indicators:
- RSI(14): {ind.get('rsi_14', 'N/A')}
- MACD line: {ind.get('macd_line', 'N/A')}, Signal: {ind.get('macd_signal', 'N/A')}, Histogram: {ind.get('macd_hist', 'N/A')}
- Bollinger %B: {ind.get('bb_pct', 'N/A')} (upper: {ind.get('bb_upper', 'N/A')}, lower: {ind.get('bb_lower', 'N/A')})
- SMA20: {ind.get('sma_20', 'N/A')}, SMA50: {ind.get('sma_50', 'N/A')}, SMA200: {ind.get('sma_200', 'N/A')}
- ATR(14): {ind.get('atr_14', 'N/A')}
- Above SMA50: {ind.get('above_sma50', 'N/A')}, Above SMA200: {ind.get('above_sma200', 'N/A')}
- Golden cross: {ind.get('golden_cross', 'N/A')}

Return a JSON object with: ticker, signal (BUY/SELL/HOLD), strength (0.0-1.0), reasons (list of strings)."""

    result = ollama_client.call_structured(
        system_prompt=SYSTEM_PROMPT,
        user_prompt=user_prompt,
        output_schema=TechnicalSignalOutput,
        max_tokens=256,
        fallback_fn=rule_based_signal,
        fallback_input=snapshot,
    )

    return {
        "ticker": result.ticker,
        "signal": result.signal,
        "strength": result.strength,
        "reasons": result.reasons,
        "source": "llm",
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def run(snapshots: list[dict]) -> list[dict]:
    """Analyze all snapshots and return technical signals."""
    signals = []
    for snapshot in snapshots:
        try:
            sig = analyze_ticker(snapshot)
            signals.append(sig)
        except Exception as e:
            logger.error("Technical agent failed for %s: %s", snapshot.get("ticker"), e)
            # Always produce a signal — fall back to rule-based
            fallback = rule_based_signal(snapshot)
            fallback["timestamp"] = datetime.now(timezone.utc).isoformat()
            signals.append(fallback)
    return signals
