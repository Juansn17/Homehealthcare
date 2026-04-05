"""Orchestrator Agent — Claude Sonnet.

Consolidates technical and fundamental signals, ranks opportunities,
and produces final BUY/SELL/HOLD decisions with confidence scores.
"""
from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from typing import Any

from langchain_core.messages import AIMessage, HumanMessage
from pydantic import BaseModel, Field, field_validator

from trading.llm import claude_client

logger = logging.getLogger(__name__)

SYSTEM_PROMPT = """You are an algorithmic trading decision agent managing a US equity portfolio.
You have received signals from two specialized analysts:
1. Technical Analyst: signals based on price action and indicators
2. Fundamental Analyst: signals based on valuation, earnings, and macro context

Your job:
1. Identify the top 3-5 highest conviction opportunities (both BUY and SELL/reduce signals)
2. For each, produce a final action (BUY/SELL/HOLD) and confidence score (0-100)
3. Weigh both signals equally; if they disagree, explain why one takes precedence
4. Flag any concentration risk or macro headwinds that should limit position sizes
5. Reject opportunities where analysts disagree strongly (one BUY, one SELL with both > 0.6 strength)

Only produce decisions for tickers NOT already in current_positions (for BUY signals).
For tickers IN current_positions, you may produce SELL signals if conditions have deteriorated.

Output ONLY the JSON array of decisions — no prose."""


class TradeDecision(BaseModel):
    ticker: str
    action: str = Field(pattern="^(BUY|SELL|HOLD)$")
    confidence: int = Field(ge=0, le=100)
    ta_signal: str = ""
    fa_signal: str = ""
    reasoning: str
    concerns: str = ""


class OrchestratorOutput(BaseModel):
    decisions: list[TradeDecision]

    @field_validator("decisions")
    @classmethod
    def limit_decisions(cls, v):
        return v[:5]  # max 5 decisions per cycle


def run(
    technical_signals: list[dict],
    fundamental_signals: list[dict],
    portfolio_value: float,
    current_positions: dict[str, Any],
    available_cash: float,
) -> tuple[list[dict], list[Any]]:
    """Consolidate signals and return (decisions, messages)."""
    # Build lookup dicts for easy access
    ta_by_ticker = {s["ticker"]: s for s in technical_signals}
    fa_by_ticker = {s["ticker"]: s for s in fundamental_signals}

    # All tickers that appear in either signal set
    all_tickers = set(ta_by_ticker.keys()) | set(fa_by_ticker.keys())

    # Build combined signal summary for prompt
    signal_lines = []
    for ticker in sorted(all_tickers):
        ta = ta_by_ticker.get(ticker)
        fa = fa_by_ticker.get(ticker)
        ta_str = f"TA:{ta['signal']}({ta['strength']:.2f})" if ta else "TA:N/A"
        fa_str = f"FA:{fa['signal']}({fa['strength']:.2f})" if fa else "FA:N/A"
        ta_reasons = "; ".join(ta.get("reasons", [])[:2]) if ta else ""
        fa_reasons = "; ".join(fa.get("reasons", [])[:2]) if fa else ""
        line = f"- {ticker}: {ta_str} | {fa_str}"
        if ta_reasons:
            line += f" | TA reasons: {ta_reasons}"
        if fa_reasons:
            line += f" | FA reasons: {fa_reasons}"
        signal_lines.append(line)

    signals_text = "\n".join(signal_lines) or "No signals available"

    positions_text = (
        json.dumps(
            {t: {"qty": p.get("qty"), "avg_entry": p.get("avg_entry_price")}
             for t, p in current_positions.items()},
            indent=2,
        )
        if current_positions
        else "No current positions"
    )

    user_prompt = f"""Current portfolio:
- Total value: ${portfolio_value:,.0f}
- Available cash: ${available_cash:,.0f}
- Open positions: {positions_text}

Combined analyst signals ({len(all_tickers)} tickers):
{signals_text}

Produce a JSON object with key "decisions" containing an array of up to 5 trade decisions.
Each decision: ticker, action (BUY/SELL/HOLD), confidence (0-100), ta_signal, fa_signal, reasoning (1-2 sentences), concerns."""

    output = claude_client.call_structured(
        system_prompt=SYSTEM_PROMPT,
        user_prompt=user_prompt,
        output_schema=OrchestratorOutput,
        max_tokens=1024,
    )

    decisions = [
        {
            "ticker": d.ticker,
            "action": d.action,
            "confidence": d.confidence,
            "ta_signal": d.ta_signal,
            "fa_signal": d.fa_signal,
            "reasoning": d.reasoning,
            "concerns": d.concerns,
            "timestamp": datetime.now(timezone.utc).isoformat(),
        }
        for d in output.decisions
        if d.action != "HOLD"  # only actionable decisions flow to risk layer
    ]

    messages = [
        HumanMessage(content=user_prompt),
        AIMessage(content=json.dumps({"decisions": [d.dict() for d in output.decisions]})),
    ]

    logger.info(
        "Orchestrator decisions: %s",
        [(d["ticker"], d["action"], d["confidence"]) for d in decisions],
    )
    return decisions, messages
