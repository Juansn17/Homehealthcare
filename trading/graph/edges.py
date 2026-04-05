"""Conditional edge routing logic for the trading graph."""
from __future__ import annotations

from trading.graph.state import TradingState


def route_after_clock(state: TradingState) -> str:
    """Route based on market phase."""
    phase = state.get("phase", "closed")
    if phase == "market_open":
        return "data_fetch"
    if phase == "pre_market":
        return "data_fetch"   # pre-market: fetch data but execution is blocked by risk
    return "sleep_node"


def route_after_risk(state: TradingState) -> str:
    """Route to execution if any trades are approved, else halt."""
    if state.get("error"):
        return "halt_node"
    approved = [r for r in state.get("risk_assessments", []) if r.get("approved")]
    if not approved:
        return "halt_node"
    # Block execution in pre-market
    if state.get("phase") == "pre_market":
        return "halt_node"
    return "execution_node"
