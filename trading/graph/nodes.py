"""LangGraph node functions — each returns a partial TradingState update."""
from __future__ import annotations

import logging
import uuid
from datetime import datetime, timezone

from trading.config import settings
from trading.config.universe import STOCK_UNIVERSE
from trading.graph.state import TradingState
from trading.risk import manager as risk_mgr
from trading.tools import market_data, alpaca_client

logger = logging.getLogger(__name__)


def market_clock_check(state: TradingState) -> dict:
    """Determine if the market is open and set the phase."""
    from trading.scheduler.market_clock import get_market_phase
    phase = get_market_phase()
    logger.info("Market clock check: phase=%s", phase)
    return {
        "run_id": state.get("run_id") or str(uuid.uuid4()),
        "phase": phase,
        "error": None,
        "technical_signals": [],
        "fundamental_signals": [],
        "consolidated_decisions": [],
        "risk_assessments": [],
        "orders_to_place": [],
        "order_results": [],
        "messages": [],
        "cycle_summary": "",
        "daily_loss_usd": state.get("daily_loss_usd", 0.0),
    }


def data_fetch(state: TradingState) -> dict:
    """Fetch market data, indicators, and account state for all tickers."""
    logger.info("Fetching market data for %d tickers", len(STOCK_UNIVERSE))
    snapshots = market_data.build_snapshots(STOCK_UNIVERSE)

    account = alpaca_client.get_account()
    positions = alpaca_client.get_positions()

    peak = risk_mgr.update_peak_portfolio(account["portfolio_value"], settings.DB_PATH)

    return {
        "market_snapshots": snapshots,
        "portfolio_value": account["portfolio_value"],
        "current_positions": positions,
        "available_cash": account["buying_power"],
        "peak_portfolio_value": peak,
    }


def technical_agent_node(state: TradingState) -> dict:
    """Run the technical analysis agent (Gemma 4 via Ollama)."""
    from trading.agents.technical_agent import run as run_ta
    signals = run_ta(state["market_snapshots"])
    logger.info("Technical agent produced %d signals", len(signals))
    return {"technical_signals": signals}


def fundamental_agent_node(state: TradingState) -> dict:
    """Run the fundamental/macro analysis agent (Claude Sonnet)."""
    from trading.agents.fundamental_agent import run as run_fa
    # Only analyze tickers that have a non-HOLD technical signal (cost optimization)
    ta_tickers = {
        s["ticker"] for s in state.get("technical_signals", [])
        if s.get("signal") != "HOLD"
    }
    snapshots = [s for s in state["market_snapshots"] if s["ticker"] in ta_tickers]
    if not snapshots:
        snapshots = state["market_snapshots"][:5]  # fallback: analyze top 5

    signals = run_fa(snapshots)
    logger.info("Fundamental agent produced %d signals", len(signals))
    return {"fundamental_signals": signals}


def orchestrator_agent_node(state: TradingState) -> dict:
    """Consolidate signals and produce final trading decisions (Claude Sonnet)."""
    from trading.agents.orchestrator_agent import run as run_orch
    decisions, messages = run_orch(
        technical_signals=state.get("technical_signals", []),
        fundamental_signals=state.get("fundamental_signals", []),
        portfolio_value=state["portfolio_value"],
        current_positions=state["current_positions"],
        available_cash=state["available_cash"],
    )
    logger.info("Orchestrator produced %d decisions", len(decisions))
    return {
        "consolidated_decisions": decisions,
        "messages": messages,
    }


def risk_check_node(state: TradingState) -> dict:
    """Apply all risk rules to proposed decisions. Pure Python — no LLM."""
    assessments = []
    positions = state.get("current_positions", {})
    portfolio_value = state["portfolio_value"]
    available_cash = state["available_cash"]
    peak = state.get("peak_portfolio_value") or portfolio_value
    daily_loss = state.get("daily_loss_usd", 0.0)

    for decision in state.get("consolidated_decisions", []):
        ticker = decision["ticker"]
        action = decision.get("action", "HOLD")
        if action == "HOLD":
            continue

        snapshot = next(
            (s for s in state["market_snapshots"] if s["ticker"] == ticker), None
        )
        if not snapshot:
            continue

        entry_price = snapshot["price"]
        atr = snapshot["indicators"].get("atr_14")
        if not atr:
            logger.warning("No ATR for %s, skipping risk check", ticker)
            continue

        assessment = risk_mgr.assess_trade(
            ticker=ticker,
            side=action.lower(),
            entry_price=entry_price,
            atr_14=atr,
            portfolio_value=portfolio_value,
            current_positions=positions,
            available_cash=available_cash,
            peak_portfolio_value=peak,
            daily_loss_today=daily_loss,
        )
        assessments.append(assessment)

    return {"risk_assessments": assessments}


def execution_node(state: TradingState) -> dict:
    """Submit approved bracket orders to Alpaca."""
    from trading.execution.order_manager import submit_approved_orders
    approved = [r for r in state.get("risk_assessments", []) if r["approved"]]
    results = submit_approved_orders(approved, state["market_snapshots"])
    logger.info("Execution: submitted %d orders", len(results))
    return {"order_results": results}


def monitoring_node(state: TradingState) -> dict:
    """Check fills, compute cycle summary, trim message history."""
    from trading.execution.order_manager import check_fills
    updated_results = check_fills(state.get("order_results", []))

    # Trim messages to last 20 to prevent state bloat
    messages = state.get("messages", [])[-20:]

    filled = [r for r in updated_results if r["status"] == "filled"]
    rejected = [r for r in updated_results if r["status"] == "rejected"]
    summary = (
        f"Cycle {state.get('run_id', '')[:8]} | "
        f"phase={state.get('phase')} | "
        f"decisions={len(state.get('consolidated_decisions', []))} | "
        f"orders_submitted={len(updated_results)} | "
        f"filled={len(filled)} | rejected={len(rejected)} | "
        f"portfolio=${state.get('portfolio_value', 0):,.0f}"
    )
    logger.info("Cycle summary: %s", summary)

    return {
        "order_results": updated_results,
        "messages": messages,
        "cycle_summary": summary,
    }


def halt_node(state: TradingState) -> dict:
    """Log halt reason and skip execution."""
    reasons = [
        r.get("rejection_reason", "")
        for r in state.get("risk_assessments", [])
        if not r.get("approved")
    ]
    reason_str = "; ".join(reasons) if reasons else "no approved trades"
    logger.warning("Trading halted this cycle: %s", reason_str)
    return {
        "cycle_summary": f"HALTED: {reason_str}",
        "orders_to_place": [],
    }


def sleep_node(state: TradingState) -> dict:
    """Market is closed — nothing to do."""
    logger.info("Market closed (phase=%s), skipping cycle", state.get("phase"))
    return {"cycle_summary": f"Market closed: {state.get('phase')}"}
