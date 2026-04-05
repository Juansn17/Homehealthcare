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
        "approved_orders": [],
        "balance_rejections": [],
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


# ── Dual-strategy nodes ───────────────────────────────────────────────────────

def st_agent_node(state: TradingState) -> dict:
    """Short-term momentum agent (Gemma 4, intraday 15m).

    Runs in parallel with lt_agent_node. Appends to consolidated_decisions.
    Only runs when market is open (intraday execution needed).
    """
    if state.get("phase") not in ("market_open",):
        logger.info("ST agent skipped (phase=%s)", state.get("phase"))
        return {"consolidated_decisions": []}

    from trading.agents.st_agent import run as run_st
    signals = run_st(state["market_snapshots"])

    # Convert ST signals directly to decisions (ST agent is opinionated enough)
    decisions = [
        {
            "ticker": s["ticker"],
            "action": s["signal"],
            "confidence": int(s["strength"] * 100),
            "strategy": s["strategy"],
            "ta_signal": s["signal"],
            "fa_signal": "",
            "reasoning": "; ".join(s.get("reasons", [])[:2]),
            "concerns": "" if s.get("trend_aligned") else "trend misalignment",
            "timestamp": s["timestamp"],
        }
        for s in signals
        if s.get("signal") != "HOLD" and s.get("trend_aligned", True)
    ]
    logger.info("ST agent: %d actionable decisions", len(decisions))
    return {"consolidated_decisions": decisions}


def lt_agent_node(state: TradingState) -> dict:
    """Long-term fundamental agent (Claude Sonnet, daily+weekly).

    Runs in parallel with st_agent_node. Appends to consolidated_decisions.
    Runs in both market_open and pre_market phases (no intraday timing needed).
    """
    from trading.agents.lt_agent import run as run_lt
    # LT uses full universe — fundamental screening is independent of intraday
    signals = run_lt(state["market_snapshots"])

    decisions = [
        {
            "ticker": s["ticker"],
            "action": s["signal"],
            "confidence": int(s["strength"] * 100),
            "strategy": s["strategy"],
            "ta_signal": "",
            "fa_signal": s["signal"],
            "reasoning": "; ".join(s.get("reasons", [])[:2]),
            "concerns": "earnings risk" if s.get("earnings_risk") else s.get("macro_stance", ""),
            "timestamp": s["timestamp"],
        }
        for s in signals
        if s.get("signal") != "HOLD" and not s.get("earnings_risk", False)
    ]
    logger.info("LT agent: %d actionable decisions", len(decisions))
    return {"consolidated_decisions": decisions}


def st_risk_check_node(state: TradingState) -> dict:
    """Risk check for short-term decisions with ST-specific parameters."""
    from trading.config.strategy_config import SHORT_TERM
    return _strategy_risk_check(state, SHORT_TERM)


def lt_risk_check_node(state: TradingState) -> dict:
    """Risk check for long-term decisions with LT-specific parameters."""
    from trading.config.strategy_config import LONG_TERM
    return _strategy_risk_check(state, LONG_TERM)


def _strategy_risk_check(state: TradingState, strategy) -> dict:
    """Shared risk check logic parameterized by StrategyConfig."""
    from trading.config.strategy_config import StrategyConfig
    assessments = []
    positions = state.get("current_positions", {})
    portfolio_value = state["portfolio_value"]
    available_cash = state["available_cash"]
    peak = state.get("peak_portfolio_value") or portfolio_value
    daily_loss = state.get("daily_loss_usd", 0.0)

    strategy_decisions = [
        d for d in state.get("consolidated_decisions", [])
        if d.get("strategy") == strategy.name and d.get("action", "HOLD") != "HOLD"
    ]

    for decision in strategy_decisions:
        ticker = decision["ticker"]
        snapshot = next(
            (s for s in state["market_snapshots"] if s["ticker"] == ticker), None
        )
        if not snapshot:
            continue

        entry_price = snapshot["price"]
        atr = snapshot["indicators"].get("atr_14")
        if not atr:
            logger.warning("No ATR for %s (%s), skipping", ticker, strategy.name)
            continue

        # Use strategy-specific ATR multiplier and R:R ratio
        stop_distance = atr * strategy.stop_atr_multiplier
        stop_distance_pct = stop_distance / entry_price

        if stop_distance_pct < 0.01 or stop_distance_pct > 0.10:
            assessments.append({
                "ticker": ticker, "approved": False,
                "position_size_usd": 0.0, "qty": 0.0,
                "stop_loss_price": 0.0, "take_profit_price": 0.0,
                "rejection_reason": f"Stop distance {stop_distance_pct:.2%} out of range for {strategy.name}",
                "strategy": strategy.name,
            })
            continue

        # Strategy-specific position sizing
        strategy_budget = available_cash * strategy.portfolio_allocation_pct
        risk_based = (portfolio_value * strategy.position_risk_pct) / stop_distance_pct
        max_size = portfolio_value * 0.05 * (strategy.portfolio_allocation_pct * 2)
        if settings.IS_LIVE:
            max_size *= settings.LIVE_POSITION_SCALE
        position_size_usd = min(risk_based, max_size, strategy_budget * 0.3)

        if position_size_usd < 10:
            assessments.append({
                "ticker": ticker, "approved": False,
                "position_size_usd": 0.0, "qty": 0.0,
                "stop_loss_price": 0.0, "take_profit_price": 0.0,
                "rejection_reason": f"Position too small ${position_size_usd:.0f}",
                "strategy": strategy.name,
            })
            continue

        stop_loss_price = round(entry_price - stop_distance, 2)
        take_profit_price = round(entry_price + stop_distance * strategy.reward_risk_ratio, 2)
        qty = round(position_size_usd / entry_price, 4)

        # Drawdown guard applies to all strategies
        if peak > 0 and (peak - portfolio_value) / peak > settings.MAX_DRAWDOWN_PCT:
            assessments.append({
                "ticker": ticker, "approved": False,
                "position_size_usd": 0.0, "qty": 0.0,
                "stop_loss_price": 0.0, "take_profit_price": 0.0,
                "rejection_reason": "Portfolio drawdown limit reached",
                "strategy": strategy.name,
            })
            continue

        # Skip if already holding
        if ticker in positions:
            assessments.append({
                "ticker": ticker, "approved": False,
                "position_size_usd": 0.0, "qty": 0.0,
                "stop_loss_price": 0.0, "take_profit_price": 0.0,
                "rejection_reason": f"Already holding {ticker}",
                "strategy": strategy.name,
            })
            continue

        logger.info(
            "Risk APPROVED [%s] %s: size=$%.0f stop=%.2f tp=%.2f",
            strategy.name, ticker, position_size_usd, stop_loss_price, take_profit_price,
        )
        assessments.append({
            "ticker": ticker,
            "approved": True,
            "position_size_usd": round(position_size_usd, 2),
            "qty": qty,
            "stop_loss_price": stop_loss_price,
            "take_profit_price": take_profit_price,
            "rejection_reason": None,
            "strategy": strategy.name,
        })

    return {"risk_assessments": assessments}


def portfolio_balancer_node(state: TradingState) -> dict:
    """Resolve conflicts between ST and LT, enforce combined risk limits.

    Runs after both st_risk_check and lt_risk_check complete (fan-in).
    """
    from trading.graph.portfolio_balancer import balance

    approved_orders, rejections = balance(
        risk_assessments=state.get("risk_assessments", []),
        portfolio_value=state["portfolio_value"],
        available_cash=state["available_cash"],
        current_positions=state.get("current_positions", {}),
    )

    if rejections:
        for r in rejections:
            logger.info("Balancer rejected: %s", r)

    logger.info(
        "Balancer: %d orders approved (%d ST, %d LT)",
        len(approved_orders),
        sum(1 for o in approved_orders if o.get("strategy") == "short_term"),
        sum(1 for o in approved_orders if o.get("strategy") == "long_term"),
    )
    return {
        "approved_orders": approved_orders,
        "balance_rejections": rejections,
    }


def dual_execution_node(state: TradingState) -> dict:
    """Submit approved orders from the portfolio balancer."""
    from trading.execution.order_manager import submit_approved_orders
    approved = state.get("approved_orders", [])
    results = submit_approved_orders(approved, state["market_snapshots"])
    logger.info("Dual execution: submitted %d orders", len(results))
    return {"order_results": results}
