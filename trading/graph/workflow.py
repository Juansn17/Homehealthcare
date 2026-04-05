"""LangGraph workflow assembly — builds and compiles the trading graph."""
from __future__ import annotations

import logging
from pathlib import Path

from langgraph.graph import StateGraph, END
from langgraph.checkpoint.sqlite import SqliteSaver

from trading.config import settings
from trading.graph.state import TradingState
from trading.graph.nodes import (
    market_clock_check,
    data_fetch,
    technical_agent_node,
    fundamental_agent_node,
    orchestrator_agent_node,
    risk_check_node,
    execution_node,
    monitoring_node,
    halt_node,
    sleep_node,
    # dual-strategy nodes
    st_agent_node,
    lt_agent_node,
    st_risk_check_node,
    lt_risk_check_node,
    portfolio_balancer_node,
    dual_execution_node,
)
from trading.graph.edges import route_after_clock, route_after_risk, route_after_balancer

logger = logging.getLogger(__name__)


def build_graph(checkpointer=None):
    """Build and compile the trading LangGraph.

    Graph topology:
        START
          └─▶ market_clock_check
                ├─ open/pre_market ──▶ data_fetch
                └─ closed ──────────▶ sleep_node ──▶ END
                                                │
                              ┌─────────────────┘
                              ▼
                           data_fetch
                         ┌────┴────┐  (parallel fan-out)
                  technical_agent  fundamental_agent
                         └────┬────┘  (fan-in join)
                      orchestrator_agent
                              │
                         risk_check
                         ┌────┴────┐
                    approved    rejected/pre-market
                         │           │
                   execution_node  halt_node
                         │           │
                    monitoring_node  END
                         │
                         END
    """
    builder = StateGraph(TradingState)

    # ── Add nodes ─────────────────────────────────────────────────────────────
    builder.add_node("market_clock_check", market_clock_check)
    builder.add_node("data_fetch", data_fetch)
    builder.add_node("technical_agent", technical_agent_node)
    builder.add_node("fundamental_agent", fundamental_agent_node)
    builder.add_node("orchestrator_agent", orchestrator_agent_node)
    builder.add_node("risk_check", risk_check_node)
    builder.add_node("execution_node", execution_node)
    builder.add_node("monitoring_node", monitoring_node)
    builder.add_node("halt_node", halt_node)
    builder.add_node("sleep_node", sleep_node)

    # ── Entry point ───────────────────────────────────────────────────────────
    builder.set_entry_point("market_clock_check")

    # ── Conditional routing after clock check ────────────────────────────────
    builder.add_conditional_edges(
        "market_clock_check",
        route_after_clock,
        {
            "data_fetch": "data_fetch",
            "sleep_node": "sleep_node",
        },
    )

    # ── Parallel fan-out: data_fetch → both agents ────────────────────────────
    builder.add_edge("data_fetch", "technical_agent")
    builder.add_edge("data_fetch", "fundamental_agent")

    # ── Fan-in: both agents → orchestrator (LangGraph waits for both) ─────────
    builder.add_edge("technical_agent", "orchestrator_agent")
    builder.add_edge("fundamental_agent", "orchestrator_agent")

    # ── Linear: orchestrator → risk ───────────────────────────────────────────
    builder.add_edge("orchestrator_agent", "risk_check")

    # ── Conditional: risk → execution or halt ────────────────────────────────
    builder.add_conditional_edges(
        "risk_check",
        route_after_risk,
        {
            "execution_node": "execution_node",
            "halt_node": "halt_node",
        },
    )

    # ── Linear: execution → monitoring → END ─────────────────────────────────
    builder.add_edge("execution_node", "monitoring_node")
    builder.add_edge("monitoring_node", END)
    builder.add_edge("halt_node", END)
    builder.add_edge("sleep_node", END)

    # ── Compile with optional checkpointer ───────────────────────────────────
    return builder.compile(checkpointer=checkpointer)


def get_compiled_graph():
    """Return a compiled single-strategy graph with SQLite persistence."""
    db_path = Path(settings.DB_PATH)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    checkpointer = SqliteSaver.from_conn_string(str(db_path))
    graph = build_graph(checkpointer=checkpointer)
    logger.info("LangGraph compiled with SQLite checkpointer at %s", db_path)
    return graph


def build_dual_strategy_graph(checkpointer=None):
    """Build the dual-strategy (short-term + long-term parallel) graph.

    Topology:
        START
          └─▶ market_clock_check
                ├─ closed ──────────────────────────────▶ sleep_node ──▶ END
                └─ open / pre_market ──▶ data_fetch
                                              │
                            ┌─────────────────┴──────────────────┐
                            ▼                                    ▼
                     st_agent_node                        lt_agent_node
                  (Gemma 4, intraday,                  (Claude Sonnet, daily,
                   runs intraday only)                  runs always)
                            │                                    │
                     st_risk_check                        lt_risk_check
                   (ATR×1.5, 1% risk,               (ATR×3.0, 2% risk,
                    5 max positions)                  3 max positions)
                            └─────────────────┬──────────────────┘
                                              ▼
                                   portfolio_balancer
                               (conflict resolution,
                                capital allocation,
                                combined risk guard)
                                              │
                            ┌─────────────────┴──────────────────┐
                         approved                            all rejected
                            │                                    │
                   dual_execution_node                      halt_node ──▶ END
                   (bracket orders ST+LT)
                            │
                    monitoring_node ──▶ END
    """
    builder = StateGraph(TradingState)

    # ── Add nodes ─────────────────────────────────────────────────────────────
    builder.add_node("market_clock_check", market_clock_check)
    builder.add_node("data_fetch", data_fetch)
    builder.add_node("st_agent_node", st_agent_node)
    builder.add_node("lt_agent_node", lt_agent_node)
    builder.add_node("st_risk_check", st_risk_check_node)
    builder.add_node("lt_risk_check", lt_risk_check_node)
    builder.add_node("portfolio_balancer", portfolio_balancer_node)
    builder.add_node("dual_execution_node", dual_execution_node)
    builder.add_node("monitoring_node", monitoring_node)
    builder.add_node("halt_node", halt_node)
    builder.add_node("sleep_node", sleep_node)

    # ── Entry ─────────────────────────────────────────────────────────────────
    builder.set_entry_point("market_clock_check")

    # ── Clock → data or sleep ────────────────────────────────────────────────
    builder.add_conditional_edges(
        "market_clock_check",
        route_after_clock,
        {"data_fetch": "data_fetch", "sleep_node": "sleep_node"},
    )

    # ── data_fetch → parallel fan-out to both agents ─────────────────────────
    builder.add_edge("data_fetch", "st_agent_node")
    builder.add_edge("data_fetch", "lt_agent_node")

    # ── Each agent feeds its own risk check (still parallel) ─────────────────
    builder.add_edge("st_agent_node", "st_risk_check")
    builder.add_edge("lt_agent_node", "lt_risk_check")

    # ── Both risk checks fan-in at portfolio_balancer ─────────────────────────
    builder.add_edge("st_risk_check", "portfolio_balancer")
    builder.add_edge("lt_risk_check", "portfolio_balancer")

    # ── Balancer → execution or halt ─────────────────────────────────────────
    builder.add_conditional_edges(
        "portfolio_balancer",
        route_after_balancer,
        {"dual_execution_node": "dual_execution_node", "halt_node": "halt_node"},
    )

    # ── Execution → monitoring → END ─────────────────────────────────────────
    builder.add_edge("dual_execution_node", "monitoring_node")
    builder.add_edge("monitoring_node", END)
    builder.add_edge("halt_node", END)
    builder.add_edge("sleep_node", END)

    return builder.compile(checkpointer=checkpointer)


def get_compiled_dual_graph():
    """Return a compiled dual-strategy graph with SQLite persistence."""
    db_path = Path(settings.DB_PATH).parent / "trading_dual_state.db"
    db_path.parent.mkdir(parents=True, exist_ok=True)
    checkpointer = SqliteSaver.from_conn_string(str(db_path))
    graph = build_dual_strategy_graph(checkpointer=checkpointer)
    logger.info("Dual-strategy LangGraph compiled at %s", db_path)
    return graph
