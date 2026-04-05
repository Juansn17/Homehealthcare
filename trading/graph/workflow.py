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
)
from trading.graph.edges import route_after_clock, route_after_risk

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
    """Return a compiled graph with SQLite persistence."""
    db_path = Path(settings.DB_PATH)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    checkpointer = SqliteSaver.from_conn_string(str(db_path))
    graph = build_graph(checkpointer=checkpointer)
    logger.info("LangGraph compiled with SQLite checkpointer at %s", db_path)
    return graph
