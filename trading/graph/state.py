"""LangGraph TradingState — the single source of truth for each trading cycle."""
from __future__ import annotations

import operator
from typing import Annotated, TypedDict

from langgraph.graph.message import add_messages


class MarketSnapshot(TypedDict):
    ticker: str
    price: float
    volume: float
    volume_ratio: float
    ohlcv: dict          # latest bar {open, high, low, close, volume}
    indicators: dict     # RSI, MACD, BB, SMA, ATR, etc.


class TechnicalSignal(TypedDict):
    ticker: str
    signal: str          # "BUY" | "SELL" | "HOLD"
    strength: float      # 0.0–1.0
    reasons: list[str]
    source: str          # "llm" | "rule_based" | "st_agent" | "lt_agent"
    strategy: str        # "short_term" | "long_term" | "" (single-strategy mode)
    timestamp: str


class FundamentalSignal(TypedDict):
    ticker: str
    signal: str          # "BUY" | "SELL" | "HOLD"
    strength: float      # 0.0–1.0
    pe_ratio: float | None
    earnings_surprise: float | None
    macro_context: str
    reasons: list[str]
    timestamp: str


class RiskAssessment(TypedDict):
    ticker: str
    approved: bool
    position_size_usd: float
    qty: float
    stop_loss_price: float
    take_profit_price: float
    rejection_reason: str | None
    strategy: str        # "short_term" | "long_term" | "" (single-strategy mode)


class OrderResult(TypedDict):
    ticker: str
    order_id: str
    side: str            # "buy" | "sell"
    qty: float
    fill_price: float | None
    status: str          # "submitted" | "filled" | "rejected"
    timestamp: str


class TradingState(TypedDict):
    # ── Control flow ──────────────────────────────────────────────────────────
    run_id: str
    phase: str           # "market_open" | "pre_market" | "after_hours" | "closed"
    error: str | None

    # ── Market inputs ─────────────────────────────────────────────────────────
    market_snapshots: list[MarketSnapshot]
    portfolio_value: float
    current_positions: dict[str, dict]   # ticker -> Alpaca position dict
    available_cash: float
    peak_portfolio_value: float | None

    # ── Agent outputs (parallel-safe: both agents append to these lists) ──────
    technical_signals: Annotated[list[TechnicalSignal], operator.add]
    fundamental_signals: Annotated[list[FundamentalSignal], operator.add]

    # ── Decision layer ────────────────────────────────────────────────────────
    # Both strategies append to these lists concurrently (fan-in safe)
    consolidated_decisions: Annotated[list[dict], operator.add]
    risk_assessments: Annotated[list[RiskAssessment], operator.add]
    # Final balanced orders after portfolio_balancer resolves conflicts
    approved_orders: list[dict]
    balance_rejections: list[str]

    # ── Execution ─────────────────────────────────────────────────────────────
    orders_to_place: list[dict]
    order_results: Annotated[list[OrderResult], operator.add]

    # ── Audit / trace ─────────────────────────────────────────────────────────
    messages: Annotated[list, add_messages]  # LLM reasoning trace (last 20 kept)
    cycle_summary: str
    daily_loss_usd: float
