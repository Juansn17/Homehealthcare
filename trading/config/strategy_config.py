"""Strategy configuration — parameters that distinguish Short-Term from Long-Term.

Usage:
    from trading.config.strategy_config import SHORT_TERM, LONG_TERM, StrategyConfig
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class StrategyConfig:
    name: str                    # "short_term" | "long_term"
    label: str                   # human-readable

    # ── Time horizon ──────────────────────────────────────────────────────────
    horizon_days_min: int        # minimum expected hold (days)
    horizon_days_max: int        # maximum expected hold (days)
    intraday_interval: str       # interval for signal bars: "15m", "1h", "1d"

    # ── Risk parameters (override settings.py defaults per strategy) ─────────
    stop_atr_multiplier: float   # ATR multiplier for stop-loss placement
    reward_risk_ratio: float     # minimum take-profit / stop-loss distance ratio
    position_risk_pct: float     # fraction of portfolio risked per trade
    max_open_positions: int      # maximum simultaneous positions for this strategy

    # ── Capital allocation ────────────────────────────────────────────────────
    portfolio_allocation_pct: float  # fraction of total portfolio allocated here

    # ── Signal emphasis ───────────────────────────────────────────────────────
    ta_weight: float             # weight given to technical signals (0–1)
    fa_weight: float             # weight given to fundamental signals (0–1)

    # ── LLM routing ───────────────────────────────────────────────────────────
    use_local_llm: bool          # True → Gemma (Ollama); False → Claude API

    # ── Scheduling ───────────────────────────────────────────────────────────
    runs_intraday: bool          # True → runs every REBALANCE_INTERVAL_MINUTES
    runs_eod: bool               # True → runs once at pre/post market


# ── Strategy definitions ─────────────────────────────────────────────────────

SHORT_TERM = StrategyConfig(
    name="short_term",
    label="Short-Term Momentum",
    # Horizon: 1–5 trading days
    horizon_days_min=1,
    horizon_days_max=5,
    intraday_interval="15m",
    # Tighter risk — smaller moves, quicker exits
    stop_atr_multiplier=1.5,
    reward_risk_ratio=1.5,
    position_risk_pct=0.01,       # 1% of portfolio per trade
    max_open_positions=5,
    portfolio_allocation_pct=0.35,
    # Heavily technical — intraday momentum, RSI extremes, MACD cross
    ta_weight=0.80,
    fa_weight=0.20,
    use_local_llm=True,            # Gemma 4 for fast intraday screening
    runs_intraday=True,
    runs_eod=False,
)

LONG_TERM = StrategyConfig(
    name="long_term",
    label="Long-Term Fundamental",
    # Horizon: 2–8 weeks (10–40 trading days)
    horizon_days_min=10,
    horizon_days_max=40,
    intraday_interval="1d",
    # Wider risk — ride larger moves, fewer false stops
    stop_atr_multiplier=3.0,
    reward_risk_ratio=3.0,
    position_risk_pct=0.02,       # 2% of portfolio per trade
    max_open_positions=3,
    portfolio_allocation_pct=0.65,
    # Fundamentals-first — earnings momentum, sector rotation, macro
    ta_weight=0.30,
    fa_weight=0.70,
    use_local_llm=False,           # Claude Sonnet for deep reasoning
    runs_intraday=False,
    runs_eod=True,
)

ALL_STRATEGIES: dict[str, StrategyConfig] = {
    SHORT_TERM.name: SHORT_TERM,
    LONG_TERM.name: LONG_TERM,
}
