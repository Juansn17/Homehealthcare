"""Risk management — pure Python, no LLM override possible.

All hard limits enforced here. This module runs as the risk_check node
in the LangGraph graph before any order is submitted.
"""
from __future__ import annotations

import logging
import sqlite3
from pathlib import Path

import numpy as np

from trading.config import settings
from trading.config.universe import SECTOR_MAP

logger = logging.getLogger(__name__)

# ── Constants (overridable via settings) ────────────────────────────────────
MIN_STOP_DISTANCE = 0.01   # 1% minimum stop-loss distance
MAX_STOP_DISTANCE = 0.08   # 8% maximum stop-loss distance
MIN_RR_RATIO = 2.0         # minimum reward/risk ratio
MAX_CORRELATION = 0.85     # reject if corr with existing positions > this
MAX_SECTOR_PCT = 0.30      # max 30% of portfolio in one sector


def assess_trade(
    ticker: str,
    side: str,
    entry_price: float,
    atr_14: float,
    portfolio_value: float,
    current_positions: dict[str, dict],  # ticker -> {qty, avg_entry_price, market_value, side}
    available_cash: float,
    peak_portfolio_value: float | None = None,
    daily_loss_today: float = 0.0,
) -> dict:
    """Evaluate a proposed trade and return a RiskAssessment dict.

    Returns:
        {
            "ticker": str,
            "approved": bool,
            "position_size_usd": float,
            "qty": float,
            "stop_loss_price": float,
            "take_profit_price": float,
            "rejection_reason": str | None,
        }
    """
    def _reject(reason: str) -> dict:
        logger.info("Risk REJECTED %s: %s", ticker, reason)
        return {
            "ticker": ticker,
            "approved": False,
            "position_size_usd": 0.0,
            "qty": 0.0,
            "stop_loss_price": 0.0,
            "take_profit_price": 0.0,
            "rejection_reason": reason,
        }

    # ── 1. Max drawdown guard ─────────────────────────────────────────────────
    if peak_portfolio_value and peak_portfolio_value > 0:
        drawdown = (peak_portfolio_value - portfolio_value) / peak_portfolio_value
        if drawdown > settings.MAX_DRAWDOWN_PCT:
            return _reject(
                f"Portfolio drawdown {drawdown:.1%} exceeds limit {settings.MAX_DRAWDOWN_PCT:.1%}"
            )

    # ── 2. Daily loss limit ───────────────────────────────────────────────────
    if settings.IS_LIVE and daily_loss_today > settings.MAX_DAILY_LOSS_USD:
        return _reject(f"Daily loss limit ${settings.MAX_DAILY_LOSS_USD:.0f} reached")

    # ── 3. Stop-loss calculation (ATR-based) ──────────────────────────────────
    stop_distance = atr_14 * 2.0
    stop_distance_pct = stop_distance / entry_price

    if stop_distance_pct < MIN_STOP_DISTANCE:
        return _reject(f"Stop-loss distance {stop_distance_pct:.2%} < minimum {MIN_STOP_DISTANCE:.2%}")
    if stop_distance_pct > MAX_STOP_DISTANCE:
        return _reject(f"Stop-loss distance {stop_distance_pct:.2%} > maximum {MAX_STOP_DISTANCE:.2%}")

    if side.lower() == "buy":
        stop_loss_price = entry_price - stop_distance
        take_profit_price = entry_price + stop_distance * MIN_RR_RATIO
    else:
        stop_loss_price = entry_price + stop_distance
        take_profit_price = entry_price - stop_distance * MIN_RR_RATIO

    # ── 4. Position sizing ────────────────────────────────────────────────────
    risk_based_size = (portfolio_value * settings.MAX_PORTFOLIO_RISK_PCT) / stop_distance_pct
    max_size = portfolio_value * settings.POSITION_SIZE_PCT
    if settings.IS_LIVE:
        max_size *= settings.LIVE_POSITION_SCALE

    position_size_usd = min(risk_based_size, max_size, available_cash * 0.95)
    if position_size_usd < 10:
        return _reject(f"Position size ${position_size_usd:.2f} too small (< $10)")

    qty = position_size_usd / entry_price

    # ── 5. Duplicate position check ───────────────────────────────────────────
    if ticker in current_positions:
        return _reject(f"Already holding position in {ticker}")

    # ── 6. Sector concentration check ────────────────────────────────────────
    sector = SECTOR_MAP.get(ticker)
    if sector:
        sector_exposure = sum(
            pos.get("market_value", 0.0)
            for t, pos in current_positions.items()
            if SECTOR_MAP.get(t) == sector
        )
        sector_pct = (sector_exposure + position_size_usd) / portfolio_value
        if sector_pct > MAX_SECTOR_PCT:
            return _reject(
                f"Sector '{sector}' concentration {sector_pct:.1%} would exceed {MAX_SECTOR_PCT:.1%}"
            )

    logger.info(
        "Risk APPROVED %s: size=$%.0f, qty=%.2f, stop=%.2f, tp=%.2f",
        ticker, position_size_usd, qty, stop_loss_price, take_profit_price,
    )
    return {
        "ticker": ticker,
        "approved": True,
        "position_size_usd": round(position_size_usd, 2),
        "qty": round(qty, 4),
        "stop_loss_price": round(stop_loss_price, 2),
        "take_profit_price": round(take_profit_price, 2),
        "rejection_reason": None,
    }


def update_peak_portfolio(portfolio_value: float, db_path: str) -> float:
    """Persist and return the peak portfolio value in SQLite."""
    Path(db_path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(db_path)
    try:
        conn.execute(
            "CREATE TABLE IF NOT EXISTS risk_state (key TEXT PRIMARY KEY, value REAL)"
        )
        row = conn.execute(
            "SELECT value FROM risk_state WHERE key='peak_portfolio_value'"
        ).fetchone()
        peak = row[0] if row else portfolio_value
        if portfolio_value > peak:
            peak = portfolio_value
            conn.execute(
                "INSERT OR REPLACE INTO risk_state (key, value) VALUES ('peak_portfolio_value', ?)",
                (peak,),
            )
            conn.commit()
        return peak
    finally:
        conn.close()


def compute_portfolio_drawdown(portfolio_value: float, db_path: str) -> float:
    """Return current drawdown from peak as a fraction (0.0 to 1.0)."""
    peak = update_peak_portfolio(portfolio_value, db_path)
    if peak <= 0:
        return 0.0
    return max(0.0, (peak - portfolio_value) / peak)
