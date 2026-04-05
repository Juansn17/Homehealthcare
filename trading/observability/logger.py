"""Structured logging, trade audit trail, and performance tracking."""
from __future__ import annotations

import csv
import logging
import os
from datetime import datetime, timezone
from pathlib import Path

import structlog

from trading.config import settings

_configured = False


def configure_logging():
    """Configure structlog for JSON structured logging. Call once at startup."""
    global _configured
    if _configured:
        return

    log_dir = Path(settings.LOG_PATH)
    log_dir.mkdir(parents=True, exist_ok=True)

    # Standard library logging to file
    log_file = log_dir / "app.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(message)s",
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(),  # also print to console
        ],
    )

    structlog.configure(
        processors=[
            structlog.stdlib.filter_by_level,
            structlog.stdlib.add_logger_name,
            structlog.stdlib.add_log_level,
            structlog.processors.TimeStamper(fmt="iso"),
            structlog.processors.StackInfoRenderer(),
            structlog.processors.format_exc_info,
            structlog.processors.JSONRenderer(),
        ],
        wrapper_class=structlog.stdlib.BoundLogger,
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
    )

    _configured = True


def log_trade(
    ticker: str,
    side: str,
    qty: float,
    entry_price: float,
    stop_loss: float,
    take_profit: float,
    order_id: str,
    status: str,
    technical_signal: str = "",
    fundamental_signal: str = "",
    orchestrator_confidence: int = 0,
    risk_approved: bool = True,
):
    """Append a trade record to the CSV audit log."""
    log_dir = Path(settings.LOG_PATH)
    log_dir.mkdir(parents=True, exist_ok=True)
    audit_file = log_dir / "trades.csv"

    fieldnames = [
        "timestamp", "ticker", "side", "qty", "entry_price", "stop_loss",
        "take_profit", "order_id", "status", "technical_signal",
        "fundamental_signal", "orchestrator_confidence", "risk_approved",
    ]
    row = {
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "ticker": ticker,
        "side": side,
        "qty": qty,
        "entry_price": entry_price,
        "stop_loss": stop_loss,
        "take_profit": take_profit,
        "order_id": order_id,
        "status": status,
        "technical_signal": technical_signal,
        "fundamental_signal": fundamental_signal,
        "orchestrator_confidence": orchestrator_confidence,
        "risk_approved": risk_approved,
    }

    write_header = not audit_file.exists()
    with open(audit_file, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)


def log_performance(
    portfolio_value: float,
    daily_pnl: float,
    num_trades: int,
    peak_value: float | None = None,
):
    """Append EOD performance record."""
    log_dir = Path(settings.LOG_PATH)
    log_dir.mkdir(parents=True, exist_ok=True)
    perf_file = log_dir / "performance.csv"

    fieldnames = ["date", "portfolio_value", "daily_pnl", "daily_pnl_pct", "num_trades", "max_drawdown_pct"]
    daily_pnl_pct = daily_pnl / (portfolio_value - daily_pnl) if (portfolio_value - daily_pnl) != 0 else 0
    drawdown = 0.0
    if peak_value and peak_value > 0:
        drawdown = max(0.0, (peak_value - portfolio_value) / peak_value)

    row = {
        "date": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "portfolio_value": round(portfolio_value, 2),
        "daily_pnl": round(daily_pnl, 2),
        "daily_pnl_pct": round(daily_pnl_pct * 100, 3),
        "num_trades": num_trades,
        "max_drawdown_pct": round(drawdown * 100, 3),
    }

    write_header = not perf_file.exists()
    with open(perf_file, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        if write_header:
            writer.writeheader()
        writer.writerow(row)
