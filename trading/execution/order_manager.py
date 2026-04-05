"""Order lifecycle management — bracket order submission and fill tracking."""
from __future__ import annotations

import logging
import time
from datetime import datetime, timezone

from trading.config import settings
from trading.tools import alpaca_client

logger = logging.getLogger(__name__)


def submit_approved_orders(
    risk_assessments: list[dict],
    market_snapshots: list[dict],
) -> list[dict]:
    """Submit bracket orders for all approved risk assessments.

    In live mode with REQUIRE_CONFIRMATION=True, prints order details
    and prompts for CLI approval before submitting.
    """
    results = []
    snapshot_map = {s["ticker"]: s for s in market_snapshots}

    for assessment in risk_assessments:
        if not assessment.get("approved"):
            continue

        ticker = assessment["ticker"]
        snapshot = snapshot_map.get(ticker)
        if not snapshot:
            logger.warning("No snapshot for %s, skipping order", ticker)
            continue

        side = "buy"  # orchestrator only sends BUY decisions through risk for now
        qty = assessment["qty"]
        stop_loss = assessment["stop_loss_price"]
        take_profit = assessment["take_profit_price"]
        entry_price = snapshot["price"]

        if settings.IS_LIVE and settings.REQUIRE_CONFIRMATION:
            print(
                f"\n[LIVE TRADE CONFIRMATION REQUIRED]\n"
                f"  Ticker:      {ticker}\n"
                f"  Side:        {side.upper()}\n"
                f"  Qty:         {qty:.2f} shares\n"
                f"  Entry:       ~${entry_price:.2f}\n"
                f"  Stop-loss:   ${stop_loss:.2f}\n"
                f"  Take-profit: ${take_profit:.2f}\n"
                f"  Size:        ${assessment['position_size_usd']:,.0f}\n"
            )
            confirm = input("Approve? [y/N]: ").strip().lower()
            if confirm != "y":
                logger.info("Order for %s manually rejected", ticker)
                results.append(_make_result(ticker, side, qty, "manually_rejected"))
                continue

        try:
            order = alpaca_client.submit_bracket_order(
                ticker=ticker,
                qty=qty,
                side=side,
                stop_loss_price=stop_loss,
                take_profit_price=take_profit,
            )
            logger.info("Order submitted: %s %s %.2f @ ~$%.2f", side, ticker, qty, entry_price)
            results.append({
                "ticker": ticker,
                "order_id": order["order_id"],
                "side": side,
                "qty": qty,
                "fill_price": None,
                "status": "submitted",
                "timestamp": datetime.now(timezone.utc).isoformat(),
            })
        except Exception as e:
            logger.error("Order submission failed for %s: %s", ticker, e)
            results.append(_make_result(ticker, side, qty, "submission_failed"))

    return results


def check_fills(order_results: list[dict]) -> list[dict]:
    """Poll Alpaca for fill status on submitted orders."""
    updated = []
    for result in order_results:
        if result.get("status") not in ("submitted",):
            updated.append(result)
            continue
        order_id = result.get("order_id")
        if not order_id:
            updated.append(result)
            continue
        try:
            status = alpaca_client.get_order(order_id)
            fill_price = status.get("filled_avg_price") or None
            new_status = _map_alpaca_status(status.get("status", ""))
            updated.append({**result, "status": new_status, "fill_price": fill_price})
        except Exception as e:
            logger.warning("Could not check fill for order %s: %s", order_id, e)
            updated.append(result)
    return updated


def _make_result(ticker: str, side: str, qty: float, status: str) -> dict:
    return {
        "ticker": ticker,
        "order_id": "",
        "side": side,
        "qty": qty,
        "fill_price": None,
        "status": status,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }


def _map_alpaca_status(alpaca_status: str) -> str:
    filled = {"filled", "partially_filled"}
    rejected = {"rejected", "canceled", "expired", "suspended"}
    if alpaca_status in filled:
        return "filled"
    if alpaca_status in rejected:
        return "rejected"
    return "submitted"
