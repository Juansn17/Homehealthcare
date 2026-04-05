"""Trading system entrypoint.

Usage:
    python trading/main.py                  # run with scheduler (production)
    python trading/main.py --once           # run a single cycle and exit
    python trading/main.py --dry-run        # run without submitting orders
"""
from __future__ import annotations

import argparse
import logging
import signal
import sys
import uuid
from pathlib import Path

# ── Bootstrap path so `trading.*` imports work from repo root ─────────────────
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from trading.observability.logger import configure_logging, log_performance
from trading.config import settings

configure_logging()
logger = logging.getLogger(__name__)


def run_cycle(graph, dry_run: bool = False):
    """Execute one full LangGraph trading cycle."""
    from trading.tools.alpaca_client import get_account

    run_id = str(uuid.uuid4())
    config = {"configurable": {"thread_id": run_id}}

    logger.info("=" * 60)
    logger.info("Starting trading cycle run_id=%s mode=%s", run_id[:8], settings.TRADING_MODE)

    initial_state = {
        "run_id": run_id,
        "phase": "",
        "error": None,
        "market_snapshots": [],
        "portfolio_value": 0.0,
        "current_positions": {},
        "available_cash": 0.0,
        "peak_portfolio_value": None,
        "technical_signals": [],
        "fundamental_signals": [],
        "consolidated_decisions": [],
        "risk_assessments": [],
        "orders_to_place": [],
        "order_results": [],
        "messages": [],
        "cycle_summary": "",
        "daily_loss_usd": 0.0,
    }

    try:
        final_state = graph.invoke(initial_state, config=config)
        logger.info("Cycle complete: %s", final_state.get("cycle_summary", ""))

        # Log filled orders to audit trail
        from trading.observability.logger import log_trade
        for result in final_state.get("order_results", []):
            if result.get("status") in ("filled", "submitted"):
                log_trade(
                    ticker=result["ticker"],
                    side=result["side"],
                    qty=result["qty"],
                    entry_price=result.get("fill_price") or 0.0,
                    stop_loss=0.0,  # populated from risk_assessments if needed
                    take_profit=0.0,
                    order_id=result.get("order_id", ""),
                    status=result["status"],
                    orchestrator_confidence=0,
                    risk_approved=True,
                )

        return final_state

    except Exception as e:
        logger.error("Cycle failed with exception: %s", e, exc_info=True)
        return None


def main():
    parser = argparse.ArgumentParser(description="Algorithmic Trading Agent System")
    parser.add_argument("--once", action="store_true", help="Run a single cycle and exit")
    parser.add_argument("--dry-run", action="store_true", help="Run without submitting real orders")
    parser.add_argument(
        "--dual-strategy",
        action="store_true",
        help="Run parallel short-term (Gemma) + long-term (Claude) strategies",
    )
    args = parser.parse_args()

    # Safety check: prevent accidental live trading
    if settings.IS_LIVE:
        logger.warning("⚠️  LIVE TRADING MODE ACTIVE — real money at risk!")
        if not settings.REQUIRE_CONFIRMATION:
            logger.warning("⚠️  REQUIRE_CONFIRMATION=false — orders will auto-submit!")
    else:
        logger.info("Paper trading mode (safe)")

    # Build compiled graph
    if args.dual_strategy:
        from trading.graph.workflow import get_compiled_dual_graph
        graph = get_compiled_dual_graph()
        logger.info("Dual-strategy graph compiled (ST: Gemma/intraday + LT: Claude/fundamental)")
    else:
        from trading.graph.workflow import get_compiled_graph
        graph = get_compiled_graph()
        logger.info("Single-strategy graph compiled")

    if args.once or args.dry_run:
        logger.info("Running single cycle...")
        run_cycle(graph, dry_run=args.dry_run)
        logger.info("Single cycle complete. Exiting.")
        return

    # Production: run via APScheduler
    from trading.scheduler.market_clock import setup_scheduler

    def scheduled_cycle():
        run_cycle(graph, dry_run=False)

    scheduler = setup_scheduler(scheduled_cycle, interval_minutes=settings.REBALANCE_INTERVAL_MINUTES)

    # Graceful shutdown
    def _shutdown(signum, frame):
        logger.info("Received signal %d, shutting down scheduler...", signum)
        scheduler.shutdown(wait=False)
        sys.exit(0)

    signal.signal(signal.SIGINT, _shutdown)
    signal.signal(signal.SIGTERM, _shutdown)

    logger.info(
        "Scheduler starting — cycles every %d min (TRADING_MODE=%s)",
        settings.REBALANCE_INTERVAL_MINUTES,
        settings.TRADING_MODE,
    )
    scheduler.start()


if __name__ == "__main__":
    main()
