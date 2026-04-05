"""NYSE market calendar and APScheduler integration."""
from __future__ import annotations

import logging
from datetime import datetime, time

import pytz

logger = logging.getLogger(__name__)

_EASTERN = pytz.timezone("America/New_York")

# NYSE hours in ET
_MARKET_OPEN = time(9, 30)
_MARKET_CLOSE = time(15, 45)   # slightly before 16:00 to avoid EOD issues
_PRE_MARKET_START = time(9, 0)


def get_market_phase() -> str:
    """Return the current NYSE market phase.

    Returns: "market_open" | "pre_market" | "after_hours" | "closed"
    """
    try:
        import pandas_market_calendars as mcal
        nyse = mcal.get_calendar("NYSE")
        now_et = datetime.now(_EASTERN)
        today_str = now_et.strftime("%Y-%m-%d")

        schedule = nyse.schedule(start_date=today_str, end_date=today_str)
        if schedule.empty:
            # Not a trading day (weekend/holiday)
            return "closed"

        current_time = now_et.time()

        if _PRE_MARKET_START <= current_time < _MARKET_OPEN:
            return "pre_market"
        elif _MARKET_OPEN <= current_time <= _MARKET_CLOSE:
            return "market_open"
        else:
            return "after_hours"

    except Exception as e:
        logger.warning("Market calendar check failed (%s), assuming open", e)
        # Safe default: assume open to not miss trades; risk layer will catch issues
        now_et = datetime.now(_EASTERN)
        t = now_et.time()
        if time(9, 30) <= t <= time(16, 0):
            return "market_open"
        return "closed"


def is_market_open() -> bool:
    return get_market_phase() == "market_open"


def setup_scheduler(run_cycle_fn, interval_minutes: int = 30):
    """Configure and return an APScheduler BlockingScheduler.

    run_cycle_fn: callable — the LangGraph cycle to run each interval.
    """
    from apscheduler.schedulers.blocking import BlockingScheduler
    from apscheduler.triggers.cron import CronTrigger
    from apscheduler.triggers.interval import IntervalTrigger

    scheduler = BlockingScheduler(timezone=_EASTERN)

    # Pre-market fundamentals fetch (09:05 ET every weekday)
    scheduler.add_job(
        run_cycle_fn,
        trigger=CronTrigger(day_of_week="mon-fri", hour=9, minute=5, timezone=_EASTERN),
        id="pre_market_cycle",
        name="Pre-market fundamental cycle",
        max_instances=1,
        coalesce=True,
    )

    # Intraday loop (every N minutes during market hours, weekdays)
    scheduler.add_job(
        run_cycle_fn,
        trigger=IntervalTrigger(minutes=interval_minutes, timezone=_EASTERN),
        id="intraday_cycle",
        name=f"Intraday cycle (every {interval_minutes}m)",
        max_instances=1,
        coalesce=True,
        misfire_grace_time=60,
    )

    logger.info(
        "Scheduler configured: pre-market @ 09:05 ET + every %d min intraday",
        interval_minutes,
    )
    return scheduler
