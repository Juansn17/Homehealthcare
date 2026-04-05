"""Earnings, fundamental, and macroeconomic data fetching."""
from __future__ import annotations

import functools
import logging
import time
from typing import Any

import yfinance as yf

logger = logging.getLogger(__name__)

# Simple in-process TTL cache for macro data (changes slowly)
_macro_cache: dict[str, tuple[float, Any]] = {}
MACRO_TTL_SECONDS = 4 * 3600  # 4 hours


def get_earnings_data(ticker: str) -> dict:
    """Fetch earnings calendar and recent EPS surprises via yfinance."""
    try:
        t = yf.Ticker(ticker)
        info = t.info or {}
        earnings = {}

        # Upcoming earnings date
        cal = t.calendar
        if cal is not None and not cal.empty:
            earnings["next_earnings_date"] = str(cal.iloc[0].get("Earnings Date", ""))

        # Recent quarterly EPS
        quarterly = t.quarterly_earnings
        if quarterly is not None and not quarterly.empty:
            recent = quarterly.head(4)
            earnings["quarterly_eps"] = recent.to_dict(orient="records")

        # Analyst estimates
        earnings["forward_eps"] = info.get("forwardEps")
        earnings["trailing_eps"] = info.get("trailingEps")
        earnings["earnings_growth"] = info.get("earningsGrowth")
        earnings["revenue_growth"] = info.get("revenueGrowth")
        return earnings
    except Exception as e:
        logger.warning("earnings_data failed for %s: %s", ticker, e)
        return {}


def get_valuation(ticker: str) -> dict:
    """Fetch valuation metrics via yfinance."""
    try:
        info = yf.Ticker(ticker).info or {}
        return {
            "pe_ratio": info.get("trailingPE"),
            "forward_pe": info.get("forwardPE"),
            "peg_ratio": info.get("pegRatio"),
            "price_to_book": info.get("priceToBook"),
            "profit_margin": info.get("profitMargins"),
            "debt_to_equity": info.get("debtToEquity"),
            "return_on_equity": info.get("returnOnEquity"),
            "enterprise_value": info.get("enterpriseValue"),
            "market_cap": info.get("marketCap"),
            "sector": info.get("sector"),
            "industry": info.get("industry"),
        }
    except Exception as e:
        logger.warning("valuation failed for %s: %s", ticker, e)
        return {}


def get_news_headlines(ticker: str, max_items: int = 5) -> list[dict]:
    """Fetch recent news headlines for a ticker."""
    try:
        news = yf.Ticker(ticker).news or []
        return [
            {
                "title": n.get("title", ""),
                "publisher": n.get("publisher", ""),
                "published_at": str(n.get("providerPublishTime", "")),
                "url": n.get("link", ""),
            }
            for n in news[:max_items]
        ]
    except Exception as e:
        logger.warning("news failed for %s: %s", ticker, e)
        return []


def get_macro_indicators() -> dict:
    """Fetch key macro indicators via FRED API (cached 4 hours)."""
    cache_key = "macro"
    if cache_key in _macro_cache:
        ts, val = _macro_cache[cache_key]
        if time.time() - ts < MACRO_TTL_SECONDS:
            return val

    result = {}
    try:
        from trading.config.settings import FRED_API_KEY
        if not FRED_API_KEY:
            raise ValueError("FRED_API_KEY not set")
        from fredapi import Fred
        fred = Fred(api_key=FRED_API_KEY)

        series = {
            "fed_funds_rate": "FEDFUNDS",
            "cpi_yoy": "CPIAUCSL",
            "unemployment_rate": "UNRATE",
            "t10y2y_spread": "T10Y2Y",
            "gdp_growth": "A191RL1Q225SBEA",
        }
        for label, series_id in series.items():
            try:
                s = fred.get_series(series_id, observation_start="2023-01-01")
                result[label] = float(s.dropna().iloc[-1]) if not s.dropna().empty else None
            except Exception:
                result[label] = None

        _macro_cache[cache_key] = (time.time(), result)
    except Exception as e:
        logger.warning("FRED macro fetch failed: %s", e)
        result["error"] = str(e)

    return result


def get_sector_performance() -> dict:
    """Return 1-month and 3-month returns for major sector ETFs."""
    cache_key = "sector_perf"
    if cache_key in _macro_cache:
        ts, val = _macro_cache[cache_key]
        if time.time() - ts < MACRO_TTL_SECONDS:
            return val

    from trading.config.universe import SECTOR_ETFS
    import yfinance as yf

    result = {}
    try:
        etfs = list(SECTOR_ETFS.keys())
        data = yf.download(etfs, period="3mo", interval="1d", progress=False, auto_adjust=True)
        close = data["Close"] if "Close" in data.columns else data["close"]
        for etf in etfs:
            if etf in close.columns:
                prices = close[etf].dropna()
                if len(prices) >= 2:
                    ret_1m = float((prices.iloc[-1] / prices.iloc[-21] - 1)) if len(prices) >= 21 else None
                    ret_3m = float((prices.iloc[-1] / prices.iloc[0] - 1))
                    result[etf] = {
                        "sector": SECTOR_ETFS[etf],
                        "return_1m": ret_1m,
                        "return_3m": ret_3m,
                    }
        _macro_cache[cache_key] = (time.time(), result)
    except Exception as e:
        logger.warning("sector_performance failed: %s", e)

    return result
