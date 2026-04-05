"""Market data fetching — yfinance primary, Alpaca bars fallback."""
from __future__ import annotations

import logging
from typing import Optional

import pandas as pd
import yfinance as yf

from trading.tools import indicators as ind

logger = logging.getLogger(__name__)


def get_ohlcv(ticker: str, period: str = "1y", interval: str = "1d") -> pd.DataFrame:
    """Fetch OHLCV data via yfinance.

    period: "6mo", "1y", "2y"
    interval: "1d", "1h", "15m"
    """
    try:
        df = yf.download(ticker, period=period, interval=interval, progress=False, auto_adjust=True)
        if df.empty:
            raise ValueError(f"Empty DataFrame for {ticker}")
        df.columns = [c.lower() for c in df.columns]
        df.index.name = "timestamp"
        return df.dropna()
    except Exception as e:
        logger.warning("yfinance failed for %s (%s), trying Alpaca fallback", ticker, e)
        return _alpaca_fallback(ticker, interval)


def _alpaca_fallback(ticker: str, interval: str) -> pd.DataFrame:
    """Fallback to Alpaca historical bars."""
    from trading.tools.alpaca_client import get_bars
    tf_map = {"1d": "1Day", "1h": "1Hour", "15m": "15Min"}
    bars = get_bars(ticker, timeframe=tf_map.get(interval, "1Day"), limit=200)
    df = pd.DataFrame(bars)
    df["timestamp"] = pd.to_datetime(df["timestamp"])
    df = df.set_index("timestamp")
    return df


def build_market_snapshot(ticker: str) -> dict:
    """Full snapshot: OHLCV + computed indicators for one ticker."""
    df = get_ohlcv(ticker, period="1y", interval="1d")
    computed = ind.compute_all(df)
    last = df.iloc[-1]
    avg_volume = df["volume"].tail(20).mean()
    return {
        "ticker": ticker,
        "price": float(last["close"]),
        "volume": float(last["volume"]),
        "volume_ratio": float(last["volume"] / avg_volume) if avg_volume > 0 else 1.0,
        "ohlcv": {
            "open": float(last["open"]),
            "high": float(last["high"]),
            "low": float(last["low"]),
            "close": float(last["close"]),
            "volume": float(last["volume"]),
        },
        "indicators": computed,
    }


def build_snapshots(tickers: list[str]) -> list[dict]:
    """Build snapshots for a list of tickers."""
    snapshots = []
    for ticker in tickers:
        try:
            snapshots.append(build_market_snapshot(ticker))
        except Exception as e:
            logger.error("Failed to build snapshot for %s: %s", ticker, e)
    return snapshots
