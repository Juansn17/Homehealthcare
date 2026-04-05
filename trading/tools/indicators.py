"""Technical indicator computation using pandas-ta."""
from __future__ import annotations

import pandas as pd
import pandas_ta as ta


def compute_all(df: pd.DataFrame) -> dict:
    """Compute all indicators on a daily OHLCV DataFrame.

    Returns a flat dict of indicator values for the most recent bar.
    """
    close = df["close"]
    high = df["high"]
    low = df["low"]
    volume = df["volume"]

    result: dict = {}

    # ── RSI (14) ─────────────────────────────────────────────────────────────
    rsi = ta.rsi(close, length=14)
    result["rsi_14"] = _last(rsi)

    # ── MACD (12, 26, 9) ─────────────────────────────────────────────────────
    macd = ta.macd(close, fast=12, slow=26, signal=9)
    if macd is not None and not macd.empty:
        result["macd_line"] = _last(macd["MACD_12_26_9"])
        result["macd_signal"] = _last(macd["MACDs_12_26_9"])
        result["macd_hist"] = _last(macd["MACDh_12_26_9"])

    # ── Bollinger Bands (20, 2) ───────────────────────────────────────────────
    bb = ta.bbands(close, length=20, std=2)
    if bb is not None and not bb.empty:
        result["bb_upper"] = _last(bb["BBU_20_2.0"])
        result["bb_middle"] = _last(bb["BBM_20_2.0"])
        result["bb_lower"] = _last(bb["BBL_20_2.0"])
        result["bb_pct"] = _last(bb["BBP_20_2.0"])  # %B: position within bands

    # ── Moving Averages ───────────────────────────────────────────────────────
    result["sma_20"] = _last(ta.sma(close, length=20))
    result["sma_50"] = _last(ta.sma(close, length=50))
    result["sma_200"] = _last(ta.sma(close, length=200))
    result["ema_9"] = _last(ta.ema(close, length=9))

    # ── ATR (14) — used for stop-loss sizing ─────────────────────────────────
    atr = ta.atr(high, low, close, length=14)
    result["atr_14"] = _last(atr)

    # ── Trend context ─────────────────────────────────────────────────────────
    last_close = float(close.iloc[-1])
    sma50 = result.get("sma_50")
    sma200 = result.get("sma_200")
    if sma50 and sma200:
        result["golden_cross"] = sma50 > sma200
        result["above_sma200"] = last_close > sma200
        result["above_sma50"] = last_close > sma50

    return result


def _last(series: pd.Series | None) -> float | None:
    if series is None or series.dropna().empty:
        return None
    return float(series.dropna().iloc[-1])


def rule_based_signal(snapshot: dict) -> dict:
    """Pure Python fallback signal when Gemma is unavailable."""
    ind = snapshot.get("indicators", {})
    rsi = ind.get("rsi_14")
    macd_hist = ind.get("macd_hist")
    above_sma50 = ind.get("above_sma50", False)
    above_sma200 = ind.get("above_sma200", False)
    bb_pct = ind.get("bb_pct")

    score = 0
    reasons = []

    # RSI signals
    if rsi is not None:
        if rsi < 30:
            score += 2
            reasons.append(f"RSI oversold ({rsi:.1f})")
        elif rsi > 70:
            score -= 2
            reasons.append(f"RSI overbought ({rsi:.1f})")
        elif rsi < 45:
            score += 1
            reasons.append(f"RSI below midpoint ({rsi:.1f})")
        elif rsi > 60:
            score -= 1
            reasons.append(f"RSI above midpoint ({rsi:.1f})")

    # MACD
    if macd_hist is not None:
        if macd_hist > 0:
            score += 1
            reasons.append("MACD histogram positive")
        else:
            score -= 1
            reasons.append("MACD histogram negative")

    # Trend
    if above_sma200:
        score += 1
        reasons.append("Price above SMA200 (uptrend)")
    else:
        score -= 1
        reasons.append("Price below SMA200 (downtrend)")

    if above_sma50:
        score += 1
        reasons.append("Price above SMA50")

    # Bollinger
    if bb_pct is not None:
        if bb_pct < 0.1:
            score += 1
            reasons.append(f"Near lower Bollinger Band (BB%={bb_pct:.2f})")
        elif bb_pct > 0.9:
            score -= 1
            reasons.append(f"Near upper Bollinger Band (BB%={bb_pct:.2f})")

    if score >= 2:
        signal = "BUY"
        strength = min(1.0, score / 6.0)
    elif score <= -2:
        signal = "SELL"
        strength = min(1.0, abs(score) / 6.0)
    else:
        signal = "HOLD"
        strength = 0.3

    return {
        "ticker": snapshot["ticker"],
        "signal": signal,
        "strength": round(strength, 2),
        "reasons": reasons,
        "source": "rule_based",
    }
