"""Alpaca Markets REST client — paper and live trading."""
from __future__ import annotations

from alpaca.trading.client import TradingClient
from alpaca.trading.requests import (
    MarketOrderRequest,
    LimitOrderRequest,
    BracketOrderRequest,
    GetOrdersRequest,
)
from alpaca.trading.enums import OrderSide, TimeInForce, QueryOrderStatus
from alpaca.data.historical import StockHistoricalDataClient
from alpaca.data.requests import StockBarsRequest
from alpaca.data.timeframe import TimeFrame, TimeFrameUnit

from trading.config import settings


def _trading_client() -> TradingClient:
    return TradingClient(
        api_key=settings.ALPACA_API_KEY,
        secret_key=settings.ALPACA_SECRET_KEY,
        paper=not settings.IS_LIVE,
    )


def _data_client() -> StockHistoricalDataClient:
    return StockHistoricalDataClient(
        api_key=settings.ALPACA_API_KEY,
        secret_key=settings.ALPACA_SECRET_KEY,
    )


def get_account() -> dict:
    """Return current account info (buying power, equity, etc.)."""
    client = _trading_client()
    account = client.get_account()
    return {
        "portfolio_value": float(account.portfolio_value),
        "buying_power": float(account.buying_power),
        "cash": float(account.cash),
        "equity": float(account.equity),
        "daytrade_count": account.daytrade_count,
    }


def get_positions() -> dict[str, dict]:
    """Return current open positions keyed by ticker."""
    client = _trading_client()
    positions = client.get_all_positions()
    return {
        p.symbol: {
            "qty": float(p.qty),
            "avg_entry_price": float(p.avg_entry_price),
            "market_value": float(p.market_value),
            "unrealized_pl": float(p.unrealized_pl),
            "unrealized_plpc": float(p.unrealized_plpc),
            "side": p.side.value,
        }
        for p in positions
    }


def submit_bracket_order(
    ticker: str,
    qty: float,
    side: str,
    stop_loss_price: float,
    take_profit_price: float,
) -> dict:
    """Submit a bracket order (entry + stop-loss + take-profit)."""
    client = _trading_client()
    order_side = OrderSide.BUY if side.lower() == "buy" else OrderSide.SELL
    request = MarketOrderRequest(
        symbol=ticker,
        qty=qty,
        side=order_side,
        time_in_force=TimeInForce.DAY,
        order_class="bracket",
        stop_loss={"stop_price": round(stop_loss_price, 2)},
        take_profit={"limit_price": round(take_profit_price, 2)},
    )
    order = client.submit_order(request)
    return {
        "order_id": str(order.id),
        "ticker": ticker,
        "side": side,
        "qty": float(order.qty),
        "status": order.status.value,
        "created_at": str(order.created_at),
    }


def get_order(order_id: str) -> dict:
    """Fetch current status of an order."""
    client = _trading_client()
    order = client.get_order_by_id(order_id)
    return {
        "order_id": str(order.id),
        "status": order.status.value,
        "filled_qty": float(order.filled_qty or 0),
        "filled_avg_price": float(order.filled_avg_price or 0),
    }


def cancel_order(order_id: str) -> bool:
    """Cancel a pending order. Returns True if successful."""
    client = _trading_client()
    try:
        client.cancel_order_by_id(order_id)
        return True
    except Exception:
        return False


def get_bars(ticker: str, timeframe: str = "1Day", limit: int = 200) -> list[dict]:
    """Fetch historical OHLCV bars from Alpaca.

    timeframe options: "1Min", "5Min", "15Min", "1Hour", "1Day"
    """
    tf_map = {
        "1Min": TimeFrame(1, TimeFrameUnit.Minute),
        "5Min": TimeFrame(5, TimeFrameUnit.Minute),
        "15Min": TimeFrame(15, TimeFrameUnit.Minute),
        "1Hour": TimeFrame(1, TimeFrameUnit.Hour),
        "1Day": TimeFrame(1, TimeFrameUnit.Day),
    }
    timeframe_obj = tf_map.get(timeframe, TimeFrame(1, TimeFrameUnit.Day))
    client = _data_client()
    request = StockBarsRequest(symbol_or_symbols=ticker, timeframe=timeframe_obj, limit=limit)
    bars = client.get_stock_bars(request)
    result = []
    for bar in bars[ticker]:
        result.append({
            "timestamp": str(bar.timestamp),
            "open": float(bar.open),
            "high": float(bar.high),
            "low": float(bar.low),
            "close": float(bar.close),
            "volume": float(bar.volume),
        })
    return result
