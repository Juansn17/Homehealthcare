"""Central configuration loaded from environment variables."""
import os
from pathlib import Path
from dotenv import load_dotenv

# Load .env from the trading/ directory
_ROOT = Path(__file__).resolve().parent.parent
load_dotenv(_ROOT / ".env")


def _get(key: str, default: str | None = None, required: bool = False) -> str:
    val = os.getenv(key, default)
    if required and not val:
        raise RuntimeError(f"Required environment variable '{key}' is not set.")
    return val


# ── Alpaca ──────────────────────────────────────────────────────────────────
ALPACA_API_KEY: str = _get("ALPACA_API_KEY", required=True)
ALPACA_SECRET_KEY: str = _get("ALPACA_SECRET_KEY", required=True)
ALPACA_BASE_URL: str = _get("ALPACA_BASE_URL", "https://paper-api.alpaca.markets")

# ── LLMs ────────────────────────────────────────────────────────────────────
ANTHROPIC_API_KEY: str = _get("ANTHROPIC_API_KEY", required=True)
CLAUDE_MODEL: str = _get("CLAUDE_MODEL", "claude-sonnet-4-6")

OLLAMA_BASE_URL: str = _get("OLLAMA_BASE_URL", "http://localhost:11434")
OLLAMA_MODEL: str = _get("OLLAMA_MODEL", "gemma3:4b")

# ── FRED ─────────────────────────────────────────────────────────────────────
FRED_API_KEY: str | None = _get("FRED_API_KEY")

# ── Trading mode ─────────────────────────────────────────────────────────────
TRADING_MODE: str = _get("TRADING_MODE", "paper")
IS_LIVE: bool = TRADING_MODE == "live"

# ── Persistence / logging ────────────────────────────────────────────────────
DB_PATH: str = _get("DB_PATH", str(_ROOT / "data" / "trading_state.db"))
LOG_PATH: str = _get("LOG_PATH", str(_ROOT / "logs"))

# ── Risk parameters ──────────────────────────────────────────────────────────
MAX_PORTFOLIO_RISK_PCT: float = float(_get("MAX_PORTFOLIO_RISK_PCT", "0.02"))
MAX_DRAWDOWN_PCT: float = float(_get("MAX_DRAWDOWN_PCT", "0.10"))
POSITION_SIZE_PCT: float = float(_get("POSITION_SIZE_PCT", "0.05"))
MAX_DAILY_LOSS_USD: float = float(_get("MAX_DAILY_LOSS_USD", "500"))

# ── Scheduling ───────────────────────────────────────────────────────────────
REBALANCE_INTERVAL_MINUTES: int = int(_get("REBALANCE_INTERVAL_MINUTES", "30"))

# ── Live safety ──────────────────────────────────────────────────────────────
REQUIRE_CONFIRMATION: bool = _get("REQUIRE_CONFIRMATION", "true").lower() == "true"
LIVE_POSITION_SCALE: float = float(_get("LIVE_POSITION_SCALE", "0.25"))

# ── Alerts ───────────────────────────────────────────────────────────────────
SLACK_WEBHOOK_URL: str | None = _get("SLACK_WEBHOOK_URL")
ALERT_EMAIL: str | None = _get("ALERT_EMAIL")
