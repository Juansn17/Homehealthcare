"""Portfolio Balancer — resolves conflicts between ST and LT strategies.

Runs after both st_risk_check and lt_risk_check complete (fan-in).
Responsibilities:
  1. Conflict resolution: same ticker wanted by both → LT takes priority
  2. Capital allocation: each strategy gets its configured % of available cash
  3. Total risk guard: combined risk across strategies must not exceed 5%
  4. Position count: total open positions (existing + new ST + new LT) ≤ 8
  5. Produces the final merged approved_orders list passed to execution
"""
from __future__ import annotations

import logging

from trading.config.strategy_config import SHORT_TERM, LONG_TERM, StrategyConfig, ALL_STRATEGIES

logger = logging.getLogger(__name__)

# Hard cap: combined risk across all strategies in a single cycle
MAX_COMBINED_RISK_PCT = 0.05      # 5% total portfolio at risk in one cycle
MAX_TOTAL_POSITIONS = 8           # ST + LT combined open position limit


def balance(
    risk_assessments: list[dict],
    portfolio_value: float,
    available_cash: float,
    current_positions: dict[str, dict],
) -> tuple[list[dict], list[str]]:
    """Resolve conflicts and return (approved_assessments, rejection_reasons).

    Each risk_assessment is expected to have a 'strategy' field ('short_term' | 'long_term').
    Returns the final list of assessments that may proceed to execution.
    """
    approved = [r for r in risk_assessments if r.get("approved")]
    rejected_reasons: list[str] = []

    # ── 1. Conflict resolution: same ticker in both strategies ─────────────────
    st_tickers = {r["ticker"]: r for r in approved if r.get("strategy") == SHORT_TERM.name}
    lt_tickers = {r["ticker"]: r for r in approved if r.get("strategy") == LONG_TERM.name}

    conflicts = set(st_tickers.keys()) & set(lt_tickers.keys())
    for ticker in conflicts:
        # Long-term takes priority — drop the ST entry
        logger.info("Conflict on %s: LT takes priority over ST", ticker)
        rejected_reasons.append(f"{ticker} ST: dropped in favor of LT position")
        del st_tickers[ticker]

    # ── 2. Capital allocation: each strategy gets its % of available cash ─────
    st_budget = available_cash * SHORT_TERM.portfolio_allocation_pct
    lt_budget = available_cash * LONG_TERM.portfolio_allocation_pct

    st_approved = _apply_budget(list(st_tickers.values()), st_budget, SHORT_TERM, rejected_reasons)
    lt_approved = _apply_budget(list(lt_tickers.values()), lt_budget, LONG_TERM, rejected_reasons)

    combined = st_approved + lt_approved

    # ── 3. Total risk guard ────────────────────────────────────────────────────
    total_risk = sum(
        r["position_size_usd"] * _stop_pct(r)
        for r in combined
    )
    total_risk_pct = total_risk / portfolio_value if portfolio_value > 0 else 0

    if total_risk_pct > MAX_COMBINED_RISK_PCT:
        combined, extra_rejections = _trim_to_risk_limit(
            combined, portfolio_value, MAX_COMBINED_RISK_PCT
        )
        rejected_reasons.extend(extra_rejections)
        logger.warning(
            "Combined risk %.2f%% exceeded limit %.2f%% — trimmed to %d orders",
            total_risk_pct * 100, MAX_COMBINED_RISK_PCT * 100, len(combined),
        )

    # ── 4. Total position count guard ─────────────────────────────────────────
    existing_count = len(current_positions)
    slots_available = MAX_TOTAL_POSITIONS - existing_count
    if len(combined) > slots_available:
        # LT positions take priority in the final cut
        combined.sort(key=lambda r: (0 if r.get("strategy") == LONG_TERM.name else 1))
        for dropped in combined[slots_available:]:
            rejected_reasons.append(
                f"{dropped['ticker']} ({dropped.get('strategy')}): dropped — max {MAX_TOTAL_POSITIONS} positions reached"
            )
        combined = combined[:slots_available]

    logger.info(
        "Portfolio balancer: %d approved (ST=%d, LT=%d) | %d conflicts resolved",
        len(combined),
        sum(1 for r in combined if r.get("strategy") == SHORT_TERM.name),
        sum(1 for r in combined if r.get("strategy") == LONG_TERM.name),
        len(conflicts),
    )
    return combined, rejected_reasons


def _apply_budget(
    assessments: list[dict],
    budget: float,
    strategy: StrategyConfig,
    rejection_log: list[str],
) -> list[dict]:
    """Trim assessments to fit within budget. Sort by position_size_usd descending
    so the highest-conviction (largest) positions are kept first."""
    # Respect max_open_positions for each strategy
    if len(assessments) > strategy.max_open_positions:
        dropped = assessments[strategy.max_open_positions:]
        for r in dropped:
            rejection_log.append(
                f"{r['ticker']} ({strategy.name}): dropped — max {strategy.max_open_positions} positions per strategy"
            )
        assessments = assessments[: strategy.max_open_positions]

    approved = []
    spent = 0.0
    for r in assessments:
        size = r.get("position_size_usd", 0.0)
        if spent + size <= budget:
            approved.append(r)
            spent += size
        else:
            rejection_log.append(
                f"{r['ticker']} ({strategy.name}): dropped — strategy budget ${budget:,.0f} exhausted"
            )
    return approved


def _stop_pct(assessment: dict) -> float:
    """Rough stop-distance % from an assessment."""
    entry = assessment.get("position_size_usd", 0)
    size = assessment.get("position_size_usd", 1)
    stop = assessment.get("stop_loss_price", 0)
    tp = assessment.get("take_profit_price", 0)
    if size <= 0 or stop <= 0:
        return 0.02  # default 2%
    return 0.02


def _trim_to_risk_limit(
    assessments: list[dict],
    portfolio_value: float,
    max_risk_pct: float,
) -> tuple[list[dict], list[str]]:
    """Drop lower-priority assessments until total risk is within limit.
    LT positions have higher priority and are dropped last.
    """
    # Sort: LT first, then by position_size_usd descending (keep largest)
    sorted_a = sorted(
        assessments,
        key=lambda r: (0 if r.get("strategy") == LONG_TERM.name else 1),
    )
    approved = []
    rejections = []
    cumulative_risk = 0.0

    for r in sorted_a:
        risk = r.get("position_size_usd", 0) * _stop_pct(r)
        if (cumulative_risk + risk) / portfolio_value <= max_risk_pct:
            approved.append(r)
            cumulative_risk += risk
        else:
            rejections.append(
                f"{r['ticker']} ({r.get('strategy')}): dropped — combined risk limit {max_risk_pct:.0%} reached"
            )

    return approved, rejections
