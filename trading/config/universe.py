"""Stock universe: liquid US large-caps + sector ETFs as macro proxies."""

# Core equity universe (~20 liquid names across sectors)
STOCK_UNIVERSE: list[str] = [
    # Technology
    "AAPL", "MSFT", "NVDA", "GOOGL", "META",
    # Consumer / E-commerce
    "AMZN", "TSLA",
    # Financials
    "JPM", "V", "MA",
    # Healthcare
    "UNH", "JNJ",
    # Energy
    "XOM", "CVX",
    # Industrials
    "CAT", "HON",
    # Communication
    "T", "VZ",
]

# Sector ETFs — used as macro/sector rotation proxies (not traded directly)
SECTOR_ETFS: dict[str, str] = {
    "XLK": "Technology",
    "XLF": "Financials",
    "XLE": "Energy",
    "XLV": "Healthcare",
    "XLY": "Consumer Discretionary",
    "XLP": "Consumer Staples",
    "XLI": "Industrials",
    "XLU": "Utilities",
    "XLRE": "Real Estate",
    "XLB": "Materials",
    "XLC": "Communication",
}

# Broad market benchmarks
BENCHMARKS: list[str] = ["SPY", "QQQ", "IWM"]

# Sector mapping for concentration checks (ticker → GICS sector)
SECTOR_MAP: dict[str, str] = {
    "AAPL": "Technology", "MSFT": "Technology", "NVDA": "Technology",
    "GOOGL": "Communication", "META": "Communication",
    "AMZN": "Consumer Discretionary", "TSLA": "Consumer Discretionary",
    "JPM": "Financials", "V": "Financials", "MA": "Financials",
    "UNH": "Healthcare", "JNJ": "Healthcare",
    "XOM": "Energy", "CVX": "Energy",
    "CAT": "Industrials", "HON": "Industrials",
    "T": "Communication", "VZ": "Communication",
}
