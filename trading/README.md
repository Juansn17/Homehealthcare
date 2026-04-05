# Trading Agents — Portfolio Automation System

Sistema multi-agente de trading algorítmico para acciones US, construido con **LangGraph**, **Claude Sonnet** y **Gemma 4 local**. Opera en **Alpaca Markets** (paper trading por defecto, con switch a live).

---

## Arquitectura general

```
┌─────────────────────────────────────────────────────────────────┐
│                        LangGraph Graph                          │
│                                                                 │
│  START → market_clock_check                                     │
│               │                                                 │
│          ┌────┴────┐                                            │
│        open      closed → sleep → END                           │
│          │                                                      │
│       data_fetch  (OHLCV + indicadores + cuenta Alpaca)         │
│          │                                                      │
│    ┌─────┴──────┐   ← fan-out paralelo                          │
│    ▼            ▼                                               │
│ technical_   fundamental_                                       │
│  agent        agent                                             │
│ (Gemma 4)   (Claude Sonnet)                                     │
│    └─────┬──────┘   ← fan-in (espera ambos)                     │
│          ▼                                                      │
│   orchestrator_agent  (Claude Sonnet)                           │
│   consolida señales → top 5 decisiones                          │
│          │                                                      │
│      risk_check  (Python puro — sin LLM)                        │
│          │                                                      │
│    ┌─────┴──────┐                                               │
│  aprobado    rechazado                                          │
│    │             │                                              │
│ execution_   halt_node → END                                    │
│  node                                                           │
│ (Alpaca bracket orders)                                         │
│    │                                                            │
│ monitoring_node → END                                           │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### Los tres agentes

| Agente | LLM | Rol |
|--------|-----|-----|
| **Technical Agent** | Gemma 4 (local, Ollama) | Screena el universo completo con RSI, MACD, Bollinger Bands, SMAs, ATR. Genera señales BUY/SELL/HOLD con convicción 0–1. Fallback a reglas Python si Ollama no está disponible. |
| **Fundamental Agent** | Claude Sonnet (API) | Analiza P/E, earnings, crecimiento de revenue, macro FRED (tasas, inflación, spread yield curve) y rotación sectorial. Solo evalúa los tickers que pasaron el filtro técnico. |
| **Orchestrator Agent** | Claude Sonnet (API) | Recibe ambas señales, rankea las top 3–5 oportunidades, asigna confianza 0–100, rechaza si los analistas se contradicen con alta convicción. Produce decisiones finales. |

---

## Estructura del proyecto

```
trading/
├── main.py                        # Entrypoint: --once, --dry-run, o scheduler
├── requirements.txt
├── .env.example                   # Plantilla de variables de entorno
│
├── config/
│   ├── settings.py                # Config central (env vars)
│   └── universe.py               # Universo de acciones (~18 large-caps US)
│
├── graph/
│   ├── state.py                   # TradingState TypedDict (estado del ciclo)
│   ├── nodes.py                   # Funciones de cada nodo del grafo
│   ├── edges.py                   # Routing condicional
│   └── workflow.py                # Ensamblado + compilación del grafo
│
├── agents/
│   ├── technical_agent.py         # Gemma 4 + fallback reglas Python
│   ├── fundamental_agent.py       # Claude Sonnet + FRED + yfinance
│   └── orchestrator_agent.py      # Claude Sonnet, structured output
│
├── tools/
│   ├── alpaca_client.py           # Alpaca REST (órdenes, posiciones, barras)
│   ├── market_data.py             # yfinance + fallback Alpaca bars
│   ├── indicators.py              # pandas-ta: RSI, MACD, BB, SMA, ATR
│   └── news_macro.py              # Earnings, valuación, FRED, sector ETFs
│
├── llm/
│   ├── claude_client.py           # Anthropic SDK, structured output, retry
│   └── ollama_client.py           # Ollama REST, Gemma 4, fallback handler
│
├── risk/
│   └── manager.py                 # Límites duros: sizing, stop-loss, drawdown
│
├── execution/
│   └── order_manager.py           # Bracket orders, confirmación, fill tracking
│
├── scheduler/
│   └── market_clock.py            # Calendario NYSE + APScheduler
│
├── observability/
│   └── logger.py                  # structlog JSON, audit CSV, performance CSV
│
├── persistence/                   # LangGraph SQLite checkpointer (auto)
├── data/                          # trading_state.db (generado)
└── logs/                          # app.log, trades.csv, performance.csv (generado)
```

---

## Setup inicial

### 1. Dependencias Python

```bash
cd /path/to/Homehealthcare

python3 -m venv trading/.venv
source trading/.venv/bin/activate   # Windows: trading\.venv\Scripts\activate

pip install -r trading/requirements.txt
```

### 2. Gemma 4 local (Ollama)

```bash
# Instalar Ollama
curl -fsSL https://ollama.ai/install.sh | sh

# Descargar Gemma 4 (4B — ~3GB en disco)
ollama pull gemma3:4b

# Verificar que corre
ollama run gemma3:4b "responde solo: OK"
```

### 3. Variables de entorno

```bash
cp trading/.env.example trading/.env
```

Edita `trading/.env` con tus credenciales:

```dotenv
# Alpaca (paper por defecto — crea cuenta gratis en alpaca.markets)
ALPACA_API_KEY=tu_api_key
ALPACA_SECRET_KEY=tu_secret_key
ALPACA_BASE_URL=https://paper-api.alpaca.markets

# Claude API (console.anthropic.com)
ANTHROPIC_API_KEY=tu_anthropic_key

# FRED API — gratis: https://fred.stlouisfed.org/docs/api/api_key.html
FRED_API_KEY=tu_fred_key   # opcional, degradará graciosamente sin él

# Ollama (local, no necesita key)
OLLAMA_BASE_URL=http://localhost:11434
OLLAMA_MODEL=gemma3:4b
```

---

## Uso

### Ciclo único (testing)

```bash
python trading/main.py --once
```

Ejecuta un ciclo completo: clock check → data → TA → FA → orquestador → risk → ejecución. Ideal para verificar que todo el pipeline funciona antes de activar el scheduler.

### Scheduler continuo (producción paper)

```bash
python trading/main.py
```

Levanta APScheduler con:
- `09:05 ET` — ciclo pre-mercado (solo fundamentales, sin ejecución)
- Cada `30 min` durante horario NYSE — ciclo completo

### Modo live (dinero real)

```bash
# En trading/.env:
ALPACA_BASE_URL=https://api.alpaca.markets
TRADING_MODE=live
REQUIRE_CONFIRMATION=true      # solicita confirmación CLI por orden
LIVE_POSITION_SCALE=0.25       # empieza al 25% del tamaño máximo
```

> **Advertencia:** `TRADING_MODE=live` ejecuta órdenes con dinero real en Alpaca. Asegúrate de haber validado el sistema en paper trading durante al menos 2 semanas antes de hacer el switch.

---

## Risk management

Todas las siguientes reglas son **Python puro** — ningún LLM puede sobrescribirlas.

| Regla | Parámetro | Descripción |
|-------|-----------|-------------|
| **Sizing por riesgo** | `MAX_PORTFOLIO_RISK_PCT=0.02` | Máximo 2% del portafolio arriesgado por trade |
| **Tamaño máximo** | `POSITION_SIZE_PCT=0.05` | Máximo 5% del portafolio por posición |
| **Stop-loss** | ATR(14) × 2.0 | Basado en volatilidad real del activo |
| **Distancia mínima** | 1% | Rechaza si el stop está muy cerca |
| **Distancia máxima** | 8% | Rechaza si el stop requiere asumir demasiado riesgo |
| **Take-profit mínimo** | Ratio 2:1 | Mínimo 2x la distancia al stop como target |
| **Max drawdown** | `MAX_DRAWDOWN_PCT=0.10` | Congela nuevas entradas si el portafolio cae >10% desde el pico |
| **Concentración sectorial** | 30% | Máximo 30% del portafolio en un mismo sector GICS |
| **Posición duplicada** | — | Rechaza si ya hay posición abierta en el ticker |
| **Pérdida diaria** | `MAX_DAILY_LOSS_USD=500` | Solo en modo live — detiene el día si se alcanza el límite |

---

## Universo de acciones

18 large-caps US distribuidas en 6 sectores:

| Sector | Tickers |
|--------|---------|
| Technology | AAPL, MSFT, NVDA |
| Communication | GOOGL, META, T, VZ |
| Consumer Disc. | AMZN, TSLA |
| Financials | JPM, V, MA |
| Healthcare | UNH, JNJ |
| Energy | XOM, CVX |
| Industrials | CAT, HON |

Edita `trading/config/universe.py` para ajustar el universo.

---

## Observabilidad

Todos los artefactos se generan automáticamente en `trading/logs/`:

| Archivo | Contenido |
|---------|-----------|
| `app.log` | Logs estructurados JSON (structlog): cada nodo, llamadas LLM, errores |
| `trades.csv` | Audit trail por orden: ticker, side, qty, entry, stop, target, status, señales |
| `performance.csv` | Registro EOD: portfolio value, P&L diario, num trades, drawdown |

### Estado persistido (LangGraph SQLite)

`trading/data/trading_state.db` almacena el estado completo del grafo después de cada nodo. Si el proceso se cae a mitad de ciclo, puede retomar desde el último checkpoint.

```bash
# Inspeccionar estado
sqlite3 trading/data/trading_state.db ".tables"

# Ver últimos cycles
sqlite3 trading/data/trading_state.db "SELECT * FROM checkpoints ORDER BY step DESC LIMIT 5;"
```

---

## Checklist paper trading → live

Completa estas validaciones antes de encender dinero real:

- [ ] Sistema corre 10 días consecutivos sin crashes
- [ ] Bracket orders (stop-loss + take-profit) se crean y cancelan correctamente en paper
- [ ] Drawdown guard se activa correctamente al forzar una pérdida del 10%
- [ ] Estado se recupera después de matar y reiniciar el proceso
- [ ] Outputs de Gemma 4 y Claude pasan validación Pydantic el 100% de los ciclos
- [ ] `trades.csv` y `performance.csv` se generan correctamente
- [ ] Win rate > 45% en 20+ trades (mínimo estadístico)
- [ ] Ningún bug de position sizing (verificar `qty` y `position_size_usd` en logs)

---

## Ajuste de parámetros

Todos los parámetros de riesgo y scheduling se controlan desde `trading/.env`:

```dotenv
MAX_PORTFOLIO_RISK_PCT=0.02       # riesgo por trade (2%)
MAX_DRAWDOWN_PCT=0.10             # halt si drawdown > 10%
POSITION_SIZE_PCT=0.05            # máximo por posición (5%)
MAX_DAILY_LOSS_USD=500            # límite de pérdida diaria (solo live)
REBALANCE_INTERVAL_MINUTES=30     # frecuencia de ciclos intraday
LIVE_POSITION_SCALE=0.25          # escala de tamaño en live (25% → 50% → 100%)
```

---

## Tecnologías

| Categoría | Tecnología |
|-----------|-----------|
| Orquestación | [LangGraph](https://langchain-ai.github.io/langgraph/) |
| LLM cloud | [Claude Sonnet](https://www.anthropic.com) via Anthropic SDK |
| LLM local | [Gemma 4](https://ollama.ai/library/gemma3) via Ollama |
| Broker | [Alpaca Markets](https://alpaca.markets) (alpaca-py v2) |
| Datos de mercado | yfinance, Alpaca Historical Bars |
| Indicadores técnicos | pandas-ta |
| Datos macro | FRED API (fredapi) |
| Scheduling | APScheduler + pandas_market_calendars |
| Persistencia | LangGraph SQLite Checkpointer |
| Logging | structlog |
