"""Ollama client for local Gemma 4 inference with rule-based fallback."""
from __future__ import annotations

import json
import logging
from typing import Any, Type

import httpx
from pydantic import BaseModel
from tenacity import retry, stop_after_attempt, wait_fixed, retry_if_exception_type

from trading.config import settings

logger = logging.getLogger(__name__)

_TIMEOUT = httpx.Timeout(120.0, connect=10.0)


def _is_available() -> bool:
    """Check if Ollama is reachable."""
    try:
        r = httpx.get(f"{settings.OLLAMA_BASE_URL}/api/tags", timeout=5.0)
        return r.status_code == 200
    except Exception:
        return False


def call(prompt: str, system_prompt: str = "", max_tokens: int = 512) -> str:
    """Call Ollama chat endpoint. Raises if unavailable."""
    messages = []
    if system_prompt:
        messages.append({"role": "system", "content": system_prompt})
    messages.append({"role": "user", "content": prompt})

    payload = {
        "model": settings.OLLAMA_MODEL,
        "messages": messages,
        "stream": False,
        "options": {"num_predict": max_tokens},
    }
    r = httpx.post(
        f"{settings.OLLAMA_BASE_URL}/api/chat",
        json=payload,
        timeout=_TIMEOUT,
    )
    r.raise_for_status()
    return r.json()["message"]["content"]


def call_structured(
    system_prompt: str,
    user_prompt: str,
    output_schema: Type[BaseModel],
    max_tokens: int = 512,
    fallback_fn=None,
    fallback_input: Any = None,
) -> BaseModel:
    """Call Gemma via Ollama and parse into a Pydantic model.

    Falls back to fallback_fn(fallback_input) if Ollama is unavailable or
    if the output fails to parse.
    """
    schema_json = json.dumps(output_schema.model_json_schema(), indent=2)
    enriched_system = (
        f"{system_prompt}\n\n"
        f"Respond ONLY with valid JSON matching this schema:\n{schema_json}\n"
        f"No explanation, no markdown fences — just the JSON."
    )

    def _try_ollama() -> BaseModel | None:
        if not _is_available():
            logger.warning("Ollama not reachable at %s — using fallback", settings.OLLAMA_BASE_URL)
            return None
        raw = call(user_prompt, system_prompt=enriched_system, max_tokens=max_tokens)
        text = raw.strip()
        if text.startswith("```"):
            text = text.split("```")[1]
            if text.startswith("json"):
                text = text[4:]
        data = json.loads(text)
        return output_schema.model_validate(data)

    try:
        result = _try_ollama()
        if result is not None:
            return result
    except Exception as e:
        logger.warning("Ollama structured call failed: %s — using fallback", e)

    # Fallback
    if fallback_fn is not None and fallback_input is not None:
        logger.info("Using rule-based fallback for %s", output_schema.__name__)
        fallback_data = fallback_fn(fallback_input)
        return output_schema.model_validate(fallback_data)

    raise RuntimeError(
        f"Ollama unavailable and no fallback provided for {output_schema.__name__}"
    )
