"""Anthropic / Claude client with structured output and retry logic."""
from __future__ import annotations

import json
import logging
from typing import Any, Type

import anthropic
from pydantic import BaseModel
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from trading.config import settings

logger = logging.getLogger(__name__)

_client: anthropic.Anthropic | None = None


def _get_client() -> anthropic.Anthropic:
    global _client
    if _client is None:
        _client = anthropic.Anthropic(api_key=settings.ANTHROPIC_API_KEY)
    return _client


@retry(
    retry=retry_if_exception_type((anthropic.RateLimitError, anthropic.APIConnectionError)),
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=2, max=30),
)
def call(
    system_prompt: str,
    user_prompt: str,
    max_tokens: int = 1024,
    model: str | None = None,
) -> str:
    """Raw text call to Claude. Returns the assistant's response text."""
    client = _get_client()
    response = client.messages.create(
        model=model or settings.CLAUDE_MODEL,
        max_tokens=max_tokens,
        system=system_prompt,
        messages=[{"role": "user", "content": user_prompt}],
    )
    return response.content[0].text


def call_structured(
    system_prompt: str,
    user_prompt: str,
    output_schema: Type[BaseModel],
    max_tokens: int = 1024,
    model: str | None = None,
    max_fix_attempts: int = 2,
) -> BaseModel:
    """Call Claude and parse the response into a Pydantic model.

    If the first response fails to parse, asks Claude to fix its own output.
    """
    schema_json = json.dumps(output_schema.model_json_schema(), indent=2)
    enriched_system = (
        f"{system_prompt}\n\n"
        f"You MUST respond with valid JSON that matches this schema exactly:\n"
        f"```json\n{schema_json}\n```\n"
        f"Return ONLY the JSON object — no prose, no markdown fences."
    )

    raw = call(enriched_system, user_prompt, max_tokens=max_tokens, model=model)

    for attempt in range(max_fix_attempts + 1):
        try:
            # Strip potential markdown fences
            text = raw.strip()
            if text.startswith("```"):
                text = text.split("```")[1]
                if text.startswith("json"):
                    text = text[4:]
            data = json.loads(text)
            return output_schema.model_validate(data)
        except Exception as e:
            if attempt == max_fix_attempts:
                logger.error("Claude structured output failed after %d attempts: %s", attempt + 1, e)
                raise ValueError(f"Could not parse Claude output into {output_schema.__name__}: {e}") from e
            # Ask Claude to fix it
            fix_prompt = (
                f"Your previous response failed JSON schema validation with error: {e}\n"
                f"Original response: {raw}\n\n"
                f"Please return ONLY the corrected JSON object."
            )
            raw = call(enriched_system, fix_prompt, max_tokens=max_tokens, model=model)
