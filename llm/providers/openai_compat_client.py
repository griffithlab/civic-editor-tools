import json
import time
from typing import Optional

import openai
from openai import OpenAI

from .base import LLMClient, LLMRefusalError, LLMResult, validate_with_retry


def _is_unsupported_response_format_error(exc: openai.BadRequestError) -> bool:
    """Narrow check for exactly the case this fallback exists for: a server (most local
    OpenAI-compatible servers today - Ollama, LM Studio, vLLM) that doesn't support strict
    JSON-schema response_format. Never treat an arbitrary 400 as a reason to fall back."""
    if getattr(exc, "status_code", None) not in (400, 422):
        return False
    message = str(exc).lower()
    return "response_format" in message or "json_schema" in message


class OpenAICompatibleClient(LLMClient):
    """LLMClient backed by the `openai` SDK's chat.completions API, pointed at a configurable
    base_url - covers OpenAI itself and any OpenAI-compatible local server (Ollama, LM Studio,
    vLLM).

    Structured output tries chat.completions.parse(response_format=<PydanticModel>) first. If the
    server rejects that (most local servers today), falls back to JSON mode
    (response_format={"type": "json_object"}) with the target JSON schema spelled out in the
    system prompt, then validates the result with validate_with_retry().
    """

    provider_name = "openai_compatible"

    def __init__(self, model: str, base_url: str, api_key: Optional[str] = None, **default_params):
        self.model = model
        self.base_url = base_url
        self.default_params = default_params
        self._client = OpenAI(base_url=base_url, api_key=api_key or "not-needed")

    def complete(self, *, system: str, messages: list, output_schema=None, **opts) -> LLMResult:
        max_tokens = opts.pop("max_tokens", self.default_params.get("max_tokens", 4096))
        full_messages = [{"role": "system", "content": system}] + list(messages)

        start = time.monotonic()

        if output_schema is None:
            response = self._client.chat.completions.create(
                model=self.model, messages=full_messages, max_tokens=max_tokens,
            )
            self._require_finish_reason(response.choices[0])
            raw_text = response.choices[0].message.content or ""
            return self._to_result(
                response, parsed=None, raw_text=raw_text,
                latency_ms=(time.monotonic() - start) * 1000, retries=0,
            )

        used_fallback = False
        try:
            response = self._client.chat.completions.parse(
                model=self.model, messages=full_messages, max_tokens=max_tokens,
                response_format=output_schema,
            )
        except openai.BadRequestError as exc:
            if not _is_unsupported_response_format_error(exc):
                raise
            used_fallback = True
            response = self._complete_via_json_mode(full_messages, output_schema, max_tokens)

        self._require_finish_reason(response.choices[0])
        raw_text = response.choices[0].message.content or ""

        state = {"raw_text": raw_text, "response": response}

        def retry_fn(error_message: str) -> str:
            corrected_messages = full_messages + [
                {"role": "assistant", "content": state["raw_text"]},
                {
                    "role": "user",
                    "content": (
                        f"That response failed schema validation with this error:\n{error_message}\n"
                        f"Reply again with ONLY corrected JSON matching the schema, nothing else."
                    ),
                },
            ]
            if used_fallback:
                new_response = self._complete_via_json_mode(corrected_messages, output_schema, max_tokens)
            else:
                new_response = self._client.chat.completions.parse(
                    model=self.model, messages=corrected_messages, max_tokens=max_tokens,
                    response_format=output_schema,
                )
            self._require_finish_reason(new_response.choices[0])
            new_raw = new_response.choices[0].message.content or ""
            state["raw_text"] = new_raw
            state["response"] = new_response
            return new_raw

        parsed, retries = validate_with_retry(output_schema, raw_text, retry_fn)

        return self._to_result(
            state["response"], parsed=parsed, raw_text=state["raw_text"],
            latency_ms=(time.monotonic() - start) * 1000, retries=retries,
        )

    def _complete_via_json_mode(self, full_messages, output_schema, max_tokens):
        schema_json = json.dumps(output_schema.model_json_schema())
        augmented_messages = list(full_messages)
        augmented_messages[0] = {
            "role": "system",
            "content": (
                full_messages[0]["content"]
                + "\n\nRespond with ONLY a single JSON object matching this JSON Schema, "
                "and no other text:\n" + schema_json
            ),
        }
        return self._client.chat.completions.create(
            model=self.model, messages=augmented_messages, max_tokens=max_tokens,
            response_format={"type": "json_object"},
        )

    @staticmethod
    def _require_finish_reason(choice) -> None:
        """Only "stop" is safe to trust as a complete, schema-compliant response. Anything else
        ("length" = truncated, "content_filter" = refused, "tool_calls" = unexpected here) is a
        hard failure, not something validate_with_retry() should try to fix."""
        if choice.finish_reason != "stop":
            raise LLMRefusalError(stop_reason=choice.finish_reason, raw_text=choice.message.content or "")

    def _to_result(self, response, *, parsed, raw_text, latency_ms, retries) -> LLMResult:
        usage = response.usage
        return LLMResult(
            parsed=parsed,
            raw_text=raw_text,
            provider=self.provider_name,
            model=self.model,
            input_tokens=usage.prompt_tokens if usage else 0,
            output_tokens=usage.completion_tokens if usage else 0,
            cache_read_input_tokens=None,
            cache_creation_input_tokens=None,
            latency_ms=latency_ms,
            stop_reason=response.choices[0].finish_reason,
            retries=retries,
        )
