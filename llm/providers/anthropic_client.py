import time
from typing import Optional

import anthropic

from .base import LLMClient, LLMRefusalError, LLMResult, validate_with_retry

_FORCED_TOOL_NAME = "emit_result"


def _build_system_param(system: str, cache_system: bool):
    """Wrap the system prompt as a single cacheable content block by default. Callers that pass
    cache_system=False get a plain string instead (e.g. for very short/one-off prompts where a
    cache write wouldn't pay for itself)."""
    if not cache_system:
        return system
    return [{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}]


def _extract_text(response) -> str:
    """Concatenate every text block in a Message/ParsedMessage's content."""
    return "".join(block.text for block in response.content if block.type == "text")


def _extract_tool_input_json(response, tool_name: str) -> str:
    import json

    for block in response.content:
        if block.type == "tool_use" and block.name == tool_name:
            return json.dumps(block.input)
    raise RuntimeError(f"Forced tool-use response did not contain a '{tool_name}' tool call")


def _is_unsupported_structured_output_error(exc: anthropic.BadRequestError) -> bool:
    """Narrow check for exactly the case this fallback exists for: an older SDK/API surface (or
    a non-standard proxy) that doesn't recognize output_format/output_config. Never treat an
    arbitrary 400 (e.g. a bad model name, an invalid message) as a reason to silently retry via a
    completely different mechanism."""
    if getattr(exc, "status_code", None) != 400:
        return False
    message = str(exc).lower()
    return "output_format" in message or "output_config" in message


class AnthropicClient(LLMClient):
    """LLMClient backed by the official `anthropic` SDK.

    Structured output uses the current (GA) mechanism, client.messages.parse(output_format=...),
    which returns response.parsed_output already validated against the given Pydantic model. If
    that raises a 400 specifically about output_format/output_config being unrecognized (an older
    pinned SDK/API version, or a non-standard proxy), this falls back to forced tool use - current
    Claude models do NOT reject forced tool_choice; structured outputs and strict tool use are
    documented as complementary, separate features, so this fallback is a defensive measure, not
    a workaround for current-model behavior.
    """

    provider_name = "anthropic"

    def __init__(self, model: str, api_key: Optional[str] = None, **default_params):
        self.model = model
        self.default_params = default_params
        self._client = anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()

    def complete(self, *, system: str, messages: list, output_schema=None, **opts) -> LLMResult:
        max_tokens = opts.pop("max_tokens", self.default_params.get("max_tokens", 4096))
        thinking = opts.pop("thinking", self.default_params.get("thinking"))
        cache_system = opts.pop("cache_system", True)

        system_param = _build_system_param(system, cache_system)
        extra = {"thinking": thinking} if thinking else {}

        start = time.monotonic()

        if output_schema is None:
            response = self._client.messages.create(
                model=self.model, max_tokens=max_tokens, system=system_param, messages=messages, **extra,
            )
            self._require_stop_reason(response, {"end_turn"})
            return self._to_result(
                response, parsed=None, raw_text=_extract_text(response),
                latency_ms=(time.monotonic() - start) * 1000, retries=0,
            )

        used_fallback = False
        try:
            response = self._client.messages.parse(
                model=self.model, max_tokens=max_tokens, system=system_param, messages=messages,
                output_format=output_schema, **extra,
            )
        except anthropic.BadRequestError as exc:
            if not _is_unsupported_structured_output_error(exc):
                raise
            used_fallback = True
            response = self._complete_via_forced_tool(system_param, messages, output_schema, max_tokens, extra)

        self._require_stop_reason(response, {"tool_use"} if used_fallback else {"end_turn"})
        raw_text = _extract_tool_input_json(response, _FORCED_TOOL_NAME) if used_fallback else _extract_text(response)

        state = {"raw_text": raw_text, "response": response}

        def retry_fn(error_message: str) -> str:
            corrected_messages = messages + [
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
                new_response = self._complete_via_forced_tool(
                    system_param, corrected_messages, output_schema, max_tokens, extra
                )
                self._require_stop_reason(new_response, {"tool_use"})
                new_raw = _extract_tool_input_json(new_response, _FORCED_TOOL_NAME)
            else:
                new_response = self._client.messages.parse(
                    model=self.model, max_tokens=max_tokens, system=system_param,
                    messages=corrected_messages, output_format=output_schema, **extra,
                )
                self._require_stop_reason(new_response, {"end_turn"})
                new_raw = _extract_text(new_response)

            state["raw_text"] = new_raw
            state["response"] = new_response
            return new_raw

        parsed, retries = validate_with_retry(output_schema, raw_text, retry_fn)

        return self._to_result(
            state["response"], parsed=parsed, raw_text=state["raw_text"],
            latency_ms=(time.monotonic() - start) * 1000, retries=retries,
        )

    def _complete_via_forced_tool(self, system_param, messages, output_schema, max_tokens, extra):
        tool = {
            "name": _FORCED_TOOL_NAME,
            "description": f"Emit the final {output_schema.__name__} result.",
            "input_schema": output_schema.model_json_schema(),
        }
        return self._client.messages.create(
            model=self.model, max_tokens=max_tokens, system=system_param, messages=messages,
            tools=[tool], tool_choice={"type": "tool", "name": _FORCED_TOOL_NAME}, **extra,
        )

    @staticmethod
    def _require_stop_reason(response, expected: set) -> None:
        """Only end_turn (or tool_use, for the forced-tool fallback) is safe to trust: on
        'refusal' the docs state the refusal message takes precedence over schema constraints,
        and on 'max_tokens' the output may be truncated/invalid. Treat any other stop_reason as
        a hard failure rather than something validate_with_retry() should try to fix."""
        if response.stop_reason not in expected:
            raise LLMRefusalError(stop_reason=response.stop_reason, raw_text=_extract_text(response))

    def _to_result(self, response, *, parsed, raw_text, latency_ms, retries) -> LLMResult:
        usage = response.usage
        return LLMResult(
            parsed=parsed,
            raw_text=raw_text,
            provider=self.provider_name,
            model=self.model,
            input_tokens=usage.input_tokens,
            output_tokens=usage.output_tokens,
            cache_read_input_tokens=usage.cache_read_input_tokens,
            cache_creation_input_tokens=usage.cache_creation_input_tokens,
            latency_ms=latency_ms,
            stop_reason=response.stop_reason,
            retries=retries,
        )
