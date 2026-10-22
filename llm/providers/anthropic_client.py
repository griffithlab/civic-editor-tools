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
                [response], parsed=None, raw_text=_extract_text(response),
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

        state = {"raw_text": raw_text, "responses": [response]}

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
            state["responses"].append(new_response)
            return new_raw

        parsed, retries = validate_with_retry(output_schema, raw_text, retry_fn)

        return self._to_result(
            state["responses"], parsed=parsed, raw_text=state["raw_text"],
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

    def _to_result(self, responses, *, parsed, raw_text, latency_ms, retries) -> LLMResult:
        """Build the LLMResult from every response this call got (more than one after a
        validation retry), summing their usage. Anthropic reports cache reads/writes separately
        from input_tokens, so they're added back in to make input_tokens the total prompt size."""
        cache_read = sum(r.usage.cache_read_input_tokens or 0 for r in responses)
        cache_write = sum(r.usage.cache_creation_input_tokens or 0 for r in responses)
        return LLMResult(
            parsed=parsed,
            raw_text=raw_text,
            provider=self.provider_name,
            model=self.model,
            input_tokens=sum(_usage_tokens(r.usage)[0] for r in responses) + cache_read + cache_write,
            output_tokens=sum(_usage_tokens(r.usage)[1] for r in responses),
            cache_read_input_tokens=cache_read,
            cache_creation_input_tokens=cache_write,
            latency_ms=latency_ms,
            stop_reason=responses[-1].stop_reason,
            retries=retries,
        )


def _usage_tokens(usage) -> tuple:
    """(input_tokens, output_tokens) from a Message's usage. Some Anthropic-format proxies (e.g.
    the WashU gateway) leave those None and report OpenAI-style prompt_tokens/completion_tokens
    instead, which the SDK keeps as extra fields."""
    input_tokens = usage.input_tokens if usage.input_tokens is not None else getattr(usage, "prompt_tokens", None)
    output_tokens = usage.output_tokens if usage.output_tokens is not None else getattr(usage, "completion_tokens", None)
    return input_tokens or 0, output_tokens or 0


def check_model_available(sdk_client, model: str, profile_name: str):
    """
    Look up `model` with the Anthropic models API (free, no tokens used) so a misspelled or
    non-Anthropic model id, or a rejected API key, is reported clearly up front rather than as a
    404/401 on every review request. Returns (error, warning): error explains why the model
    can't be used, warning reports that the check itself couldn't be completed; both may be None.
    """
    import difflib

    try:
        sdk_client.models.retrieve(model)
        return None, None
    except anthropic.NotFoundError:
        pass
    except anthropic.AuthenticationError:
        return f"Anthropic rejected the API key for profile '{profile_name}'", None
    except Exception as exc:
        return None, f"could not check model '{model}' with Anthropic ({exc}); continuing"

    try:
        available = [m.id for m in list_models(sdk_client)]
    except Exception:
        available = []
    suggestions = difflib.get_close_matches(model, available, n=3, cutoff=0.6)
    hint = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
    if available:
        hint += f"\n  Models available to your API key: {', '.join(available)}"
    return (f"model '{model}' is not available from Anthropic (profile '{profile_name}' talks "
            f"directly to the Anthropic API, so --llm-model must be an Anthropic model id).{hint}"), None


def list_models(sdk_client) -> list:
    """Every model available to this Anthropic API key, newest first (the API's own order)."""
    return list(sdk_client.models.list(limit=1000))
