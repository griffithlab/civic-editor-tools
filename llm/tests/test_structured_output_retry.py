"""End-to-end (mocked HTTP) checks that a reply failing schema validation is retried once by
validate_with_retry() on the primary structured-output path of the openai- and anthropic-SDK
clients - those SDKs' own .parse() helpers validate internally and would raise first."""

import json

import httpx
import pydantic
import pytest

from llm.providers.base import LLMOutputValidationError
from llm.providers.entra_gateway_client import EntraGatewayAnthropicClient, EntraGatewayClient, EntraTokenProvider


class Answer(pydantic.BaseModel):
    answer: str


def _token_provider() -> EntraTokenProvider:
    provider = EntraTokenProvider(token_url="https://login.example/token", client_id="id",
                                  client_secret="secret", scope="api://example/.default")
    provider.get_token = lambda: "token-x"
    return provider


def _openai_client(replies, seen_bodies):
    replies = iter(replies)

    def handler(request: httpx.Request) -> httpx.Response:
        seen_bodies.append(json.loads(request.content))
        content, prompt_tokens, completion_tokens = next(replies)
        return httpx.Response(200, json={
            "id": "c", "object": "chat.completion", "created": 0, "model": "m",
            "choices": [{"index": 0, "finish_reason": "stop", "message": {"role": "assistant", "content": content}}],
            "usage": {"prompt_tokens": prompt_tokens, "completion_tokens": completion_tokens,
                      "total_tokens": prompt_tokens + completion_tokens},
        })

    return EntraGatewayClient(
        model="m", base_url="https://gateway.example/models/v2", token_provider=_token_provider(),
        api_key="k", http_client=httpx.Client(transport=httpx.MockTransport(handler)),
    )


def test_openai_format_schema_mismatch_is_retried_once():
    seen = []
    client = _openai_client([('{"wrong": 1}', 100, 10), ('{"answer": "ok"}', 150, 20)], seen)

    result = client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)

    assert result.parsed.answer == "ok"
    assert result.retries == 1
    assert (result.input_tokens, result.output_tokens) == (250, 30)
    assert seen[0]["response_format"]["type"] == "json_schema"
    assert seen[0]["response_format"]["json_schema"]["strict"] is True
    # the retry shows the model its rejected reply and the validation error
    assert seen[1]["messages"][-2] == {"role": "assistant", "content": '{"wrong": 1}'}
    assert "failed schema validation" in seen[1]["messages"][-1]["content"]


def test_openai_format_schema_mismatch_twice_raises_validation_error():
    client = _openai_client([('{"wrong": 1}', 1, 1), ("not json", 1, 1)], [])

    with pytest.raises(LLMOutputValidationError) as excinfo:
        client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)

    assert excinfo.value.raw_output == "not json"


def _anthropic_http_module():
    """The HTTP library the installed anthropic SDK is built on (newer releases use httpx2)."""
    try:
        import httpx2
        return httpx2
    except ImportError:
        return httpx


def _anthropic_client(replies, seen_bodies):
    http = _anthropic_http_module()
    replies = iter(replies)

    def handler(request):
        seen_bodies.append(json.loads(request.content))
        text = next(replies)
        return http.Response(200, json={
            "id": "msg", "type": "message", "role": "assistant", "model": "claude-sonnet-5",
            "content": [{"type": "text", "text": text}], "stop_reason": "end_turn", "stop_sequence": None,
            "usage": {"input_tokens": 10, "output_tokens": 5,
                      "cache_read_input_tokens": 0, "cache_creation_input_tokens": 0},
        })

    return EntraGatewayAnthropicClient(
        model="claude-sonnet-5", base_url="https://gateway.example/models/v2", token_provider=_token_provider(),
        api_key="k", http_client=http.Client(transport=http.MockTransport(handler)),
    )


def test_anthropic_format_schema_mismatch_is_retried_once():
    seen = []
    client = _anthropic_client(['{"wrong": 1}', '{"answer": "ok"}'], seen)

    result = client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)

    assert result.parsed.answer == "ok"
    assert result.retries == 1
    assert (result.input_tokens, result.output_tokens) == (20, 10)
    output_format = seen[0]["output_config"]["format"]
    assert output_format["type"] == "json_schema"
    assert "answer" in output_format["schema"]["properties"]
    assert seen[1]["messages"][-2] == {"role": "assistant", "content": '{"wrong": 1}'}


def test_anthropic_format_schema_mismatch_twice_raises_validation_error():
    client = _anthropic_client(['{"wrong": 1}', '{"also": "wrong"}'], [])

    with pytest.raises(LLMOutputValidationError):
        client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)


#--- the forced-tool fallback, used only when a server rejects output_config ---

_OUTPUT_CONFIG_UNSUPPORTED = {"type": "error", "error": {"type": "invalid_request_error",
                                                          "message": "output_config: Extra inputs are not permitted"}}
_FORCED_TOOL_CHOICE_UNSUPPORTED = {"type": "error", "error": {
    "type": "invalid_request_error",
    "message": 'tool_choice: type "tool" and "any" are not supported for this model.'}}


def _message(content, stop_reason):
    return {"id": "msg", "type": "message", "role": "assistant", "model": "claude-opus-5-5",
            "content": content, "stop_reason": stop_reason, "stop_sequence": None,
            "usage": {"input_tokens": 10, "output_tokens": 5,
                      "cache_read_input_tokens": 0, "cache_creation_input_tokens": 0}}


def _tool_call(answer):
    return _message([{"type": "tool_use", "id": "t1", "name": "emit_result", "input": {"answer": answer}}], "tool_use")


def _fallback_client(server, seen_bodies):
    """server(body) -> (status, json); an anthropic-format gateway client over a mock transport."""
    http = _anthropic_http_module()

    def handler(request):
        body = json.loads(request.content)
        seen_bodies.append(body)
        status, payload = server(body)
        return http.Response(status, json=payload)

    return EntraGatewayAnthropicClient(
        model="claude-opus-5-5", base_url="https://gateway.example/models/v2", token_provider=_token_provider(),
        api_key="k", http_client=http.Client(transport=http.MockTransport(handler)),
    )


def _no_output_config_server(forced_tool_choice_ok, auto_reply=None):
    def server(body):
        if "output_config" in body:
            return 400, _OUTPUT_CONFIG_UNSUPPORTED
        if body["tool_choice"]["type"] == "tool":
            return (200, _tool_call("forced")) if forced_tool_choice_ok else (400, _FORCED_TOOL_CHOICE_UNSUPPORTED)
        return 200, auto_reply or _tool_call("auto")
    return server


def _complete(client):
    return client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)


def test_tool_fallback_forces_the_tool_where_the_model_allows_it():
    seen = []
    result = _complete(_fallback_client(_no_output_config_server(forced_tool_choice_ok=True), seen))

    assert result.parsed.answer == "forced"
    assert [b["tool_choice"]["type"] for b in seen if "tools" in b] == ["tool"]


def test_tool_fallback_switches_to_auto_tool_choice_when_forcing_is_rejected():
    seen = []
    client = _fallback_client(_no_output_config_server(forced_tool_choice_ok=False), seen)

    assert _complete(client).parsed.answer == "auto"
    tool_requests = [b for b in seen if "tools" in b]
    assert [b["tool_choice"]["type"] for b in tool_requests] == ["tool", "auto"]
    # the instruction to call the tool goes after the cached system block, leaving it intact
    system = tool_requests[1]["system"]
    assert system[0]["text"] == "s" and "cache_control" in system[0]
    assert "emit_result" in system[-1]["text"]

    # remembered: the next call goes straight to auto
    seen.clear()
    assert _complete(client).parsed.answer == "auto"
    assert [b["tool_choice"]["type"] for b in seen if "tools" in b] == ["auto"]


def test_tool_fallback_auto_mode_accepts_a_plain_text_json_reply():
    seen = []
    text_reply = _message([{"type": "text", "text": '{"answer": "text"}'}], "end_turn")
    client = _fallback_client(_no_output_config_server(forced_tool_choice_ok=False, auto_reply=text_reply), seen)

    assert _complete(client).parsed.answer == "text"


def test_tool_fallback_other_bad_requests_are_not_swallowed():
    import anthropic

    def server(body):
        if "output_config" in body:
            return 400, _OUTPUT_CONFIG_UNSUPPORTED
        return 400, {"type": "error", "error": {"type": "invalid_request_error", "message": "messages: invalid"}}

    with pytest.raises(anthropic.BadRequestError):
        _complete(_fallback_client(server, []))
