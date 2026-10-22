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
