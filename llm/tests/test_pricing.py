import json
from types import SimpleNamespace

import httpx
import pydantic
import pytest

from llm import pricing, registry, run_task
from llm.providers import entra_gateway_client
from llm.providers import profiles as profiles_module
from llm.providers.base import LLMResult
from llm.providers.entra_gateway_client import EntraGatewayClient, EntraGatewayGeminiClient, EntraTokenProvider
from llm.testing import FakeClient

GATEWAY_CREDENTIALS = {
    "llm_gateway_client_id": "test-client-id",
    "llm_gateway_client_secret": "test-client-secret",
    "llm_gateway_api_key": "test-gateway-api-key",
}


class Answer(pydantic.BaseModel):
    answer: str


def _result(input_tokens, output_tokens, cache_read=None, cache_write=None, **fields) -> LLMResult:
    return LLMResult(
        parsed=None, raw_text="", provider="fake", model="fake-model",
        input_tokens=input_tokens, output_tokens=output_tokens,
        cache_read_input_tokens=cache_read, cache_creation_input_tokens=cache_write,
        latency_ms=1.0, stop_reason="end_turn", **fields,
    )


def test_estimate_cost_plain_input_and_output():
    prices = pricing.ModelPricing(input=2.0, output=10.0)
    assert pricing.estimate_cost(_result(1_000_000, 100_000), prices) == pytest.approx(2.0 + 1.0)


def test_estimate_cost_prices_cache_reads_and_writes_separately():
    prices = pricing.ModelPricing(input=4.0, output=20.0, cache_write=5.0, cache_read=0.20)
    # 1M prompt tokens: 600k read from cache, 100k written to it, 300k uncached
    result = _result(1_000_000, 0, cache_read=600_000, cache_write=100_000)
    assert pricing.estimate_cost(result, prices) == pytest.approx(0.3 * 4.0 + 0.6 * 0.20 + 0.1 * 5.0)


def test_estimate_cost_cache_rates_fall_back_to_input_rate():
    prices = pricing.ModelPricing(input=1.0, output=1.0)
    result = _result(1_000_000, 0, cache_read=500_000, cache_write=500_000)
    assert pricing.estimate_cost(result, prices) == pytest.approx(1.0)


def test_estimate_cost_unknown_pricing():
    assert pricing.estimate_cost(_result(10, 10), None) is None
    assert pricing.format_cost(None) == "unknown"
    assert pricing.format_cost(0.031987) == "$0.0320"


def test_gateway_cost_converts_per_1k_to_per_million():
    prices = pricing.ModelPricing.from_gateway_cost(
        {"currency": "USD", "input_per_1k_tokens": 0.002, "output_per_1k_tokens": 0.006,
         "cached_input_per_1k_tokens": 0.0005}
    )
    assert (prices.input, prices.output, prices.cache_read) == pytest.approx((2.0, 6.0, 0.5))
    assert prices.cache_write is None
    assert prices.source == "gateway discovery"


@pytest.mark.parametrize("cost", [None, {}, {"currency": "EUR", "input_per_1k_tokens": 1, "output_per_1k_tokens": 1},
                                  {"currency": "USD", "input_per_1k_tokens": 1}])
def test_gateway_cost_unusable_gives_none(cost):
    assert pricing.ModelPricing.from_gateway_cost(cost) is None


def test_claude_profile_pricing_table_parses():
    profile = profiles_module.get_profile("claude")
    assert "pricing" in profile.pricing_source
    prices = profiles_module._profile_pricing(profile, profile.model)
    assert prices is not None and prices.input > 0 and prices.output > prices.input
    assert profiles_module._profile_pricing(profile, "not-a-real-model") is None


def test_get_client_attaches_gateway_discovery_pricing(monkeypatch):
    discovery = {"models": [{
        "id": "gpt-5.4-mini", "api_format": "openai", "capabilities": ["chat"], "app_has_access": True,
        "cost": {"currency": "USD", "input_per_1k_tokens": 0.001, "output_per_1k_tokens": 0.004},
    }]}
    monkeypatch.setattr(EntraTokenProvider, "get_token", lambda self: "token-x")
    monkeypatch.setattr(entra_gateway_client, "fetch_discovery", lambda *args, **kwargs: discovery)

    client = profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS)

    assert (client.pricing.input, client.pricing.output) == pytest.approx((1.0, 4.0))


def test_get_client_gateway_discovery_failure_leaves_pricing_unknown(monkeypatch):
    def fail(*args, **kwargs):
        raise RuntimeError("discovery down")

    monkeypatch.setattr(EntraTokenProvider, "get_token", lambda self: "token-x")
    monkeypatch.setattr(entra_gateway_client, "fetch_discovery", fail)

    client = profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS)

    assert client.pricing is None


def test_run_task_sets_and_logs_estimated_cost(review_evidence_item_input, tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    parsed = task.schema(summary="fine", revision_assessments=[])
    canned = _result(1000, 500)
    canned.parsed, canned.raw_text = parsed, parsed.model_dump_json()
    client = FakeClient(responses=[canned])
    client.pricing = pricing.ModelPricing(input=2.0, output=10.0)
    log_path = tmp_path / "provenance.jsonl"

    result = run_task("review_evidence_item", inputs=review_evidence_item_input, profile="claude",
                      log_path=log_path, client=client)

    assert result.estimated_cost_usd == pytest.approx((1000 * 2.0 + 500 * 10.0) / 1_000_000)
    record = json.loads(log_path.read_text().strip())
    assert record["estimated_cost_usd"] == pytest.approx(result.estimated_cost_usd)


def _token_provider() -> EntraTokenProvider:
    provider = EntraTokenProvider(token_url="https://login.example/token", client_id="id",
                                  client_secret="secret", scope="api://example/.default")
    provider.get_token = lambda: "token-x"
    return provider


def test_openai_format_usage_is_summed_over_every_response_and_counts_cached_tokens():
    # (the SDK's .parse() validates responses itself, so a retry is exercised via _to_result directly)
    def response(prompt, completion, cached):
        return SimpleNamespace(
            usage=SimpleNamespace(prompt_tokens=prompt, completion_tokens=completion,
                                  prompt_tokens_details=SimpleNamespace(cached_tokens=cached)),
            choices=[SimpleNamespace(finish_reason="stop")],
        )

    client = EntraGatewayClient(model="m", base_url="https://gateway.example/models/v2",
                                token_provider=_token_provider(), api_key="k")
    result = client._to_result([response(100, 10, 40), response(150, 20, 100)],
                               parsed=None, raw_text="", latency_ms=1.0, retries=1)

    assert (result.input_tokens, result.output_tokens, result.cache_read_input_tokens) == (250, 30, 140)


def test_openai_format_without_cached_token_details_reports_none():
    def handler(request: httpx.Request) -> httpx.Response:
        return httpx.Response(200, json={
            "id": "c", "object": "chat.completion", "created": 0, "model": "m",
            "choices": [{"index": 0, "finish_reason": "stop", "message": {"role": "assistant", "content": "hi"}}],
            "usage": {"prompt_tokens": 3, "completion_tokens": 2, "total_tokens": 5},
        })

    client = EntraGatewayClient(
        model="m", base_url="https://gateway.example/models/v2", token_provider=_token_provider(),
        api_key="k", http_client=httpx.Client(transport=httpx.MockTransport(handler)),
    )
    result = client.complete(system="s", messages=[{"role": "user", "content": "hi"}])

    assert result.cache_read_input_tokens is None


def test_anthropic_format_input_tokens_include_cache_and_are_summed_over_every_response():
    from llm.providers.anthropic_client import AnthropicClient

    def response(input_tokens, output_tokens, cache_read, cache_write):
        return SimpleNamespace(
            usage=SimpleNamespace(input_tokens=input_tokens, output_tokens=output_tokens,
                                  cache_read_input_tokens=cache_read, cache_creation_input_tokens=cache_write),
            stop_reason="end_turn",
        )

    client = AnthropicClient(model="claude-sonnet-5", api_key="k")
    result = client._to_result([response(10, 5, 0, 1000), response(30, 7, 1000, None)],
                               parsed=None, raw_text="", latency_ms=1.0, retries=1)

    assert result.input_tokens == 10 + 1000 + 30 + 1000
    assert (result.cache_read_input_tokens, result.cache_creation_input_tokens) == (1000, 1000)
    assert result.output_tokens == 12


def test_anthropic_format_falls_back_to_openai_style_usage_fields():
    # the WashU gateway's Anthropic-format endpoint reports prompt_tokens/completion_tokens instead
    from llm.providers.anthropic_client import AnthropicClient

    response = SimpleNamespace(
        usage=SimpleNamespace(input_tokens=None, output_tokens=None, cache_read_input_tokens=None,
                              cache_creation_input_tokens=None, prompt_tokens=9, completion_tokens=31),
        stop_reason="end_turn",
    )
    client = AnthropicClient(model="claude-sonnet-5", api_key="k")
    result = client._to_result([response], parsed=None, raw_text="", latency_ms=1.0, retries=0)

    assert (result.input_tokens, result.output_tokens) == (9, 31)


def test_gemini_usage_is_summed_over_a_validation_retry():
    class Response:
        status_code = 200

        def __init__(self, text, cached):
            usage = {"promptTokenCount": 50, "candidatesTokenCount": 4, "thoughtsTokenCount": 6}
            if cached is not None:
                usage["cachedContentTokenCount"] = cached
            self._body = {"candidates": [{"content": {"parts": [{"text": text}]}, "finishReason": "STOP"}],
                          "usageMetadata": usage}

        def json(self):
            return self._body

    class Session:
        def __init__(self):
            self.responses = [Response('{"wrong": 1}', None), Response('{"answer": "ok"}', 20)]

        def post(self, *args, **kwargs):
            return self.responses.pop(0)

    client = EntraGatewayGeminiClient(model="gemini-2.5-flash", base_url="https://gateway.example/models/v2",
                                      token_provider=_token_provider(), api_key="k", session=Session())
    result = client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)

    assert result.retries == 1
    assert (result.input_tokens, result.output_tokens, result.cache_read_input_tokens) == (100, 20, 20)
