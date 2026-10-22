import json

import httpx
import pytest

from llm.providers import entra_gateway_client
from llm.providers import profiles as profiles_module
from llm.providers.entra_gateway_client import EntraGatewayClient, EntraTokenProvider, fetch_discovery

GATEWAY_CREDENTIALS = {
    "llm_gateway_client_id": "test-client-id",
    "llm_gateway_client_secret": "test-client-secret",
    "llm_gateway_api_key": "test-gateway-api-key",
}


class FakeTokenResponse:
    def __init__(self, status_code: int, body: dict):
        self.status_code = status_code
        self._body = body
        self.text = json.dumps(body)

    def json(self):
        return self._body

    def raise_for_status(self):
        if self.status_code >= 400:
            raise RuntimeError(f"HTTP {self.status_code}")


class FakeTokenSession:
    """Stands in for `requests`: hands out token-1, token-2, ... and records each POST."""

    def __init__(self, expires_in: int = 3600, status_code: int = 200, error_body: dict = None):
        self.expires_in = expires_in
        self.status_code = status_code
        self.error_body = error_body
        self.posts = []
        self.gets = []

    def post(self, url, data=None, timeout=None):
        self.posts.append({"url": url, "data": data})
        if self.status_code != 200:
            return FakeTokenResponse(self.status_code, self.error_body or {})
        return FakeTokenResponse(200, {"access_token": f"token-{len(self.posts)}", "expires_in": self.expires_in})

    def get(self, url, headers=None, timeout=None):
        self.gets.append({"url": url, "headers": headers})
        return FakeTokenResponse(200, {"models": [{"id": "gpt-5.4-mini"}]})


def _token_provider(session) -> EntraTokenProvider:
    return EntraTokenProvider(
        token_url="https://login.example/token", client_id="test-client-id",
        client_secret="test-client-secret", scope="api://example/.default", session=session,
    )


def _chat_completion(content: str) -> dict:
    return {
        "id": "chatcmpl-1", "object": "chat.completion", "created": 0, "model": "gpt-5.4-mini",
        "choices": [{"index": 0, "finish_reason": "stop", "message": {"role": "assistant", "content": content}}],
        "usage": {"prompt_tokens": 3, "completion_tokens": 2, "total_tokens": 5},
    }


def test_token_provider_sends_client_credentials_and_caches_token():
    session = FakeTokenSession()
    provider = _token_provider(session)

    assert provider.get_token() == "token-1"
    assert provider.get_token() == "token-1"
    assert len(session.posts) == 1
    assert session.posts[0]["data"] == {
        "grant_type": "client_credentials",
        "client_id": "test-client-id",
        "client_secret": "test-client-secret",
        "scope": "api://example/.default",
    }


def test_token_provider_refreshes_near_expiry(monkeypatch):
    clock = {"now": 1000.0}
    monkeypatch.setattr(entra_gateway_client.time, "monotonic", lambda: clock["now"])
    session = FakeTokenSession(expires_in=3600)
    provider = _token_provider(session)

    assert provider.get_token() == "token-1"
    clock["now"] += 3600 - entra_gateway_client.TOKEN_REFRESH_MARGIN_SECONDS - 1
    assert provider.get_token() == "token-1"
    clock["now"] += 2
    assert provider.get_token() == "token-2"


def test_token_provider_error_surfaces_entra_description_not_secret():
    session = FakeTokenSession(
        status_code=401,
        error_body={"error": "invalid_client", "error_description": "AADSTS7000215: Invalid client secret provided."},
    )
    with pytest.raises(RuntimeError) as exc_info:
        _token_provider(session).get_token()

    assert "AADSTS7000215" in str(exc_info.value)
    assert "test-client-secret" not in str(exc_info.value)


def test_gateway_client_sends_bearer_token_and_api_key():
    seen_requests = []

    def handler(request: httpx.Request) -> httpx.Response:
        seen_requests.append(request)
        return httpx.Response(200, json=_chat_completion("hello"))

    client = EntraGatewayClient(
        model="gpt-5.4-mini", base_url="https://gateway.example/models/v2",
        token_provider=_token_provider(FakeTokenSession()), api_key="test-gateway-api-key",
        http_client=httpx.Client(transport=httpx.MockTransport(handler)),
    )
    result = client.complete(system="be brief", messages=[{"role": "user", "content": "hi"}])

    assert result.raw_text == "hello"
    assert result.provider == "entra_gateway"
    request = seen_requests[0]
    assert str(request.url) == "https://gateway.example/models/v2/chat/completions"
    assert request.headers["Authorization"] == "Bearer token-1"
    assert request.headers["X-Api-Key"] == "test-gateway-api-key"
    assert json.loads(request.content)["model"] == "gpt-5.4-mini"


def test_gateway_client_retries_once_with_fresh_token_on_401():
    seen_auth = []

    def handler(request: httpx.Request) -> httpx.Response:
        seen_auth.append(request.headers["Authorization"])
        if request.headers["Authorization"] == "Bearer token-1":
            return httpx.Response(401, json={"error": {"message": "token expired"}})
        return httpx.Response(200, json=_chat_completion("ok"))

    session = FakeTokenSession()
    client = EntraGatewayClient(
        model="gpt-5.4-mini", base_url="https://gateway.example/models/v2",
        token_provider=_token_provider(session), api_key="test-gateway-api-key",
        http_client=httpx.Client(transport=httpx.MockTransport(handler)),
    )
    result = client.complete(system="s", messages=[{"role": "user", "content": "hi"}])

    assert result.raw_text == "ok"
    assert seen_auth == ["Bearer token-1", "Bearer token-2"]
    assert len(session.posts) == 2


def test_fetch_discovery_sends_both_credentials():
    session = FakeTokenSession()
    body = fetch_discovery(
        "https://gateway.example/discovery", _token_provider(session), "test-gateway-api-key", session=session,
    )

    assert body == {"models": [{"id": "gpt-5.4-mini"}]}
    assert session.gets[0]["headers"] == {"Authorization": "Bearer token-1", "X-Api-Key": "test-gateway-api-key"}


def test_get_client_wustl_profile_requires_conf_credentials():
    with pytest.raises(ValueError) as exc_info:
        profiles_module.get_client("wustl", credentials={"llm_gateway_client_id": "x"})

    message = str(exc_info.value)
    assert "llm_gateway_client_secret" in message
    assert "llm_gateway_api_key" in message
    assert "*.civic.conf" in message


SAMPLE_DISCOVERY = {
    "models": [
        {"id": "gpt-5.4-mini", "api_format": "openai", "capabilities": ["chat", "tools"], "app_has_access": True},
        {"id": "gpt-5.4-nano", "api_format": "openai", "capabilities": ["chat"], "app_has_access": True},
        {"id": "claude-sonnet-4-6", "api_format": "anthropic", "capabilities": ["chat"], "app_has_access": True},
        {"id": "claude-opus-5-5-direct", "api_format": "anthropic", "capabilities": ["chat"], "app_has_access": False},
        {"id": "gemini-2.5-flash", "api_format": "gemini", "capabilities": ["chat"], "app_has_access": True},
        {"id": "text-embedding-3-small", "api_format": "openai", "capabilities": ["embeddings"], "app_has_access": True},
    ],
}


def _patch_gateway_network(monkeypatch, discovery=SAMPLE_DISCOVERY):
    monkeypatch.setattr(EntraTokenProvider, "get_token", lambda self: "token-x")
    monkeypatch.setattr(entra_gateway_client, "fetch_discovery", lambda *args, **kwargs: discovery)


def test_get_client_wustl_profile_with_model_override(monkeypatch):
    _patch_gateway_network(monkeypatch)

    client = profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="gpt-5.4-nano")

    assert isinstance(client, EntraGatewayClient)
    assert client.model == "gpt-5.4-nano"
    assert client.base_url == "https://aiapi.wustl.edu/models/v2"


def test_get_client_wustl_misspelled_model_suggests_close_match(monkeypatch):
    _patch_gateway_network(monkeypatch)

    with pytest.raises(ValueError) as exc_info:
        profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="gpt-5.4-min")

    message = str(exc_info.value)
    assert "not offered by the gateway" in message
    assert "gpt-5.4-mini" in message


def test_get_client_wustl_model_without_access_says_so(monkeypatch):
    _patch_gateway_network(monkeypatch)

    with pytest.raises(ValueError) as exc_info:
        profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="claude-opus-5-5-direct")

    assert "does not have access" in str(exc_info.value)


def test_get_client_wustl_rejects_non_chat_model(monkeypatch):
    _patch_gateway_network(monkeypatch)

    with pytest.raises(ValueError) as exc_info:
        profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="text-embedding-3-small")

    assert "not a chat model" in str(exc_info.value)


def test_get_client_wustl_routes_anthropic_format_models_to_messages_client(monkeypatch):
    from llm.providers.entra_gateway_client import EntraGatewayAnthropicClient

    _patch_gateway_network(monkeypatch)

    client = profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="claude-sonnet-4-6")

    assert isinstance(client, EntraGatewayAnthropicClient)
    assert client.model == "claude-sonnet-4-6"
    assert client.base_url == "https://aiapi.wustl.edu/models/v2"


def test_get_client_wustl_routes_gemini_format_models_to_generate_content_client(monkeypatch):
    from llm.providers.entra_gateway_client import EntraGatewayGeminiClient

    _patch_gateway_network(monkeypatch)

    client = profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="gemini-2.5-flash")

    assert isinstance(client, EntraGatewayGeminiClient)
    assert client.model == "gemini-2.5-flash"
    assert client._url == "https://aiapi.wustl.edu/models/v2/generateContent"


def test_get_client_wustl_rejects_unknown_api_format(monkeypatch):
    discovery = {"models": [{"id": "odd-model", "api_format": "something-new", "capabilities": ["chat"], "app_has_access": True}]}
    _patch_gateway_network(monkeypatch, discovery=discovery)

    with pytest.raises(ValueError) as exc_info:
        profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="odd-model")

    assert "'something-new' API format" in str(exc_info.value)


def _anthropic_http_module():
    """The HTTP library the installed anthropic SDK is built on (newer releases use httpx2, older
    ones httpx) - its http_client argument must come from the same package."""
    try:
        import httpx2
        return httpx2
    except ImportError:
        return httpx


def test_anthropic_format_gateway_client_sends_both_credentials_to_messages_endpoint():
    import pydantic
    from llm.providers.entra_gateway_client import EntraGatewayAnthropicClient

    http = _anthropic_http_module()

    class Answer(pydantic.BaseModel):
        answer: str

    seen_requests = []

    def handler(request: http.Request) -> http.Response:
        seen_requests.append(request)
        return http.Response(200, json={
            "id": "msg_1", "type": "message", "role": "assistant", "model": "claude-sonnet-4-6",
            "content": [{"type": "text", "text": '{"answer": "ok"}'}],
            "stop_reason": "end_turn", "stop_sequence": None,
            "usage": {"input_tokens": 3, "output_tokens": 2,
                      "cache_read_input_tokens": 0, "cache_creation_input_tokens": 0},
        })

    client = EntraGatewayAnthropicClient(
        model="claude-sonnet-4-6", base_url="https://gateway.example/models/v2",
        token_provider=_token_provider(FakeTokenSession()), api_key="test-gateway-api-key",
        http_client=http.Client(transport=http.MockTransport(handler)),
    )
    result = client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)

    assert result.parsed.answer == "ok"
    assert result.provider == "entra_gateway"
    request = seen_requests[0]
    assert str(request.url) == "https://gateway.example/models/v2/messages"
    assert request.headers["X-Api-Key"] == "test-gateway-api-key"
    assert request.headers["Authorization"] == "Bearer token-1"
    assert json.loads(request.content)["model"] == "claude-sonnet-4-6"


def test_anthropic_format_gateway_client_retries_once_with_fresh_token_on_401():
    from llm.providers.entra_gateway_client import EntraGatewayAnthropicClient

    http = _anthropic_http_module()

    seen_auth = []

    def handler(request: http.Request) -> http.Response:
        seen_auth.append(request.headers["Authorization"])
        if request.headers["Authorization"] == "Bearer token-1":
            return http.Response(401, json={"type": "error", "error": {"type": "authentication_error", "message": "expired"}})
        return http.Response(200, json={
            "id": "msg_1", "type": "message", "role": "assistant", "model": "claude-sonnet-4-6",
            "content": [{"type": "text", "text": "hello"}], "stop_reason": "end_turn", "stop_sequence": None,
            "usage": {"input_tokens": 1, "output_tokens": 1},
        })

    session = FakeTokenSession()
    client = EntraGatewayAnthropicClient(
        model="claude-sonnet-4-6", base_url="https://gateway.example/models/v2",
        token_provider=_token_provider(session), api_key="test-gateway-api-key",
        http_client=http.Client(transport=http.MockTransport(handler)),
    )
    result = client.complete(system="s", messages=[{"role": "user", "content": "hi"}])

    assert result.raw_text == "hello"
    assert seen_auth == ["Bearer token-1", "Bearer token-2"]


def test_get_client_wustl_discovery_failure_only_warns(monkeypatch, capsys):
    monkeypatch.setattr(EntraTokenProvider, "get_token", lambda self: "token-x")

    def failing_discovery(*args, **kwargs):
        raise RuntimeError("discovery down")

    monkeypatch.setattr(entra_gateway_client, "fetch_discovery", failing_discovery)

    client = profiles_module.get_client("wustl", credentials=GATEWAY_CREDENTIALS, model="anything")

    assert client.model == "anything"
    assert "discovery down" in capsys.readouterr().out


def test_get_client_anthropic_key_from_conf_takes_precedence_over_env(monkeypatch):
    captured = {}

    class RecordingAnthropicClient:
        def __init__(self, model, api_key=None, **default_params):
            captured["model"] = model
            captured["api_key"] = api_key
            self._client = None

    import llm.providers.anthropic_client as anthropic_client_module
    monkeypatch.setattr(anthropic_client_module, "AnthropicClient", RecordingAnthropicClient)
    monkeypatch.setattr(anthropic_client_module, "check_model_available", lambda *args: (None, None))
    monkeypatch.setenv("ANTHROPIC_API_KEY", "from-env")

    profiles_module.get_client("claude", credentials={"anthropic_api_key": "from-conf"})
    assert captured["api_key"] == "from-conf"

    profiles_module.get_client("claude", credentials={}, model="claude-other")
    assert captured["api_key"] == "from-env"
    assert captured["model"] == "claude-other"


def test_discovery_summary_lists_only_accessible_chat_models_by_default():
    from llm.__main__ import _format_discovery_summary

    discovery = {
        "models": [
            {"id": "gpt-5.4-mini", "api_format": "openai", "context_window": 400000,
             "capabilities": ["chat"], "cost": {"input_per_1k_tokens": 0.000825, "output_per_1k_tokens": 0.00495},
             "tags": ["recommended"], "app_has_access": True},
            {"id": "claude-opus-5-5-direct", "api_format": "anthropic", "context_window": 1000000,
             "capabilities": ["chat"], "cost": {}, "tags": [], "app_has_access": False},
            {"id": "text-embedding-3-small", "api_format": "openai", "context_window": 15000000,
             "capabilities": ["embeddings"], "cost": {}, "tags": [], "app_has_access": True},
        ],
    }

    summary = _format_discovery_summary(discovery)
    assert "gpt-5.4-mini" in summary
    assert "0.82" in summary and "4.95" in summary  # per-1k costs shown per 1M tokens
    assert "claude-opus-5-5-direct" not in summary
    assert "text-embedding-3-small" not in summary
    assert "2 hidden" in summary

    summary_all = _format_discovery_summary(discovery, show_all=True)
    assert "claude-opus-5-5-direct" in summary_all
    assert "text-embedding-3-small" in summary_all


def test_gateway_client_switches_to_max_completion_tokens_when_required():
    import pydantic

    class Answer(pydantic.BaseModel):
        answer: str

    sent_bodies = []

    def handler(request: httpx.Request) -> httpx.Response:
        body = json.loads(request.content)
        sent_bodies.append(body)
        if "max_tokens" in body:
            return httpx.Response(400, json={"error": {
                "message": "Unsupported parameter: 'max_tokens' is not supported with this model. Use 'max_completion_tokens' instead.",
                "type": "invalid_request_error", "param": "max_tokens", "code": "unsupported_parameter",
            }})
        return httpx.Response(200, json=_chat_completion('{"answer": "ok"}'))

    client = EntraGatewayClient(
        model="gpt-5.4-mini", base_url="https://gateway.example/models/v2",
        token_provider=_token_provider(FakeTokenSession()), api_key="test-gateway-api-key",
        http_client=httpx.Client(transport=httpx.MockTransport(handler)), max_tokens=16384,
    )

    first = client.complete(system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer)
    assert first.parsed.answer == "ok"
    assert "max_tokens" in sent_bodies[0]
    assert sent_bodies[1]["max_completion_tokens"] == 16384
    assert "max_tokens" not in sent_bodies[1]

    #remembered: the next call goes straight to max_completion_tokens, no failed attempt first
    client.complete(system="s", messages=[{"role": "user", "content": "again"}], output_schema=Answer)
    assert len(sent_bodies) == 3
    assert "max_completion_tokens" in sent_bodies[2]


def _anthropic_error(error_class, status_code: int):
    request = httpx.Request("GET", "https://api.anthropic.com/v1/models/x")
    return error_class("error", response=httpx.Response(status_code, request=request), body=None)


class FakeAnthropicModels:
    def __init__(self, retrieve_error=None, available=()):
        self.retrieve_error = retrieve_error
        self.available = available

    def retrieve(self, model):
        if self.retrieve_error:
            raise self.retrieve_error

    def list(self, limit=None):
        return [type("Model", (), {"id": model_id})() for model_id in self.available]


class FakeAnthropicSDK:
    def __init__(self, **kwargs):
        self.models = FakeAnthropicModels(**kwargs)


def test_anthropic_check_model_available_ok():
    from llm.providers.anthropic_client import check_model_available

    assert check_model_available(FakeAnthropicSDK(), "claude-opus-5", "claude") == (None, None)


def test_anthropic_check_model_unknown_model_explains_and_suggests():
    import anthropic
    from llm.providers.anthropic_client import check_model_available

    sdk = FakeAnthropicSDK(
        retrieve_error=_anthropic_error(anthropic.NotFoundError, 404),
        available=["claude-opus-5", "claude-sonnet-5"],
    )

    error, warning = check_model_available(sdk, "gpt-5.4-mini", "claude")
    assert "not available from Anthropic" in error
    assert "Did you mean" not in error
    assert warning is None

    assert "Models available to your API key: claude-opus-5, claude-sonnet-5" in error

    error, _ = check_model_available(sdk, "claude-opus5", "claude")
    assert "Did you mean: claude-opus-5" in error


def test_format_anthropic_models_table():
    from datetime import datetime, timezone
    from llm.__main__ import _format_anthropic_models

    models = [
        type("Model", (), {"id": "claude-opus-5", "display_name": "Claude Opus 5",
                           "created_at": datetime(2026, 5, 1, tzinfo=timezone.utc)})(),
        type("Model", (), {"id": "claude-sonnet-5", "display_name": "Claude Sonnet 5", "created_at": None})(),
    ]

    table = _format_anthropic_models(models)
    assert "claude-opus-5" in table and "Claude Opus 5" in table and "2026-05-01" in table
    assert "claude-sonnet-5" in table
    assert "2 model(s) available" in table


def test_anthropic_check_model_bad_api_key():
    import anthropic
    from llm.providers.anthropic_client import check_model_available

    sdk = FakeAnthropicSDK(retrieve_error=_anthropic_error(anthropic.AuthenticationError, 401))
    error, _ = check_model_available(sdk, "claude-opus-5", "claude")
    assert "rejected the API key" in error


def test_anthropic_check_model_other_failure_only_warns():
    from llm.providers.anthropic_client import check_model_available

    sdk = FakeAnthropicSDK(retrieve_error=RuntimeError("network down"))
    error, warning = check_model_available(sdk, "claude-opus-5", "claude")
    assert error is None
    assert "network down" in warning


def test_rewrite_sdk_messages_path_only_touches_messages():
    from llm.providers.entra_gateway_client import _rewrite_sdk_messages_path

    http = _anthropic_http_module()

    messages = http.Request("POST", "https://gateway.example/models/v2/v1/messages")
    _rewrite_sdk_messages_path(messages)
    assert str(messages.url) == "https://gateway.example/models/v2/messages"

    other = http.Request("GET", "https://gateway.example/models/v2/v1/models")
    _rewrite_sdk_messages_path(other)
    assert str(other.url) == "https://gateway.example/models/v2/v1/models"


class FakeGeminiResponse:
    def __init__(self, status_code: int, body: dict):
        self.status_code = status_code
        self._body = body
        self.text = json.dumps(body)

    def json(self):
        return self._body


class FakeGeminiSession:
    """Stands in for `requests` for EntraGatewayGeminiClient: returns queued responses in order
    and records each POST."""

    def __init__(self, responses):
        self.responses = list(responses)
        self.posts = []

    def post(self, url, headers=None, json=None, timeout=None):
        self.posts.append({"url": url, "headers": headers, "json": json})
        return self.responses.pop(0)


def _gemini_ok(text: str, finish_reason: str = "STOP") -> FakeGeminiResponse:
    return FakeGeminiResponse(200, {
        "candidates": [{"content": {"role": "model", "parts": [{"text": text}]}, "finishReason": finish_reason}],
        "usageMetadata": {"promptTokenCount": 10, "candidatesTokenCount": 4, "thoughtsTokenCount": 6},
    })


def _gemini_client(session, token_session=None):
    from llm.providers.entra_gateway_client import EntraGatewayGeminiClient

    return EntraGatewayGeminiClient(
        model="gemini-2.5-flash", base_url="https://gateway.example/models/v2",
        token_provider=_token_provider(token_session or FakeTokenSession()), api_key="test-gateway-api-key",
        session=session, max_tokens=16384,
    )


def test_gemini_client_request_shape_and_structured_output():
    import pydantic

    class Answer(pydantic.BaseModel):
        answer: str

    session = FakeGeminiSession([_gemini_ok('{"answer": "ok"}')])
    result = _gemini_client(session).complete(
        system="be careful", messages=[{"role": "user", "content": "hi"}], output_schema=Answer,
    )

    assert result.parsed.answer == "ok"
    assert result.provider == "entra_gateway"
    assert result.input_tokens == 10
    assert result.output_tokens == 10  # candidates + thoughts

    post = session.posts[0]
    assert post["url"] == "https://gateway.example/models/v2/generateContent"
    assert post["headers"]["Authorization"] == "Bearer token-1"
    assert post["headers"]["X-Api-Key"] == "test-gateway-api-key"
    body = post["json"]
    assert body["model"] == "gemini-2.5-flash"
    assert body["contents"] == [{"role": "user", "parts": [{"text": "hi"}]}]
    assert body["generationConfig"] == {"maxOutputTokens": 16384, "responseMimeType": "application/json"}
    system_text = body["systemInstruction"]["parts"][0]["text"]
    assert system_text.startswith("be careful")
    assert '"answer"' in system_text  # the JSON schema is spelled out


def test_gemini_client_schema_retry_sends_prior_answer_as_model_turn():
    import pydantic

    class Answer(pydantic.BaseModel):
        answer: str

    session = FakeGeminiSession([_gemini_ok('{"wrong": 1}'), _gemini_ok('{"answer": "fixed"}')])
    result = _gemini_client(session).complete(
        system="s", messages=[{"role": "user", "content": "hi"}], output_schema=Answer,
    )

    assert result.parsed.answer == "fixed"
    assert result.retries == 1
    retry_contents = session.posts[1]["json"]["contents"]
    assert [c["role"] for c in retry_contents] == ["user", "model", "user"]
    assert retry_contents[1]["parts"][0]["text"] == '{"wrong": 1}'


def test_gemini_client_refreshes_token_once_on_401():
    token_session = FakeTokenSession()
    session = FakeGeminiSession([FakeGeminiResponse(401, {"error": "expired"}), _gemini_ok("hello")])
    result = _gemini_client(session, token_session).complete(system="s", messages=[{"role": "user", "content": "hi"}])

    assert result.raw_text == "hello"
    assert [p["headers"]["Authorization"] for p in session.posts] == ["Bearer token-1", "Bearer token-2"]


def test_gemini_client_retries_transient_errors(monkeypatch):
    monkeypatch.setattr(entra_gateway_client.time, "sleep", lambda seconds: None)
    session = FakeGeminiSession([FakeGeminiResponse(503, {"error": "busy"}), _gemini_ok("hello")])
    result = _gemini_client(session).complete(system="s", messages=[{"role": "user", "content": "hi"}])

    assert result.raw_text == "hello"
    assert len(session.posts) == 2


def test_gemini_client_http_error_raises_with_status_and_body():
    from llm.providers.entra_gateway_client import GatewayHTTPError

    session = FakeGeminiSession([FakeGeminiResponse(400, {"error": {"message": "bad request"}})])
    with pytest.raises(GatewayHTTPError) as exc_info:
        _gemini_client(session).complete(system="s", messages=[{"role": "user", "content": "hi"}])

    assert "Error code: 400" in str(exc_info.value)
    assert "bad request" in str(exc_info.value)


def test_gemini_client_truncation_and_blocks_are_refusals():
    from llm.providers.base import LLMRefusalError

    session = FakeGeminiSession([_gemini_ok('{"answer": "cut', finish_reason="MAX_TOKENS")])
    with pytest.raises(LLMRefusalError) as exc_info:
        _gemini_client(session).complete(system="s", messages=[{"role": "user", "content": "hi"}])
    assert exc_info.value.stop_reason == "MAX_TOKENS"

    session = FakeGeminiSession([FakeGeminiResponse(200, {"promptFeedback": {"blockReason": "SAFETY"}})])
    with pytest.raises(LLMRefusalError) as exc_info:
        _gemini_client(session).complete(system="s", messages=[{"role": "user", "content": "hi"}])
    assert exc_info.value.stop_reason == "SAFETY"
