import json
import time
from typing import Optional

import anthropic
import openai
import requests

from .anthropic_client import AnthropicClient
from .base import LLMClient, LLMRefusalError, LLMResult, validate_with_retry
from .openai_compat_client import OpenAICompatibleClient

#refresh the bearer token this many seconds before Entra says it expires, so a request is never
#sent with a token that expires mid-flight
TOKEN_REFRESH_MARGIN_SECONDS = 300

#names of the *.civic.conf settings an entra_gateway profile reads its credentials from
GATEWAY_CREDENTIAL_KEYS = ("llm_gateway_client_id", "llm_gateway_client_secret", "llm_gateway_api_key")


class EntraTokenProvider:
    """Fetches and caches a Microsoft Entra client-credentials bearer token, requesting a new one
    shortly before the cached one expires (Entra tokens last ~1 hour, so a long review run needs
    several). Never prints or logs the client secret or the token."""

    def __init__(self, token_url: str, client_id: str, client_secret: str, scope: str, session=None):
        self.token_url = token_url
        self.client_id = client_id
        self.client_secret = client_secret
        self.scope = scope
        self._session = session or requests
        self._token = None
        self._expires_at = 0.0

    def get_token(self) -> str:
        if self._token is None or time.monotonic() >= self._expires_at - TOKEN_REFRESH_MARGIN_SECONDS:
            self._fetch_token()
        return self._token

    def invalidate(self) -> None:
        """Forget the cached token so the next get_token() call fetches a fresh one."""
        self._token = None
        self._expires_at = 0.0

    def _fetch_token(self) -> None:
        resp = self._session.post(
            self.token_url,
            data={
                "grant_type": "client_credentials",
                "client_id": self.client_id,
                "client_secret": self.client_secret,
                "scope": self.scope,
            },
            timeout=(20, 60),
        )
        if resp.status_code != 200:
            #Entra's error body (e.g. AADSTS7000215 'invalid client secret') explains the problem
            #without echoing the secret back, so it's safe and useful to surface
            try:
                body = resp.json()
                detail = body.get("error_description") or body.get("error") or resp.text
            except ValueError:
                detail = resp.text
            raise RuntimeError(f"Failed to obtain LLM gateway bearer token (HTTP {resp.status_code}): {detail}")

        body = resp.json()
        self._token = body["access_token"]
        self._expires_at = time.monotonic() + float(body.get("expires_in", 3600))


class EntraGatewayClient(OpenAICompatibleClient):
    """LLMClient for an OpenAI-compatible LLM API gateway that requires two credentials on every
    request: a short-lived Microsoft Entra bearer token (proves who is calling, sent as
    'Authorization: Bearer') and a long-lived API key (ties the call to the account's model access
    and quota, sent as 'X-Api-Key').

    Everything else (structured output, JSON-mode fallback, finish_reason checks) is inherited
    from OpenAICompatibleClient. Before each complete() call the bearer token is refreshed if it is
    near expiry, and a 401 triggers one retry with a freshly fetched token.
    """

    provider_name = "entra_gateway"

    def __init__(self, model: str, base_url: str, token_provider: EntraTokenProvider, api_key: str,
                 http_client=None, **default_params):
        self.model = model
        self.base_url = base_url
        self.default_params = default_params
        self._token_provider = token_provider
        self._client = openai.OpenAI(
            base_url=base_url,
            api_key=token_provider.get_token(),
            default_headers={"X-Api-Key": api_key},
            http_client=http_client,
        )

    def complete(self, *, system: str, messages: list, output_schema=None, **opts) -> LLMResult:
        self._client.api_key = self._token_provider.get_token()
        try:
            return super().complete(system=system, messages=messages, output_schema=output_schema, **opts)
        except openai.AuthenticationError:
            #the token may have been revoked or expired early - retry once with a fresh one
            self._token_provider.invalidate()
            self._client.api_key = self._token_provider.get_token()
            return super().complete(system=system, messages=messages, output_schema=output_schema, **opts)


_SDK_MESSAGES_PATH_SUFFIX = "/v1/messages"


def _rewrite_sdk_messages_path(request) -> None:
    """httpx request hook: <base_url>/v1/messages (what the anthropic SDK sends) -> the gateway's
    <base_url>/messages. Leaves every other path alone."""
    path = request.url.path
    if path.endswith(_SDK_MESSAGES_PATH_SUFFIX):
        request.url = request.url.copy_with(path=path[: -len(_SDK_MESSAGES_PATH_SUFFIX)] + "/messages")


class EntraGatewayAnthropicClient(AnthropicClient):
    """LLMClient for the gateway's Anthropic-format models (api_format 'anthropic' in discovery),
    which the gateway only serves via the Anthropic Messages API, not chat completions.

    The Entra token is passed as the SDK's auth_token (sent as 'Authorization: Bearer') and the
    gateway API key as an explicit 'X-Api-Key' default header - not as api_key, since anthropic
    SDK 2.x refuses a client given both api_key and auth_token. The gateway serves Messages at
    <base_url>/messages (alongside <base_url>/chat/completions for OpenAI-format models), but the
    SDK always posts to <base_url>/v1/messages, so a request hook rewrites that path just before
    each request is sent. Everything else (structured output, forced-tool fallback, stop_reason checks) is
    inherited from AnthropicClient. Token refresh and the one-retry-on-401 match
    EntraGatewayClient.
    """

    provider_name = "entra_gateway"

    def __init__(self, model: str, base_url: str, token_provider: EntraTokenProvider, api_key: str,
                 http_client=None, **default_params):
        self.model = model
        self.base_url = base_url
        self.default_params = default_params
        self._token_provider = token_provider
        http_client = http_client or anthropic.DefaultHttpxClient()
        hooks = http_client.event_hooks
        http_client.event_hooks = {**hooks, "request": list(hooks.get("request", [])) + [_rewrite_sdk_messages_path]}
        self._client = anthropic.Anthropic(
            base_url=base_url,
            auth_token=token_provider.get_token(),
            default_headers={"X-Api-Key": api_key},
            http_client=http_client,
        )

    def complete(self, *, system: str, messages: list, output_schema=None, **opts) -> LLMResult:
        self._client.auth_token = self._token_provider.get_token()
        try:
            return super().complete(system=system, messages=messages, output_schema=output_schema, **opts)
        except anthropic.AuthenticationError:
            self._token_provider.invalidate()
            self._client.auth_token = self._token_provider.get_token()
            return super().complete(system=system, messages=messages, output_schema=output_schema, **opts)


#HTTP statuses worth retrying for the Gemini client (the openai/anthropic SDKs retry these themselves)
_RETRYABLE_STATUS_CODES = {429, 500, 502, 503, 504}
_MAX_RETRIES = 2


class GatewayHTTPError(Exception):
    """A non-success response from the gateway, with its body (which never echoes credentials)."""

    def __init__(self, status_code: int, body: str):
        self.status_code = status_code
        self.body = body
        super().__init__(f"Error code: {status_code} - {body[:500]}")


class EntraGatewayGeminiClient(LLMClient):
    """LLMClient for the gateway's Gemini-format models (api_format 'gemini' in discovery), served
    at <base_url>/generateContent with a Gemini generateContent request body - except that the
    model id goes in the body, not the URL, so Google's SDK doesn't fit and plain HTTP is used.

    Structured output uses Gemini's JSON mode (responseMimeType 'application/json') with the
    target JSON schema spelled out in the system instruction, then validate_with_retry() - the
    same approach as OpenAICompatibleClient's JSON-mode fallback. Like the other gateway clients,
    the bearer token is refreshed near expiry and a 401 triggers one retry with a fresh token;
    429/5xx responses are retried with backoff, as the openai/anthropic SDKs do for theirs.
    """

    provider_name = "entra_gateway"

    def __init__(self, model: str, base_url: str, token_provider: EntraTokenProvider, api_key: str,
                 session=None, **default_params):
        self.model = model
        self.base_url = base_url
        self.default_params = default_params
        self._token_provider = token_provider
        self._api_key = api_key
        self._session = session or requests
        self._url = base_url.rstrip("/") + "/generateContent"
        token_provider.get_token()  # fail fast on bad credentials, like the other gateway clients

    def complete(self, *, system: str, messages: list, output_schema=None, **opts) -> LLMResult:
        max_tokens = opts.pop("max_tokens", self.default_params.get("max_tokens", 4096))

        system_text = system
        if output_schema is not None:
            system_text += (
                "\n\nRespond with ONLY a single JSON object matching this JSON Schema, and no other text:\n"
                + json.dumps(output_schema.model_json_schema())
            )

        start = time.monotonic()
        body = self._post(self._request_body(system_text, messages, max_tokens, json_mode=output_schema is not None))
        raw_text = self._extract_text(body)

        if output_schema is None:
            return self._to_result([body], parsed=None, raw_text=raw_text,
                                   latency_ms=(time.monotonic() - start) * 1000, retries=0)

        state = {"raw_text": raw_text, "bodies": [body]}

        def retry_fn(error_message: str) -> str:
            corrected_messages = list(messages) + [
                {"role": "assistant", "content": state["raw_text"]},
                {
                    "role": "user",
                    "content": (
                        f"That response failed schema validation with this error:\n{error_message}\n"
                        f"Reply again with ONLY corrected JSON matching the schema, nothing else."
                    ),
                },
            ]
            new_body = self._post(self._request_body(system_text, corrected_messages, max_tokens, json_mode=True))
            state["raw_text"] = self._extract_text(new_body)
            state["bodies"].append(new_body)
            return state["raw_text"]

        parsed, retries = validate_with_retry(output_schema, raw_text, retry_fn)

        return self._to_result(state["bodies"], parsed=parsed, raw_text=state["raw_text"],
                               latency_ms=(time.monotonic() - start) * 1000, retries=retries)

    def _request_body(self, system_text: str, messages: list, max_tokens: int, json_mode: bool) -> dict:
        generation_config = {"maxOutputTokens": max_tokens}
        if json_mode:
            generation_config["responseMimeType"] = "application/json"
        return {
            "model": self.model,
            "systemInstruction": {"parts": [{"text": system_text}]},
            "contents": [
                {"role": "model" if m["role"] == "assistant" else "user", "parts": [{"text": m["content"]}]}
                for m in messages
            ],
            "generationConfig": generation_config,
        }

    def _post(self, request_body: dict) -> dict:
        refreshed_after_401 = False
        attempt = 0
        while True:
            resp = self._session.post(
                self._url,
                headers={
                    "Authorization": f"Bearer {self._token_provider.get_token()}",
                    "X-Api-Key": self._api_key,
                    "Content-Type": "application/json",
                },
                json=request_body,
                timeout=(20, 600),
            )
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 401 and not refreshed_after_401:
                refreshed_after_401 = True
                self._token_provider.invalidate()
                continue
            if resp.status_code in _RETRYABLE_STATUS_CODES and attempt < _MAX_RETRIES:
                attempt += 1
                time.sleep(2 ** attempt)
                continue
            raise GatewayHTTPError(resp.status_code, resp.text)

    @staticmethod
    def _extract_text(body: dict) -> str:
        """Concatenate the first candidate's text parts (skipping any thought-summary parts).
        Anything but a normal STOP - a blocked prompt, a safety stop, truncation at
        maxOutputTokens - is a hard failure, not something validate_with_retry() should fix."""
        candidates = body.get("candidates") or []
        if not candidates:
            block_reason = (body.get("promptFeedback") or {}).get("blockReason", "NO_CANDIDATES")
            raise LLMRefusalError(stop_reason=block_reason, raw_text="")

        candidate = candidates[0]
        parts = (candidate.get("content") or {}).get("parts") or []
        text = "".join(part.get("text", "") for part in parts if not part.get("thought"))

        finish_reason = candidate.get("finishReason", "STOP")
        if finish_reason != "STOP":
            raise LLMRefusalError(stop_reason=finish_reason, raw_text=text)
        return text

    def _to_result(self, bodies: list, *, parsed, raw_text, latency_ms, retries) -> LLMResult:
        """Build the LLMResult from every response body this call got (more than one after a
        validation retry), summing their usage. promptTokenCount already includes any cached
        tokens."""
        usages = [body.get("usageMetadata") or {} for body in bodies]
        cached = [u["cachedContentTokenCount"] for u in usages if u.get("cachedContentTokenCount") is not None]
        return LLMResult(
            parsed=parsed,
            raw_text=raw_text,
            provider=self.provider_name,
            model=self.model,
            input_tokens=sum(u.get("promptTokenCount", 0) for u in usages),
            #thinking tokens are billed as output, so they're counted here too
            output_tokens=sum(u.get("candidatesTokenCount", 0) + u.get("thoughtsTokenCount", 0) for u in usages),
            cache_read_input_tokens=sum(cached) if cached else None,
            cache_creation_input_tokens=None,
            latency_ms=latency_ms,
            stop_reason=bodies[-1]["candidates"][0].get("finishReason", "STOP"),
            retries=retries,
        )


def fetch_discovery(discovery_url: str, token_provider: EntraTokenProvider, api_key: str, session=None) -> dict:
    """Call the gateway's discovery endpoint, which describes the models this account can
    currently use (capabilities, costs, etc.), and return the parsed JSON."""
    session = session or requests
    resp = session.get(
        discovery_url,
        headers={"Authorization": f"Bearer {token_provider.get_token()}", "X-Api-Key": api_key},
        timeout=(20, 60),
    )
    resp.raise_for_status()
    return resp.json()


#discovery api_format values this package has a client for (see profiles.get_client)
SUPPORTED_GATEWAY_API_FORMATS = ("openai", "anthropic", "gemini")


def check_model_against_discovery(model: str, discovery: dict) -> tuple:
    """
    Check a requested model id against a gateway discovery response, since the gateway itself
    answers both a misspelled model id and a real access problem with the same 403 'does not
    have access' message. Returns (error, info): error is a message explaining why the model
    can't be used (unknown id - with close-match suggestions - no access, not a chat model, or
    an api_format this package can't speak yet), or None; info is the model's discovery entry.
    """
    import difflib

    models = {m.get("id"): m for m in discovery.get("models", []) if m.get("id")}

    if model not in models:
        accessible = [model_id for model_id, m in models.items() if m.get("app_has_access")]
        suggestions = difflib.get_close_matches(model, accessible, n=3, cutoff=0.6)
        hint = f" Did you mean: {', '.join(suggestions)}?" if suggestions else ""
        return (f"model '{model}' is not offered by the gateway.{hint} "
                f"List available models with: python -m llm discover --profile <profile> --civic-config <file>"), None

    info = models[model]
    if not info.get("app_has_access"):
        return (f"your gateway account does not have access to model '{model}' yet - "
                f"request it in the gateway dashboard"), None

    if "chat" not in info.get("capabilities", []):
        return f"model '{model}' is not a chat model (capabilities: {', '.join(info.get('capabilities', []))})", info

    api_format = info.get("api_format", "openai")
    if api_format not in SUPPORTED_GATEWAY_API_FORMATS:
        return (f"model '{model}' uses the gateway's '{api_format}' API format, which isn't supported "
                f"yet (supported: {', '.join(SUPPORTED_GATEWAY_API_FORMATS)})"), info

    return None, info
