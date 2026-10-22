import dataclasses
import os
import tomllib
from pathlib import Path
from typing import Optional

from .base import LLMClient

YELLOW = "\033[33m"
RESET = "\033[0m"

_PROVIDERS_DIR = Path(__file__).resolve().parent
DEFAULT_PROFILES_PATH = _PROVIDERS_DIR.parent / "profiles.toml"


@dataclasses.dataclass
class Profile:
    name: str
    provider: str
    model: str
    base_url: Optional[str]
    api_key_env: Optional[str]
    supports: list
    default_params: dict
    api_key_config: Optional[str] = None  # name of the optional *.civic.conf setting holding the API key
    token_url: Optional[str] = None  # entra_gateway only
    scope: Optional[str] = None  # entra_gateway only
    discovery_url: Optional[str] = None  # entra_gateway only


def load_profiles(path=None) -> dict:
    """
    Parse the provider-profiles TOML file (stdlib tomllib - the file holds no secrets, API keys
    are only ever read from the user's *.civic.conf file or the environment) into
    {profile_name: Profile}.

    path defaults to DEFAULT_PROFILES_PATH, overridable by the LLM_PROFILES_PATH env var,
    or by passing path explicitly.
    """
    if path is None:
        path = os.environ.get("LLM_PROFILES_PATH", DEFAULT_PROFILES_PATH)
    path = Path(path)

    if not path.exists():
        raise FileNotFoundError(f"Profiles file not found: {path}")

    with path.open("rb") as f:
        raw = tomllib.load(f)

    profiles = {}
    for name, values in raw.items():
        profiles[name] = Profile(
            name=name,
            provider=values["provider"],
            model=values["model"],
            base_url=values.get("base_url"),
            api_key_env=values.get("api_key_env") or None,
            supports=values.get("supports", []),
            default_params=values.get("default_params", {}),
            api_key_config=values.get("api_key_config") or None,
            token_url=values.get("token_url"),
            scope=values.get("scope"),
            discovery_url=values.get("discovery_url"),
        )

    return profiles


def get_profile(profile_name: str, profiles_path=None) -> Profile:
    """Load profiles and return the one named profile_name, or raise KeyError."""
    profiles = load_profiles(profiles_path)
    if profile_name not in profiles:
        raise KeyError(f"Unknown profile '{profile_name}'. Known profiles: {sorted(profiles)}")
    return profiles[profile_name]


def _require_gateway_credentials(profile: Profile, credentials: dict) -> dict:
    """Return the entra_gateway credentials from `credentials`, or raise naming what's missing."""
    from .entra_gateway_client import GATEWAY_CREDENTIAL_KEYS

    missing = [key for key in GATEWAY_CREDENTIAL_KEYS if not credentials.get(key)]
    if missing:
        raise ValueError(
            f"Profile '{profile.name}' needs {', '.join(missing)} set in your *.civic.conf file "
            f"(see example.civic.conf.template)"
        )
    return {key: credentials[key] for key in GATEWAY_CREDENTIAL_KEYS}


def build_gateway_token_provider(profile: Profile, credentials: dict):
    """Return (EntraTokenProvider, gateway_api_key) for an entra_gateway profile."""
    from .entra_gateway_client import EntraTokenProvider

    creds = _require_gateway_credentials(profile, credentials)
    token_provider = EntraTokenProvider(
        token_url=profile.token_url,
        client_id=creds["llm_gateway_client_id"],
        client_secret=creds["llm_gateway_client_secret"],
        scope=profile.scope,
    )
    return token_provider, creds["llm_gateway_api_key"]


def _gateway_model_api_format(profile: Profile, model: str, token_provider, gateway_api_key: str) -> str:
    """
    Look `model` up in the gateway's discovery list and return its api_format ('openai' or
    'anthropic'), which decides which request format - and so which client - to use. Fails fast,
    with a useful message, if the model is misspelled or not available to this account (the
    gateway's own error for both is an ambiguous 403). If discovery is unavailable, warns and
    guesses the format from the model id rather than blocking on it.
    """
    from .entra_gateway_client import check_model_against_discovery, fetch_discovery

    guessed_format = "anthropic" if model.startswith("claude") else "gemini" if model.startswith("gemini") else "openai"
    if not profile.discovery_url:
        return guessed_format

    try:
        discovery = fetch_discovery(profile.discovery_url, token_provider, gateway_api_key)
    except Exception as exc:
        print(f"{YELLOW}WARNING: could not check model '{model}' against the gateway's model list ({exc}); "
              f"continuing, assuming the '{guessed_format}' API format.{RESET}")
        return guessed_format

    error, info = check_model_against_discovery(model, discovery)
    if error:
        raise ValueError(error)
    return info.get("api_format", "openai")


def resolve_api_key(profile: Profile, credentials: Optional[dict] = None) -> Optional[str]:
    """A profile's API key: its api_key_config setting from the *.civic.conf credentials if set
    there, else its api_key_env environment variable, else None."""
    credentials = credentials or {}
    api_key = credentials.get(profile.api_key_config) if profile.api_key_config else None
    if not api_key and profile.api_key_env:
        api_key = os.environ.get(profile.api_key_env)
    return api_key


def get_client(profile_name: str, profiles_path=None, credentials: Optional[dict] = None,
               model: Optional[str] = None) -> LLMClient:
    """
    Look up profile_name and construct + return the matching LLMClient. Never reads secrets from
    the profiles file itself.

    credentials: optional {setting_name: value} of the LLM settings from the user's *.civic.conf
    file (see utils.civic_graphql_utils.load_civic_config). A profile's api_key_config setting,
    when present there, takes precedence over its api_key_env environment variable. entra_gateway
    profiles read their client id/secret and gateway API key only from here.

    model: optional override of the profile's default model, e.g. to try a different model
    through the same provider.
    """
    profile = get_profile(profile_name, profiles_path)
    credentials = credentials or {}
    model = model or profile.model

    api_key = resolve_api_key(profile, credentials)

    if profile.provider == "anthropic":
        from .anthropic_client import AnthropicClient, check_model_available
        client = AnthropicClient(model=model, api_key=api_key, **profile.default_params)
        error, warning = check_model_available(client._client, model, profile.name)
        if error:
            raise ValueError(error)
        if warning:
            print(f"{YELLOW}WARNING: {warning}{RESET}")
        return client

    if profile.provider == "openai_compatible":
        from .openai_compat_client import OpenAICompatibleClient
        return OpenAICompatibleClient(
            model=model,
            base_url=profile.base_url,
            api_key=api_key,
            **profile.default_params,
        )

    if profile.provider == "entra_gateway":
        from .entra_gateway_client import EntraGatewayAnthropicClient, EntraGatewayClient, EntraGatewayGeminiClient
        token_provider, gateway_api_key = build_gateway_token_provider(profile, credentials)
        api_format = _gateway_model_api_format(profile, model, token_provider, gateway_api_key)
        client_class = {
            "anthropic": EntraGatewayAnthropicClient,
            "gemini": EntraGatewayGeminiClient,
        }.get(api_format, EntraGatewayClient)
        return client_class(
            model=model,
            base_url=profile.base_url,
            token_provider=token_provider,
            api_key=gateway_api_key,
            **profile.default_params,
        )

    raise ValueError(f"Unknown provider '{profile.provider}' for profile '{profile_name}'")
