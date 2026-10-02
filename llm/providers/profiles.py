import dataclasses
import os
import tomllib
from pathlib import Path
from typing import Optional

from .base import LLMClient

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


def load_profiles(path=None) -> dict:
    """
    Parse the provider-profiles TOML file (stdlib tomllib - the file holds no secrets, API keys
    are only ever read from the environment) into {profile_name: Profile}.

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
        )

    return profiles


def get_profile(profile_name: str, profiles_path=None) -> Profile:
    """Load profiles and return the one named profile_name, or raise KeyError."""
    profiles = load_profiles(profiles_path)
    if profile_name not in profiles:
        raise KeyError(f"Unknown profile '{profile_name}'. Known profiles: {sorted(profiles)}")
    return profiles[profile_name]


def get_client(profile_name: str, profiles_path=None) -> LLMClient:
    """
    Look up profile_name, read its API key from the environment variable it names (never from
    the profiles file itself), and construct + return the matching LLMClient.
    """
    profile = get_profile(profile_name, profiles_path)
    api_key = os.environ.get(profile.api_key_env) if profile.api_key_env else None

    if profile.provider == "anthropic":
        from .anthropic_client import AnthropicClient
        return AnthropicClient(model=profile.model, api_key=api_key, **profile.default_params)

    if profile.provider == "openai_compatible":
        from .openai_compat_client import OpenAICompatibleClient
        return OpenAICompatibleClient(
            model=profile.model,
            base_url=profile.base_url,
            api_key=api_key,
            **profile.default_params,
        )

    raise ValueError(f"Unknown provider '{profile.provider}' for profile '{profile_name}'")
