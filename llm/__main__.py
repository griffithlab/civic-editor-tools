"""
Smoke-test CLI for llm. Library usage is `from llm import run_task`; this CLI exists only to
make plumbing (task discovery, prompt assembly, a provider call, provenance logging) easy to
exercise by hand without writing a script.

    python -m llm run review_evidence_item --input <file.json> --profile local [--fake]
    python -m llm run review_evidence_item --input <file.json> --profile wustl --civic-config <file> [--model <name>]
    python -m llm discover --profile wustl --civic-config <file>
    python -m llm discover --profile claude [--civic-config <file>]
    python -m llm list-tasks
"""

import argparse
import json
import sys
from pathlib import Path

from . import registry
from .providers import profiles as profiles_module
from .providers.base import LLMOutputValidationError, LLMRefusalError, LLMResult
from .runner import run_task
from .testing import FakeClient


def parse_args():
    parser = argparse.ArgumentParser(
        prog="python -m llm",
        description="Smoke-test CLI for llm.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run a task against one input and print the result")
    run_parser.add_argument("task_name", type=str, help="Task name, e.g. review_evidence_item")
    run_parser.add_argument(
        "--input", dest="input_path", type=str, required=True,
        help="Path to a JSON file matching the task's `inputs` shape",
    )
    run_parser.add_argument(
        "--profile", dest="profile", type=str, default="local",
        help="Provider profile name from llm/profiles.toml (default: local)",
    )
    run_parser.add_argument(
        "--fake", dest="fake", action="store_true",
        help="Use a canned FakeClient instead of a real provider - zero network, zero config",
    )
    run_parser.add_argument(
        "--fake-output", dest="fake_output_path", type=str, default=None,
        help="With --fake, a JSON file for the canned output (default: a generic placeholder "
             "generated from the task's schema)",
    )
    run_parser.add_argument(
        "--model", dest="model", type=str, default=None,
        help="Override the profile's default model",
    )
    run_parser.add_argument(
        "--civic-config", dest="civic_config", type=str, default=None,
        help="Path to a *.civic.conf file to read LLM credentials from (required for entra_gateway "
             "profiles such as wustl; optional for claude, which otherwise uses ANTHROPIC_API_KEY)",
    )
    run_parser.add_argument("--content-root", dest="content_root", type=str, default=None)
    run_parser.add_argument("--profiles-file", dest="profiles_path", type=str, default=None)
    run_parser.add_argument("--log-file", dest="log_path", type=str, default=None)

    discover_parser = subparsers.add_parser(
        "discover", help="List the models a profile's credentials can use (anthropic or entra_gateway profiles)",
    )
    discover_parser.add_argument(
        "--profile", dest="profile", type=str, default="wustl",
        help="An anthropic (e.g. claude) or entra_gateway (e.g. wustl) profile name from llm/profiles.toml (default: wustl)",
    )
    discover_parser.add_argument(
        "--civic-config", dest="civic_config", type=str, default=None,
        help="Path to a *.civic.conf file holding the profile's credentials (required for entra_gateway "
             "profiles; optional for claude, which otherwise uses ANTHROPIC_API_KEY)",
    )
    discover_parser.add_argument(
        "--json", dest="raw_json", action="store_true",
        help="Print the provider's full JSON response instead of the summary table",
    )
    discover_parser.add_argument(
        "--all", dest="show_all", action="store_true",
        help="entra_gateway only: include models your account can't use yet and non-chat (e.g. embedding) models",
    )
    discover_parser.add_argument("--profiles-file", dest="profiles_path", type=str, default=None)

    list_parser = subparsers.add_parser("list-tasks", help="List discovered tasks")
    list_parser.add_argument("--content-root", dest="content_root", type=str, default=None)

    return parser.parse_args()


def _placeholder_from_json_schema(schema_dict: dict, defs: dict):
    """Best-effort generic example value for a JSON Schema fragment - just enough to satisfy
    pydantic validation for a --fake smoke test, not a realistic example."""
    if "$ref" in schema_dict:
        return _placeholder_from_json_schema(defs[schema_dict["$ref"].split("/")[-1]], defs)

    if "enum" in schema_dict:
        return schema_dict["enum"][0]

    if "anyOf" in schema_dict:
        for option in schema_dict["anyOf"]:
            if option.get("type") != "null":
                return _placeholder_from_json_schema(option, defs)
        return None

    json_type = schema_dict.get("type")
    if json_type == "string":
        return "example"
    if json_type == "integer":
        return 0
    if json_type == "number":
        return 0.0
    if json_type == "boolean":
        return True
    if json_type == "array":
        return []
    if json_type == "object":
        required = set(schema_dict.get("required", []))
        props = schema_dict.get("properties", {})
        return {key: _placeholder_from_json_schema(value, defs) for key, value in props.items() if key in required}

    return None


def _build_fake_response(schema: type, fake_output_path) -> LLMResult:
    if fake_output_path:
        output_dict = json.loads(Path(fake_output_path).read_text())
    else:
        full_schema = schema.model_json_schema()
        output_dict = _placeholder_from_json_schema(full_schema, full_schema.get("$defs", {}))

    parsed = schema.model_validate(output_dict)

    return LLMResult(
        parsed=parsed,
        raw_text=json.dumps(output_dict),
        provider="fake",
        model="fake-model",
        input_tokens=0,
        output_tokens=0,
        cache_read_input_tokens=None,
        cache_creation_input_tokens=None,
        latency_ms=0.0,
        stop_reason="end_turn",
        retries=0,
    )


def _load_llm_credentials(civic_config) -> dict:
    """Read just the LLM credentials from a *.civic.conf file (with all of load_civic_config()'s
    filename/permission checks). Imported lazily: utils/ lives at the repo root, alongside llm/."""
    if not civic_config:
        return {}
    from utils.civic_graphql_utils import load_civic_config
    return load_civic_config(civic_config)["llm_credentials"]


def _format_discovery_summary(discovery: dict, show_all: bool = False) -> str:
    """Render a gateway discovery response as a compact table: one row per model, costs in USD
    per million tokens. By default only chat models this account can use (app_has_access) are
    listed. Only the 'openai' api_format is currently supported by the entra_gateway client."""
    models = discovery.get("models", [])
    shown = [
        m for m in models
        if show_all or (m.get("app_has_access") and "chat" in m.get("capabilities", []))
    ]

    def per_million(value):
        return "-" if value is None else f"{value * 1000:.2f}"

    header = ("model", "format", "context", "in $/1M", "out $/1M", "reasoning", "access", "tags")
    rows = [header]
    for m in shown:
        cost = m.get("cost", {})
        rows.append((
            m.get("id", "?"),
            m.get("api_format", "?"),
            f"{m.get('context_window', 0) // 1000}k",
            per_million(cost.get("input_per_1k_tokens")),
            per_million(cost.get("output_per_1k_tokens")),
            "yes" if "reasoning" in m.get("capabilities", []) else "",
            "yes" if m.get("app_has_access") else "no",
            ",".join(m.get("tags", [])),
        ))

    widths = [max(len(str(row[i])) for row in rows) for i in range(len(header))]
    lines = ["  ".join(str(value).ljust(width) for value, width in zip(row, widths)).rstrip() for row in rows]
    lines.insert(1, "  ".join("-" * width for width in widths))

    hidden = len(models) - len(shown)
    lines.append("")
    lines.append(f"{len(shown)} model(s) shown" + (f", {hidden} hidden (use --all to include them)" if hidden else ""))

    defaults = discovery.get("defaults")
    if defaults:
        lines.append("Gateway defaults: " + ", ".join(f"{key}={value}" for key, value in defaults.items()))

    limits = discovery.get("limits")
    if limits:
        lines.append("Gateway limits: " + ", ".join(f"{key}={value}" for key, value in limits.items()))

    return "\n".join(lines)


def _format_anthropic_models(models: list) -> str:
    """Render Anthropic models API results (newest first) as a compact table."""
    rows = [("model", "name", "released")]
    for m in models:
        created = getattr(m, "created_at", None)
        rows.append((m.id, getattr(m, "display_name", "") or "", created.date().isoformat() if created else ""))

    widths = [max(len(str(row[i])) for row in rows) for i in range(len(rows[0]))]
    lines = ["  ".join(str(value).ljust(width) for value, width in zip(row, widths)).rstrip() for row in rows]
    lines.insert(1, "  ".join("-" * width for width in widths))
    lines.append("")
    lines.append(f"{len(models)} model(s) available to this API key. Use the 'model' column with --llm-model.")
    return "\n".join(lines)


def _discover_anthropic(profile, civic_config, raw_json: bool) -> None:
    import anthropic

    from .providers.anthropic_client import list_models

    api_key = profiles_module.resolve_api_key(profile, _load_llm_credentials(civic_config))
    try:
        sdk_client = anthropic.Anthropic(api_key=api_key) if api_key else anthropic.Anthropic()
        models = list_models(sdk_client)
    except Exception as exc:
        sys.exit(f"Discovery failed: {exc}")
    if raw_json:
        print(json.dumps([m.model_dump(mode="json") for m in models], indent=2))
    else:
        print(_format_anthropic_models(models))


def _discover(profile_name: str, civic_config: str, profiles_path=None, raw_json: bool = False,
              show_all: bool = False) -> None:
    from .providers.entra_gateway_client import fetch_discovery

    profile = profiles_module.get_profile(profile_name, profiles_path)
    if profile.provider == "anthropic":
        _discover_anthropic(profile, civic_config, raw_json)
        return
    if profile.provider != "entra_gateway" or not profile.discovery_url:
        sys.exit(f"Profile '{profile_name}' doesn't support discovery (only anthropic and entra_gateway profiles do)")

    try:
        token_provider, gateway_api_key = profiles_module.build_gateway_token_provider(
            profile, _load_llm_credentials(civic_config)
        )
        discovery = fetch_discovery(profile.discovery_url, token_provider, gateway_api_key)
    except Exception as exc:
        sys.exit(f"Discovery failed: {exc}")
    if raw_json:
        print(json.dumps(discovery, indent=2))
    else:
        print(_format_discovery_summary(discovery, show_all=show_all))


def main(
    command: str,
    task_name=None,
    input_path=None,
    profile=None,
    fake: bool = False,
    fake_output_path=None,
    content_root=None,
    profiles_path=None,
    log_path=None,
    model=None,
    civic_config=None,
    raw_json: bool = False,
    show_all: bool = False,
):
    if command == "discover":
        _discover(profile, civic_config, profiles_path, raw_json=raw_json, show_all=show_all)
        return

    if command == "list-tasks":
        root = Path(content_root) if content_root else registry.DEFAULT_CONTENT_ROOT
        tasks = registry.discover_tasks(root)
        if not tasks:
            print(f"No tasks found under {root}")
            return
        for name, task_dir in tasks.items():
            print(f"{name}\t{task_dir}")
        return

    inputs = json.loads(Path(input_path).read_text())

    client = None
    if fake:
        root = Path(content_root) if content_root else registry.DEFAULT_CONTENT_ROOT
        task = registry.load_task(task_name, root)
        client = FakeClient(responses=[_build_fake_response(task.schema, fake_output_path)])

    try:
        result = run_task(
            task_name, inputs=inputs, profile=profile,
            content_root=content_root, profiles_path=profiles_path, log_path=log_path,
            client=client, credentials=_load_llm_credentials(civic_config), model=model,
        )
    except LLMOutputValidationError as e:
        print(f"Model output failed schema validation:\n{e.raw_output}\n\nError: {e.validation_error}", file=sys.stderr)
        sys.exit(1)
    except LLMRefusalError as e:
        print(f"Model did not complete normally (stop_reason={e.stop_reason!r}):\n{e.raw_text}", file=sys.stderr)
        sys.exit(1)

    print(result.parsed.model_dump_json(indent=2))
    print(
        f"\n--- provider={result.provider} model={result.model} stop_reason={result.stop_reason} "
        f"latency_ms={result.latency_ms:.1f} retries={result.retries} "
        f"input_tokens={result.input_tokens} output_tokens={result.output_tokens} "
        f"cache_read_input_tokens={result.cache_read_input_tokens} ---"
    )


if __name__ == "__main__":
    args = parse_args()
    main(
        command=args.command,
        task_name=getattr(args, "task_name", None),
        input_path=getattr(args, "input_path", None),
        profile=getattr(args, "profile", None),
        fake=getattr(args, "fake", False),
        fake_output_path=getattr(args, "fake_output_path", None),
        content_root=getattr(args, "content_root", None),
        profiles_path=getattr(args, "profiles_path", None),
        log_path=getattr(args, "log_path", None),
        model=getattr(args, "model", None),
        civic_config=getattr(args, "civic_config", None),
        raw_json=getattr(args, "raw_json", False),
        show_all=getattr(args, "show_all", False),
    )
