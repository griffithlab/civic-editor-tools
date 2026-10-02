"""
Smoke-test CLI for llm. Library usage is `from llm import run_task`; this CLI exists only to
make plumbing (task discovery, prompt assembly, a provider call, provenance logging) easy to
exercise by hand without writing a script.

    python -m llm run review_evidence_item --input <file.json> --profile local [--fake]
    python -m llm list-tasks
"""

import argparse
import json
import sys
from pathlib import Path

from . import registry
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
    run_parser.add_argument("--content-root", dest="content_root", type=str, default=None)
    run_parser.add_argument("--profiles-file", dest="profiles_path", type=str, default=None)
    run_parser.add_argument("--log-file", dest="log_path", type=str, default=None)

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
):
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
            client=client,
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
        content_root=args.content_root,
        profiles_path=getattr(args, "profiles_path", None),
        log_path=getattr(args, "log_path", None),
    )
