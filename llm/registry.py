import dataclasses
import importlib.util
import json
from pathlib import Path

import yaml

_MODULE_DIR = Path(__file__).resolve().parent
DEFAULT_CONTENT_ROOT = _MODULE_DIR / "civic-editorial-review"

_REQUIRED_FRONTMATTER_KEYS = {"name", "description", "version", "knowledge", "requires", "recommended_profiles"}


@dataclasses.dataclass
class TaskDefinition:
    name: str
    description: str
    version: str
    knowledge: list
    requires: list
    recommended_profiles: list
    instructions_source: str  # raw Jinja2 template text, frontmatter already stripped
    schema: type
    examples_dir: Path
    task_dir: Path


def discover_tasks(content_root) -> dict:
    """Return {task_name: task_dir} for every tasks/*/instructions.md under content_root."""
    tasks_dir = Path(content_root) / "tasks"
    if not tasks_dir.exists():
        return {}

    discovered = {}
    for child in sorted(tasks_dir.iterdir()):
        if child.is_dir() and (child / "instructions.md").exists():
            discovered[child.name] = child

    return discovered


def load_task(name: str, content_root) -> TaskDefinition:
    """Parse instructions.md frontmatter, dynamically import <task_dir>/schema.py, and return a
    fully assembled TaskDefinition."""
    task_dir = Path(content_root) / "tasks" / name
    instructions_path = task_dir / "instructions.md"
    schema_path = task_dir / "schema.py"
    examples_dir = task_dir / "examples"

    if not instructions_path.exists():
        raise FileNotFoundError(f"No such task '{name}': missing {instructions_path}")
    if not schema_path.exists():
        raise FileNotFoundError(f"Task '{name}' is missing {schema_path}")

    frontmatter, instructions_source = _split_frontmatter(instructions_path.read_text())

    missing = _REQUIRED_FRONTMATTER_KEYS - frontmatter.keys()
    if missing:
        raise ValueError(f"{instructions_path} frontmatter is missing required key(s): {sorted(missing)}")

    schema = _import_schema(schema_path)

    return TaskDefinition(
        name=frontmatter["name"],
        description=frontmatter["description"],
        version=str(frontmatter["version"]),
        knowledge=frontmatter["knowledge"] or [],
        requires=frontmatter["requires"] or [],
        recommended_profiles=frontmatter["recommended_profiles"] or [],
        instructions_source=instructions_source,
        schema=schema,
        examples_dir=examples_dir,
        task_dir=task_dir,
    )


def load_knowledge(content_root, filenames: list) -> str:
    """Concatenate knowledge/<filename> contents, in the order given, each under a heading.
    Raises if a listed file is missing - fail fast rather than silently dropping stable context."""
    knowledge_dir = Path(content_root) / "knowledge"
    sections = []

    for filename in filenames:
        path = knowledge_dir / filename
        if not path.exists():
            raise FileNotFoundError(f"Task references a knowledge file that doesn't exist: {path}")
        sections.append(f"## {filename}\n\n{path.read_text().strip()}")

    return "\n\n".join(sections)


def load_examples(examples_dir) -> list:
    """Load NNN_input.json/NNN_output.json pairs from examples_dir in sorted (numeric) order.
    Returns [] if the directory doesn't exist or has no examples yet."""
    examples_dir = Path(examples_dir)
    if not examples_dir.exists():
        return []

    pairs = []
    for input_path in sorted(examples_dir.glob("*_input.json")):
        prefix = input_path.name[: -len("_input.json")]
        output_path = examples_dir / f"{prefix}_output.json"
        if not output_path.exists():
            raise FileNotFoundError(f"Example {input_path} has no matching {output_path}")
        pairs.append((json.loads(input_path.read_text()), json.loads(output_path.read_text())))

    return pairs


def _split_frontmatter(text: str) -> tuple:
    """Split a '---\\n<yaml>\\n---\\n<body>' file into (frontmatter_dict, body_str)."""
    if not text.startswith("---"):
        raise ValueError("instructions.md must start with a '---' YAML frontmatter block")

    parts = text.split("---", 2)
    if len(parts) < 3:
        raise ValueError("instructions.md must have a closing '---' after the frontmatter block")

    _, frontmatter_text, body = parts
    frontmatter = yaml.safe_load(frontmatter_text) or {}

    return frontmatter, body.lstrip("\n")


_schema_cache = {}


def _import_schema(schema_path: Path) -> type:
    """
    Dynamically import a task's schema.py and return its module-level SCHEMA attribute. Cached by
    resolved path so repeated load_task() calls for the same task (run_task() calls load_task()
    on every invocation) return the SAME class object rather than re-executing schema.py and
    minting a new, structurally-identical-but-distinct class each time - important since callers
    may reasonably import a task's schema directly (e.g. `from ... import schema; schema.SCHEMA`)
    and expect `isinstance`/`is` checks against a run_task() result to hold.

    Not registered in sys.modules - each schema.py is a standalone file with no need to be
    importable from elsewhere, and this avoids any risk of module-name collisions across
    tasks/content roots.
    """
    resolved_path = schema_path.resolve()
    if resolved_path in _schema_cache:
        return _schema_cache[resolved_path]

    spec = importlib.util.spec_from_file_location(f"_llm_schema_{schema_path.parent.name}", schema_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "SCHEMA"):
        raise ValueError(f"{schema_path} must define a module-level `SCHEMA = <PydanticModel>`")

    _schema_cache[resolved_path] = module.SCHEMA
    return module.SCHEMA
