import pytest

from llm import registry


def test_discover_tasks_finds_sample_task(fixture_content_root):
    tasks = registry.discover_tasks(fixture_content_root)
    assert tasks == {"sample_task": fixture_content_root / "tasks" / "sample_task"}


def test_discover_tasks_empty_when_no_tasks_dir(tmp_path):
    assert registry.discover_tasks(tmp_path) == {}


def test_load_task_parses_frontmatter(fixture_content_root):
    task = registry.load_task("sample_task", fixture_content_root)
    assert task.name == "sample_task"
    assert task.version == "0.1.0"
    assert task.knowledge == ["sample_knowledge.md"]
    assert task.requires == []
    assert task.recommended_profiles == ["claude"]


def test_load_task_strips_frontmatter_from_instructions_source(fixture_content_root):
    task = registry.load_task("sample_task", fixture_content_root)
    assert "---" not in task.instructions_source
    assert "sample_task instructions body" in task.instructions_source


def test_load_task_imports_schema(fixture_content_root):
    task = registry.load_task("sample_task", fixture_content_root)
    instance = task.schema(answer="hi")
    assert instance.answer == "hi"


def test_import_schema_is_cached_across_load_task_calls(fixture_content_root):
    """Regression test: load_task() re-imports schema.py on every call (run_task() calls it on
    every invocation), so without caching, two calls for the same task would return
    structurally-identical but distinct classes, breaking isinstance/`is` checks."""
    task1 = registry.load_task("sample_task", fixture_content_root)
    task2 = registry.load_task("sample_task", fixture_content_root)
    assert task1.schema is task2.schema


def test_load_task_missing_required_frontmatter_key_raises(tmp_path):
    task_dir = tmp_path / "tasks" / "broken_task"
    task_dir.mkdir(parents=True)
    (task_dir / "instructions.md").write_text("---\nname: broken_task\n---\nbody")
    (task_dir / "schema.py").write_text("from pydantic import BaseModel\nclass X(BaseModel):\n    pass\nSCHEMA = X\n")

    with pytest.raises(ValueError, match="missing required key"):
        registry.load_task("broken_task", tmp_path)


def test_load_task_missing_schema_file_raises(tmp_path):
    task_dir = tmp_path / "tasks" / "no_schema"
    task_dir.mkdir(parents=True)
    (task_dir / "instructions.md").write_text(
        "---\nname: no_schema\ndescription: x\nversion: '0.1'\n"
        "knowledge: []\nrequires: []\nrecommended_profiles: []\n---\nbody"
    )

    with pytest.raises(FileNotFoundError):
        registry.load_task("no_schema", tmp_path)


def test_load_task_missing_instructions_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        registry.load_task("does_not_exist", tmp_path)


def test_import_schema_missing_schema_attribute_raises(tmp_path):
    task_dir = tmp_path / "tasks" / "no_schema_attr"
    task_dir.mkdir(parents=True)
    (task_dir / "instructions.md").write_text(
        "---\nname: no_schema_attr\ndescription: x\nversion: '0.1'\n"
        "knowledge: []\nrequires: []\nrecommended_profiles: []\n---\nbody"
    )
    (task_dir / "schema.py").write_text("x = 1\n")  # no module-level SCHEMA

    with pytest.raises(ValueError, match="SCHEMA"):
        registry.load_task("no_schema_attr", tmp_path)


def test_load_knowledge_concatenates_in_order(fixture_content_root):
    text = registry.load_knowledge(fixture_content_root, ["sample_knowledge.md"])
    assert "## sample_knowledge.md" in text
    assert "Sample knowledge content" in text


def test_load_knowledge_missing_file_raises(fixture_content_root):
    with pytest.raises(FileNotFoundError):
        registry.load_knowledge(fixture_content_root, ["does_not_exist.md"])


def test_load_examples_returns_pairs_in_order(fixture_content_root):
    examples_dir = fixture_content_root / "tasks" / "sample_task" / "examples"
    examples = registry.load_examples(examples_dir)
    assert len(examples) == 1
    example_input, example_output = examples[0]
    assert example_input == {"question": "What is 2+2?"}
    assert example_output == {"answer": "4"}


def test_load_examples_missing_dir_returns_empty_list(tmp_path):
    assert registry.load_examples(tmp_path / "does_not_exist") == []


def test_load_examples_unmatched_input_raises(tmp_path):
    (tmp_path / "001_input.json").write_text("{}")
    with pytest.raises(FileNotFoundError):
        registry.load_examples(tmp_path)
