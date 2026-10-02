"""
Regression guard for the REAL, human-authored civic-editorial-review/tasks/review_evidence_item/
content (as opposed to the synthetic fixtures/content_root/ used by the other test files).

Re-run this whenever instructions.md, schema.py, or knowledge/evidence_item_basics.md is
hand-edited, to catch frontmatter/schema breakage introduced while writing real SOP content.
"""

from llm import registry, prompt


def test_review_evidence_item_task_loads():
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    assert task.name == "review_evidence_item"
    assert task.description
    assert task.version


def test_review_evidence_item_requires_structured_output():
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    assert "structured_output" in task.requires


def test_review_evidence_item_references_evidence_item_basics_knowledge():
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    assert "evidence_item_basics.md" in task.knowledge
    # and that referenced file must actually exist and load cleanly
    text = registry.load_knowledge(registry.DEFAULT_CONTENT_ROOT, task.knowledge)
    assert len(text) > 0


def test_review_evidence_item_schema_field_names():
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    assert set(task.schema.model_fields.keys()) == {
        "summary", "revision_assessments", "issues", "suggested_revisions",
    }


def test_review_evidence_item_prompt_assembles_without_error(review_evidence_item_input):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    system_text, user_text = prompt.assemble_prompt(task, review_evidence_item_input)

    assert len(system_text) > 0
    assert review_evidence_item_input["evidence_item"]["statement"] in user_text


def test_review_evidence_item_examples_load_and_validate():
    """Guards against a future examples/ addition breaking the NNN_input/NNN_output pairing
    convention, or an example output drifting out of sync with the schema."""
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    examples = registry.load_examples(task.examples_dir)
    assert len(examples) >= 1
    for example_input, example_output in examples:
        assert "evidence_item" in example_input
        assert "revisions" in example_input
        validated = task.schema.model_validate(example_output)  # raises if invalid

        # every input revision must get exactly one matching assessment, by revision_id
        input_ids = {revision["revision_id"] for revision in example_input["revisions"]}
        output_ids = {assessment.revision_id for assessment in validated.revision_assessments}
        assert input_ids == output_ids
