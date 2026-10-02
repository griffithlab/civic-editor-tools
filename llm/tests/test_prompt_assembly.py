from llm import prompt, registry


def test_system_text_is_stable_across_different_inputs(fixture_content_root):
    """The whole point of the knowledge->instructions->examples ordering: system_text must be
    byte-identical regardless of the per-item inputs, so it's actually a usable cache prefix."""
    task = registry.load_task("sample_task", fixture_content_root)

    system1, user1 = prompt.assemble_prompt(task, {"question": "A"})
    system2, user2 = prompt.assemble_prompt(task, {"question": "B"})

    assert system1 == system2
    assert user1 != user2


def test_assemble_prompt_order_knowledge_before_instructions(fixture_content_root):
    task = registry.load_task("sample_task", fixture_content_root)
    system_text, _ = prompt.assemble_prompt(task, {"question": "A"})

    knowledge_pos = system_text.index("Sample knowledge content")
    instructions_pos = system_text.index("sample_task v0.1.0")
    assert knowledge_pos < instructions_pos


def test_assemble_prompt_includes_examples_after_instructions(fixture_content_root):
    task = registry.load_task("sample_task", fixture_content_root)
    system_text, _ = prompt.assemble_prompt(task, {"question": "A"})

    instructions_pos = system_text.index("sample_task v0.1.0")
    examples_pos = system_text.index("## Examples")
    example_content_pos = system_text.index("What is 2+2?")

    assert instructions_pos < examples_pos < example_content_pos


def test_render_instructions_uses_task_metadata_only(fixture_content_root):
    task = registry.load_task("sample_task", fixture_content_root)
    rendered = prompt.render_instructions(task)
    assert "sample_task v0.1.0" in rendered


def test_render_input_block_generic_nested_rendering():
    inputs = {
        "evidence_item": {"id": 1, "statement": "x"},
        "therapies": ["Foo", "Bar"],
        "source_full_text": None,
    }
    block = prompt.render_input_block(inputs)

    assert "## Evidence Item" in block
    assert "**Id**: 1" in block
    assert "**Statement**: x" in block
    assert "## Therapies" in block
    assert "- Foo" in block and "- Bar" in block
    assert "## Source Full Text" in block
    assert "(none)" in block
