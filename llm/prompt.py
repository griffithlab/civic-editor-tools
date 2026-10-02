import json

import jinja2

from . import registry


def render_instructions(task: "registry.TaskDefinition") -> str:
    """
    Render a task's instructions.md body as a Jinja2 template.

    Context is task metadata ONLY (name, description, version, requires, knowledge,
    recommended_profiles) - deliberately NOT the per-item `inputs` dict. This keeps the rendered
    instructions byte-identical across every run_task() call for a given task/version, which is
    what makes the whole knowledge+instructions+examples block a stable, cacheable prefix.
    """
    template = jinja2.Template(task.instructions_source)
    return template.render(
        name=task.name,
        description=task.description,
        version=task.version,
        requires=task.requires,
        knowledge=task.knowledge,
        recommended_profiles=task.recommended_profiles,
    )


def render_input_block(inputs: dict) -> str:
    """
    Render the per-item `inputs` dict into a fixed, generic (not task-specific) labelled Markdown
    block. This is the ONLY part of the prompt that varies call-to-call for a given task.
    """
    sections = [f"## {_title(key)}\n\n{_render_value(value)}" for key, value in inputs.items()]
    return "\n\n".join(sections)


def assemble_prompt(task: "registry.TaskDefinition", inputs: dict) -> tuple:
    """
    Returns (system_text, user_text):
      system_text = knowledge (stable) + rendered instructions (stable) + examples (stable)
      user_text   = render_input_block(inputs)  (volatile)
    Order matches the task-framework spec exactly: knowledge -> instructions -> examples ->
    per-item input. This is exactly what makes the Anthropic cache breakpoint (placed at the end
    of system_text, see providers/anthropic_client.py) meaningful.
    """
    content_root = task.task_dir.parent.parent
    sections = []

    if task.knowledge:
        sections.append(registry.load_knowledge(content_root, task.knowledge))

    sections.append(render_instructions(task))

    examples = registry.load_examples(task.examples_dir)
    if examples:
        sections.append(_render_examples(examples))

    system_text = "\n\n".join(sections)
    user_text = render_input_block(inputs)

    return system_text, user_text


def _title(key: str) -> str:
    return key.replace("_", " ").title()


def _render_value(value, indent: int = 0) -> str:
    prefix = "  " * indent

    if value is None:
        return f"{prefix}(none)"

    if isinstance(value, dict):
        lines = []
        for key, sub_value in value.items():
            if isinstance(sub_value, (dict, list)):
                lines.append(f"{prefix}- **{_title(key)}**:")
                lines.append(_render_value(sub_value, indent + 1))
            else:
                lines.append(f"{prefix}- **{_title(key)}**: {sub_value}")
        return "\n".join(lines)

    if isinstance(value, list):
        if not value:
            return f"{prefix}(none)"
        lines = []
        for item in value:
            if isinstance(item, (dict, list)):
                lines.append(_render_value(item, indent))
            else:
                lines.append(f"{prefix}- {item}")
        return "\n".join(lines)

    return f"{prefix}{value}"


def _render_examples(examples: list) -> str:
    blocks = []
    for i, (example_input, example_output) in enumerate(examples, start=1):
        blocks.append(
            f"### Example {i}\n\n"
            f"Input:\n```json\n{json.dumps(example_input, indent=2)}\n```\n\n"
            f"Output:\n```json\n{json.dumps(example_output, indent=2)}\n```"
        )
    return "## Examples\n\n" + "\n\n".join(blocks)
