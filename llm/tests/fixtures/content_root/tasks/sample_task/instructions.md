---
name: sample_task
description: A tiny synthetic task used only by llm's own test suite.
version: "0.1.0"
knowledge:
  - sample_knowledge.md
requires: []
recommended_profiles:
  - claude
---
This is the sample_task instructions body ({{ name }} v{{ version }}). It exists only to test
frontmatter parsing, Jinja2 rendering of task metadata, and prompt assembly ordering.
