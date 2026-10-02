---
name: civic-editorial-review
description: >
  Editorial review assistance for content that already exists in the CIViC knowledgebase
  (evidence items, assertions, and their pending revisions) — evaluating whether a pending
  revision should be accepted, identifying errors or unsupported claims, and suggesting
  corrections. This is NOT for curating new evidence items from the literature, and it is
  NOT a substitute for human moderation: every recommendation is advisory, and a human
  editor makes the final accept/reject decision.
---

# civic-editorial-review

This directory holds the shared knowledge and per-task instructions used by the `llm` package's
provider layer (see `llm/` at the repo root) when it asks an LLM to help review pending CIViC
revisions. It is structured so it can also be used directly as a Claude Code skill.

## Scope

In scope:
- Reviewing all currently open pending revisions on an existing evidence item or assertion
  together, as a collection, against the item's current fields and any supplied source material
  (ideally the full text of the cited publication — CIViC curation is expected to be based on
  full articles, not abstracts alone). Revisions on the same item routinely interact (e.g. two
  revisions together forming one coherent change), so they're reviewed together rather than one
  at a time.
- Flagging factual inconsistencies, unsupported claims, or formatting/style issues in a
  revision, or across the set of pending revisions.
- Suggesting alternative wording or values for a revision, for a human to consider.

Out of scope:
- Curating brand-new evidence items from primary literature.
- Making or submitting any change to CIViC. Every task here only ever returns a structured,
  advisory recommendation; the calling script (`review_evidence_items.py`) is the only thing
  that can submit a mutation, and only after an interactive human confirmation.

## Layout

- `knowledge/` — background reference material shared across tasks (e.g. what fields an
  evidence item has, how CIViC's evidence levels/significance values are defined). Loaded as a
  stable prefix ahead of any task's own instructions, so it's a good place for anything that
  doesn't change per-task.
- `tasks/<task_name>/` — one directory per review task:
  - `instructions.md` — YAML frontmatter (name, description, version, knowledge, requires,
    recommended_profiles) plus a Jinja2 prompt template body.
  - `schema.py` — a Pydantic model (exposed as module-level `SCHEMA`) describing the task's
    structured output.
  - `examples/` — optional `NNN_input.json`/`NNN_output.json` pairs used as few-shot examples.

## Tasks

| Task | Path | Purpose |
|---|---|---|
| `review_evidence_item` | `tasks/review_evidence_item/` | Review all currently open pending revisions on an evidence item together and return a per-revision accept/reject/accept-with-changes/needs-human-review recommendation, plus any cross-revision issues. |
