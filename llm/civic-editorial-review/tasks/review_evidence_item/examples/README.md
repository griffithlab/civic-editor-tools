# Examples for `review_evidence_item`

Follow this convention so `llm.registry.load_examples()` can discover and order them
automatically:

- Name files in matched, numbered pairs: `001_input.json`, `001_output.json`, `002_input.json`,
  `002_output.json`, and so on.
- `NNN_input.json` must match the `inputs` dict shape documented in
  `civic-editorial-review/tasks/review_evidence_item/instructions.md`: an `evidence_item` object,
  a `revisions` list (each with `revision_id`, `field_name`, `current_value`, `proposed_value`,
  `comment`, `revisor_id`), and an optional `source_full_text` string — ideally the full text of
  the cited publication, since CIViC curation is expected to rely on full articles rather than
  abstracts alone.
- `NNN_output.json` must be a JSON object that validates against `EvidenceItemReview` in
  `../schema.py`, with one `revision_assessments` entry per revision in the input, matched by
  `revision_id`.
- Examples are loaded in ascending numeric order and inserted into the prompt between the task
  instructions and the per-item input being reviewed, as few-shot demonstrations.

Example `001` uses real, public data from CIViC evidence item 8123
(https://civicdb.org/evidence/8123/revisions) — both of its real pending revisions, reviewed
together, since it's a good illustration of exactly the kind of cross-revision consistency
reasoning this task exists for. Real examples are fine to use (this one was a deliberate choice,
not a default), but be mindful that an example bakes in one specific judgment call as an implicit
"gold standard" for as long as it stays in this directory — if you add another real example, prefer
one where the correct call is genuinely well-supported by the evidence item's own stated fields
(as PP4 is here), not a borderline or disputed case.
