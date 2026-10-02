import json

from llm import registry, run_task
from llm.providers.base import LLMClient, LLMOutputValidationError, LLMResult
from llm.testing import FakeClient


def _make_result(schema, **fields) -> LLMResult:
    parsed = schema(**fields)
    return LLMResult(
        parsed=parsed,
        raw_text=parsed.model_dump_json(),
        provider="fake",
        model="fake-model",
        input_tokens=1,
        output_tokens=1,
        cache_read_input_tokens=None,
        cache_creation_input_tokens=None,
        latency_ms=1.0,
        stop_reason="end_turn",
        retries=0,
    )


def _canned_review(schema, revision_ids) -> LLMResult:
    return _make_result(
        schema,
        summary="fine",
        revision_assessments=[
            {"revision_id": rid, "field_name": "x", "recommendation": "accept", "confidence": "high", "rationale": "ok"}
            for rid in revision_ids
        ],
    )


def test_run_task_returns_fake_client_result(review_evidence_item_input, tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    client = FakeClient(responses=[_canned_review(task.schema, [90001, 90002])])

    result = run_task(
        "review_evidence_item", inputs=review_evidence_item_input, profile="claude",
        log_path=tmp_path / "provenance.jsonl", client=client,
    )

    assert result.parsed.summary == "fine"
    assert len(result.parsed.revision_assessments) == 2
    assert len(client.calls) == 1


def test_run_task_passes_assembled_prompt_to_client(review_evidence_item_input, tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    client = FakeClient(responses=[_canned_review(task.schema, [90001, 90002])])

    run_task(
        "review_evidence_item", inputs=review_evidence_item_input, profile="claude",
        log_path=tmp_path / "provenance.jsonl", client=client,
    )

    call = client.calls[0]
    assert "CIViC Evidence Item Basics" in call["system"]
    assert review_evidence_item_input["evidence_item"]["statement"] in call["messages"][0]["content"]
    assert call["output_schema"] is task.schema


def test_run_task_logs_provenance_on_success(review_evidence_item_input, tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    client = FakeClient(responses=[_canned_review(task.schema, [90001, 90002])])
    log_path = tmp_path / "provenance.jsonl"

    run_task(
        "review_evidence_item", inputs=review_evidence_item_input, profile="claude",
        log_path=log_path, client=client,
    )

    lines = log_path.read_text().strip().split("\n")
    assert len(lines) == 1
    record = json.loads(lines[0])
    assert record["status"] == "ok"
    assert record["task_name"] == "review_evidence_item"
    assert record["profile"] == "claude"


def test_run_task_propagates_and_logs_validation_error(review_evidence_item_input, tmp_path):
    class AlwaysFailsClient(LLMClient):
        provider_name = "fake"
        model = "fake-model"

        def complete(self, **kwargs):
            raise LLMOutputValidationError(raw_output="garbage", validation_error="nope")

    log_path = tmp_path / "provenance.jsonl"

    try:
        run_task(
            "review_evidence_item", inputs=review_evidence_item_input, profile="claude",
            log_path=log_path, client=AlwaysFailsClient(),
        )
        assert False, "expected LLMOutputValidationError to propagate"
    except LLMOutputValidationError:
        pass

    lines = log_path.read_text().strip().split("\n")
    assert len(lines) == 1
    record = json.loads(lines[0])
    assert record["status"] == "error"
    assert record["raw_output"] == "garbage"
    assert record["parsed_output"] is None


def test_run_task_uses_custom_content_root(fixture_content_root, tmp_path):
    task = registry.load_task("sample_task", fixture_content_root)
    client = FakeClient(responses=[_make_result(task.schema, answer="42")])

    result = run_task(
        "sample_task", inputs={"question": "What is 6*7?"}, profile="claude",
        content_root=fixture_content_root, log_path=tmp_path / "p.jsonl", client=client,
    )

    assert result.parsed.answer == "42"
