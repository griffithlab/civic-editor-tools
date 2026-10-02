import contextlib
import io

from llm import registry, run_task
from llm.providers.base import LLMResult
from llm.testing import FakeClient


def _make_result(schema) -> LLMResult:
    parsed = schema(
        summary="fine",
        revision_assessments=[
            {"revision_id": 90001, "field_name": "x", "recommendation": "accept", "confidence": "high", "rationale": "ok"}
        ],
    )
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


def test_warning_printed_when_profile_lacks_required_support(review_evidence_item_input, tmp_path):
    """review_evidence_item requires structured_output; the committed 'local' profile deliberately
    doesn't list it, since most local models don't reliably honor strict JSON schema mode."""
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    client = FakeClient(responses=[_make_result(task.schema)])

    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        run_task(
            "review_evidence_item", inputs=review_evidence_item_input, profile="local",
            log_path=tmp_path / "p.jsonl", client=client,
        )

    output = captured.getvalue()
    assert "WARNING" in output
    assert "structured_output" in output
    assert "local" in output


def test_no_warning_when_profile_supports_requirements(review_evidence_item_input, tmp_path):
    """The committed 'claude' profile lists structured_output, so no warning should print."""
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    client = FakeClient(responses=[_make_result(task.schema)])

    captured = io.StringIO()
    with contextlib.redirect_stdout(captured):
        run_task(
            "review_evidence_item", inputs=review_evidence_item_input, profile="claude",
            log_path=tmp_path / "p.jsonl", client=client,
        )

    assert "WARNING" not in captured.getvalue()


def test_no_crash_with_unknown_profile_name_and_explicit_client(review_evidence_item_input, tmp_path):
    """An unregistered profile name combined with an explicit `client` (the FakeClient test seam)
    must not raise - there's simply nothing to warn about."""
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    client = FakeClient(responses=[_make_result(task.schema)])

    result = run_task(
        "review_evidence_item", inputs=review_evidence_item_input, profile="totally-made-up",
        log_path=tmp_path / "p.jsonl", client=client,
    )

    assert result.parsed.summary == "fine"
