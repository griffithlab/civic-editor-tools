import json

from llm import registry
from llm.providers.base import LLMOutputValidationError, LLMResult
from llm.provenance import ProvenanceLogger


class DummyClient:
    provider_name = "fake"
    model = "fake-model"


def _result() -> LLMResult:
    return LLMResult(
        parsed=None,
        raw_text="raw",
        provider="fake",
        model="fake-model",
        input_tokens=10,
        output_tokens=5,
        cache_read_input_tokens=3,
        cache_creation_input_tokens=0,
        latency_ms=42.0,
        stop_reason="end_turn",
        retries=0,
    )


def test_log_run_writes_one_jsonl_line(tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    logger = ProvenanceLogger(log_path=tmp_path / "p.jsonl")

    logger.log_run(task, "claude", DummyClient(), _result(), "SYSTEM\n\nUSER")

    lines = (tmp_path / "p.jsonl").read_text().strip().split("\n")
    assert len(lines) == 1
    record = json.loads(lines[0])
    assert record["status"] == "ok"
    assert record["task_name"] == "review_evidence_item"
    assert record["provider"] == "fake"
    assert record["usage"]["cache_read_input_tokens"] == 3


def test_log_run_writes_prompt_sidecar_file(tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    logger = ProvenanceLogger(log_path=tmp_path / "p.jsonl")

    logger.log_run(task, "claude", DummyClient(), _result(), "SYSTEM\n\nUSER")

    record = json.loads((tmp_path / "p.jsonl").read_text().strip())
    prompt_file = tmp_path / record["prompt_path"]
    assert prompt_file.exists()
    assert prompt_file.read_text() == "SYSTEM\n\nUSER"
    assert record["prompt_sha256"] in record["prompt_path"]


def test_log_run_does_not_duplicate_sidecar_for_identical_prompt(tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    logger = ProvenanceLogger(log_path=tmp_path / "p.jsonl")

    logger.log_run(task, "claude", DummyClient(), _result(), "SAME TEXT")
    logger.log_run(task, "claude", DummyClient(), _result(), "SAME TEXT")

    prompts_dir = tmp_path / "prompts"
    assert len(list(prompts_dir.glob("*.txt"))) == 1
    assert len((tmp_path / "p.jsonl").read_text().strip().split("\n")) == 2


def test_log_run_records_failure(tmp_path):
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    logger = ProvenanceLogger(log_path=tmp_path / "p.jsonl")
    error = LLMOutputValidationError(raw_output="bad", validation_error="nope")

    logger.log_run(task, "claude", DummyClient(), None, "SYSTEM\n\nUSER", error=error)

    record = json.loads((tmp_path / "p.jsonl").read_text().strip())
    assert record["status"] == "error"
    assert record["raw_output"] == "bad"
    assert record["parsed_output"] is None
    assert "nope" in record["error"]
