import pytest
from pydantic import BaseModel, ConfigDict

from llm.providers.base import LLMOutputValidationError, validate_with_retry


class Foo(BaseModel):
    x: int


def test_validate_with_retry_success_first_try():
    parsed, retries = validate_with_retry(Foo, '{"x": 1}', lambda err: "unused")
    assert parsed.x == 1
    assert retries == 0


def test_validate_with_retry_success_on_retry():
    calls = []

    def retry_fn(err):
        calls.append(err)
        return '{"x": 2}'

    parsed, retries = validate_with_retry(Foo, "not json", retry_fn)
    assert parsed.x == 2
    assert retries == 1
    assert len(calls) == 1


def test_validate_with_retry_raises_after_two_failures():
    with pytest.raises(LLMOutputValidationError) as exc_info:
        validate_with_retry(Foo, "not json", lambda err: "still not json")
    assert exc_info.value.raw_output == "still not json"
    assert exc_info.value.validation_error


def test_validate_with_retry_retry_fn_receives_error_message():
    received = {}

    def retry_fn(err):
        received["err"] = err
        return '{"x": 3}'

    validate_with_retry(Foo, "not json at all", retry_fn)
    assert received["err"]


def test_validate_with_retry_catches_schema_violation_not_just_bad_json():
    """A syntactically valid JSON document that violates the schema (wrong type) must also be
    treated as a validation failure eligible for retry, not just a JSON parse error."""
    parsed, retries = validate_with_retry(Foo, '{"x": "not an int"}', lambda err: '{"x": 5}')
    assert parsed.x == 5
    assert retries == 1


def test_validate_with_retry_rejects_extra_fields_with_forbid_config():
    class Strict(BaseModel):
        model_config = ConfigDict(extra="forbid")
        x: int

    parsed, retries = validate_with_retry(Strict, '{"x": 1, "y": 2}', lambda err: '{"x": 1}')
    assert retries == 1
    assert parsed.x == 1
