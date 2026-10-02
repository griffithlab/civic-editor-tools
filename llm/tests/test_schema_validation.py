import pytest
from pydantic import ValidationError

from llm import registry


@pytest.fixture
def evidence_item_review_schema():
    task = registry.load_task("review_evidence_item", registry.DEFAULT_CONTENT_ROOT)
    return task.schema


def test_valid_minimal_payload(evidence_item_review_schema):
    instance = evidence_item_review_schema(
        summary="ok",
        revision_assessments=[
            {"revision_id": 1, "field_name": "description", "recommendation": "accept", "confidence": "high", "rationale": "fine"}
        ],
    )
    assert instance.issues == []
    assert instance.suggested_revisions == []


def test_valid_full_payload(evidence_item_review_schema):
    instance = evidence_item_review_schema(
        summary="mostly fine, one cross-field issue",
        revision_assessments=[
            {"revision_id": 1, "field_name": "significance", "recommendation": "accept", "confidence": "high", "rationale": "supported"},
            {"revision_id": 2, "field_name": "evidence_direction", "recommendation": "accept_with_changes", "confidence": "medium", "rationale": "mostly right"},
        ],
        issues=[{"field": "description", "severity": "minor", "description": "typo", "supporting_evidence": None}],
        suggested_revisions=[
            {"field": "description", "current_value": "a", "suggested_value": "b", "rationale": "clarity"}
        ],
    )
    assert len(instance.revision_assessments) == 2
    assert instance.revision_assessments[0].recommendation.value == "accept"
    assert instance.issues[0].severity.value == "minor"
    assert instance.suggested_revisions[0].suggested_value == "b"


def test_invalid_recommendation_enum_rejected(evidence_item_review_schema):
    with pytest.raises(ValidationError):
        evidence_item_review_schema(
            summary="ok",
            revision_assessments=[
                {"revision_id": 1, "field_name": "x", "recommendation": "maybe", "confidence": "high", "rationale": "y"}
            ],
        )


def test_missing_required_field_rejected(evidence_item_review_schema):
    with pytest.raises(ValidationError):
        evidence_item_review_schema(revision_assessments=[])  # missing summary


def test_revision_assessments_required(evidence_item_review_schema):
    """Unlike issues/suggested_revisions, revision_assessments has no default - every reviewed
    revision must get an assessment, so an empty/omitted list should still validate (zero
    revisions is a valid, if unusual, state) but the field itself must be present."""
    with pytest.raises(ValidationError):
        evidence_item_review_schema(summary="ok")


def test_extra_field_rejected(evidence_item_review_schema):
    with pytest.raises(ValidationError):
        evidence_item_review_schema(
            summary="ok",
            revision_assessments=[],
            made_up_field="x",
        )


def test_nested_revision_assessment_extra_field_rejected(evidence_item_review_schema):
    with pytest.raises(ValidationError):
        evidence_item_review_schema(
            summary="ok",
            revision_assessments=[
                {"revision_id": 1, "field_name": "x", "recommendation": "accept", "confidence": "high", "rationale": "y", "made_up": "z"}
            ],
        )


def test_nested_issue_extra_field_rejected(evidence_item_review_schema):
    with pytest.raises(ValidationError):
        evidence_item_review_schema(
            summary="ok",
            revision_assessments=[],
            issues=[{"field": "x", "severity": "minor", "description": "y", "made_up": "z"}],
        )


def test_model_validate_json_round_trip(evidence_item_review_schema):
    payload = (
        '{"summary": "no", "revision_assessments": ['
        '{"revision_id": 1, "field_name": "x", "recommendation": "reject", "confidence": "low", "rationale": "y"}'
        '], "issues": [], "suggested_revisions": []}'
    )
    instance = evidence_item_review_schema.model_validate_json(payload)
    assert instance.revision_assessments[0].recommendation.value == "reject"


def test_json_schema_has_additional_properties_false(evidence_item_review_schema):
    """extra='forbid' must produce additionalProperties: false, which both Anthropic strict mode
    and OpenAI-compatible json_schema mode rely on."""
    schema_dict = evidence_item_review_schema.model_json_schema()
    assert schema_dict.get("additionalProperties") is False
