"""Pydantic output schema for the review_evidence_item task.

Loaded dynamically by llm.registry.load_task(), which looks for a module-level `SCHEMA`
attribute (not a fixed class name) so this class can be renamed freely.

Reviews ALL currently open pending revisions on one evidence item together, as a collection,
rather than one revision in isolation - matching how a human editor actually reviews (and matching
civic_graphql_utils.accept_revisions(), which already accepts a list of revision ids as one
action). The output is a per-revision assessment list (so a human can still accept any subset)
plus cross-revision issues that don't belong to any single revision.
"""

from enum import Enum
from typing import Optional

from pydantic import BaseModel, ConfigDict, Field


class Recommendation(str, Enum):
    ACCEPT = "accept"
    REJECT = "reject"
    ACCEPT_WITH_CHANGES = "accept_with_changes"
    NEEDS_HUMAN_REVIEW = "needs_human_review"


class Confidence(str, Enum):
    LOW = "low"
    MEDIUM = "medium"
    HIGH = "high"


class IssueSeverity(str, Enum):
    MINOR = "minor"
    MODERATE = "moderate"
    MAJOR = "major"


class Issue(BaseModel):
    model_config = ConfigDict(extra="forbid")

    field: str
    severity: IssueSeverity
    description: str
    supporting_evidence: Optional[str] = None


class SuggestedRevision(BaseModel):
    model_config = ConfigDict(extra="forbid")

    field: str
    current_value: str
    suggested_value: str
    rationale: str


class RevisionAssessment(BaseModel):
    model_config = ConfigDict(extra="forbid")

    # revision_id/field_name are echoed back from the input so a caller can match each assessment
    # to the revision it's about without relying on list order or field_name uniqueness alone.
    revision_id: int
    field_name: str
    recommendation: Recommendation
    confidence: Confidence
    rationale: str


class EvidenceItemReview(BaseModel):
    model_config = ConfigDict(extra="forbid")

    summary: str
    revision_assessments: list[RevisionAssessment]
    issues: list[Issue] = Field(default_factory=list)
    suggested_revisions: list[SuggestedRevision] = Field(default_factory=list)


# Fixed convention the registry loader looks for.
SCHEMA = EvidenceItemReview
