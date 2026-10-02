import dataclasses
import json
from abc import ABC, abstractmethod
from typing import Callable, Optional

import pydantic
from pydantic import BaseModel


@dataclasses.dataclass
class LLMResult:
    """The result of one LLMClient.complete() call."""

    parsed: Optional[BaseModel]
    raw_text: str
    provider: str
    model: str
    input_tokens: int
    output_tokens: int
    cache_read_input_tokens: Optional[int]
    cache_creation_input_tokens: Optional[int]
    latency_ms: float
    stop_reason: str
    retries: int = 0


class LLMOutputValidationError(Exception):
    """Raised when a model's output still fails schema validation after one retry."""

    def __init__(self, raw_output: str, validation_error: str):
        self.raw_output = raw_output
        self.validation_error = validation_error
        super().__init__(f"Model output failed schema validation twice: {validation_error}")


class LLMRefusalError(Exception):
    """Raised when the model's stop_reason indicates it did not produce a usable completion
    (e.g. a safety refusal, or being cut off at max_tokens) - these are not schema validation
    failures, so retrying validate_with_retry() against them would be meaningless."""

    def __init__(self, stop_reason: str, raw_text: str):
        self.stop_reason = stop_reason
        self.raw_text = raw_text
        super().__init__(f"Model did not complete normally (stop_reason={stop_reason!r})")


class LLMClient(ABC):
    """Common interface every provider-specific client implements."""

    provider_name: str

    @abstractmethod
    def complete(
        self,
        *,
        system: str,
        messages: list,
        output_schema: Optional[type] = None,
        **opts,
    ) -> LLMResult:
        """Send one request. If output_schema is given, the returned LLMResult.parsed is a
        validated instance of it; otherwise .parsed is None and .raw_text holds the plain text
        response."""
        raise NotImplementedError


def validate_with_retry(
    schema: type[BaseModel],
    raw_text: str,
    retry_fn: Callable[[str], str],
) -> tuple[BaseModel, int]:
    """
    Validate raw_text against schema. On failure, call retry_fn(error_message) once to get a
    corrected raw_text, then validate again.

    Returns (parsed, retries_used). Raises LLMOutputValidationError if the second attempt also
    fails validation.
    """
    try:
        return schema.model_validate_json(raw_text), 0
    except (json.JSONDecodeError, pydantic.ValidationError) as first_error:
        second_raw = retry_fn(str(first_error))
        try:
            return schema.model_validate_json(second_raw), 1
        except (json.JSONDecodeError, pydantic.ValidationError) as second_error:
            raise LLMOutputValidationError(raw_output=second_raw, validation_error=str(second_error)) from second_error
