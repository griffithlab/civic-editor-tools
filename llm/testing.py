from .providers.base import LLMClient, LLMResult


class FakeClient(LLMClient):
    """
    No-network stand-in for LLMClient, for testing the registry/prompt/runner/provenance layers
    without hitting any real provider. NOT used to test per-provider retry-on-invalid-JSON logic
    (that lives in providers.base.validate_with_retry(), unit-tested directly with a stub
    retry_fn - see the provider client test files).

    Takes a list of canned LLMResult objects and returns them in order, one per complete() call.
    Every call is recorded in .calls for test assertions.
    """

    provider_name = "fake"
    model = "fake-model"

    def __init__(self, responses: list):
        self._responses = list(responses)
        self.calls = []

    def complete(self, *, system: str, messages: list, output_schema=None, **opts) -> LLMResult:
        self.calls.append({"system": system, "messages": messages, "output_schema": output_schema, "opts": opts})

        if not self._responses:
            raise RuntimeError("FakeClient has no queued responses left")

        return self._responses.pop(0)
