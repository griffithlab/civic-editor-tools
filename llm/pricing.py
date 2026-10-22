import dataclasses
from typing import Optional


@dataclasses.dataclass
class ModelPricing:
    """
    A model's token prices in USD per million tokens. cache_write/cache_read are optional: when a
    provider doesn't price them separately they're billed at the plain input rate. `source` says
    where the prices came from (shown alongside estimates so a stale price table is noticeable).
    """

    input: float
    output: float
    cache_write: Optional[float] = None
    cache_read: Optional[float] = None
    source: str = ""

    @classmethod
    def from_profile_entry(cls, entry: dict, source: str) -> "ModelPricing":
        """From one model's entry in a profile's [<profile>.pricing] table in profiles.toml."""
        return cls(
            input=float(entry["input"]),
            output=float(entry["output"]),
            cache_write=float(entry["cache_write"]) if "cache_write" in entry else None,
            cache_read=float(entry["cache_read"]) if "cache_read" in entry else None,
            source=source,
        )

    @classmethod
    def from_gateway_cost(cls, cost: dict) -> Optional["ModelPricing"]:
        """From a model's 'cost' object in a gateway discovery response (prices per 1k tokens),
        or None if it doesn't give both an input and an output price in USD."""
        if not cost or cost.get("currency", "USD") != "USD":
            return None
        if cost.get("input_per_1k_tokens") is None or cost.get("output_per_1k_tokens") is None:
            return None
        cached = cost.get("cached_input_per_1k_tokens")
        return cls(
            input=cost["input_per_1k_tokens"] * 1000,
            output=cost["output_per_1k_tokens"] * 1000,
            cache_read=cached * 1000 if cached is not None else None,
            source="gateway discovery",
        )


def estimate_cost(result, pricing: Optional[ModelPricing]) -> Optional[float]:
    """
    Estimated USD cost of an LLMResult, or None if pricing is unknown. result.input_tokens is
    the total prompt size; its cache_read_input_tokens/cache_creation_input_tokens subsets are
    priced at the cache rates (falling back to the input rate) and the rest at the input rate.
    Output tokens include any reasoning/thinking tokens, which are billed as output.
    """
    if pricing is None:
        return None
    cache_read = result.cache_read_input_tokens or 0
    cache_write = result.cache_creation_input_tokens or 0
    uncached = max(result.input_tokens - cache_read - cache_write, 0)
    read_rate = pricing.cache_read if pricing.cache_read is not None else pricing.input
    write_rate = pricing.cache_write if pricing.cache_write is not None else pricing.input
    return (
        uncached * pricing.input
        + cache_read * read_rate
        + cache_write * write_rate
        + result.output_tokens * pricing.output
    ) / 1_000_000


def format_cost(cost: Optional[float]) -> str:
    """'$0.0321' style (4 decimals, enough for single cheap calls), or 'unknown'."""
    return "unknown" if cost is None else f"${cost:.4f}"
