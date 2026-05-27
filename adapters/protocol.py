"""
Adapter protocol.

An Adapter is anything that knows how to populate one or more of the
seven fixture sections (orthology, chemistry, genetics, expression,
reproducibility, selectivity, structure) for a target claim.

This file is the *contract*. Adapters do not need to inherit from
anything - they just need to implement `get(section, claim) -> dict`.

If you write a new adapter for, say, Open Targets or AlphaFold or a
proprietary cell-painting dataset, conforming to this protocol means
it drops straight into CompositeAdapter alongside the built-in ones.

Example minimal adapter:

    class MyOpenTargetsAdapter:
        def __init__(self, api_key=None):
            self._api_key = api_key

        def get(self, section, claim):
            if section != "genetics":
                return {}
            # ... fetch from Open Targets ...
            return {
                "gwas_hits": 12,
                "mendelian_evidence": True,
                "open_targets_score": 0.87,
            }

Plugging it in:

    from adapters import CompositeAdapter, FixtureAdapter, UniProtAdapter

    composite = CompositeAdapter([
        MyOpenTargetsAdapter(api_key=os.environ["OT_API_KEY"]),  # highest priority
        UniProtAdapter(),
        FixtureAdapter(fixture),                                  # lowest priority
    ])

CompositeAdapter walks adapters in priority order, taking the first
non-null value for each KEY (not each section). So a partial response
from your adapter gets augmented by the fixture for missing fields.
"""
from __future__ import annotations

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class Adapter(Protocol):
    """The single-method protocol all adapters must implement.

    Args:
        section: One of "orthology", "chemistry", "genetics",
                 "expression", "reproducibility", "selectivity",
                 "structure". See VALID_FIXTURE_SECTIONS in
                 validate_claim.py for the authoritative list.
        claim:   A TargetClaim. The adapter may use claim.uniprot_id,
                 claim.target_symbol, claim.indication, etc.

    Returns:
        A dict of {field_name: value}. Empty dict means "I have no
        data for this section." Missing keys mean "I have no data
        for these specific fields"; CompositeAdapter will fall through
        to lower-priority adapters for the missing keys.

    Conventions:
        - NEVER raise on missing data. Return {}.
        - It IS OK to raise on a programming error (bad section name,
          malformed claim) - that's how the engine surfaces bugs.
        - Adapters that hit the network should cache to .ae_cache/
          (or whatever AE_CACHE_DIR is set to) so audits are
          deterministically reproducible.
        - Adapters should respect AE_OFFLINE=1 by serving from cache
          only and returning {} on cache miss.
    """

    def get(self, section: str, claim: Any) -> dict[str, Any]: ...


# Documented section names. Adapters can serve any subset.
SECTIONS = (
    "orthology",
    "chemistry",
    "genetics",
    "expression",
    "reproducibility",
    "selectivity",
    "structure",
)
