# Open Targets adapter

`adapters/opentargets.py` provides `OpenTargetsAdapter`, a live-evidence
adapter for the **`genetics`** section. It conforms to the adapter protocol
(`adapters/protocol.py`): a single `get(section, claim) -> dict`. It drops
into `CompositeAdapter` above `FixtureAdapter`, so live values win and the
fixture backfills any field the adapter can't resolve.

## Evidence, not verdict

The adapter **supplies evidence only**. Direction-of-effect *concordance with
the claim's mechanism* is a rule decision and stays in R3's SHA-locked logic,
not here. The adapter returns the gene's genetic direction, its provenance,
and a quality tier; R3 compares that to `claim.mechanism_direction` and emits
the verdict. Keeping the split here is what lets the audit logic stay hashable
and the adapter stay swappable — the ruleset SHA does not change when you swap
or upgrade this adapter.

## What `get("genetics", claim)` returns

Backward-compatible keys R3 already reads (unchanged):

- `gwas_hits` — count of `gwas_credible_sets` rows
- `mendelian_evidence` — bool, ClinVar (`eva`) or `gene_burden` present
- `somatic_driver_evidence` — bool
- `open_targets_score` — overall target–disease association score (or `None`)

Direction-of-effect evidence (consumed by the R3 v1.5.0 reads):

- `doe_gene_direction` — `GoF_increases_risk` | `LoF_increases_risk` |
  `conflicting` | `unknown`
- `doe_tier` — `strong` | `moderate` | `weak` | `absent` (or `conflicting`)
- `doe_sources` — which datasources contributed (provenance)
- `doe_n_assessments` — number of evidence rows considered
- `_provenance` — `{ensembl_id, efo_id, endpoint}`

Direction is **derived, not read off a single field** (Open Targets folded
Genetics into the unified Platform; there is no clean `directionOfEffect`
scalar in the GraphQL API). Cleanliness order is `gene_burden` (strong) >
ClinVar/`eva` (moderate) > `gwas_credible_sets` + molQTL coloc (weak). The
first source with an unambiguous call sets the tier; a later source only
downgrades the result to `conflicting` if it disagrees. Disagreement is
**surfaced, never averaged**.

If identifiers can't be resolved, or there's no evidence, `get` returns `{}`
so the composite falls through to lower-priority adapters — it never raises on
missing data.

## Wiring it in

```python
from adapters import CompositeAdapter, FixtureAdapter, OpenTargetsAdapter

composite = CompositeAdapter([
    OpenTargetsAdapter(),   # highest priority: live genetics
    FixtureAdapter(fixture),  # lowest priority: backfill
])
```

It is intentionally **not** added to `default_composite()` — opt in
explicitly, so offline/deterministic runs stay the default.

## Caching and offline use

The adapter caches every GraphQL response to disk for deterministic replay:

- `AE_CACHE_DIR` — cache location (defaults to `.ae_cache/`).
- `AE_OFFLINE=1` — serve from cache only; return `{}` on a cache miss instead
  of hitting the network.

Network calls retry with exponential backoff on HTTP 429/502/503.

Requires the optional `live` dependency group: `pip install -e ".[live]"`
(installs `requests`).

## Known limitation (flagged in source)

The three per-datasource directional extractors — `_burden_direction`,
`_clinvar_direction`, `_gwas_coloc_direction` — are deliberately thin. The
exact GraphQL row fields are datasource-specific and **must be confirmed
against the live Open Targets schema** before the derived direction is relied
on in production. Until confirmed, each extractor returns `None` ("no
directional signal from this row"), so unconfirmed rows degrade safely to
`unknown` / `absent` rather than producing a wrong call.
