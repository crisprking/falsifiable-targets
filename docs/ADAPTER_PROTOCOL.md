# Adapter protocol

How to write a new data adapter for falsifiable-targets.

---

## What an adapter does

An **adapter** answers questions about a target claim using a specific
data source. Built-in adapters:

| Adapter | Sections served | Source |
|---|---|---|
| `FixtureAdapter` | all 7 | hand-encoded YAML data in the claim file |
| `UniProtAdapter` | `structure`, `expression`, `orthology` | live UniProt REST API |
| `ChEMBLAdapter` | `chemistry` (and paralog counts) | live ChEMBL REST API |

Multiple adapters get composed via `CompositeAdapter`, which walks them
in priority order and merges per-key. This is how a partial live response
(UniProt returns `pdb_count` but not `alphafold_plddt`) gets augmented
from a fixture (which supplies `alphafold_plddt`).

---

## The interface

Just one method:

```python
def get(self, section: str, claim: TargetClaim) -> dict:
    ...
```

`section` is one of seven names. `claim` is the target claim object
(with `target_symbol`, `uniprot_id`, `indication`, `mechanism`,
`claim_type`). Return: a dict of fields the adapter knows.

That's it. No base class to inherit. The `Adapter` Protocol in
`adapters/protocol.py` is `runtime_checkable`, so `isinstance(my, Adapter)`
works to verify conformance.

---

## Section names and expected fields

This table is the canonical contract. If your adapter populates a field
not listed here, document it; if your adapter populates one that *is*
listed here, follow the type convention.

### `orthology`

| Field | Type | Used by | Meaning |
|---|---|---|---|
| `sources_agreeing` | int | R1 | How many ortholog DBs return a human ortholog (strict) |
| `sources_total` | int | R1 | How many ortholog DBs were consulted |
| `sources_agreeing_uniprot_lower_bound` | int | R1 (fallback) | Live UniProt cross-references — lower bound on agreement |
| `sources_total_uniprot_lower_bound` | int | R1 (fallback) | Live UniProt cross-references — total consulted |

### `chemistry`

| Field | Type | Used by | Meaning |
|---|---|---|---|
| `chembl_distinct_compounds` | int | R2, R6 | Distinct compounds with activity vs this target |
| `chembl_target_id` | str | R6 | ChEMBL target ID (e.g., "CHEMBL2929") |
| `chembl_pfam_class_collapse_fraction` | float \[0,1\] | R6 | Fraction of binders that collapse onto a paralog class |
| `chembl_paralog_compound_counts` | dict[str, int] | R6 (heuristic) | Per-paralog compound count (live ChEMBL adapter) |
| `max_phase` | int \[0,4\] | R2 | Highest clinical phase reached |

### `genetics`

| Field | Type | Used by | Meaning |
|---|---|---|---|
| `gwas_hits` | int | R3 | Number of genome-wide-significant GWAS associations |
| `mendelian_evidence` | bool | R3 | Loss-of-function in humans alters indication phenotype |
| `somatic_driver` | bool | R3 | Recurrent somatic driver in indication (oncology) |
| `open_targets_score` | float \[0,1\] | R3 | Open Targets composite (if integrated) |

### `expression`

| Field | Type | Used by | Meaning |
|---|---|---|---|
| `target_tissue_expressed` | bool | R4 | Detectably expressed in indication-relevant tissue |
| `uniprot_tissue_text` | str | R4 | Free-text tissue annotation from UniProt |
| `hpa_tissue_tpm_max` | float | R4 | Max TPM across HPA tissue panel (if HPA adapter) |

### `reproducibility`

| Field | Type | Used by | Meaning |
|---|---|---|---|
| `independent_replications` | int | R5 | Number of independent labs reproducing the result |
| `retracted` | bool | R5 | Primary claim source is retracted |
| `pubpeer_serious_concerns` | bool | R5 | Substantive PubPeer concerns on primary source |

### `selectivity`

| Field | Type | Used by | Meaning |
|---|---|---|---|
| `selectivity_data` | bool | R7 | Any selectivity data exists (vs paralogs / off-targets) |
| `selectivity_index_log` | float | R7 | log10(IC50_offtarget / IC50_target) |
| `off_targets_in_indication_relevant_tissue` | bool | R7 | Off-target hits express in indication tissue |
| `closest_paralog_offtarget_log` | float | R7 | Same as selectivity_index_log but for nearest paralog |

### `structure`

| Field | Type | Used by | Meaning |
|---|---|---|---|
| `structure_resolved` | bool | (informational) | Any structural data available |
| `pdb_count` | int | (informational) | Number of PDB structures |
| `alphafold_plddt` | float | (informational) | AlphaFold confidence |

---

## Caching

Live adapters MUST cache responses to `.ae_cache/` (or whatever
`AE_CACHE_DIR` is set to) so audits are byte-deterministic across runs.

Pattern in the built-in adapters:

```python
def _fetch_json(self, url):
    cache_path = CACHE_DIR / quote(url, safe='') + ".json"
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    if os.environ.get("AE_OFFLINE"):
        return None  # cache miss in offline mode
    response = urllib.request.urlopen(url, timeout=15)
    data = json.loads(response.read())
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps(data))
    return data
```

Two important properties this gives you:

1. **Reproducibility.** Once a target has been audited, its
   `.ae_cache/` contents are sufficient to re-run the audit
   offline indefinitely.
2. **CI compatibility.** Setting `AE_OFFLINE=1` in CI prevents
   network calls; tests must use mock data.

---

## Testing your adapter

Hermetic tests only (no live network in CI). Pattern:

```python
def test_my_adapter_returns_genetics():
    mock_data = {"my-api://endpoint/PCSK9": {"gwas": 47}}
    adapter = MyAdapter(mock_data=mock_data)
    result = adapter.get("genetics", pcsk9_claim)
    assert result["gwas_hits"] == 47
```

The `mock_data` constructor parameter is the convention all live
adapters expose. When set, the adapter serves from it instead of
hitting the network. See `tests/test_adapters.py` for canonical
examples.

For live integration tests (manual, not in CI):

```python
import pytest

@pytest.mark.live
def test_my_adapter_live():
    adapter = MyAdapter()  # no mock_data => live mode
    result = adapter.get("genetics", pcsk9_claim)
    assert result["gwas_hits"] > 0
```

Run with `pytest -m live`. CI skips these via `pytest -m "not live"`.

---

## Worked example: a (hypothetical) Open Targets adapter

```python
# adapters/open_targets.py
import json
import os
import urllib.parse
import urllib.request
from pathlib import Path


CACHE_DIR = Path(os.environ.get("AE_CACHE_DIR", ".ae_cache"))


class OpenTargetsAdapter:
    """Live Open Targets adapter. Serves the 'genetics' section."""

    BASE_URL = "https://api.platform.opentargets.org/api/v4/graphql"

    def __init__(self, mock_data=None):
        self._mock = mock_data

    def get(self, section, claim):
        if section != "genetics":
            return {}
        if claim.uniprot_id is None:
            return {}

        target_data = self._fetch_target(claim.uniprot_id)
        if not target_data:
            return {}

        return {
            "open_targets_score": target_data.get("associatedDiseaseScore"),
            "gwas_hits": target_data.get("gwasHits", 0),
            "mendelian_evidence": target_data.get("hasMendelianEvidence", False),
        }

    def _fetch_target(self, uniprot_id):
        if self._mock is not None:
            return self._mock.get(uniprot_id)
        # ... live fetch with cache ...
```

That's a complete adapter. Plug it into the composite:

```python
from adapters import CompositeAdapter, FixtureAdapter, UniProtAdapter, ChEMBLAdapter
from adapters.open_targets import OpenTargetsAdapter

composite = CompositeAdapter([
    OpenTargetsAdapter(),  # highest priority for genetics
    UniProtAdapter(),
    ChEMBLAdapter(),
    FixtureAdapter(claim_fixture),
])
```

R3 will now see live Open Targets data for `gwas_hits` and friends.

---

## Backwards compatibility

When you add fields to an existing section, **never remove** existing
fields. R3 might still need `gwas_hits` even after you add
`open_targets_score`. Audits run under older claim YAMLs need their
existing field semantics preserved.
