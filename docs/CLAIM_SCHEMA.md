# Claim YAML schema

The canonical specification for a claim file. Use this when writing your
own claims; validate against it with `ft-validate claim.yaml`.

---

## Minimal example

```yaml
claim:
  target_symbol: "PCSK9"
  uniprot_id: "Q8NBP7"
  indication: "Hypercholesterolemia"
  mechanism: |
    PCSK9 loss-of-function variants lower LDL-C in humans;
    antibody inhibition mimics the LoF phenotype.
  claim_type: validated_mechanism
```

The minimal claim has four required fields and nothing else. The audit
will run, but rules with no data abstain or use defaults.

---

## Full example with fixture

A fixture supplies data the live adapters can't fetch (or that you want
to override). See `claims/tyk2_psoriasis.yaml` for a real-world example.

```yaml
claim:
  target_symbol: "TYK2"
  uniprot_id: "P29597"
  indication: "Psoriasis"
  mechanism: |
    TYK2 mediates IL-12/IL-23/Type-I IFN signaling.
  claim_type: validated_mechanism
  source_url: "https://www.fda.gov/..."
  citation_doi: "10.xxxx/xxxxx"

fixture:
  orthology:
    sources_agreeing: 4
    sources_total: 4
  genetics:
    gwas_hits: 23
    mendelian_evidence: true
  expression:
    target_tissue_expressed: true
  reproducibility:
    independent_replications: 12
    retracted: false
  selectivity:
    selectivity_data_exists: true
    selectivity_index_log: 2.1
```

---

## Field reference

### `claim` (required)

| Field | Type | Required | Description |
|---|---|---|---|
| `target_symbol` | string | yes | Gene/protein symbol (e.g. "TYK2"). 1-64 chars. |
| `indication` | string | yes | Therapeutic indication. 1-256 chars. |
| `mechanism` | string (multi-line) | yes | Mechanism of action narrative. 1-2048 chars. |
| `claim_type` | enum | yes | One of `validated_mechanism`, `novel_target`, `chemistry_series`, `extraordinary_claim` |
| `uniprot_id` | string | no | UniProt accession. 6-12 alphanumeric chars. |
| `ensembl_id` | string | no | Ensembl gene ID. |
| `chembl_id` | string | no | ChEMBL target ID. |
| `source_url` | string | no | URL of primary claim source (paper, FDA notice, etc.) |
| `citation_doi` | string | no | DOI of primary claim source. |
| `citation_pmid` | string | no | PubMed ID. |
| `author` | string | no | Who authored this claim YAML (not the science). |
| `date` | string | no | YYYY-MM-DD when claim was filed. |
| `notes` | string | no | Free-text notes. |

Any field starting with `x_` is reserved for user extensions and won't
be flagged by the validator.

### `fixture` (optional)

A dict of section names → data dicts. See
[`docs/ADAPTER_PROTOCOL.md`](ADAPTER_PROTOCOL.md) for the canonical
field list per section.

Valid section names:
- `orthology` — R1 data
- `chemistry` — R2 and R6 data
- `genetics` — R3 data
- `expression` — R4 data
- `reproducibility` — R5 data
- `selectivity` — R7 data
- `structure` — informational

### `expected_verdict` and `id` (sentinel files only)

Sentinel files have two extra top-level fields used by the calibration
suite. Regular claim files don't need these.

```yaml
id: PCSK9_LDL
expected_verdict: SURVIVED
claim:
  ...
fixture:
  ...
```

---

## Honesty discipline

When you hand-set a fixture value, it's your job to make sure it
reflects published evidence. Conventions:

1. **Cite the source in a YAML comment** above the value:
   ```yaml
   genetics:
     # POETYK PSO-1 Phase 3: deucravacitinib > placebo.
     # PMID 35290654, doi:10.1056/NEJMoa2104953
     gwas_hits: 23
     mendelian_evidence: true
   ```
2. **Don't massage values to engineer a verdict.** The framework's
   value depends on the inputs being honest.
3. **If you don't know, leave it out.** A missing field → rule
   abstains. A made-up field → rule passes/fails on bad evidence.
4. **Use the `notes` field for caveats** that don't fit into a
   structured field.

---

## Common mistakes

- **Forgetting `claim_type`.** Defaults don't exist. Pick one.
- **Wrong `claim_type` for the audit.** A novel-target claim audited
  as `validated_mechanism` will skip R1 (orthology) — usually wrong.
- **UniProt ID with non-alphanumeric chars.** Hyphen is OK (e.g.
  isoform IDs); spaces and slashes are not.
- **Putting fixture data under the wrong section.** `gwas_hits` goes
  under `genetics`, not `chemistry`.

The validator catches all of these:

```bash
$ ft-validate claims/your_claim.yaml
FAIL: claims/your_claim.yaml
  - claim.claim_type: 'foo' not in ['chemistry_series', 'extraordinary_claim', 'novel_target', 'validated_mechanism']
  - fixture has unknown section(s) ['molecular']. Valid: ['chemistry', 'expression', 'genetics', 'orthology', 'reproducibility', 'selectivity', 'structure']
```

---

## Versioning

This schema is at v1.0 (defined alongside ruleset v1.2.0 and tool v1.4.0).
Schema changes are minor unless they remove or rename required fields.
Backwards-compat: a v1.0 claim YAML will validate against any future v1.x
schema; only the tool's v2.0 will allow breaking changes.
