# Audit limitations under v1.3.x

**Read this before drawing conclusions from any audit report produced
by this framework.** Updates [`AUDIT_LIMITATIONS_v1_2.md`](AUDIT_LIMITATIONS_v1_2.md)
to reflect changes in v1.3.0. Where this document differs from the v1.2.x
version, this one applies.

Document covers ruleset **v1.2.0** (SHA
`35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221`),
active in releases v1.3.0 and v1.3.1.

## 1. What changed between v1.2.x and v1.3.x

**R6_chemistry_class_collapse was upgraded** (rule version 1.0.0 → 1.2.0).
Two changes:

1. **Scope expanded.** R6 now applies to `validated_mechanism` claims
   in addition to `novel_target` and `chemistry_series`. The
   class-collapse risk for kinase-family validated mechanisms (canonical
   example: pan-JAK overlap for TYK2) is no longer silently out of scope.

2. **New heuristic substantive-caveat path.** When the chemistry section
   contains `chembl_paralog_compound_counts` (populated by the live
   ChEMBLAdapter when a paralog map is configured), R6 checks whether any
   paralog has ≥ 2× the primary's distinct-compound count. If so, R6
   emits a substantive caveat. The existing exact-fraction falsification
   path (≥ 0.80) is preserved.

The v1.2.x R6 scope gap named in §2 of `AUDIT_LIMITATIONS_v1_2.md`
is closed.

## 2. Per-rule data sources under v1.3.x

| Rule | Applies to | Live data populated? | Fixture-supplied fields |
|---|---|---|---|
| R1_orthology | `novel_target`, `extraordinary_claim` | Lower-bound only (UniProt cross-references count, not interpretation) | `sources_agreeing`, `sources_total`, `ortholog_count_human_pathogen` |
| R2_chemistry_support | `chemistry_series`, `validated_mechanism` | **Yes** — `chembl_distinct_compounds` from live ChEMBL | (none required) |
| R3_genetics_support | `validated_mechanism` | No live adapter yet | `gwas_hits`, `mendelian_evidence`, `somatic_driver_evidence`, `loss_of_function_phenotype`, `clinical_outcome_contested` |
| R4_expression | `validated_mechanism`, `novel_target` | Partial — `target_tissue_expressed` set when UniProt has any tissue annotation | `target_tissue_expressed` (fallback) |
| R5_replication | all | No live adapter yet | `retracted`, `rebuttals_count`, `independent_replications` |
| R6_chemistry_class_collapse | `novel_target`, `chemistry_series`, **`validated_mechanism`** (new in v1.2.0) | **Yes** — `chembl_paralog_compound_counts` from live ChEMBL when paralog map has entry | `chembl_pfam_class_collapse_fraction` (for exact path), `chembl_pfam_class_collapse_target_symbol` |
| R7_selectivity_counterscreen | `chemistry_series`, `novel_target`, `validated_mechanism` | No live adapter yet | `selectivity_data`, `off_targets_in_indication_relevant_tissue`, `selectivity_index_log` |

**Summary across the three completed audits:**

| Audit | Live-data rules | Fixture-data rules | Not-applicable rules |
|---|---|---|---|
| Ipi1 (novel target) | R5 | R1, R7 | R2, R3, R4, R6 |
| TYK2 v1.1.0 (validated mechanism) | R2, R4 | R3, R5, R7 | R1, R6 |
| **TYK2 v1.3.0** (validated mechanism, **R6 now in scope**) | R2, R4, **R6** | R3, R5, R7 | R1 |

The v1.3.0 TYK2 audit is the framework's *deepest* live-data audit to
date: 3 of 6 applicable rules driven by ChEMBL+UniProt fetches, 3
driven by fixture-encoded literature values, R1 not applicable.

## 3. The R6 heuristic — what it does and does not measure

R6 v1.2.0's new substantive-caveat path measures *pool-size overshadow*:
does a paralog have meaningfully more ChEMBL chemistry than the primary?
This is NOT the same as *compound-level overlap* (do the same compounds
hit primary and paralog with comparable potency?).

For TYK2 specifically (see [`AUDIT_TYK2_v1_3.md`](AUDIT_TYK2_v1_3.md)):
the JAK family pool sizes are comparable to TYK2's (max ratio 1.55×),
so R6 does not flag overshadow. The pan-JAK *binding-mode* overlap
that deucravacitinib's JH2-domain mechanism resolves is a different
question that R6 v1.2.0 does not audit. **A SURVIVED verdict at R6
v1.2.0 should not be read as "the framework confirmed selectivity"; it
should be read as "the framework checked one specific axis (pool size)
and found no overshadow."**

Compound-level overlap is a v1.4+ feature. Its specification is in
[`AUDIT_TYK2_v1_3.md`](AUDIT_TYK2_v1_3.md) §"Honest v1.4 questions".

## 4. Remaining v1.3.x gaps

After the v1.3.0 R6 upgrade, the remaining gaps from v1.2.x are:

- **R3 has no live adapter.** Open Targets / GWAS Catalog integration
  is queued for v1.4+. The genetics verdict still depends on
  fixture-encoded literature values for validated_mechanism claims.
- **R5 has no live adapter.** Europe PMC + Retraction Watch integration
  is queued for v1.4+. A paper retracted between fixture encoding
  and audit date would not be caught.
- **R7 has no live adapter.** Live kinome-selectivity computation is
  queued for v1.4+ (likely shares infrastructure with the R6 compound-
  level overlap computation).
- **R4 is structurally weak.** Detecting "UniProt has any tissue
  annotation" is not the same as "target expressed in the indication-
  relevant tissue specifically." GTEx-by-indication is v1.5+.
- **ChEMBL adapter `limit=1000` per fetch.** All chemistry counts are
  bounded above by this. For most rules this is sufficient (R2 is a
  "> 0" check, R6 ratios are stable under similar capping). Pagination
  beyond 1000 is v1.4 if any rule ever needs precise enumeration.

## 5. Reproducibility guarantees

Unchanged from v1.2.x. SHA-locked audits produce byte-identical
verdicts given identical inputs. Live-mode runs on different days may
return slightly different absolute counts as ChEMBL adds records, but
the cached `.ae_cache/` from a given audit run freezes those numbers
for re-runs in `AE_OFFLINE=1` mode.

## 6. v1.4.0 plan

Roadmap in priority order:

1. **Compound-level overlap computation for R6** — the structural
   complement to v1.3.0's pool-size heuristic. Per primary-target
   compound, count how many hit any paralog at comparable potency.
2. **Live Open Targets adapter** for R3 (genetics).
3. **Live Europe PMC + Retraction Watch adapter** for R5 (replication).
4. **Cache archival** in `run_audit.py` (`--archive` flag bundles
   `.ae_cache/` snapshot with the JSON report for exact-reproducibility
   shipping).
5. **More external audits** to calibrate R6's 2.0× threshold against
   the distribution of real-world paralog ratios.

v1.5+: Live kinome-selectivity computation for R7; GTEx-by-indication
for R4; pagination in ChEMBL adapter for absolute compound enumeration.
