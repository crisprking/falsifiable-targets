# Audit limitations under v1.2.x

**Read this before drawing conclusions from any audit report produced by
this framework.** The purpose of this document is to state explicitly,
for every rule, what evidence the rule actually uses and where that
evidence comes from. A verdict is only as strong as its inputs; this
document is the audit trail for the inputs.

Document covers ruleset **v1.1.0** (SHA
`2f9aab7d0ebc209f62c16eb35be31bc5b65fa2eb09adc02bea5bff5176269b32`),
which is the active ruleset in releases v1.2.0 and v1.2.1.

## 1. Per-rule data sources

For each rule, three columns: which claim types invoke it, what fields
the rule consumes, and where those fields actually come from in the
current adapter stack.

| Rule | Applies to | Required inputs | Source under v1.2.x |
|---|---|---|---|
| R1_orthology | `novel_target`, `extraordinary_claim` | `sources_agreeing`, `sources_total` | **Fixture only.** Live UniProt reports a *lower bound* (`sources_agreeing_uniprot_lower_bound`) but does not yet determine whether the cross-referenced DBs actually agree on the orthology call. The fixture's `sources_agreeing` value is the load-bearing input. |
| R2_chemistry_support | `chemistry_series`, `validated_mechanism` | `chembl_distinct_compounds` | **Live ChEMBL.** Counts distinct `molecule_chembl_id` values across up to 1000 bioactivity records. (Pagination is a v1.3+ improvement; 1000 records is enough for the rule's "n < 1 → falsify" logic but not for precise enumeration.) |
| R3_genetics_support | `validated_mechanism` | `gwas_hits`, `mendelian_evidence`, `somatic_driver_evidence`, `loss_of_function_phenotype`, `clinical_outcome_contested` | **Fixture only.** No live GWAS/Open Targets adapter exists yet. The TYK2 audit's R3 verdict reflects hand-encoded values paraphrased from primary literature (cited in the claim YAML). Live Open Targets integration is queued for v1.3.0. |
| R4_expression | `validated_mechanism`, `novel_target` | `target_tissue_expressed` (bool), `uniprot_tissue_text` (informational) | **Mixed.** Live UniProt populates `target_tissue_expressed = True` whenever the entry has any tissue-specificity comment. **This is structurally weak**: it tests "does UniProt have any expression annotation," not "is the target expressed in the indication-relevant tissue." A more semantic check (GTEx, Human Protein Atlas by tissue) is a v1.3+ improvement. |
| R5_replication | all claim types | `retracted`, `rebuttals_count`, `independent_replications` | **Fixture only.** No live retraction-check adapter exists yet. Europe PMC + Retraction Watch integration is queued for v1.3.0. |
| R6_chemistry_class_collapse | `novel_target`, `chemistry_series` | `chembl_pfam_class_collapse_fraction`, `chembl_pfam_class_collapse_target_symbol` | **Fixture only when applicable.** Crucially: R6 does **not** apply to `validated_mechanism` claims under v1.1.0. The TYK2 audit's `R6_NOT_APPLICABLE` result reflects this scope restriction, not a clean pass on the class-collapse axis. See §2 below. |
| R7_selectivity_counterscreen | `chemistry_series`, `novel_target`, `validated_mechanism` | `selectivity_data` (bool), `off_targets_in_indication_relevant_tissue` (bool), `selectivity_index_log` (float) | **Fixture only.** No live kinome-selectivity adapter exists yet (the join across paralog activities is complex). Queued for v1.4+. |

**Summary across the two completed audits:**

| Audit | Live-data rules | Fixture-data rules | Not-applicable rules |
|---|---|---|---|
| Ipi1 (novel target) | R5 (replication, fixture-only, no retraction) | R1, R7 | R2, R3, R4, R6 |
| TYK2 (validated mechanism) | R2 (chemistry count), R4 (expression presence) | R3, R5, R7 | R1, R6 |

A reader of the TYK2 SURVIVED verdict should weight it accordingly:
**two rules genuinely audited against live public data, three audited
from hand-encoded literature values, two not applicable.** The verdict
is correct under v1.1.0 logic but the depth of independent evidence
is thinner than a 5-out-of-7 pass-rate alone implies.

## 2. R6 scope: the gap the TYK2 audit surfaced

R6 (chemistry_class_collapse) detects the "phantom evidence" failure
mode: the ChEMBL compounds claimed as evidence for a target are, on
inspection, mostly inhibitors of a different protein in the same Pfam
class. The canonical example encoded in the sentinel suite is the
SAT/HDAC4 case from the Madurella audit.

Under v1.1.0, R6 only applies to `novel_target` and `chemistry_series`
claims. The pan-JAK selectivity question for TYK2 - **the textbook
class-collapse risk for kinase-family validated mechanisms** - is
therefore not audited. A TYK2 ChEMBL ligand could overlap heavily with
JAK1/2/3 actives and R6 would not flag it because the claim type is
`validated_mechanism`.

This is a deliberate v1.1.0 scope choice (R6 was designed around novel-
target failure modes) that the TYK2 audit revealed to be too narrow for
real-world kinase chemistry. **Expanding R6 to validated_mechanism is
the headline change planned for v1.3.0.** Implementation will require:

1. R6.applies_to() expanded to include `validated_mechanism`
2. A live class-collapse computation in ChEMBLAdapter that fetches
   paralog activities (e.g. JAK1 CHEMBL2835, JAK2 CHEMBL2971, JAK3
   CHEMBL2148 for the TYK2 case) and computes per-compound overlap
3. A new sentinel that exercises validated_mechanism class-collapse
4. Re-run of the TYK2 audit with R6 active; verdict reported honestly
   whether SURVIVED or FALSIFIED_WITH_CAVEATS

Until v1.3.0 ships, the v1.2.x TYK2 audit's `SURVIVED` verdict should
be read as "passed all rules that applied; the most relevant kinase-
specific risk was out of scope."

## 3. Reproducibility guarantees

The SHA lock guarantees that two runs of the same audit under the same
ruleset version produce the same verdict **if the same inputs were
provided**. It does *not* guarantee that two live-mode runs on different
days produce identical outputs - ChEMBL adds records, UniProt updates
annotations, and the cached responses in `.ae_cache/` are timestamped
informally by directory mtime, not by the report itself.

For exact reproducibility of an external audit:

1. Pin the ruleset SHA (already in every report)
2. Pin the claim SHA (already in every report)
3. Archive the `.ae_cache/` directory from the audit run
4. Re-run with `AE_OFFLINE=1` to force cache-only mode

Cache archival is currently a manual step. v1.3+ will add a `--archive`
flag to `run_audit.py` that bundles the cache snapshot with the JSON
report.

## 4. What v1.2.x is suitable for

- **Internal calibration of in-flight target hypotheses** (the Ipi1 use
  case): apply the framework to your own pre-publication targets to
  surface gaps the headline narrative might be papering over.
- **Public auditing of well-validated mechanisms with rich ChEMBL data**
  (the TYK2 use case): the live UniProt and ChEMBL adapters carry
  enough weight to make R2 and R4 verdicts meaningful.
- **Comparative auditing across same-class claims**, where the audit's
  internal asymmetries (which rules fire substantively) carry as much
  information as the final verdict.

## 5. What v1.2.x is NOT suitable for

- **Auditing kinase-class validated mechanisms with confidence that the
  class-collapse axis was checked.** R6 doesn't apply; explicitly
  caveat any such audit until v1.3.0.
- **Audits where the genetics evidence is the load-bearing claim and
  primary literature has not been read into the fixture.** R3 has no
  live data source; an audit of a GWAS-driven claim will pass R3 on
  whatever the fixture asserts.
- **Audits of recently-retracted papers without manual retraction-check
  data in the fixture.** R5 has no live Retraction Watch lookup; a
  retracted-but-not-yet-encoded paper will not fail R5.
- **Anything where the audit's *output verdict* is the only artifact
  consulted.** Always read the per-rule input data in the JSON report
  and confirm the load-bearing rule actually ran on live data.

## 6. v1.3.0 plan

The above limitations define the v1.3.0 roadmap. In priority order:

1. **R6 scope expansion** + **live class-collapse adapter** for the
   kinase-class blind spot
2. **Live Open Targets adapter** for R3 (genetics)
3. **Live Europe PMC + Retraction Watch adapter** for R5 (replication)
4. **Cache archival** in `run_audit.py` for exact reproducibility
5. **TYK2 re-audit** under the expanded ruleset; verdict reported
   honestly whether the result changes

v1.4+ deferred: live kinome-selectivity computation for R7; semantic
tissue-expression lookup for R4 (GTEx by indication); pagination in
ChEMBL adapter beyond the 1000-record limit.
