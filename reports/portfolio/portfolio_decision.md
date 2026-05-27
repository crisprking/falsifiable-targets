# Portfolio audit decision sheet

_Generated 2026-05-27T17:06:13.854053+00:00_

## Summary by risk tier

- **TIER_1_DROP**: 2
- **TIER_0_FIX_INPUT**: 4
- **TIER_2_INVESTIGATE**: 4
- **TIER_3_PROCEED**: 2

## Per-claim detail (highest risk first)

| Tier | Target | Verdict | Caveats | Fired | Action items |
|---|---|---|---|---|---|
| `TIER_1_DROP` | **OCT4** | FALSIFIED | 0 | `R5_replication` | R5_replication FALSIFIED @ public_data_lookup: PubMed / Retraction Watch lookup. Primary claim is retracted (2014). |
| `TIER_1_DROP` | **SAT_PATHOGEN_HOMOLOG_X** | FALSIFIED | 1 | `R6_chemistry_class_collapse` | R6_chemistry_class_collapse FALSIFIED @ public_data_lookup: Group ChEMBL hits by Pfam class. 91% are HDAC4 inhibitors - the chemistry support for SAT_PATHOGEN_HOMOLOG_X is class-collapsed phantom evid |
| `TIER_0_FIX_INPUT` | **(unparsed)** | ERROR | 0 | `` | ERROR: invalid claim_type 'validatedmechanism' in /kaggle/working/falsifiable-targets-main/tests/hardtest_claims/s6a_malformed.yaml. Valid: ['validated_mechanism', 'novel_target', 'chemistry_series',  |
| `TIER_0_FIX_INPUT` | **(unparsed)** | ERROR | 0 | `` | ERROR: /kaggle/working/falsifiable-targets-main/tests/hardtest_claims/s6b_missing_claim_type.yaml missing required claim fields: ['claim_type'] |
| `TIER_0_FIX_INPUT` | **(unparsed)** | ERROR | 0 | `` | ERROR: claim file is empty: /kaggle/working/falsifiable-targets-main/tests/hardtest_claims/s6c_empty.yaml |
| `TIER_0_FIX_INPUT` | **(unparsed)** | ERROR | 0 | `` | ERROR: YAML parse error in /kaggle/working/falsifiable-targets-main/tests/hardtest_claims/s6d_yaml_parse_error.yaml: |
| `TIER_2_INVESTIGATE` | **BACE1** | FALSIFIED_WITH_CAVEATS | 1 | `-` | R3_genetics_support: clinical outcome contested; framework does not predict trial outcomes - mechanism validity does not entail clinical efficacy |
| `TIER_2_INVESTIGATE` | **FUNGAL_KINASE_FKx1** | FALSIFIED_WITH_CAVEATS | 2 | `-` | R1_orthology: orthology sources disagree (1/4, strict majority requires 3); selectivity claim weakened until resolved structurally; R7_selectivity_counterscreen: no selectivity data for a novel-target |
| `TIER_2_INVESTIGATE` | **KRAS** | FALSIFIED_WITH_CAVEATS | 1 | `-` | R3_genetics_support: clinical outcome contested; framework does not predict trial outcomes - mechanism validity does not entail clinical efficacy |
| `TIER_2_INVESTIGATE` | **ZIKV_NS5** | FALSIFIED_WITH_CAVEATS | 1 | `-` | R1_orthology: orthology sources disagree (2/4, strict majority requires 3); selectivity claim weakened until resolved structurally |
| `TIER_3_PROCEED` | **PCSK9** | SURVIVED | 0 | `-` | data gap: R2_chemistry_support abstained (no input) |
| `TIER_3_PROCEED` | **CCR5** | SURVIVED | 0 | `-` | all rules cleared |
