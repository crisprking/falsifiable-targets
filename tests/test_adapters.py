"""Adapter tests - hermetic. No live network. Mock payloads cover the
key UniProt and ChEMBL response shapes the adapters parse."""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from adapters import (
    ChEMBLAdapter,
    CompositeAdapter,
    FixtureAdapter,
    UniProtAdapter,
    default_composite,
)
from smoke_test import ClaimType, TargetClaim

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def pcsk9_claim():
    return TargetClaim(
        target_symbol="PCSK9",
        indication="hypercholesterolemia",
        mechanism="LoF reduces LDL-C",
        claim_type=ClaimType.VALIDATED_MECHANISM,
        uniprot_id="Q8NBP7",
    )


@pytest.fixture
def ipi1_claim():
    """Novel target, no UniProt ID. Adapters should return {} for everything."""
    return TargetClaim(
        target_symbol="Ipi1 (KXX81897.1)",
        indication="Madurella mycetoma",
        mechanism="Pre-60S ribosome biogenesis",
        claim_type=ClaimType.NOVEL_TARGET,
        uniprot_id=None,
    )


# Mock UniProt response: PCSK9 with PDB structures and tissue expression
PCSK9_UNIPROT_MOCK = {
    "Q8NBP7": {
        "primaryAccession": "Q8NBP7",
        "uniProtKBCrossReferences": [
            {"database": "PDB", "id": "2P4E"},
            {"database": "PDB", "id": "2PMW"},
            {"database": "PDB", "id": "3BPS"},
            {"database": "AlphaFoldDB", "id": "AF-Q8NBP7-F1"},
            {"database": "OrthoDB", "id": "5479at2759"},
            {"database": "eggNOG", "id": "ENOG502QPMK"},
            {"database": "InParanoid", "id": "Q8NBP7"},
        ],
        "comments": [
            {
                "commentType": "TISSUE SPECIFICITY",
                "texts": [{"value": "Expressed in liver and small intestine."}],
            }
        ],
    }
}

# Mock ChEMBL responses
PCSK9_CHEMBL_MOCK = {
    "target_by_uniprot:Q8NBP7": {
        "targets": [{"target_chembl_id": "CHEMBL3243"}]
    },
    "activities:CHEMBL3243": {
        "activities": [
            {"molecule_chembl_id": f"CHEMBL_{i}", "max_phase": (4 if i < 3 else 2)}
            for i in range(100)
        ]
    },
}


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestUniProtAdapter:
    def test_returns_empty_when_no_uniprot_id(self, ipi1_claim):
        adapter = UniProtAdapter(mock_data=PCSK9_UNIPROT_MOCK)
        assert adapter.get("structure", ipi1_claim) == {}
        assert adapter.get("orthology", ipi1_claim) == {}
        assert adapter.get("expression", ipi1_claim) == {}

    def test_structure_section_extracts_pdb_count(self, pcsk9_claim):
        adapter = UniProtAdapter(mock_data=PCSK9_UNIPROT_MOCK)
        out = adapter.get("structure", pcsk9_claim)
        assert out["pdb_count"] == 3
        assert out["structure_resolved"] is True

    def test_expression_section_extracts_tissue(self, pcsk9_claim):
        adapter = UniProtAdapter(mock_data=PCSK9_UNIPROT_MOCK)
        out = adapter.get("expression", pcsk9_claim)
        assert out["target_tissue_expressed"] is True
        assert "liver" in out["uniprot_tissue_text"].lower()

    def test_unknown_section_returns_empty(self, pcsk9_claim):
        adapter = UniProtAdapter(mock_data=PCSK9_UNIPROT_MOCK)
        assert adapter.get("genetics", pcsk9_claim) == {}

    def test_mock_miss_returns_empty(self, pcsk9_claim):
        adapter = UniProtAdapter(mock_data={})
        assert adapter.get("structure", pcsk9_claim) == {}


class TestChEMBLAdapter:
    def test_returns_empty_when_no_uniprot_id(self, ipi1_claim):
        adapter = ChEMBLAdapter(mock_data=PCSK9_CHEMBL_MOCK)
        assert adapter.get("chemistry", ipi1_claim) == {}

    def test_chemistry_section_counts_distinct_compounds(self, pcsk9_claim):
        adapter = ChEMBLAdapter(mock_data=PCSK9_CHEMBL_MOCK)
        out = adapter.get("chemistry", pcsk9_claim)
        assert out["chembl_distinct_compounds"] == 100
        assert out["max_phase"] == 4
        assert out["chembl_target_id"] == "CHEMBL3243"

    def test_non_chemistry_section_returns_empty(self, pcsk9_claim):
        adapter = ChEMBLAdapter(mock_data=PCSK9_CHEMBL_MOCK)
        assert adapter.get("structure", pcsk9_claim) == {}

    def test_mock_miss_returns_empty_so_composite_falls_through(self, pcsk9_claim):
        """When the mock has no entry for the target, simulates a failed
        live fetch. Adapter returns {} so the composite can fall back to
        fixture data, rather than asserting a misleading zero-count."""
        adapter = ChEMBLAdapter(mock_data={})
        out = adapter.get("chemistry", pcsk9_claim)
        assert out == {}, (
            f"Expected empty dict on mock miss (fall-through behavior), got {out}"
        )

    def test_api_affirmative_zero_targets_returns_zero_compounds(self, pcsk9_claim):
        """When the API answers but reports zero targets (e.g. truly novel
        target with no ChEMBL entry), report 0 compounds - this is real data."""
        adapter = ChEMBLAdapter(mock_data={
            "target_by_uniprot:Q8NBP7": {"targets": []}
        })
        out = adapter.get("chemistry", pcsk9_claim)
        assert out == {"chembl_distinct_compounds": 0}


class TestCompositeAdapter:
    def test_higher_priority_wins_per_key(self, pcsk9_claim):
        a1 = FixtureAdapter({"structure": {"pdb_count": 99}})
        a2 = FixtureAdapter({"structure": {"pdb_count": 1, "alphafold_plddt": 95.0}})
        composite = CompositeAdapter([a1, a2])
        out = composite.get("structure", pcsk9_claim)
        # a1 wins pdb_count, a2 fills alphafold_plddt
        assert out == {"pdb_count": 99, "alphafold_plddt": 95.0}

    def test_empty_higher_priority_falls_through(self, pcsk9_claim):
        a1 = FixtureAdapter({})
        a2 = FixtureAdapter({"chemistry": {"chembl_distinct_compounds": 42}})
        composite = CompositeAdapter([a1, a2])
        assert composite.get("chemistry", pcsk9_claim) == {"chembl_distinct_compounds": 42}

    def test_all_empty_returns_empty(self, pcsk9_claim):
        composite = CompositeAdapter([FixtureAdapter({}), FixtureAdapter({})])
        assert composite.get("chemistry", pcsk9_claim) == {}


class TestDefaultComposite:
    def test_offline_returns_fixture_only(self, pcsk9_claim):
        fixture = {"structure": {"alphafold_plddt": 93.2}}
        adapter = default_composite(fixture, use_live=False)
        assert isinstance(adapter, FixtureAdapter)
        assert adapter.get("structure", pcsk9_claim) == {"alphafold_plddt": 93.2}

    def test_live_with_mocks_layers_correctly(self, pcsk9_claim):
        """The full stack: UniProt mock -> ChEMBL mock -> fixture."""
        fixture = {
            "structure": {"alphafold_plddt": 93.2, "fpocket_score": 0.749},
            "genetics": {"gwas_hits": 18},
        }
        adapter = default_composite(
            fixture, use_live=True,
            uniprot_mock=PCSK9_UNIPROT_MOCK,
            chembl_mock=PCSK9_CHEMBL_MOCK,
        )
        # Structure: UniProt provides pdb_count, fixture provides plddt/fpocket
        struct = adapter.get("structure", pcsk9_claim)
        assert struct["pdb_count"] == 3
        assert struct["structure_resolved"] is True
        assert struct["alphafold_plddt"] == 93.2
        assert struct["fpocket_score"] == 0.749
        # Chemistry: ChEMBL provides the count, fixture is silent
        chem = adapter.get("chemistry", pcsk9_claim)
        assert chem["chembl_distinct_compounds"] == 100
        # Genetics: fixture only
        gen = adapter.get("genetics", pcsk9_claim)
        assert gen["gwas_hits"] == 18


class TestEngineWithLiveAdapter:
    """End-to-end: feed the engine a composite adapter and confirm the
    rules see the merged data correctly. PCSK9 sentinel should still survive."""

    def test_pcsk9_survives_with_live_mock_adapter(self, pcsk9_claim):
        from smoke_test import RULES, RuleStatus, aggregate
        fixture = {
            "orthology": {"sources_agreeing": 4, "sources_total": 4},
            "genetics": {"gwas_hits": 18, "loss_of_function_phenotype": "lifelong low LDL", "mendelian_evidence": True},
            "reproducibility": {"retracted": False, "rebuttals_count": 0, "independent_replications": 47},
            "selectivity": {"selectivity_data": True, "off_targets_in_indication_relevant_tissue": False},
        }
        adapter = default_composite(
            fixture, use_live=True,
            uniprot_mock=PCSK9_UNIPROT_MOCK,
            chembl_mock=PCSK9_CHEMBL_MOCK,
        )
        results = []
        for rule in RULES:
            if not rule.applies_to(pcsk9_claim):
                results.append(type("R", (), {"status": RuleStatus.NOT_APPLICABLE, "caveats": [], "rule_id": rule.rule_id, "falsification_tier": None})())
                continue
            results.append(rule.evaluate(pcsk9_claim, adapter))
        from smoke_test import Verdict
        verdict, _, _, _, _ = aggregate(results)
        assert verdict == Verdict.SURVIVED, (
            f"PCSK9 expected SURVIVED with live mock adapter, got {verdict.value}. "
            f"Rules: {[(r.rule_id, r.status.value) for r in results]}"
        )


class TestParalogMap:
    """v1.2.0: ChEMBLAdapter fetches paralog compound counts when a
    paralog_map is supplied. Hermetic - uses mock data."""

    def test_paralog_counts_populated_when_map_present(self, pcsk9_claim):
        # Mock: primary Q8NBP7 has 100 compounds, paralog has 50
        mock = {
            "target_by_uniprot:Q8NBP7": {"targets": [{"target_chembl_id": "CHEMBL3243"}]},
            "activities:CHEMBL3243": {"activities": [
                {"molecule_chembl_id": f"CHEMBL_PRIMARY_{i}"} for i in range(100)
            ]},
            "target_by_uniprot:P11111": {"targets": [{"target_chembl_id": "CHEMBL_P11111"}]},
            "activities:CHEMBL_P11111": {"activities": [
                {"molecule_chembl_id": f"CHEMBL_PARALOG_{i}"} for i in range(50)
            ]},
        }
        paralog_map = {
            "Q8NBP7": {
                "primary_symbol": "PCSK9",
                "paralogs": [{"symbol": "TestParalog", "uniprot_id": "P11111"}],
            }
        }
        adapter = ChEMBLAdapter(mock_data=mock, paralog_map=paralog_map)
        out = adapter.get("chemistry", pcsk9_claim)
        assert out["chembl_distinct_compounds"] == 100
        assert out["chembl_paralog_compound_counts"] == {"TestParalog": 50}

    def test_paralog_counts_absent_when_no_map(self, pcsk9_claim):
        mock = {
            "target_by_uniprot:Q8NBP7": {"targets": [{"target_chembl_id": "CHEMBL3243"}]},
            "activities:CHEMBL3243": {"activities": [
                {"molecule_chembl_id": f"CHEMBL_PRIMARY_{i}"} for i in range(100)
            ]},
        }
        adapter = ChEMBLAdapter(mock_data=mock, paralog_map={})
        out = adapter.get("chemistry", pcsk9_claim)
        assert "chembl_paralog_compound_counts" not in out

    def test_paralog_map_loader(self, tmp_path):
        """load_paralog_map reads the YAML and returns the right shape."""
        from adapters import load_paralog_map
        p = tmp_path / "paralog_map.yaml"
        p.write_text(
            "paralog_map:\n"
            "  P29597:\n"
            "    primary_symbol: TYK2\n"
            "    paralogs:\n"
            "      - {symbol: JAK1, uniprot_id: P23458}\n"
        )
        m = load_paralog_map(p)
        assert "P29597" in m
        assert m["P29597"]["primary_symbol"] == "TYK2"
        assert m["P29597"]["paralogs"][0]["symbol"] == "JAK1"

    def test_paralog_map_loader_missing_file_returns_empty(self, tmp_path):
        from adapters import load_paralog_map
        m = load_paralog_map(tmp_path / "does_not_exist.yaml")
        assert m == {}
