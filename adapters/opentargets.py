"""
OpenTargetsAdapter — live genetics evidence for falsifiable-targets R3.

Conforms to the Adapter protocol (adapters/protocol.py): a single
`get(section, claim) -> dict`. Drops into CompositeAdapter ABOVE
FixtureAdapter so live values win and the fixture backfills any field
this adapter can't resolve.

DESIGN DECISIONS (read before extending)
-----------------------------------------
1. This adapter SUPPLIES EVIDENCE. It does NOT emit a verdict.
   Direction-of-effect *concordance with the claim's mechanism* is a
   RULE decision and belongs in R3's SHA-locked logic, not here. The
   adapter returns the gene's genetic direction + provenance + a quality
   tier; R3 compares that to claim.mechanism_direction and emits the
   verdict. This keeps the audit logic hashable and the adapter swappable.

2. Direction is DERIVED, not read off a single field. Open Targets
   Genetics was folded into the unified Platform; there is no clean
   `directionOfEffect` scalar in the GraphQL API. We derive direction
   from three sources, in descending cleanliness:
       (a) gene_burden  — rare LoF-collapsing burden vs the trait
                          (cleanest: the HMGCR/PCSK9 archetype)
       (b) eva (ClinVar) — pathogenic LoF vs GoF clinical variants
       (c) gwas_credible_sets + molQTL coloc — sign of CredibleSet.beta
                          combined with eQTL direction (noisiest; the
                          OR=1.61 "coloc barely predicts approval" caveat
                          lives here)
   The tier returned reflects which source carried the call.

3. OT's own DOE assessment is treated as ONE source among several, never
   as ground truth — that is the entire point of the flagship audit.
   If/when you wire OT's published DOE field in, add it as a fourth
   contributor with its own provenance flag; do not let it silently
   override (a).

CONFIRMED against api.platform.opentargets.org/api/v4/graphql (June 2026):
  - endpoint, `target(ensemblId)`, `search`, `disease(efoId)`,
    `evidences(efoIds, datasourceIds, size)`, `CredibleSet.beta`.
TODO-CONFIRM via live schema introspection before production:
  - exact datasourceId strings post-25.12 "exploded evidence" refactor
    (gene_burden / eva / gwas_credible_sets are the expected IDs).
  - per-datasource row fields that encode direction (beta / oddsRatio /
    variantEffect) differ by datasource after the refactor.
  Use the official OT MCP server (github.com/opentargets/open-targets-platform-mcp,
  `get_open_targets_graphql_schema`) to confirm field paths rather than guessing.
"""
from __future__ import annotations

import hashlib
import json
import os
import time
from pathlib import Path
from typing import Any

import requests

OT_GRAPHQL = "https://api.platform.opentargets.org/api/v4/graphql"

# datasourceIds expected after the 25.12 evidence refactor. CONFIRM via schema.
DS_BURDEN = "gene_burden"
DS_CLINVAR = "eva"
DS_GWAS = "gwas_credible_sets"
DIRECTION_DATASOURCES = (DS_BURDEN, DS_CLINVAR, DS_GWAS)


class OpenTargetsAdapter:
    """Live genetics adapter. NEVER raises on missing data — returns {}."""

    def __init__(
        self,
        cache_dir: str | None = None,
        offline: bool | None = None,
        timeout: float = 30.0,
        max_retries: int = 3,
    ):
        self._cache = Path(cache_dir or os.environ.get("AE_CACHE_DIR", ".ae_cache"))
        self._cache.mkdir(parents=True, exist_ok=True)
        # Respect AE_OFFLINE=1: serve from cache only, {} on miss.
        env_offline = os.environ.get("AE_OFFLINE", "") == "1"
        self._offline = env_offline if offline is None else offline
        self._timeout = timeout
        self._max_retries = max_retries

    # ---- protocol entry point -------------------------------------------
    def get(self, section: str, claim: Any) -> dict[str, Any]:
        if section != "genetics":
            return {}

        ensembl = self._resolve_target(claim)
        efo = self._resolve_disease(claim)
        if not ensembl or not efo:
            # Can't resolve identifiers -> let the fixture backfill.
            return {}

        ev = self._query_evidence(ensembl, efo)
        if ev is None:
            return {}

        rows = ev.get("rows", [])
        by_ds: dict[str, list[dict]] = {}
        for r in rows:
            by_ds.setdefault(r.get("datasourceId", ""), []).append(r)

        gene_direction, tier, sources = self._derive_direction(by_ds, ensembl, efo)

        # Backward-compatible keys R3 already reads (do not rename), plus
        # the new direction-of-effect evidence the ruleset bump will consume.
        return {
            "gwas_hits": len(by_ds.get(DS_GWAS, [])),
            "mendelian_evidence": bool(by_ds.get(DS_CLINVAR) or by_ds.get(DS_BURDEN)),
            "somatic_driver_evidence": self._has_somatic(by_ds),
            "open_targets_score": self._association_score(ensembl, efo),
            # --- new DOE evidence (R3 v1.5.0 reads these) ---
            "doe_gene_direction": gene_direction,   # GoF_increases_risk | LoF_increases_risk | conflicting | unknown
            "doe_tier": tier,                        # strong | moderate | weak | absent
            "doe_sources": sources,                  # provenance: which datasources contributed
            "doe_n_assessments": sum(len(v) for v in by_ds.values()),
            "_provenance": {"ensembl_id": ensembl, "efo_id": efo, "endpoint": OT_GRAPHQL},
        }

    # ---- direction derivation -------------------------------------------
    def _derive_direction(self, by_ds, ensembl, efo) -> tuple[str, str, list[str]]:
        """Return (gene_direction, tier, contributing_sources).

        Cleanliness order: burden > clinvar > gwas+coloc. The first source
        that yields an unambiguous call sets the tier; later sources only
        downgrade to 'conflicting' if they disagree.
        """
        calls: list[tuple[str, str, str]] = []  # (direction, tier, source)

        # (a) gene_burden — rare LoF collapsing. CONFIRM the directional field.
        for r in by_ds.get(DS_BURDEN, []):
            d = self._burden_direction(r)
            if d:
                calls.append((d, "strong", DS_BURDEN))

        # (b) ClinVar / eva — clinical LoF vs GoF.
        for r in by_ds.get(DS_CLINVAR, []):
            d = self._clinvar_direction(r)
            if d:
                calls.append((d, "moderate", DS_CLINVAR))

        # (c) GWAS credible set + molQTL coloc — sign(beta) + eQTL direction.
        for r in by_ds.get(DS_GWAS, []):
            d = self._gwas_coloc_direction(r)
            if d:
                calls.append((d, "weak", DS_GWAS))

        if not calls:
            return "unknown", "absent", []

        directions = {c[0] for c in calls}
        sources = sorted({c[2] for c in calls})
        if len(directions) > 1:
            # Genuine disagreement across sources — surface it, never average.
            return "conflicting", "conflicting", sources
        # Single agreed direction: tier = best (cleanest) contributing source.
        best_tier = next(t for t in ("strong", "moderate", "weak") if any(c[1] == t for c in calls))
        return directions.pop(), best_tier, sources

    # The three extractors below are deliberately thin and flagged: the exact
    # row fields are datasource-specific and MUST be confirmed against the live
    # schema. Returning None means "no directional signal from this row."
    def _burden_direction(self, row: dict) -> str | None:
        beta = row.get("beta")
        if beta is None:
            return None  # TODO-CONFIRM burden direction field name
        # LoF-collapsing burden: beta>0 vs risk => LoF increases risk.
        return "LoF_increases_risk" if beta > 0 else "GoF_increases_risk"

    def _clinvar_direction(self, row: dict) -> str | None:
        # TODO-CONFIRM: map variantFunctionalConsequenceId / clinicalSignificances
        # (e.g. SO loss_of_function terms) -> LoF; activating terms -> GoF.
        return None

    def _gwas_coloc_direction(self, row: dict) -> str | None:
        # TODO-CONFIRM: combine sign(CredibleSet.beta) for the disease locus with
        # the colocalising eQTL direction to infer whether risk allele raises or
        # lowers gene activity. Noisiest source (see OR=1.61 coloc caveat).
        return None

    @staticmethod
    def _has_somatic(by_ds) -> bool:
        return any(
            (r.get("datatypeId") == "somatic_mutation")
            for rows in by_ds.values() for r in rows
        )

    # ---- GraphQL calls (cached, offline-aware) --------------------------
    def _query_evidence(self, ensembl: str, efo: str) -> dict | None:
        q = """
        query Ev($ensemblId: String!, $efoId: String!, $size: Int!) {
          target(ensemblId: $ensemblId) {
            id approvedSymbol
            evidences(efoIds: [$efoId], datasourceIds: %s, size: $size) {
              count
              rows {
                datasourceId datatypeId score
                beta oddsRatio
                variantFunctionalConsequenceId
                clinicalSignificances
                studyLocusId
              }
            }
          }
        }
        """ % json.dumps(list(DIRECTION_DATASOURCES))
        data = self._post(q, {"ensemblId": ensembl, "efoId": efo, "size": 500})
        try:
            return data["data"]["target"]["evidences"]
        except (TypeError, KeyError):
            return None

    def _association_score(self, ensembl: str, efo: str) -> float | None:
        q = """
        query Assoc($ensemblId: String!, $efoId: String!) {
          disease(efoId: $efoId) {
            associatedTargets(Bs: [$ensemblId]) {
              rows { score datatypeScores { id score } }
            }
          }
        }
        """
        data = self._post(q, {"ensemblId": ensembl, "efoId": efo})
        try:
            rows = data["data"]["disease"]["associatedTargets"]["rows"]
            return rows[0]["score"] if rows else None
        except (TypeError, KeyError, IndexError):
            return None

    def _resolve_target(self, claim: Any) -> str | None:
        # Prefer an explicit ensembl_id on the claim; else resolve by symbol.
        ens = getattr(claim, "ensembl_id", None)
        if ens:
            return ens
        sym = getattr(claim, "target_symbol", None)
        if not sym:
            return None
        return self._search_id(sym, "target")

    def _resolve_disease(self, claim: Any) -> str | None:
        efo = getattr(claim, "efo_id", None)
        if efo:
            return efo
        # Free-text indication -> EFO is lossy; resolve best hit and record it.
        # RECOMMENDATION: add an explicit efo_id to the claim schema so audits
        # are deterministic rather than dependent on search ranking.
        ind = getattr(claim, "indication", None)
        if not ind:
            return None
        return self._search_id(ind, "disease")

    def _search_id(self, query_string: str, entity: str) -> str | None:
        q = """
        query S($q: String!, $entityNames: [String!]) {
          search(queryString: $q, entityNames: $entityNames) {
            hits { id entity }
          }
        }
        """
        data = self._post(q, {"q": query_string, "entityNames": [entity]})
        try:
            for hit in data["data"]["search"]["hits"]:
                if hit.get("entity") == entity:
                    return hit["id"]
        except (TypeError, KeyError):
            return None
        return None

    def _post(self, query: str, variables: dict) -> dict | None:
        key = hashlib.sha256(
            (query + json.dumps(variables, sort_keys=True)).encode()
        ).hexdigest()
        cache_file = self._cache / f"ot_{key}.json"
        if cache_file.exists():
            return json.loads(cache_file.read_text())
        if self._offline:
            return None  # AE_OFFLINE: cache-only, miss -> None
        for attempt in range(self._max_retries):
            try:
                resp = requests.post(
                    OT_GRAPHQL,
                    json={"query": query, "variables": variables},
                    timeout=self._timeout,
                )
                if resp.status_code == 200:
                    payload = resp.json()
                    cache_file.write_text(json.dumps(payload))  # deterministic replay
                    return payload
                if resp.status_code in (429, 502, 503):
                    time.sleep(2 ** attempt)
                    continue
                return None
            except requests.RequestException:
                time.sleep(2 ** attempt)
        return None
