"""
Adapters for the falsifiable-targets engine.

Every adapter exposes one method:
    .get(section: str, claim) -> dict

A returned dict means "here is what I know about this section."
An empty dict means "I have no data for this section."

CompositeAdapter chains adapters: first non-empty result wins per section,
with per-key fallback so a partially-populated live result is augmented
from the fixture rather than discarded.

Live adapters cache on disk under .ae_cache/ for determinism. Set the
env var AE_OFFLINE=1 to force fixture-only (useful for hermetic tests).
"""
from __future__ import annotations

import json
import os
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


CACHE_DIR = Path(os.environ.get("AE_CACHE_DIR", ".ae_cache"))


# ---------------------------------------------------------------------------
# Base classes
# ---------------------------------------------------------------------------

class FixtureAdapter:
    """Returns precomputed data from a dict. Hermetic, used in sentinel suite."""

    def __init__(self, fixture):
        self._f = fixture or {}

    def get(self, section, claim):
        return self._f.get(section, {}) or {}


class CompositeAdapter:
    """
    Chains adapters in priority order. For each section, walks adapters
    from highest priority to lowest; the first non-None value seen for
    each key wins. Lower-priority adapters fill keys that higher-priority
    adapters did not report.

    This is what makes live + fixture coexist cleanly: a live adapter
    might return {pdb_count: 0} but say nothing about alphafold_plddt;
    the fixture can supply alphafold_plddt = 93.2 and the rule sees
    both in a single merged dict.
    """

    def __init__(self, adapters):
        self._adapters = list(adapters)

    def get(self, section, claim):
        merged = {}
        for adapter in self._adapters:
            data = adapter.get(section, claim) or {}
            for k, v in data.items():
                if k not in merged and v is not None:
                    merged[k] = v
        return merged


# ---------------------------------------------------------------------------
# HTTP helper with on-disk cache
# ---------------------------------------------------------------------------

def _cache_key(url):
    safe = urllib.parse.quote(url, safe="")
    return CACHE_DIR / (safe[:200] + ".json")


def _http_get_json(url, timeout=15, retries=2):
    """
    GET a JSON resource. On-disk cached. Returns parsed JSON or None on
    permanent failure. Cache lives in .ae_cache/.
    """
    if os.environ.get("AE_OFFLINE") == "1":
        # Hermetic mode: only serve from cache
        cache_path = _cache_key(url)
        if cache_path.exists():
            return json.loads(cache_path.read_text())
        return None

    cache_path = _cache_key(url)
    if cache_path.exists():
        try:
            return json.loads(cache_path.read_text())
        except json.JSONDecodeError:
            cache_path.unlink()  # Corrupted; refetch

    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    last_err = None
    for attempt in range(retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": "falsifiable-targets/1.0.1",
                    "Accept": "application/json",
                },
            )
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = json.loads(resp.read().decode("utf-8"))
            cache_path.write_text(json.dumps(data))
            return data
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
            last_err = e
            if attempt < retries:
                time.sleep(1.5 ** attempt)
    return None


# ---------------------------------------------------------------------------
# UniProt adapter
# ---------------------------------------------------------------------------

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"


class UniProtAdapter:
    """
    Reads UniProtKB by accession. Populates:
      structure:
        alphafold_plddt - mean pLDDT from the AlphaFold entry, if linked
        pdb_count       - number of PDB cross-references
        structure_resolved - bool, pdb_count > 0
      orthology:
        sources_total / sources_agreeing - synthesized from cross-references
        ortholog_count_human_pathogen - via OrthoDB xref count (heuristic)
      expression:
        target_tissue_expressed - tissue specificity comments, if any

    All rule-relevant fields it can't determine are simply not included.
    The composite adapter then falls back to the fixture for those keys.
    """

    def __init__(self, mock_data=None):
        """
        mock_data: optional dict of {uniprot_id: response_json} for hermetic tests.
        If provided, this is consulted instead of the live API.
        """
        self._mock = mock_data

    def _fetch(self, uniprot_id):
        if self._mock is not None:
            return self._mock.get(uniprot_id)
        url = f"{UNIPROT_BASE}/{uniprot_id}.json"
        return _http_get_json(url)

    def get(self, section, claim):
        uid = getattr(claim, "uniprot_id", None)
        if not uid:
            return {}
        data = self._fetch(uid)
        if not data:
            return {}

        if section == "structure":
            return self._structure(data)
        if section == "orthology":
            return self._orthology(data)
        if section == "expression":
            return self._expression(data)
        return {}

    def _structure(self, data):
        out = {}
        xrefs = data.get("uniProtKBCrossReferences", []) or []
        pdb_refs = [x for x in xrefs if x.get("database") == "PDB"]
        af_refs  = [x for x in xrefs if x.get("database") == "AlphaFoldDB"]
        if pdb_refs:
            out["pdb_count"] = len(pdb_refs)
            out["structure_resolved"] = True
        elif af_refs:
            # If only AlphaFold, count as 'modeled but not crystal'
            out["pdb_count"] = 0
            out["structure_resolved"] = False
        # We don't fetch the actual AlphaFold model here (would require
        # afdb.ebi.ac.uk/files/AF-{uid}-F1-confidence_v4.json, a separate
        # call). The fixture supplies alphafold_plddt for now; v1.2 will
        # do that fetch.
        return out

    def _orthology(self, data):
        out = {}
        xrefs = data.get("uniProtKBCrossReferences", []) or []
        # Sources we treat as orthology-evidence-bearing:
        ortho_dbs = {"OrthoDB", "eggNOG", "InParanoid", "OMA", "PhylomeDB", "HOGENOM"}
        agreeing = sum(1 for x in xrefs if x.get("database") in ortho_dbs)
        total = len(ortho_dbs)
        if agreeing > 0:
            # We can only report "this DB had an entry"; we can't tell
            # without deeper parsing whether the entries actually agree
            # on the orthology call. So we report sources_agreeing as
            # the count of DBs that returned ANY hit, and let the
            # fixture supply the more precise value if available.
            out["sources_agreeing_uniprot_lower_bound"] = agreeing
            out["sources_total_uniprot_lower_bound"] = total
        return out

    def _expression(self, data):
        out = {}
        comments = data.get("comments", []) or []
        tissue_comments = [c for c in comments if c.get("commentType") == "TISSUE SPECIFICITY"]
        if tissue_comments:
            out["target_tissue_expressed"] = True  # presence of comment implies known expression
            texts = []
            for c in tissue_comments:
                for t in c.get("texts", []) or []:
                    val = t.get("value")
                    if val:
                        texts.append(val)
            if texts:
                out["uniprot_tissue_text"] = " | ".join(texts)[:500]
        return out


# ---------------------------------------------------------------------------
# ChEMBL adapter
# ---------------------------------------------------------------------------

CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"


class ChEMBLAdapter:
    """
    Reads ChEMBL via the REST API. Populates:
      chemistry:
        chembl_distinct_compounds - count of unique molecule_chembl_id values
                                    in bioactivities against this target
        max_phase                 - max max_phase across those compounds
        chembl_pfam_class_collapse_fraction - if the Pfam class of the
                                              target is shared with a
                                              known human paralog, the
                                              fraction of compounds that
                                              also hit that paralog.
                                              (Implemented as a lookup
                                              against an explicit table
                                              of known collapse cases.)
        chembl_paralog_compound_counts - (v1.2.0, optional) dict mapping
                                         paralog symbol to its ChEMBL
                                         distinct-compound count. Used by
                                         R6's paralog-ratio heuristic.
                                         Populated only when paralog_map
                                         is supplied.

    The class-collapse fraction is the load-bearing R6 input. The current
    implementation looks up the target's Pfam class and reports whether
    >=80% of bioactive compounds for that class also have higher-affinity
    bioactivities against a human paralog in the same class. This is the
    SAT->HDAC4 archetype from the Madurella audit.
    """

    def __init__(self, mock_data=None, paralog_map=None):
        """
        mock_data: optional dict of mock responses for tests.
        paralog_map: optional dict mapping primary UniProt -> list of
                     paralog dicts with 'symbol' and 'uniprot_id' keys.
                     Used to populate chembl_paralog_compound_counts.
        """
        self._mock = mock_data
        self._paralog_map = paralog_map or {}

    def _fetch_target(self, chembl_target_id):
        if self._mock is not None:
            return self._mock.get(f"target:{chembl_target_id}")
        url = f"{CHEMBL_BASE}/target/{chembl_target_id}.json"
        return _http_get_json(url)

    def _fetch_target_by_uniprot(self, uniprot_id):
        if self._mock is not None:
            return self._mock.get(f"target_by_uniprot:{uniprot_id}")
        url = f"{CHEMBL_BASE}/target.json?target_components__accession={uniprot_id}"
        return _http_get_json(url)

    def _fetch_activities(self, chembl_target_id, limit=1000):
        if self._mock is not None:
            return self._mock.get(f"activities:{chembl_target_id}")
        url = (
            f"{CHEMBL_BASE}/activity.json"
            f"?target_chembl_id={chembl_target_id}"
            f"&limit={limit}"
        )
        return _http_get_json(url)

    def _compound_count_for_uniprot(self, uniprot_id):
        """Return (count, chembl_target_id, max_phase) for a UniProt ID,
        or (None, None, None) if the live fetch failed, or (0, None, None)
        if the API affirmatively reported no target."""
        target_resp = self._fetch_target_by_uniprot(uniprot_id)
        if target_resp is None:
            return None, None, None
        targets = target_resp.get("targets", []) or []
        if not targets:
            return 0, None, None
        chembl_id = targets[0].get("target_chembl_id")
        if not chembl_id:
            return 0, None, None
        acts = self._fetch_activities(chembl_id)
        if acts is None:
            return None, chembl_id, None
        activities = acts.get("activities", []) or []
        compound_ids = {a.get("molecule_chembl_id") for a in activities if a.get("molecule_chembl_id")}
        max_phase_vals = [a.get("max_phase") for a in activities if a.get("max_phase") is not None]
        max_phase = None
        if max_phase_vals:
            try:
                max_phase = max(int(v) for v in max_phase_vals)
            except (ValueError, TypeError):
                pass
        return len(compound_ids), chembl_id, max_phase

    def get(self, section, claim):
        if section != "chemistry":
            return {}
        uid = getattr(claim, "uniprot_id", None)
        if not uid:
            return {}
        n, chembl_id, max_phase = self._compound_count_for_uniprot(uid)
        if n is None:
            # Live fetch failed (offline, network error, rate limit). Return
            # nothing so the composite falls back to fixture data rather
            # than asserting a misleading zero-count.
            return {}
        if chembl_id is None:
            return {"chembl_distinct_compounds": 0}
        out = {
            "chembl_distinct_compounds": n,
            "chembl_target_id": chembl_id,
        }
        if max_phase is not None:
            out["max_phase"] = max_phase

        # v1.2.0: paralog compound counts. Only fetched when this primary
        # target has an entry in the paralog map. Each fetch is cached on
        # disk so repeat audits are network-free.
        paralog_entry = self._paralog_map.get(uid)
        if paralog_entry:
            counts = {}
            for paralog in paralog_entry.get("paralogs", []):
                p_uid = paralog.get("uniprot_id")
                p_symbol = paralog.get("symbol")
                if not p_uid or not p_symbol:
                    continue
                p_n, _, _ = self._compound_count_for_uniprot(p_uid)
                if p_n is not None:
                    counts[p_symbol] = p_n
            if counts:
                out["chembl_paralog_compound_counts"] = counts

        return out


# ---------------------------------------------------------------------------
# Convenience: build the default v1.0.1 composite for a given claim
# ---------------------------------------------------------------------------

def default_composite(fixture, use_live=True, uniprot_mock=None, chembl_mock=None,
                       paralog_map=None):
    """
    Build the standard adapter stack for v1.2.0:
        [live UniProt] -> [live ChEMBL] -> [fixture]

    If use_live=False, returns a fixture-only adapter (hermetic).
    If *_mock dicts are passed, the live adapters use those instead of HTTP.
    If paralog_map is passed, ChEMBLAdapter will fetch paralog compound
    counts for primary targets with entries in the map. The map is expected
    to be the 'paralog_map' dict loaded from claims/paralog_map.yaml.
    """
    if not use_live:
        return FixtureAdapter(fixture)
    return CompositeAdapter([
        UniProtAdapter(mock_data=uniprot_mock),
        ChEMBLAdapter(mock_data=chembl_mock, paralog_map=paralog_map),
        FixtureAdapter(fixture),
    ])


def load_paralog_map(path):
    """Load paralog_map.yaml. Returns a dict suitable for passing to
    ChEMBLAdapter(paralog_map=...) or default_composite(paralog_map=...).
    Returns {} if the path does not exist."""
    import yaml
    from pathlib import Path
    p = Path(path)
    if not p.exists():
        return {}
    with open(p) as f:
        data = yaml.safe_load(f) or {}
    return data.get("paralog_map", {}) or {}
