"""
verify_fixes.py — proves the six direction_audit fixes offline (no network).
Mocks only the HTTP boundary (requests.post) and runs the REAL parsing/transport.
Run:  python verify_fixes.py   ->  expect 6/6 and exit 0.
"""
import glob
import os
import sys
from unittest import mock

REPO = "https://github.com/crisprking/falsifiable-targets.git"
_FTREPO = "/kaggle/working/ftrepo"  # allowed: guarded Kaggle fallback, not a hardcoded project root

def _ensure_hardened_engine():
    """Return a dir whose direction_audit.py IS the hardened engine (cis_qtl_dirs +
    OpenTargetsError). Prefers a local hardened copy; else clones the repo. Skips stale
    copies such as an attached dataset snapshot (the cause of phantom fix failures)."""
    import shutil
    import subprocess
    seen, cands = [], [os.getcwd(), _FTREPO]
    cands += [os.path.dirname(p) for p in glob.glob("/kaggle/working/**/direction_audit.py", recursive=True)]  # allowed: guarded fallback
    for c in cands:
        f = os.path.join(c, "direction_audit.py")
        if c in seen or not os.path.exists(f): continue
        seen.append(c)
        src = open(f, encoding="utf-8", errors="ignore").read()
        if "def cis_qtl_dirs" in src and "class OpenTargetsError" in src: return c
    if os.path.isdir(_FTREPO): shutil.rmtree(_FTREPO)
    subprocess.run(["git", "clone", "--depth", "1", REPO, _FTREPO], check=True, capture_output=True, text=True)
    return _FTREPO

_d = _ensure_hardened_engine()
sys.path.insert(0, _d)
sys.modules.pop("direction_audit", None)        # never trust a stale cached module
import requests

import direction_audit as da

print("direction_audit fix verification")
print("  using:", getattr(da, "__file__", "?"))
HARDENED = hasattr(da, "OpenTargetsError")
print("  hardened engine:", "YES (OpenTargetsError present)" if HARDENED else ">>> NO — this is the ORIGINAL engine; fixes will FAIL <<<")
print()

class Resp:
    def __init__(self, data=None, status=200, errors=None):
        self.status_code = status; self._d = {}
        if data is not None: self._d["data"] = data
        if errors is not None: self._d["errors"] = errors
    def json(self): return self._d

def coloc_resp(rows):   # shape a COLOC_Q response
    return Resp(data={"credibleSet": {"studyLocusId": "L", "colocalisation": {"rows": rows}}})
def cis_row(sign, h4, ens="ENSG1"):
    return {"betaRatioSignAverage": sign, "h4": h4, "rightStudyType": "eqtl",
            "otherStudyLocus": {"qtlGeneId": ens, "study": {"target": {"id": ens, "approvedSymbol": "G"}}}}

PASS, FAIL = [], []
def ok(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  [{'PASS' if cond else 'FAIL'}] {name}" + (f"  - {detail}" if detail else ""))

# A: null-h4 coloc must be DROPPED, a real-h4 coloc KEPT
try:
    with mock.patch.object(requests, "post", lambda *a, **k: coloc_resp([cis_row(1.0, None), cis_row(1.0, 0.97)])):
        got = da.cis_signs("L", "ENSG1", "G", 0.8)
    ok("A null-h4 excluded, real-h4 kept", got == [(1.0, 0.97, "eqtl")], f"got {got}")
except Exception as e:
    ok("A null-h4 excluded, real-h4 kept", False, f"raised {type(e).__name__}: {e}")

# B: biomarker rescue sha must DIFFER for concordant vs discordant on the same gene/biomarker
def b_dispatch(*a, **k):
    body = k.get("json") or (a[1] if len(a) > 1 else {})
    q = body.get("query", ""); v = body.get("variables", {})
    if "gwas_credible_sets" in q:
        if v.get("f") == "EFO_0004611":   # LDL biomarker -> 3 clean loci
            return Resp(data={"target": {"approvedSymbol": "G", "evidences": {"rows":
                   [{"credibleSet": {"studyLocusId": f"L{i}"}} for i in range(3)]}}})
        return Resp(data={"target": {"approvedSymbol": "G", "evidences": {"rows": []}}})  # disease -> no loci
    if "size:1000" in q:               # labelled_direction -> no label
        return Resp(data={"target": {"evidences": {"rows": []}}})
    if "colocalisation" in q:           # every LDL locus -> +1 cis -> inhibitor
        return Resp(data={"credibleSet": {"studyLocusId": v.get("id"), "colocalisation": {"rows": [cis_row(1.0, 0.95)]}}})
    return Resp(data={})
try:
    with mock.patch.object(requests, "post", b_dispatch):
        vc = da.verdict_with_fallback("ENSG1", "EFO_0001645", "inhibitor")   # concordant
        vd = da.verdict_with_fallback("ENSG1", "EFO_0001645", "activator")   # discordant
    same_path = vc["verdict"].endswith("VIA_BIOMARKER") and vd["verdict"].endswith("VIA_BIOMARKER")
    ok("B rescue sha differs across verdicts", same_path and vc["sha256"] != vd["sha256"],
       f"{vc['verdict']}={vc['sha256'][:12]} vs {vd['verdict']}={vd['sha256'][:12]}")
except Exception as e:
    ok("B rescue sha differs across verdicts", False, f"raised {type(e).__name__}: {e}")

# C: GraphQL errors with null data must RAISE, not return silently
try:
    with mock.patch.object(requests, "post", lambda *a, **k: Resp(errors=[{"message": "boom"}], data={"target": None})):
        try:
            da.post("q", label="t"); ok("C error+null-data raises", False, "returned instead of raising")
        except da.OpenTargetsError:
            ok("C error+null-data raises OpenTargetsError", True)
except Exception as e:
    ok("C error+null-data raises OpenTargetsError", False, f"{type(e).__name__}: {e}")

# D: a hard HTTP (400) must RAISE
try:
    with mock.patch.object(requests, "post", lambda *a, **k: Resp(status=400)):
        try:
            da.post("q", label="t"); ok("D HTTP 400 raises", False, "returned instead of raising")
        except da.OpenTargetsError:
            ok("D HTTP 400 raises OpenTargetsError", True)
except Exception as e:
    ok("D HTTP 400 raises OpenTargetsError", False, f"{type(e).__name__}: {e}")

# E: transient network errors must RETRY (not zero out on the first blip)
calls = {"n": 0}
def flaky(*a, **k):
    calls["n"] += 1
    if calls["n"] < 3: raise requests.exceptions.ConnectionError("blip")
    return Resp(data={"ok": 1})
try:
    with mock.patch.object(requests, "post", flaky), mock.patch.object(da.time, "sleep", lambda *_: None):
        d = da.post("q", retries=3, label="t")
    ok("E retried then succeeded", isinstance(d, dict) and d.get("data", {}).get("ok") == 1 and calls["n"] == 3, f"calls={calls['n']}")
except Exception as e:
    ok("E retried then succeeded", False, f"{type(e).__name__}: {e}")

# F: the locus subset must be order-independent (sorted before cap)
def gwas_resp(order):
    return Resp(data={"target": {"approvedSymbol": "G", "evidences": {"rows":
           [{"credibleSet": {"studyLocusId": x}} for x in order]}}})
try:
    with mock.patch.object(requests, "post", lambda *a, **k: gwas_resp(["L3", "L1", "L2"])):
        _, ids1 = da.gwas_locus_ids("ENSG1", "EFO_x")
    with mock.patch.object(requests, "post", lambda *a, **k: gwas_resp(["L2", "L3", "L1"])):
        _, ids2 = da.gwas_locus_ids("ENSG1", "EFO_x")
    ok("F locus subset is order-independent (sorted)", ids1 == ids2 == ["L1", "L2", "L3"], f"{ids1} vs {ids2}")
except Exception as e:
    ok("F locus subset is order-independent (sorted)", False, f"{type(e).__name__}: {e}")

# behavior preservation: a normal h4=0.97 coloc still counts exactly as before
try:
    with mock.patch.object(requests, "post", lambda *a, **k: coloc_resp([cis_row(0.93, 0.97)])):
        got = da.cis_signs("L", "ENSG1", "G", 0.8)
    ok("(preserved) real h4=0.97 coloc still counted", got == [(0.93, 0.97, "eqtl")], f"got {got}")
except Exception as e:
    ok("(preserved) real h4=0.97 coloc still counted", False, f"{type(e).__name__}: {e}")

print("\n" + "=" * 64)
print(f"{len(PASS)}/{len(PASS)+len(FAIL)} checks passed")
if FAIL:
    print("FAILED:", FAIL); sys.exit(1)
print("ALL SIX FIXES VERIFIED — and the h4=0.97 control confirms validated behavior is preserved.")
