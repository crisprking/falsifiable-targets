"""
test_direction_logic.py — locks the direction-of-effect logic that every headline
number depends on. Pure logic only (network functions are stubbed), so it runs
offline and deterministically. If any of these fail, a reported number is wrong.

Run:  python tests/test_direction_logic.py   (or: pytest tests/test_direction_logic.py)
"""
import glob
import os
import sys


def _locate_direction_audit():
    """Find direction_audit.py whether run as a repo test, a script, or pasted in a notebook."""
    try:
        return None
    except Exception:
        pass
    cands = []
    try:                       # repo layout: tests/ sits next to direction_audit.py
        here = os.path.dirname(os.path.abspath(__file__))
        cands += [os.path.join(here, ".."), here]
    except NameError:          # __file__ undefined in a notebook cell
        pass
    cands.append(os.getcwd())
    for base in ("/kaggle/working", "/kaggle/input"):  # allowed: guarded Kaggle/dataset fallback
        if os.path.isdir(base):
            cands += [os.path.dirname(p) for p in glob.glob(base + "/**/direction_audit.py", recursive=True)]
    for c in cands:
        if c and os.path.exists(os.path.join(c, "direction_audit.py")):
            return c
    raise FileNotFoundError("direction_audit.py not found — run from the repo, or attach it to the Kaggle session.")

_d = _locate_direction_audit()
if _d:
    sys.path.insert(0, _d)
import direction_audit as da

# ----------------------------------------------------------------------------
# stubs: control the network-touching layer so we test pure logic deterministically
#   LOCI[efo]  -> list of locus ids        (what gwas_locus_ids returns)
#   SIGNS[lid] -> list of (sign, h4, qtl)  (what cis_signs returns for that locus)
#   LABELS[(ensembl, efo)] -> (label, n)   (what labelled_direction returns)
# ----------------------------------------------------------------------------
LOCI, SIGNS, LABELS = {}, {}, {}
_ORIG = {n: getattr(da, n) for n in ("gwas_locus_ids", "cis_signs", "labelled_direction")}
da.gwas_locus_ids   = lambda ens, efo, cap=40: ("SYM", list(LOCI.get(efo, [])))
da.cis_signs        = lambda lid, ens, sym, h4: list(SIGNS.get(lid, []))
da.labelled_direction = lambda ens, efo: LABELS.get((ens, efo), (None, 0))

def _reset():
    LOCI.clear(); SIGNS.clear(); LABELS.clear()

def _loci(efo, per_locus_signs):
    """per_locus_signs: list of lists of bare sign floats; builds LOCI/SIGNS."""
    ids = []
    for i, signs in enumerate(per_locus_signs):
        lid = f"{efo}:L{i}"; ids.append(lid)
        SIGNS[lid] = [(s, 0.95, "eqtl") for s in signs]
    LOCI[efo] = ids

PASS = []; FAIL = []
def check(name, cond):
    (PASS if cond else FAIL).append(name)
    print(f"  {'ok ' if cond else 'FAIL'}  {name}")

# ============================================================================
print("\n[1] desired_from_label — the (directionOnTarget x directionOnTrait) truth table")
# inhibitor: losing function protects, OR gaining function raises risk
for dot in ("loss_of_function", "lof", "decrease"):
    check(f"({dot}, protect) -> inhibitor", da.desired_from_label(dot, "protect") == "inhibitor")
for dot in ("gain_of_function", "gof", "increase"):
    check(f"({dot}, risk) -> inhibitor", da.desired_from_label(dot, "risk") == "inhibitor")
# activator: losing function raises risk, OR gaining function protects
for dot in ("loss_of_function", "lof", "decrease"):
    check(f"({dot}, risk) -> activator", da.desired_from_label(dot, "risk") == "activator")
for dot in ("gain_of_function", "gof", "increase"):
    check(f"({dot}, protect) -> activator", da.desired_from_label(dot, "protect") == "activator")
# case-insensitive
check("case-insensitive ('LoF','Risk') -> activator", da.desired_from_label("LoF", "Risk") == "activator")
# invalid inputs must NOT silently produce a direction
check("unknown trait -> None", da.desired_from_label("loss_of_function", "huh") is None)
check("missing target -> None", da.desired_from_label(None, "risk") is None)
check("unknown target token -> None", da.desired_from_label("weird", "risk") is None)
check("missing trait -> None", da.desired_from_label("loss_of_function", None) is None)

print("\n[2] _flip — biomarker sign inversion")
check("flip(inhibitor)=activator", da._flip("inhibitor") == "activator")
check("flip(activator)=inhibitor", da._flip("activator") == "inhibitor")
check("flip(None)=None", da._flip(None) is None)

print("\n[3] recover — locus-consensus thresholds (MIN_LOCI=3, FRAC=0.8, MAG=0.8) + sign orientation")
_reset(); _loci("e1", [[1.0], [1.0]])                       # 2 decided < MIN_LOCI
check("2 decided loci -> recovered None (insufficient)", da.recover("ENS", "e1")["recovered"] is None)
_reset(); _loci("e2", [[1.0], [0.9], [1.0]])                # 3 positive
r = da.recover("ENS", "e2")
check("3 positive loci -> inhibitor @ 100%", r["recovered"] == "inhibitor" and r["consensus"] == 1.0)
_reset(); _loci("e3", [[-1.0], [-0.9], [-1.0]])             # 3 negative
check("3 negative loci -> activator (sign<0)", da.recover("ENS", "e3")["recovered"] == "activator")
_reset(); _loci("e4", [[1.0], [1.0], [1.0], [1.0], [-1.0]]) # 4:1 -> 0.8 == FRAC
r = da.recover("ENS", "e4")
check("4:1 split -> inhibitor @ 80% (meets FRAC)", r["recovered"] == "inhibitor" and r["consensus"] == 0.8)
_reset(); _loci("e5", [[1.0], [1.0], [1.0], [-1.0], [-1.0]])# 3:2 -> 0.6 < FRAC
check("3:2 split -> recovered None (below FRAC)", da.recover("ENS", "e5")["recovered"] is None)
_reset(); _loci("e6", [[0.5], [0.5], [0.5]])                # all |sign|<MAG
r = da.recover("ENS", "e6")
check("weak signs (|s|<MAG) -> 0 decided loci", r["n_loci_decided"] == 0 and r["recovered"] is None)
_reset(); _loci("e7", [[0.9, -0.85], [1.0], [1.0], [1.0]])  # locus0 is an internal tie -> undecided
r = da.recover("ENS", "e7")
check("intra-locus tie is undecided; 3 clean -> inhibitor", r["recovered"] == "inhibitor" and r["n_loci_decided"] == 3)

print("\n[4] verdict — branch routing, vouches flag, deterministic sha")
_reset(); LABELS[("ENS", "d")] = ("inhibitor", 5)
v = da.verdict("ENS", "d", "inhibitor")
check("label matches claim -> LABEL_CONCORDANT, vouches", v["verdict"] == "LABEL_CONCORDANT" and v["vouches"])
_reset(); LABELS[("ENS", "d")] = ("activator", 5)
v = da.verdict("ENS", "d", "inhibitor")
check("label opposes claim -> LABEL_DISCORDANT, no vouch", v["verdict"] == "LABEL_DISCORDANT" and not v["vouches"])
_reset(); _loci("d", [[1.0], [1.0], [1.0]])                 # no label, recovers inhibitor
v = da.verdict("ENS", "d", "inhibitor")
check("no label + recovered==claim -> RECOVERED_CONCORDANT", v["verdict"] == "RECOVERED_CONCORDANT" and v["vouches"])
v2 = da.verdict("ENS", "d", "activator")
check("no label + recovered!=claim -> RECOVERED_DISCORDANT", v2["verdict"] == "RECOVERED_DISCORDANT" and not v2["vouches"])
_reset(); _loci("d", [[1.0], [1.0]])
check("no label + <MIN_LOCI -> INSUFFICIENT_DIRECTION", da.verdict("ENS", "d", "inhibitor")["verdict"] == "INSUFFICIENT_DIRECTION")
_reset(); _loci("d", [[1.0], [1.0], [1.0], [-1.0], [-1.0]])
check("no label + split -> RECOVERY_CONFLICTED", da.verdict("ENS", "d", "inhibitor")["verdict"] == "RECOVERY_CONFLICTED")
try:
    da.verdict("ENS", "d", "not-a-mechanism"); check("bad mechanism raises", False)
except ValueError:
    check("bad mechanism raises ValueError", True)
_reset(); _loci("d", [[1.0], [1.0], [1.0]])
a = da.verdict("ENS", "d", "inhibitor")["sha256"]; b = da.verdict("ENS", "d", "inhibitor")["sha256"]
check("sha is deterministic", a == b and len(a) == 64)
c = da.verdict("ENS", "d", "activator")["sha256"]
check("sha changes when verdict changes", a != c)

print("\n[5] verdict_with_fallback — fires ONLY on INSUFFICIENT, never on CONFLICTED; sign-mapped; honest on disagreement")
# disease endpoint INSUFFICIENT, biomarker (LDL, +1) recovers inhibitor -> rescued concordant
da.BIOMARKER_BRIDGES = {"DIS": [("BMK", "LDL cholesterol", +1)]}
_reset(); _loci("DIS", [[1.0], [1.0]]); _loci("BMK", [[1.0], [1.0], [1.0]])
v = da.verdict_with_fallback("ENS", "DIS", "inhibitor")
check("INSUFFICIENT + bridge inhibitor -> *_CONCORDANT_VIA_BIOMARKER", v["verdict"] == "RECOVERED_CONCORDANT_VIA_BIOMARKER" and v["vouches"])
check("rescued verdict carries a sha", len(v.get("sha256", "")) == 64)
# the documented rule: CONFLICTED is NOT rescued (a split is information, not a gap)
_reset(); _loci("DIS", [[1.0], [1.0], [1.0], [-1.0], [-1.0]]); _loci("BMK", [[1.0], [1.0], [1.0]])
v = da.verdict_with_fallback("ENS", "DIS", "inhibitor")
check("RECOVERY_CONFLICTED is NOT overridden by the bridge", v["verdict"] == "RECOVERY_CONFLICTED")
# a direct call is never touched
_reset(); _loci("DIS", [[1.0], [1.0], [1.0]]); _loci("BMK", [[-1.0], [-1.0], [-1.0]])
v = da.verdict_with_fallback("ENS", "DIS", "inhibitor")
check("direct RECOVERED call left untouched", v["verdict"] == "RECOVERED_CONCORDANT" and v["source"] == "coloc-recovered")
# sign -1 bridge flips: biomarker recovers inhibitor -> mapped activator
da.BIOMARKER_BRIDGES = {"DIS": [("BMK", "HDL cholesterol", -1)]}
_reset(); _loci("DIS", [[1.0], [1.0]]); _loci("BMK", [[1.0], [1.0], [1.0]])  # recovers inhibitor
v = da.verdict_with_fallback("ENS", "DIS", "activator")
check("sign -1 bridge flips inhibitor->activator (mapped)", v["direction"] == "activator" and v["vouches"])
# two bridges that disagree after mapping -> keep the honest refusal
da.BIOMARKER_BRIDGES = {"DIS": [("B1", "LDL", +1), ("B2", "ApoB", +1)]}
_reset(); _loci("DIS", [[1.0], [1.0]]); _loci("B1", [[1.0], [1.0], [1.0]]); _loci("B2", [[-1.0], [-1.0], [-1.0]])
v = da.verdict_with_fallback("ENS", "DIS", "inhibitor")
check("bridges disagree -> stays INSUFFICIENT (no rescue)", v["verdict"] == "INSUFFICIENT_DIRECTION")

# restore the real engine functions so this never pollutes other notebook cells
for _n, _f in _ORIG.items():
    setattr(da, _n, _f)

print("\n" + "=" * 70)
print(f"PASSED {len(PASS)}   FAILED {len(FAIL)}")
if FAIL:
    print("FAILURES:", FAIL); sys.exit(1)
print("ALL DIRECTION-LOGIC TESTS PASSED — the truth table, sign convention, consensus")
print("thresholds, and fallback rules behave exactly as documented.")
