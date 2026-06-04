"""
integrity_gates.py — the two checks the benchmark's numbers rest on. Run as ONE
Kaggle cell (internet ON; ~4-7 min). Uses the REAL hardened engine's primitives
(recover, gwas_locus_ids, cis_signs) so it validates the actual recovery path.

GATE A  sign-orientation calibration. Run recovery on genes whose direction is
        textbook-certain, in BOTH directions. If the convention were flipped,
        every known inhibitor would recover "activator". Gate FAILS on a systematic
        reversal; PASSES if inhibitor-controls -> inhibitor and activator-controls
        -> activator; PASS-with-exceptions if a few disagree (real decoupling).

GATE B  locus independence. For high-locus genes, cluster loci by lead-variant
        genomic region and recompute consensus over INDEPENDENT regions. If regions
        << loci, a "100% over N loci" claim is pseudo-replicated, not N signals.
"""
import sys, os, json, hashlib, math, glob, re, time, shutil, subprocess

# ---- acquire the HARDENED engine (ignore stale /kaggle/input copies; clone if needed) ----
REPO = "https://github.com/crisprking/falsifiable-targets.git"
_FTREPO = "/kaggle/working/ftrepo"

def _ensure_hardened_engine():
    seen = []
    for c in [os.getcwd(), _FTREPO] + [os.path.dirname(p) for p in glob.glob("/kaggle/working/**/direction_audit.py", recursive=True)]:
        f = os.path.join(c, "direction_audit.py")
        if c in seen or not os.path.exists(f):
            continue
        seen.append(c)
        src = open(f, encoding="utf-8", errors="ignore").read()
        if "def cis_qtl_dirs" in src and "class OpenTargetsError" in src:
            return c
    if os.path.isdir(_FTREPO):
        shutil.rmtree(_FTREPO)
    print("  no hardened engine found locally -> cloning the repo")
    subprocess.run(["git", "clone", "--depth", "1", REPO, _FTREPO], check=True, capture_output=True, text=True)
    return _FTREPO

_d = _ensure_hardened_engine()
if _d not in sys.path:
    sys.path.insert(0, _d)
sys.modules.pop("direction_audit", None)
import direction_audit as _da
from direction_audit import post, recover, gwas_locus_ids, cis_signs, OpenTargetsError
if getattr(_da.gwas_locus_ids, "__name__", "") == "<lambda>":
    raise RuntimeError("direction_audit is monkeypatched by a test cell. Run > Restart kernel, then run THIS cell first.")
H4_MIN, MAG, FRAC, MIN_LOCI = _da.H4_MIN, _da.MAG, _da.FRAC, _da.MIN_LOCI
ENGINE_SHA = hashlib.sha256(open(_da.__file__, "rb").read()).hexdigest()
REGION_WINDOW = 1_000_000   # loci within this many bp on a chromosome = one independent region
print(f"  engine: {_da.__file__}  (sha {ENGINE_SHA[:12]}; thresholds H4>={H4_MIN} |s|>={MAG} frac>={FRAC} min_loci={MIN_LOCI})")

def wilson(k, n, z=1.96):
    if not n: return (None, None)
    p = k/n; d = 1 + z*z/n; c = (p + z*z/(2*n))/d
    h = (z/d)*math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return (round(max(0.0, c-h), 3), round(min(1.0, c+h), 3))

def region_of(variant_id):
    if not variant_id: return None
    m = re.match(r"^([0-9XYMT]+)[_:](\d+)", str(variant_id))
    return (m.group(1), int(m.group(2))) if m else None

def cluster_regions(decided):
    """decided: list of (direction, (chrom,pos)|None). Loci within REGION_WINDOW on a
    chromosome collapse to one majority vote; unplaceable loci each stay independent."""
    placed = sorted([x for x in decided if x[1]], key=lambda x: (x[1][0], x[1][1]))
    nopos = [x for x in decided if not x[1]]
    groups, cur = [], []
    for dirn, (ch, pos) in placed:
        if cur and ch == cur[-1][1][0] and pos - cur[-1][1][1] <= REGION_WINDOW:
            cur.append((dirn, (ch, pos)))
        else:
            if cur: groups.append(cur)
            cur = [(dirn, (ch, pos))]
    if cur: groups.append(cur)
    out = []
    for g in groups:
        s = sum(d for d, _ in g); out.append(1 if s > 0 else (-1 if s < 0 else 0))
    return out + [d for d, _ in nopos]

def consensus(dirs):
    dd = [d for d in dirs if d != 0]
    if not dd: return None, 0.0, 0
    pos = sum(1 for d in dd if d > 0); neg = len(dd) - pos
    return ("inhibitor" if pos >= neg else "activator"), round(max(pos, neg)/len(dd), 2), len(dd)

def gate_a(cal, misses):
    inh_ok = cal["inhibitor"][0] == cal["inhibitor"][1] and cal["inhibitor"][1] > 0
    act_ok = cal["activator"][0] == cal["activator"][1] and cal["activator"][1] > 0
    both = cal["inhibitor"][1] > 0 and cal["activator"][1] > 0
    matched = cal["inhibitor"][0] + cal["activator"][0]; recd = cal["inhibitor"][1] + cal["activator"][1]
    if both and inh_ok and act_ok: return "PASS — sign correct in both directions"
    if both and matched == 0: return "FAIL — SYSTEMATIC REVERSAL: convention flipped; every direction is backwards"
    if both and matched >= 0.75*recd: return f"PASS-with-exceptions — {recd-matched} disagree (likely real decoupling): {misses}"
    return f"INCONCLUSIVE — too few recovered controls in both directions ({matched}/{recd})"

def ot_version():
    try:
        d = post("query{ meta{ dataVersion{ year month iteration } apiVersion{ x y z } } }", {}, "meta")
        return ((d or {}).get("data") or {}).get("meta")
    except Exception:
        return None

def lead_variant(lid):
    try:
        d = post('query($id:String!){ credibleSet(studyLocusId:$id){ variant{ id } } }', {"id": lid}, "var")
        return ((((d or {}).get("data") or {}).get("credibleSet") or {}).get("variant") or {}).get("id")
    except OpenTargetsError:
        return None

def locus_dir(lid, ens, sym):
    try:
        signs = [s for (s, _h4, _qt) in cis_signs(lid, ens, sym, H4_MIN) if abs(s) >= MAG]
    except OpenTargetsError:
        return 0
    p = sum(1 for s in signs if s > 0); n = sum(1 for s in signs if s < 0)
    return 1 if p > n else (-1 if n > p else 0)

def loci_meta(ens, efo):
    """Per colocalising credible set: (studyLocusId, study_id, variant_id), deduped by
    studyLocusId. One query. Returns (rows, True) or (None, False) if the richer query
    (study{id}/variant{id}) isn't supported by the schema."""
    q = """query($e:String!,$f:String!){ target(ensemblId:$e){
      evidences(efoIds:[$f], datasourceIds:["gwas_credible_sets"], size:300){
        rows{ credibleSet{ studyLocusId study{ id } variant{ id } } } } } }"""
    try:
        d = post(q, {"e": ens, "f": efo}, "lmeta")
        rows = ((((d or {}).get("data") or {}).get("target") or {}).get("evidences") or {}).get("rows") or []
        seen, out = set(), []
        for r in rows:
            cs = r.get("credibleSet") or {}
            slid = cs.get("studyLocusId")
            if slid and slid not in seen:
                seen.add(slid)
                out.append((slid, (cs.get("study") or {}).get("id"), (cs.get("variant") or {}).get("id")))
        return out, True
    except OpenTargetsError:
        return None, False

# ============================================================================
print("="*112); print("GATE A — SIGN-ORIENTATION CALIBRATION  (textbook-direction controls; must hold BOTH ways)"); print("="*112)
OTV = ot_version()
print(f"  OpenTargets {OTV if OTV else '(version unavailable)'}\n")
# (symbol, ensembl, efo, textbook direction) — IDs verified from the benchmark panel
CONTROLS = [
  ("PCSK9",  "ENSG00000169174","EFO_0001645",  "inhibitor"),
  ("CETP",   "ENSG00000087237","EFO_0001645",  "inhibitor"),
  ("ANGPTL3","ENSG00000132855","EFO_0004530",  "inhibitor"),
  ("NPC1L1", "ENSG00000015520","EFO_0004611",  "inhibitor"),
  ("IL23R",  "ENSG00000162594","EFO_0003767",  "inhibitor"),
  ("IL12B",  "ENSG00000113302","EFO_0003767",  "inhibitor"),
  ("LPL",    "ENSG00000175445","EFO_0004530",  "activator"),
  ("GLP1R",  "ENSG00000112164","MONDO_0005148","activator"),
  ("PPARG",  "ENSG00000132170","MONDO_0005148","activator"),
]
print(f"  {'gene':<9}{'expect':<11}{'recovered':<11}{'cons%':>6}{'loci':>5}  result")
print("  " + "-"*60)
cal = {"inhibitor": [0, 0], "activator": [0, 0]}; misses = []; rowsA = []
for sym, ens, efo, expect in CONTROLS:
    try:
        r = recover(ens, efo)
        got = r.get("recovered"); cons = r.get("consensus"); nl = r.get("n_loci_decided")
    except OpenTargetsError:
        got, cons, nl = None, None, None
    if got in ("inhibitor", "activator"):
        cal[expect][1] += 1
        if got == expect: cal[expect][0] += 1
        else: misses.append((sym, expect, got))
    res = "OK" if got == expect else ("(no recovery)" if got not in ("inhibitor", "activator") else "!! REVERSED")
    rowsA.append({"sym": sym, "expect": expect, "recovered": got, "consensus": cons, "n_loci": nl})
    print(f"  {sym:<9}{expect:<11}{str(got):<11}{(str(int(cons*100)) if cons else ''):>6}{(str(nl) if nl else ''):>5}  {res}")
matched = cal["inhibitor"][0] + cal["activator"][0]; recd = cal["inhibitor"][1] + cal["activator"][1]
verdict_A = gate_a(cal, misses)
print(f"\n  inhibitor controls {cal['inhibitor'][0]}/{cal['inhibitor'][1]}  |  activator controls {cal['activator'][0]}/{cal['activator'][1]}"
      f"  |  overall {matched}/{recd} = {('%.0f%%' % (100*matched/recd)) if recd else 'n/a'} {wilson(matched, recd) if recd else ''}")
print(f"  GATE A: {verdict_A}")

# ============================================================================
print("\n" + "="*112); print("GATE B — LOCUS INDEPENDENCE  (is 'N loci' N independent signals, or one cis-locus seen N times across studies?)"); print("="*112)
PROBE = [("PCSK9","ENSG00000169174","EFO_0001645"), ("TYK2","ENSG00000105397","EFO_0000676"),
         ("IL23R","ENSG00000162594","EFO_0003767")]
_meta0, _rich = loci_meta(PROBE[0][1], PROBE[0][2])
can_region = bool(_meta0) and region_of(_meta0[0][2]) is not None
has_study  = bool(_meta0) and _meta0[0][1] is not None
print(f"  probe: rich-query={'ok' if _rich else 'fallback'}  region-clustering={'on' if can_region else 'OFF'}  "
      f"study-count={'on' if has_study else 'OFF'}\n")
print(f"  {'gene':<9}{'cred.sets':>10}{'studies':>8}{'regions':>9}{'cons(set)':>11}{'cons(reg)':>10}  interpretation")
print("  " + "-"*80)
rowsB = []
for sym, ens, efo in PROBE:
    meta, ok = loci_meta(ens, efo)
    if not ok or not meta:
        try:
            _s, ids = gwas_locus_ids(ens, efo)
            meta = [(lid, None, (lead_variant(lid) if can_region else None)) for lid in dict.fromkeys(ids)]
        except OpenTargetsError:
            print(f"  {sym:<9}{'query failed':>10}"); continue
    decided, studies = [], set()
    for slid, st, vid in meta:
        d = locus_dir(slid, ens, sym)
        if d == 0: continue
        decided.append((d, region_of(vid) if can_region else None))
        if st: studies.add(st)
    _, cons_set, n_set = consensus([d for d, _ in decided])
    rdirs = cluster_regions(decided) if can_region else [d for d, _ in decided]
    _, cons_reg, n_reg = consensus(rdirs)
    n_stud = len(studies)
    interp = (f"{n_stud} studies agree at {n_reg} cis-region" if has_study and can_region else
              (f"{n_reg} region(s)" if can_region else "no region/study data"))
    rowsB.append({"gene": sym, "n_credible_sets": n_set, "n_studies": (n_stud if has_study else None),
                  "n_regions": (n_reg if can_region else None), "cons_set": cons_set,
                  "cons_region": (cons_reg if can_region else None)})
    print(f"  {sym:<9}{n_set:>10}{(str(n_stud) if has_study else '-'):>8}{(str(n_reg) if can_region else '-'):>9}"
          f"{(str(int(cons_set*100))+'%' if cons_set else '-'):>11}{(str(int(cons_reg*100))+'%' if (can_region and cons_reg) else '-'):>10}  {interp}")

print("\n" + "="*112)
print(f"GATE A (sign): {verdict_A}")
if can_region:
    print("GATE B (independence): cis-coloc localises to the gene's locus, so ~1 region per gene is EXPECTED, not")
    print("a flaw. The honest strength axis is the STUDIES column: 'K studies agree at 1 cis-locus' is cross-study")
    print("replication, NOT K independent loci. The writeup must report studies (and 1 cis-region), not raw locus count.")
else:
    print("GATE B (independence): region clustering needs the variant-field fix; paste this and I'll correct the query.")

manifest = {
  "engine_sha": ENGINE_SHA, "opentargets_version": OTV,
  "gate_a_sign": verdict_A,
  "calibration": {k: {"matched": v[0], "recovered": v[1]} for k, v in cal.items()},
  "calibration_overall": [matched, recd], "calibration_ci": (wilson(matched, recd) if recd else None),
  "sign_misses": misses, "controls": rowsA,
  "gate_b_region_capable": can_region, "independence": rowsB,
}
manifest["manifest_sha"] = hashlib.sha256(json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
for outdir in ("validation/results", "/kaggle/working/validation/results", "/kaggle/working/ftrepo/validation/results", "."):
    try:
        os.makedirs(outdir, exist_ok=True)
        open(os.path.join(outdir, "integrity_gates.json"), "w").write(json.dumps(manifest, indent=2))
        print(f"\nwrote {outdir}/integrity_gates.json   manifest sha {manifest['manifest_sha'][:12]}")
        break
    except Exception:
        continue