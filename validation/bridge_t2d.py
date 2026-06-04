import hashlib
import json
import time

import requests

OT = "https://api.platform.opentargets.org/api/v4/graphql"
H4_MIN  = 0.8     # colocalisation posterior floor
MAG     = 0.8     # a coloc votes only if |betaRatioSignAverage| >= MAG (drops the ~0 ambiguous ones)
FRAC    = 0.8     # >= this fraction of DECIDED loci must agree to call a direction
MIN_LOCI = 3      # need at least this many decided loci to vouch
GWAS = "gwas_credible_sets"
GENETIC_SOURCES = {"gwas_credible_sets","gene_burden","eva","eva_somatic","clingen",
                   "genomics_england","orphanet","gene2phenotype","uniprot_variants",
                   "uniprot_literature","ot_genetics_portal"}   # clinical_precedence excluded (circular)
MECH = {"inhibitor":"inhibitor","antagonist":"inhibitor","blocker":"inhibitor","degrader":"inhibitor",
        "activator":"activator","agonist":"activator","potentiator":"activator"}

def post(q, variables=None, label="", retries=3):
    for a in range(retries):
        try:
            r = requests.post(OT, json={"query": q, "variables": variables or {}}, timeout=60)
        except requests.RequestException as e:
            print(f"  !! net {label}: {e}"); return None
        if r.status_code == 200:
            d = r.json()
            if "errors" in d: print(f"  !! gql {label}: {json.dumps(d['errors'])[:240]}")
            return d
        if r.status_code in (429, 502, 503): time.sleep(2 ** a); continue
        print(f"  !! HTTP {r.status_code} {label}"); return None
    return None

def gwas_locus_ids(ensembl, efo, cap=40):   # cap for the validation sweep
    q = """query($e:String!,$f:String!){ target(ensemblId:$e){ approvedSymbol
      evidences(efoIds:[$f], datasourceIds:["gwas_credible_sets"], size:300){
        rows{ credibleSet { studyLocusId } } } } }"""
    d = post(q, {"e": ensembl, "f": efo}, label="gwas-loci")
    tgt = ((d or {}).get("data") or {}).get("target") or {}
    rows = ((tgt.get("evidences") or {}).get("rows")) or []
    ids = [(r.get("credibleSet") or {}).get("studyLocusId") for r in rows]
    return tgt.get("approvedSymbol"), [i for i in dict.fromkeys(ids) if i][:cap]

COLOC_Q = """query($id:String!){ credibleSet(studyLocusId:$id){ studyLocusId
  colocalisation(studyTypes:[eqtl, pqtl]){ rows {
    betaRatioSignAverage h4 rightStudyType
    otherStudyLocus { qtlGeneId study { target { id approvedSymbol } } } } } } }"""

def cis_signs(studyLocusId, ensembl, symbol, h4_min):
    """Strong directional cis colocs as [(sign,h4,qtltype)]: gene IS the target (qtlGeneId or
    study.target, Ensembl/symbol), h4>=floor, sign present. Excludes trans (PLA2G7-type)."""
    d = post(COLOC_Q, {"id": studyLocusId}, label=f"coloc:{studyLocusId[:10]}")
    cs = ((d or {}).get("data") or {}).get("credibleSet") or {}
    rows = ((cs.get("colocalisation") or {}).get("rows")) or []
    sym_u = symbol.upper() if symbol else None
    out = []
    for r in rows:
        o = r.get("otherStudyLocus") or {}
        tgt = ((o.get("study") or {}).get("target") or {})
        is_cis = (o.get("qtlGeneId") == ensembl) or (tgt.get("id") == ensembl) \
                 or (sym_u and (tgt.get("approvedSymbol") or "").upper() == sym_u)
        h4 = r.get("h4"); sign = r.get("betaRatioSignAverage")
        if not is_cis or sign is None: continue
        if isinstance(h4, (int, float)) and h4 < h4_min: continue
        qt = "pqtl" if "pqtl" in str(r.get("rightStudyType")).lower() else \
             ("eqtl" if "eqtl" in str(r.get("rightStudyType")).lower() else "other")
        out.append((sign, h4, qt))
    return out

def recover(ensembl, efo, h4_min=H4_MIN):
    """LOCUS-CONSENSUS recovery. Each GWAS locus -> +1/-1 from its strong cis colocs
    (|sign|>=MAG); a direction is called iff >=FRAC of >=MIN_LOCI decided loci agree."""
    sym, ids = gwas_locus_ids(ensembl, efo)
    loci_dir = {}; cpos = cneg = 0
    for lid in ids:
        strong = [s for (s, _, _) in cis_signs(lid, ensembl, sym, h4_min) if abs(s) >= MAG]
        p = sum(1 for s in strong if s > 0); n = sum(1 for s in strong if s < 0)
        cpos += p; cneg += n
        loci_dir[lid] = (1 if p > n else (-1 if n > p else 0))
    decided = [d for d in loci_dir.values() if d != 0]
    lpos = sum(1 for d in decided if d > 0); lneg = len(decided) - lpos
    frac = (max(lpos, lneg) / len(decided)) if decided else 0.0
    drug = None
    if len(decided) >= MIN_LOCI and frac >= FRAC:
        drug = "inhibitor" if lpos >= lneg else "activator"
    # strong-coloc-majority as a cross-check (the simpler, pseudo-replicating rule)
    cx = None
    if cpos + cneg:
        cf = max(cpos, cneg) / (cpos + cneg)
        cx = ("inhibitor" if cpos >= cneg else "activator") if cf >= FRAC else "conflicted"
    return {"symbol": sym, "recovered": drug, "consensus": round(frac, 2),
            "n_loci_total": len(loci_dir), "n_loci_decided": len(decided),
            "loci_pos": lpos, "loci_neg": lneg, "coloc_pos": cpos, "coloc_neg": cneg,
            "crosscheck_strong_majority": cx, "loci": sorted([l for l, d in loci_dir.items() if d != 0])}

# ---- curated label (genetic only) so labelled targets bypass recovery ------
DT = {"lof":"down","loss_of_function":"down","gof":"up","gain_of_function":"up","decrease":"down","increase":"up"}
def desired_from_label(dot, dotr):
    t = DT.get(str(dot).strip().lower()) if dot else None
    tr = str(dotr).strip().lower() if dotr else None
    if not t or tr not in ("risk", "protect"): return None
    if (t, tr) in (("up","risk"),("down","protect")): return "inhibitor"
    if (t, tr) in (("down","risk"),("up","protect")): return "activator"
    return None

def labelled_direction(ensembl, efo):
    q = """query($e:String!,$f:String!){ target(ensemblId:$e){
      evidences(efoIds:[$f], size:1000){ rows{ datasourceId directionOnTarget directionOnTrait } } } }"""
    d = post(q, {"e": ensembl, "f": efo}, label="labels")
    rows = (((((d or {}).get("data") or {}).get("target") or {}).get("evidences") or {}).get("rows")) or []
    votes = {}
    for r in rows:
        ds = r.get("datasourceId")
        if ds not in GENETIC_SOURCES or ds == GWAS: continue
        des = desired_from_label(r.get("directionOnTarget"), r.get("directionOnTrait"))
        if des: votes[des] = votes.get(des, 0) + 1
    if not votes: return None, 0
    return max(votes, key=votes.get), sum(votes.values())

def verdict(ensembl, efo, mechanism, h4_min=H4_MIN):
    claim = MECH.get(str(mechanism).strip().lower())
    if claim is None: raise ValueError(f"bad mechanism {mechanism!r}")
    lab, nlab = labelled_direction(ensembl, efo)
    if lab:
        res = {"verdict": "LABEL_CONCORDANT" if lab == claim else "LABEL_DISCORDANT",
               "vouches": lab == claim, "source": "curated-label", "confidence": "high",
               "direction": lab, "detail": f"{nlab} curated genetic label rows -> {lab}"}
    else:
        rec = recover(ensembl, efo, h4_min)
        pct = int(round(rec["consensus"] * 100))
        if rec["n_loci_decided"] < MIN_LOCI:
            res = {"verdict": "INSUFFICIENT_DIRECTION", "vouches": False, "source": "none",
                   "confidence": "n/a", "direction": None, "recovery": rec,
                   "detail": f"only {rec['n_loci_decided']} decided cis-coloc locus/loci (<{MIN_LOCI})",
                   "falsifier": "too few cis molQTL coloc loci over H4>=%.2f/|sign|>=%.1f; an MR or rare-variant burden test is the entry ticket." % (h4_min, MAG)}
        elif rec["recovered"] is None:
            res = {"verdict": "RECOVERY_CONFLICTED", "vouches": False, "source": "coloc",
                   "confidence": "n/a", "direction": None, "recovery": rec,
                   "detail": f"locus consensus {pct}% over {rec['n_loci_decided']} loci is below {int(FRAC*100)}%",
                   "falsifier": "loci disagree on direction; a trait-specific MR would adjudicate."}
        else:
            concord = rec["recovered"] == claim
            res = {"verdict": "RECOVERED_CONCORDANT" if concord else "RECOVERED_DISCORDANT",
                   "vouches": concord, "source": "coloc-recovered",
                   "confidence": f"moderate (coloc-derived; {pct}% locus consensus over {rec['n_loci_decided']} loci)",
                   "direction": rec["recovered"], "recovery": rec,
                   "falsifier": "coloc-inferred direction; a rare-variant burden test or MR would upgrade it to direct."}
    payload = {"rule": "R3+recovery/locus-consensus", "ensembl": ensembl, "efo": efo, "mechanism": claim,
               "verdict": res["verdict"], "vouches": res["vouches"], "direction": res.get("direction"),
               "loci": res.get("recovery", {}).get("loci", [])}
    res["sha256"] = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    res["ensembl"], res["efo"], res["claim"] = ensembl, efo, claim
    return res


# ============================================================================
# GLYCEMIC -> T2D BRIDGE  -  same dual-endpoint test that validated LDL->CAD
#   (1) CONCORDANCE: targets resolving on BOTH type-2-diabetes and a glycemic
#       biomarker (HbA1c / fasting glucose) -- do the directions agree?
#   (2) RESCUE: T2D refusals re-run on the glycemic endpoint. DPP4 is the test.
# Sign is identity (+1): higher glucose/HbA1c -> higher T2D risk (definitional).
# Run on Kaggle, internet ON. ~6-10 min.
# ============================================================================
T2D = "MONDO_0005148"
GLYCEMIC = [("EFO_0004541","HbA1c",+1), ("EFO_0004468","fasting glucose",+1)]   # both: higher -> higher T2D risk
def flip(d): return {"inhibitor":"activator","activator":"inhibitor"}.get(d, d)

def glyc_endorsement(ens, h4_min=H4_MIN):
    """recovery on glycemic biomarkers, mapped to a T2D endorsement; require agreement if >1 speaks."""
    endos=[]
    for bm_efo,bm_name,sign in GLYCEMIC:
        r = recover(ens, bm_efo, h4_min); raw = r.get("recovered")
        if raw is None: continue
        endos.append((bm_name, raw if sign>0 else flip(raw), r))
    if not endos: return None, None, {}
    dirs={e for _,e,_ in endos}
    if len(dirs)>1: return "conflict", None, endos          # biomarkers disagree
    endo=dirs.pop(); nm,_,r = endos[0]
    return endo, nm, r

# symbol, ensembl, mechanism, note
PANEL = [
  ("GLP1R","ENSG00000112164","activator","semaglutide (agonist)"),
  ("PPARG","ENSG00000132170","activator","thiazolidinediones (agonist)"),
  ("DPP4", "ENSG00000197635","inhibitor","sitagliptin; more DPP4->more glucose"),
  ("SLC5A2","ENSG00000140675","inhibitor","empagliflozin (SGLT2); more->more glucose"),
  ("INSR", "ENSG00000171105","activator","insulin signalling; more->less glucose"),
  ("GIPR", "ENSG00000010310","contested","agonist-vs-antagonist debate (control)"),
]
REFUSED = {"RECOVERY_CONFLICTED","INSUFFICIENT_DIRECTION"}
def run_t2d(panel):
    rows=[]
    print(f"{'target':<8}{'mech':<10}{'T2D verdict':<24}{'T2Ddir':<10}{'glycdir':<10}{'loci':>5}  status")
    print("-"*92)
    for sym,ens,mech,note in panel:
        vt = verdict(ens, T2D, mech if mech in ("inhibitor","activator") else "inhibitor")
        t_dir = vt.get("direction")
        g_dir, bm, r = glyc_endorsement(ens)
        nloci = r.get("n_loci_decided")
        if g_dir=="conflict": g_dir=None
        if t_dir and g_dir:
            status = "CONCORDANT" if t_dir==g_dir else "DISCORDANT(bridge!)"
        elif vt["verdict"] in REFUSED and g_dir:
            status = f"RESCUED->{g_dir} (via {bm})"
        elif vt["verdict"] in REFUSED:
            status = "refused both"
        else:
            status = "direct only"
        rsym=r.get("symbol"); sflag="" if (rsym is None or rsym==sym) else f" !{rsym}"
        rows.append({"sym":sym,"ens":ens,"mech":mech,"t_verdict":vt["verdict"],"t_dir":t_dir,
                     "g_dir":g_dir,"bm":bm,"nloci":nloci,"status":status,
                     "sha":hashlib.sha256(f"{ens}|{T2D}|glyc|{bm}|{g_dir}".encode()).hexdigest()[:12]})
        print(f"{sym:<8}{mech:<10}{vt['verdict']:<24}{str(t_dir):<10}{str(g_dir):<10}{str(nloci):>5}  {status}{sflag}")
    return rows

print("="*92); print("GLYCEMIC -> T2D BRIDGE  -  T2D endpoint vs HbA1c/glucose"); print("="*92)
R = run_t2d(PANEL)
both=[r for r in R if r["t_dir"] and r["g_dir"]]; conc=[r for r in both if r["t_dir"]==r["g_dir"]]
disc=[r for r in both if r["t_dir"]!=r["g_dir"]]
rescue=[r for r in R if r["t_verdict"] in REFUSED and r["g_dir"] and r["mech"] in ("inhibitor","activator")]
resc_ok=[r for r in rescue if r["g_dir"]==r["mech"]]
def pct(a,b): return f"{100*a/b:.0f}%" if b else "n/a"
print("\n"+"="*92); print("HEADLINE"); print("="*92)
print(f"  (1) BRIDGE CONCORDANCE (both endpoints speak): {len(conc)}/{len(both)} = {pct(len(conc),len(both))}")
if disc: print(f"      DISCORDANT (inspect): {[(d['sym'],d['t_dir'],d['g_dir']) for d in disc]}")
print(f"  (2) T2D refusals rescued via glycemic: {len(rescue)}   matching mechanism: {len(resc_ok)}/{len(rescue)} = {pct(len(resc_ok),len(rescue))}")
print(f"      DPP4 specifically: {[r['status'] for r in R if r['sym']=='DPP4'][0]}")
print("\nPROVENANCE (glycemic-route verdicts):")
for r in R:
    if r["g_dir"]: print(f"  {r['sym']:<8}{r['status']:<30} glyc-dir={r['g_dir']:<10} loci={r['nloci']} sha={r['sha']}")
print("\nREAD: a clean concordance + a DPP4 rescue means glycemic->T2D is a second validated bridge,")
print("ready to fold into the panel (DPP4 would move it 58%->67%). If HbA1c/glucose EFOs miss, rows")
print("show 0 loci -> soft refusal, never a wrong rescue. The sign is definitional, so it cannot mislead.")
