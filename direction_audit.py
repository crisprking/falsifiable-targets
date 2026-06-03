"""
direction_audit — recover drug-target direction-of-effect from Open Targets.

Open Targets harmonises direction-of-effect (directionOnTarget / directionOnTrait)
ONLY for intrinsically-directional genetic evidence (rare-variant burden, curated
clinical). Common-variant GWAS — the bulk of genetic support — carries no direction.
This module recovers a direction for GWAS-only claims from the COLOCALISATION layer:
where a disease credible set colocalises with a cis eQTL/pQTL for the target,
betaRatioSignAverage gives sign>0 -> inhibitor, sign<0 -> activator. A locus-consensus
rule vouches only when >=FRAC of >=MIN_LOCI decided loci agree; else it refuses.

Public data only (Open Targets GraphQL). Every verdict is content-addressed (sha256).

Key entry point:  verdict(ensembl, efo, mechanism) -> dict
See README.md for the validated benchmark and the abundance!=function blind spot.
"""
import json, time, hashlib
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

def gwas_locus_ids(ensembl, efo, cap=80):
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


# ---------------------------------------------------------------------------
# v1.6  biomarker-bridge fallback
# When the disease endpoint is direction-mute (INSUFFICIENT_DIRECTION), recover
# direction from a causally-signed upstream biomarker that carries dense GWAS,
# and map the sign through. Each bridge is added ONLY after passing a dual-
# endpoint concordance check (see validation/). Fires ONLY on INSUFFICIENT --
# it NEVER overrides a RECOVERY_CONFLICTED refusal, because a genuine split is
# information, not a gap (the contested control GIPR established this rule).
# ---------------------------------------------------------------------------
def _flip(d): return {"inhibitor": "activator", "activator": "inhibitor"}.get(d, d)

# disease EFO -> [(biomarker EFO, name, risk_sign)].  +1 = higher biomarker -> higher disease risk.
BIOMARKER_BRIDGES = {
    "EFO_0001645": [("EFO_0004611", "LDL cholesterol", +1)],   # coronary artery disease (validated 4/4 concordant)
}

def verdict_with_fallback(ensembl, efo, mechanism, h4_min=H4_MIN):
    """verdict(), plus a biomarker-bridge rescue when the disease endpoint is INSUFFICIENT.
    Never overrides RECOVERY_CONFLICTED. Returns a confidence-graded *_VIA_BIOMARKER verdict
    with its own content-addressed sha when a bridge fires; otherwise returns verdict() unchanged."""
    v = verdict(ensembl, efo, mechanism, h4_min)
    if v["verdict"] != "INSUFFICIENT_DIRECTION" or efo not in BIOMARKER_BRIDGES:
        return v
    claim = MECH.get(str(mechanism).strip().lower())
    endos = []
    for bm_efo, bm_name, sign in BIOMARKER_BRIDGES[efo]:
        r = recover(ensembl, bm_efo, h4_min); raw = r.get("recovered")
        if raw is None: continue
        endos.append((bm_name, raw if sign > 0 else _flip(raw), r))
    if not endos: return v
    dirs = {e for _, e, _ in endos}
    if len(dirs) > 1: return v                       # biomarkers disagree -> keep the honest refusal
    endo = dirs.pop(); bm_name, _, r = endos[0]; concord = (endo == claim)
    out = dict(v)
    out.update(verdict="RECOVERED_CONCORDANT_VIA_BIOMARKER" if concord else "RECOVERED_DISCORDANT_VIA_BIOMARKER",
               vouches=concord, source=f"biomarker:{bm_name}", direction=endo,
               confidence=f"moderate-minus (biomarker-mediated via {bm_name}; {int(round(r['consensus']*100))}% over {r['n_loci_decided']} loci)",
               recovery=r,
               falsifier=f"direction recovered via {bm_name} (causal proxy); an MR of {bm_name}->disease, or a disease-endpoint coloc, upgrades it to direct.")
    out["sha256"] = hashlib.sha256(f"{ensembl}|{efo}|biomarker|{bm_name}|{endo}".encode()).hexdigest()
    return out
