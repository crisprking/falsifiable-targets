# === FIGURE · the boundary surface (reads YOUR saved scope_rule_grid.json) =================
# Self-contained. Reads the JSON cell 76 already saved (no re-query, no upstream cell needed).
# Auto-computes every number from the file — nothing hardcoded. Fails loud if the file is absent.
# Renders the centerpiece: coverage climbs with trait scope only by driving accuracy below the
# always-inhibitor baseline, and the pooled scorer collapses to an "always-activator" caller.
#
# Schema it expects (exactly what cell 76 writes):
#   { scope: { rule: {cov, commit, N, acc, wilson:[pt,lo,hi], bal_acc, mcc, maj,
#                      rec:{inhibitor,activator}} } }   scopes: exact/descendants/pooled

import json, os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---- locate the file (Kaggle working dir first, then cwd) -----------------------------------
CANDIDATES = ["/kaggle/working/scope_rule_grid.json", "scope_rule_grid.json",
              os.path.join(os.getcwd(), "scope_rule_grid.json")]
PATH = next((p for p in CANDIDATES if os.path.exists(p)), None)
if PATH is None:
    raise FileNotFoundError(
        "scope_rule_grid.json not found. Run the AUDIT · TRAIT-SCOPE × DECISION-RULE cell "
        "first (it saves the file), then re-run this figure cell.")
GRID = json.load(open(PATH))
print(f"loaded {PATH}")

SCOPES = [s for s in ["exact", "descendants", "pooled"] if s in GRID]
RULE   = "tier"  # the tool's own rule; rules differ by <2pts (state this in the caption)
SCOPE_LABEL = {"exact": "exact\nindication", "descendants": "indication\n+ descendants",
               "pooled": "pooled\n(all traits)"}

def cell(scope, rule=RULE):
    d = GRID[scope]
    return d.get(rule) or next(iter(d.values()))  # fall back if 'tier' absent

cov  = [cell(s)["cov"] * 100 for s in SCOPES]
acc  = [cell(s)["acc"] * 100 for s in SCOPES]
maj  = [cell(s)["maj"] * 100 for s in SCOPES]
mcc  = [cell(s)["mcc"] for s in SCOPES]
ncom = [cell(s)["commit"] for s in SCOPES]
ntot = [cell(s)["N"] for s in SCOPES]
# Wilson half-widths for accuracy error bars (lo/hi are absolute bounds in the file)
wlo  = [max(0.0, acc[i] - cell(s)["wilson"][1] * 100) for i, s in enumerate(SCOPES)]
whi  = [max(0.0, cell(s)["wilson"][2] * 100 - acc[i]) for i, s in enumerate(SCOPES)]

# pooled per-class recall for panel (b) — read the flattering corner (pooled x majority),
# the one the paper critiques; that is where {inhibitor:0.09, activator:1.0} lives.
_pdict  = GRID.get("pooled", GRID[SCOPES[-1]])
pooled  = _pdict.get("majority") or _pdict.get(RULE) or next(iter(_pdict.values()))
rec     = pooled.get("rec", {})
rec_inh = (rec.get("inhibitor") or 0.0) * 100
rec_act = (rec.get("activator") or 0.0) * 100

# ---- palette --------------------------------------------------------------------------------
C_COV, C_ACC, C_BASE = "#4C72B0", "#C44E52", "#8C8C8C"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.dpi": 140})
fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.6),
                               gridspec_kw={"width_ratios": [2.05, 1]})

# ---- Panel A: coverage bars + accuracy line + per-scope baseline ----------------------------
x = list(range(len(SCOPES)))
axA.bar(x, cov, width=0.56, color=C_COV, alpha=0.85, label="coverage (% committed)", zorder=2)
for xi, c, n, t in zip(x, cov, ncom, ntot):
    axA.text(xi, c + 1.2, f"{c:.1f}%\n({n}/{t})", ha="center", va="bottom",
             fontsize=9, color=C_COV)

axA.errorbar(x, acc, yerr=[wlo, whi], fmt="o-", color=C_ACC, lw=2, ms=8, capsize=4,
             label="accuracy on committed (Wilson 95%)", zorder=4)
axA.plot(x, maj, "s--", color=C_BASE, lw=1.6, ms=7, label="always-inhibitor baseline", zorder=3)

# shade the "below baseline" gap at the pooled scope
if "pooled" in SCOPES:
    pi = SCOPES.index("pooled")
    axA.annotate("", xy=(pi + 0.30, acc[pi]), xytext=(pi + 0.30, maj[pi]),
                 arrowprops=dict(arrowstyle="<->", color="#444", lw=1.3))
    axA.text(pi + 0.36, (acc[pi] + maj[pi]) / 2,
             f"{maj[pi]-acc[pi]:.0f} pts\nbelow chance", fontsize=8.5, va="center", color="#444")

for xi, m in zip(x, mcc):
    axA.text(xi, -7.5, f"MCC {m:+.2f}", ha="center", va="top", fontsize=8.5,
             color=("#2a7" if m > 0.20 else "#999"))

axA.set_xticks(x); axA.set_xticklabels([SCOPE_LABEL[s] for s in SCOPES])
axA.set_ylabel("percent"); axA.set_ylim(-12, 100)
axA.set_title("Coverage rises with trait scope only by inverting the sign", fontsize=11.5, pad=10)
axA.legend(loc="upper left", fontsize=8.6, framealpha=0.9)
axA.text(0.5, -0.255, "trait scope (decision rule = tier; strict/majority differ by <2 pts)",
         transform=axA.transAxes, ha="center", fontsize=8.3, color="#666")

# ---- Panel B: pooled per-class recall — the "always activator" signature --------------------
bars = axB.bar(["true\ninhibitors", "true\nactivators"], [rec_inh, rec_act],
               color=[C_ACC, "#55A868"], width=0.6, zorder=2)
for b, v in zip(bars, [rec_inh, rec_act]):
    axB.text(b.get_x() + b.get_width() / 2, v + 2, f"{v:.0f}%", ha="center", va="bottom",
             fontsize=10, fontweight="bold")
axB.axhline(50, color=C_BASE, ls=":", lw=1.2)
axB.text(1.45, 52, "constant-predictor\nfloor (bal-acc 50%)", fontsize=7.6, color="#666",
         ha="right", va="bottom")
axB.set_ylim(0, 112); axB.set_ylabel("recall when committed (%)")
axB.set_title("Pooled scorer is\n'always say activator'", fontsize=11, pad=10)

fig.suptitle("Direction-of-effect from human genetics is a biomarker-proximal signal",
             fontsize=13, fontweight="bold", y=1.02)
fig.tight_layout(rect=[0, 0.02, 1, 0.98])

OUT = "/kaggle/working/fig_boundary_surface.png" if os.path.isdir("/kaggle/working") \
      else "fig_boundary_surface.png"
fig.savefig(OUT, bbox_inches="tight", dpi=200)
print(f"saved {OUT}")
print(f"  panel A scopes={SCOPES}  coverage={[round(c,1) for c in cov]}  "
      f"accuracy={[round(a,1) for a in acc]}  baseline={[round(m,1) for m in maj]}")
print(f"  panel B pooled recall: inhibitor={rec_inh:.0f}%  activator={rec_act:.0f}%")
