#!/usr/bin/env python3
"""Supplementary tables and figures.

Tables are copies of the authoritative files in results/, renamed to their S-numbers so
a reader can find them from the text, with a provenance line recording which analysis
script produced each. Nothing is recomputed for a table: a supplementary table that
disagreed with the result file it came from would be worse than no table.

Figures follow the same style rules as the main ones (see make_figures.py): Liberation
Sans at 6-7 pt, no titles or explanatory text inside the panel, panel letters outside.
"""
import json
import shutil
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
P = REPO_ROOT
R, S = P / "results", P / "supplementary"
S.mkdir(exist_ok=True)

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#111111", "#4a4a4a", "#8a8a8a"
GRID, SPINE, SURF = "#e8e8e4", "#4a4a4a", "#ffffff"
MM = 1 / 25.4
W1, W2 = 89 * MM, 183 * MM
FAM = "Liberation Sans" if any(f.name == "Liberation Sans" for f in fm.fontManager.ttflist) \
      else "DejaVu Sans"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FAM],
    "font.size": 6.5, "axes.labelsize": 7, "xtick.labelsize": 6, "ytick.labelsize": 6,
    "legend.fontsize": 6, "axes.edgecolor": SPINE, "axes.linewidth": 0.6,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": SPINE, "ytick.color": SPINE,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.45,
    "axes.axisbelow": True, "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF, "savefig.bbox": "tight", "savefig.pad_inches": 0.015,
    "lines.linewidth": 0.9, "lines.markersize": 2.8, "legend.frameon": False,
    "legend.handlelength": 1.2, "legend.handletextpad": 0.45, "legend.borderpad": 0.0,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, letter, dx=-0.17, dy=1.03):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8, fontweight="bold",
            va="bottom", ha="left", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(S / f"{name}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"  {name}.pdf / .png")


# ================================================================ TABLES
TABLES = [
    ("S1", "primary_order.tsv",
     "FunAqua primary result: all 100 tested orders", "analysis_primary.py"),
    ("S2", "draw_stability.tsv",
     "Per-order significance across 200 rarefaction draws", "analysis_draw_stability.py"),
    ("S3", "matched_fungi.tsv",
     "Fungal orders on the 204 matched pairs", "analysis_euk_matched.py"),
    ("S4", "matched_euk.tsv",
     "Non-fungal eukaryote orders on the same 204 pairs", "analysis_euk_matched.py"),
    ("S5", "kingdom_consistent.tsv",
     "Per-kingdom metrics on one consistent footing", "analysis_kingdom_consistent.py"),
    ("S6", "levins_per_order.tsv",
     "Levins' B and the preference index, per order", "analysis_levins.py"),
    ("S7", "neon_order.tsv",
     "NEON per-order result: all 103 tested orders", "analysis_neon_validation.py"),
    ("S8", "neon_vs_funaqua.tsv",
     "The 27 orders comparable between FunAqua and NEON", "analysis_neon_validation.py"),
    ("S9", "neon_pairs.tsv",
     "The 149 NEON water/benthic pairs", "analysis_neon_validation.py"),
    ("S10", "zoosporic.tsv",
     "Flagellum trait assignment and preference, both datasets", "analysis_zoosporic.py"),
    ("S11", "neon_depth_probe.tsv",
     "Agreement under capping of NEON depth", "exploratory_neon_depth.py"),
    ("S12", "euk_robustness.tsv",
     "Per-order values behind the five robustness checks", "analysis_euk_robustness.py"),
]

print("=== supplementary tables")
prov = []
for sid, src, title, script in TABLES:
    p = R / src
    if not p.exists():
        print(f"  MISSING {src} -> {sid} skipped")
        continue
    d = pd.read_csv(p, sep="\t")
    out = S / f"Table_{sid}.tsv"
    shutil.copyfile(p, out)
    prov.append(dict(table=sid, title=title, rows=len(d), cols=d.shape[1],
                     source=f"results/{src}", script=script))
    print(f"  Table_{sid}.tsv  {len(d):4d} x {d.shape[1]:2d}  {title}")

# The robustness checks are a summary, not a per-order table, so they get their own.
rb = json.load((R / "euk_robustness.json").open())
rows = [
    dict(check="1 draw stability", statistic="reversals of the ordering",
         value=f"{rb['check1']['reversals']}/{rb['check1']['n_draws']}",
         threshold="<= 5 % of draws", passed=rb["check1"]["pass"]),
    dict(check="2 matched on power", statistic="stratified rank test p",
         value=f"{rb['check2']['stratified_p']:.4g}",
         threshold=f"fungi lower in all bins ({rb['check2']['bins_fungi_lower']}/{rb['check2']['bins']})",
         passed=rb["check2"]["pass"]),
    dict(check="3 each kingdom", statistic="kingdoms with fungi lower",
         value=f"{rb['check3']['n_lower']}/{len(rb['check3']['kingdoms'])}",
         threshold="all, significant against >= 2", passed=rb["check3"]["pass"]),
    dict(check="4 neighbouring ranks", statistic="ranks with fungi lower",
         value=f"{sum(r['fungi_lower'] for r in rb['check4']['ranks'])}/{len(rb['check4']['ranks'])}",
         threshold="majority of ranks", passed=rb["check4"]["pass"]),
    dict(check="5 OTUs per order", statistic="partial Spearman p",
         value=f"{rb['check5']['partial_p']:.4g}",
         threshold="gap survives adjustment", passed=rb["check5"]["pass"]),
]
pd.DataFrame(rows).to_csv(S / "Table_S13.tsv", sep="\t", index=False)
prov.append(dict(table="S13", title="The five pre-registered robustness checks",
                 rows=len(rows), cols=5, source="results/euk_robustness.json",
                 script="analysis_euk_robustness.py"))
print(f"  Table_S13.tsv  {len(rows):4d} x  5  The five pre-registered robustness checks")
pd.DataFrame(prov).to_csv(S / "table_provenance.tsv", sep="\t", index=False)

# ================================================================ FIGURE S1
# The two metrics agree. Currently this is two correlation coefficients in the text;
# as a panel it is the evidence that the result is not an artefact of our metric.
lp = pd.read_csv(R / "levins_per_order.tsv", sep="\t")
lv = json.load((R / "levins.json").open())
fig, axes = plt.subplots(1, 2, figsize=(W2 * 0.72, 52 * MM), sharey=True)
fig.subplots_adjust(wspace=0.08)
for ax, (grp, col, key) in zip(axes, (("fungi", BLUE, "agreement_fungi"),
                                      ("other eukaryotes", ORANGE, "agreement_other"))):
    g = lp[lp.group == grp]
    ax.scatter(g.B, g.effect, s=11, facecolor=col, edgecolor="white", linewidth=0.4)
    a = lv[key]
    ax.text(0.97, 0.95, f"$\\rho$ = {a['rho']:+.2f}\n$n$ = {a['n']}",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.3, color=INK2)
    ax.set_xlabel("Levins' $B$ (2 compartments)")
    ax.set_xlim(1.0, 2.02)
axes[0].set_ylabel("|preference − 0.5|")
axes[0].set_ylim(0, 0.52)
panel(axes[0], "a", dx=-0.19); panel(axes[1], "b", dx=-0.06)
save(fig, "Figure_S1_metric_agreement")

# ================================================================ FIGURE S2
# The flagellum trait test, which replicated in NEON and failed in FunAqua.
z = pd.read_csv(R / "zoosporic.tsv", sep="\t")
zj = json.load((R / "zoosporic.json").open())
fig, axes = plt.subplots(1, 2, figsize=(W2 * 0.62, 52 * MM), sharey=True)
fig.subplots_adjust(wspace=0.10)
for ax, (ds, lab) in zip(axes, (("funaqua", "FunAqua"), ("neon", "NEON"))):
    sub = z[(z.dataset == ds) & (z.trait.isin(["flagellated", "aflagellate"]))]
    groups = [("flagellated", BLUE), ("aflagellate", MUTED)]
    for i, (tr, col) in enumerate(groups):
        v = sub.loc[sub.trait == tr, "water_fraction"].dropna()
        xj = np.random.default_rng(7).normal(i, 0.055, len(v))
        ax.scatter(xj, v, s=9, facecolor=col, edgecolor="white", linewidth=0.35, zorder=3)
        ax.plot([i - 0.26, i + 0.26], [v.median()] * 2, color=INK, linewidth=1.1, zorder=4)
    ax.axhline(0.5, color=MUTED, linewidth=0.5, dashes=(2.5, 2))
    ax.set_xticks([0, 1]); ax.set_xticklabels(["flagellated", "aflagellate"])
    ax.set_xlim(-0.5, 1.5); ax.set_ylim(0, 1)
    p = zj[ds].get("mw_p")
    ax.text(0.5, 0.97, f"{lab}\n$p$ = {p:.2g}", transform=ax.transAxes,
            ha="center", va="top", fontsize=6.3, color=INK2)
axes[0].set_ylabel("water preference index")
panel(axes[0], "a", dx=-0.26); panel(axes[1], "b", dx=-0.08)
save(fig, "Figure_S2_flagellum_trait")

# ================================================================ FIGURE S3
# Rank comparability: why metaPR2 could not be used. Recomputed here from the two
# taxonomies so the panel cannot drift from the claim it supports.
print("\n=== rank granularity audit for Figure S3")
UN = {".", "", "unidentified", "unassigned", "NA"}
FD = FUNAQUA_DIR
ft = pd.read_csv(FD / "FunAqua_fungi_tax.csv", index_col=0, dtype=str, keep_default_na=False)
et = pd.read_csv(FD / "FunAqua_euk_tax.csv", index_col=0, dtype=str, keep_default_na=False)
et = et[et.Kingdom.str.strip() != "Fungi"]
fa_rows = []
for rank in ("Phylum", "Class", "Order", "Family"):
    fv = ft[rank].str.strip(); fv = fv[~fv.isin(UN)]
    ev = et[rank].str.strip(); ev = ev[~ev.isin(UN)]
    fa_rows.append(dict(dataset="FunAqua (ITS)", rank=rank,
                        fungi=len(ft) / max(fv.nunique(), 1),
                        other=len(et) / max(ev.nunique(), 1)))

pr2_path = METAPR2_DIR / "asv.xlsx"
pr_rows = []
if pr2_path.exists():
    asv = pd.read_excel(pr2_path)
    mask = pd.Series(False, index=asv.index)
    for c in ("domain", "supergroup", "division", "subdivision", "class", "order",
              "family", "genus", "species"):
        mask |= asv[c].astype(str).str.contains("Fungi", na=False)
    nf, nff = (~mask).sum(), mask.sum()
    for rank, lab in (("division", "Phylum"), ("class", "Class"),
                      ("order", "Order"), ("family", "Family")):
        pr_rows.append(dict(dataset="metaPR2 (18S)", rank=lab,
                            fungi=nff / max(asv.loc[mask, rank].nunique(), 1),
                            other=nf / max(asv.loc[~mask, rank].nunique(), 1)))
else:
    print("  metaPR2 asv.xlsx not present; Figure S3 will show FunAqua only")

gr = pd.DataFrame(fa_rows + pr_rows)
gr["ratio"] = gr.fungi / gr.other
gr.to_csv(S / "Table_S14.tsv", sep="\t", index=False)
print(gr.to_string(index=False, float_format=lambda v: f"{v:.2f}"))

fig, ax = plt.subplots(figsize=(W1, 46 * MM))
order = ["Phylum", "Class", "Order", "Family"]
w = 0.36
for i, (ds, col) in enumerate((("FunAqua (ITS)", BLUE), ("metaPR2 (18S)", ORANGE))):
    sub = gr[gr.dataset == ds].set_index("rank").reindex(order)
    if sub.ratio.isna().all():
        continue
    ax.bar(np.arange(len(order)) + (i - 0.5) * w, sub.ratio, width=w,
           color=col, linewidth=0, label=ds, zorder=3)
ax.axhline(1.0, color=INK, linewidth=0.6, zorder=4)
ax.set_xticks(range(len(order))); ax.set_xticklabels(order)
ax.set_ylabel("sequences per label,\nfungi ÷ non-fungal")
ax.yaxis.grid(True, linewidth=0.45, color=GRID)
ax.legend(loc="upper right")           # upper left sat on top of the tallest bar
ax.set_ylim(0, max(2.6, gr.ratio.max() * 1.28))
save(fig, "Figure_S3_rank_comparability")

print(f"\nsupplementary written to {S}")
