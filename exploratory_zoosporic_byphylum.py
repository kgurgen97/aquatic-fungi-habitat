#!/usr/bin/env python3
"""POST HOC, added 2026-10-21 after the pre-registered flagellum test split:
NEON supported it (p 5.1e-08), FunAqua did not (p 0.179, direction reversed).

ZOOSPORIC_PLAN.md already warned that phylum is a coarse proxy for the trait. This
asks which phylum drives the split, and nothing here changes the pre-registered answer
(NOT REPLICATED). Any claim arising from this needs a new dataset, not these two.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
P = REPO_ROOT
d = pd.read_csv(P / "results/zoosporic.tsv", sep="\t")

print("=== median water_fraction by phylum, both datasets")
hdr = f"{'phylum':30s} {'trait':12s} " + " ".join(f"{k:>22s}" for k in ("funaqua", "neon"))
print(hdr)
rows = []
for ph in sorted(d.phylum.unique()):
    line = f"{ph:30s} {d[d.phylum==ph].trait.iat[0]:12s} "
    rec = {"phylum": ph}
    for ds in ("funaqua", "neon"):
        s = d[(d.phylum == ph) & (d.dataset == ds)].water_fraction
        line += f"{(f'{s.median():.3f} (n={len(s)})' if len(s) else '-'):>22s} "
        rec[ds] = float(s.median()) if len(s) else None
        rec[f"n_{ds}"] = int(len(s))
    print(line)
    rows.append(rec)

print("\n=== the two flagellated phyla with enough orders, compared directly")
for ds in ("funaqua", "neon"):
    sub = d[d.dataset == ds]
    ch = sub[sub.phylum == "Chytridiomycota"].water_fraction.to_numpy()
    ro = sub[sub.phylum == "Rozellomycota"].water_fraction.to_numpy()
    if len(ch) >= 3 and len(ro) >= 3:
        u = stats.mannwhitneyu(ch, ro, alternative="two-sided")
        print(f"  {ds:9s} Chytridiomycota n={len(ch):2d} median {np.median(ch):.3f}   "
              f"Rozellomycota n={len(ro):2d} median {np.median(ro):.3f}   "
              f"Mann-Whitney p={u.pvalue:.4g}")

print("\n=== Chytridiomycota alone vs the aflagellate group")
for ds in ("funaqua", "neon"):
    sub = d[d.dataset == ds]
    ch = sub[sub.phylum == "Chytridiomycota"].water_fraction.to_numpy()
    af = sub[sub.trait == "aflagellate"].water_fraction.to_numpy()
    u = stats.mannwhitneyu(ch, af, alternative="two-sided")
    rbc = 2 * u.statistic / (len(ch) * len(af)) - 1
    print(f"  {ds:9s} chytrids n={len(ch):2d} median {np.median(ch):.3f}   "
          f"aflagellate n={len(af):2d} median {np.median(af):.3f}   "
          f"p={u.pvalue:.4g}  rbc={rbc:+.3f}")

print("\n=== Rozellomycota alone vs the aflagellate group")
for ds in ("funaqua", "neon"):
    sub = d[d.dataset == ds]
    ro = sub[sub.phylum == "Rozellomycota"].water_fraction.to_numpy()
    af = sub[sub.trait == "aflagellate"].water_fraction.to_numpy()
    if len(ro) < 3:
        print(f"  {ds:9s} too few Rozellomycota ({len(ro)})"); continue
    u = stats.mannwhitneyu(ro, af, alternative="two-sided")
    print(f"  {ds:9s} Rozellomycota n={len(ro):2d} median {np.median(ro):.3f}   "
          f"aflagellate n={len(af):2d} median {np.median(af):.3f}   p={u.pvalue:.4g}")

print("\n=== how many of these orders are formally described vs environmental clades?")
env = d[d.order.str.contains(r"^(GS|Branch|Agar|Capnodiales$)|incertae|Incertae", regex=True)]
for ds in ("funaqua", "neon"):
    sub = d[d.dataset == ds]
    for ph in ("Chytridiomycota", "Rozellomycota"):
        s = sub[sub.phylum == ph]
        e = s[s.order.str.contains(r"^(GS|Branch)|incertae|Incertae", regex=True)]
        print(f"  {ds:9s} {ph:18s} {len(s):2d} orders, "
              f"{len(e):2d} are environmental clade labels ({len(e)/len(s):.0%})")

json.dump(rows, open(P / "results/zoosporic_byphylum.json", "w"), indent=1)
print("\nwrote results/zoosporic_byphylum.json")
