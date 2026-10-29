#!/usr/bin/env python3
"""Is there ANY rank at which PR2 treats fungi and protists with comparable granularity?

If yes, a rank-matched test is defensible even though the `order` column is not.
Granularity is measured two ways: ASVs per label, and the share of labels that are
genuine order/class-rank names rather than higher ranks or placeholders.
"""
import collections, re
from pathlib import Path
import pandas as pd

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = METAPR2_DIR
asv = pd.read_excel(D / "asv.xlsx")
RANKS = ["supergroup", "division", "subdivision", "class", "order", "family", "genus"]

mask = pd.Series(False, index=asv.index)
for c in RANKS + ["domain", "species"]:
    mask |= asv[c].astype(str).str.contains("Fungi", na=False)
nf, f = (~mask).sum(), mask.sum()
print(f"fungal ASVs {f}, non-fungal {nf}\n")

PLACE = re.compile(r"_X+$|Incertae|incertae|_sp\.?$|_[0-9A-Z]$", re.I)

print(f"{'rank':12s} {'fungal labels':>14s} {'ASVs/label':>11s} | "
      f"{'other labels':>13s} {'ASVs/label':>11s} | {'granularity ratio':>18s} {'placeholder f/o':>16s}")
rows = []
for r in RANKS:
    fl = asv.loc[mask, r].astype(str); nl = asv.loc[~mask, r].astype(str)
    fu, nu = fl.nunique(), nl.nunique()
    fa, na = f / max(fu, 1), nf / max(nu, 1)
    fp = sum(bool(PLACE.search(x)) for x in set(fl)) / max(fu, 1)
    np_ = sum(bool(PLACE.search(x)) for x in set(nl)) / max(nu, 1)
    ratio = fa / na if na else float("nan")
    print(f"{r:12s} {fu:14d} {fa:11.0f} | {nu:13d} {na:11.0f} | "
          f"{ratio:18.2f} {fp:7.0%}/{np_:<8.0%}")
    rows.append((r, fu, fa, nu, na, ratio, fp, np_))

print("\nA ratio near 1 means a fungal label aggregates about as much diversity as a")
print("non-fungal one. Far above 1 means fungal labels are coarser, which inflates")
print("fungal saturation and would push fungi toward looking generalist by construction.")

best = min(rows, key=lambda x: abs(x[5] - 1) if x[5] == x[5] else 9e9)
print(f"\nclosest to parity: {best[0]}  (ratio {best[5]:.2f}, "
       f"{best[1]} fungal vs {best[3]} non-fungal labels)")

print("\n=== what FunAqua offers at the comparable ranks, for reference")
FD = FUNAQUA_DIR
ft = pd.read_csv(FD / "FunAqua_fungi_tax.csv", index_col=0, dtype=str, keep_default_na=False)
et = pd.read_csv(FD / "FunAqua_euk_tax.csv", index_col=0, dtype=str, keep_default_na=False)
et = et[et.Kingdom.str.strip() != "Fungi"]
UN = {".", "", "unidentified", "unassigned", "NA"}
for r in ("Phylum", "Class", "Order", "Family"):
    fv = ft[r].str.strip(); fv = fv[~fv.isin(UN)]
    ev = et[r].str.strip(); ev = ev[~ev.isin(UN)]
    print(f"  {r:8s} fungi {fv.nunique():4d} labels ({len(ft)/max(fv.nunique(),1):6.0f} OTUs each)"
          f"   other euk {ev.nunique():4d} ({len(et)/max(ev.nunique(),1):6.0f} each)"
          f"   ratio {(len(ft)/max(fv.nunique(),1))/((len(et)/max(ev.nunique(),1))):.2f}")
