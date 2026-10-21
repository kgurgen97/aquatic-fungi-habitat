#!/usr/bin/env python3
"""Does PR2's `order` column hold the same taxonomic RANK for fungi and for protists?

If fungal 'orders' are mostly subphyla while protist 'orders' are real orders, then
aggregating at that column lumps far more diversity for fungi - which makes fungal
presence/absence saturate and pushes fungi toward looking generalist. That is exactly
the artefact the replication is meant to rule out, so the dataset would manufacture the
result instead of testing it. This audits the labels by their nomenclatural suffix.
"""
import collections, re
from pathlib import Path
import pandas as pd

D = Path.home() / "research/datasets/metapr2/extract"
asv = pd.read_excel(D / "asv.xlsx")

mask = pd.Series(False, index=asv.index)
for c in ("domain", "supergroup", "division", "subdivision", "class", "order",
          "family", "genus", "species"):
    mask |= asv[c].astype(str).str.contains("Fungi", na=False)

RANK = [
    ("order (-ales)",        re.compile(r"ales$")),
    ("class (-mycetes/-ea)", re.compile(r"mycetes$|phyceae$|ea$")),
    ("subphylum (-mycotina/-incertae)", re.compile(r"mycotina$")),
    ("phylum (-mycota/-phyta)", re.compile(r"mycota$|phyta$|zoa$")),
    ("family (-aceae/-idae)", re.compile(r"aceae$|idae$")),
    ("placeholder (_X / _XX)", re.compile(r"_X+$|_sp\.?$|Incertae|incertae", re.I)),
]


def audit(labels, title):
    labels = [l for l in labels if l and l not in ("nan", "")]
    uniq = sorted(set(labels))
    print(f"\n=== {title}: {len(uniq)} distinct labels")
    hit = collections.Counter()
    unmatched = []
    for l in uniq:
        for name, pat in RANK:
            if pat.search(l):
                hit[name] += 1
                break
        else:
            unmatched.append(l)
    tot = len(uniq)
    for name, _ in RANK:
        if hit[name]:
            print(f"  {name:34s} {hit[name]:4d}  {hit[name]/tot:6.1%}")
    print(f"  {'unclassifiable by suffix':34s} {len(unmatched):4d}  {len(unmatched)/tot:6.1%}")
    true_order = hit["order (-ales)"] / tot
    print(f"  -> share of labels that are genuinely ORDER rank: {true_order:.1%}")
    if unmatched:
        print(f"  examples not matched: {unmatched[:10]}")
    return true_order, uniq


f_share, f_lab = audit(asv.loc[mask, "order"].astype(str).tolist(), "FUNGAL 'orders'")
n_share, n_lab = audit(asv.loc[~mask, "order"].astype(str).tolist(), "NON-FUNGAL 'orders'")

print("\n" + "=" * 70)
print(f"genuinely order-rank labels: fungi {f_share:.1%}  vs  non-fungal {n_share:.1%}")
print(f"distinct labels:             fungi {len(f_lab)}      vs  non-fungal {len(n_lab)}")
print(f"ASVs per label:              fungi {int(mask.sum()/max(len(f_lab),1))}"
      f"        vs  non-fungal {int((~mask).sum()/max(len(n_lab),1))}")

print("\nall 43 fungal labels, for the record:")
for l in f_lab:
    print(f"  {l}")
