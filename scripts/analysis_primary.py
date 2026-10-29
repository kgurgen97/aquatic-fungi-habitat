#!/usr/bin/env python3
"""FunAqua primary endpoint, exactly as ANALYSIS_PLAN.md specifies.
Within-pair rarefaction to the pair minimum, then an exact paired test per order."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
pd.set_option("display.width", 210)

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = FUNAQUA_DIR
P = REPO_ROOT
(P / "results").mkdir(exist_ok=True)
RANK = "Order"
MIN_READS = 100
MIN_DISCORDANT = 20
SEED = 20261021
UNASSIGNED = {".", "", "unidentified", "unassigned", "NA"}

with open(D / "FunAqua_fungi_otu.csv") as fh:
    cols = fh.readline().rstrip("\n").split(",")
dt = {c: np.int32 for c in cols[1:]}; dt[cols[0]] = str
otu = pd.read_csv(D / "FunAqua_fungi_otu.csv", index_col=0, dtype=dt)
tax = pd.read_csv(D / "FunAqua_fungi_tax.csv", index_col=0, dtype=str,
                  keep_default_na=False).reindex(otu.index)
md = pd.read_csv(P / "metadata_clean.tsv", sep="\t", dtype=str,
                 keep_default_na=False).set_index("sample_code").loc[otu.columns]
md["water"] = md.sample_source.str.lower().str.contains("water")

keep = ~tax[RANK].str.strip().isin(UNASSIGNED)
otu, lab = otu[keep.values], tax.loc[keep, RANK].str.strip().values
print(f"{keep.sum()} of {len(keep)} OTUs have a named {RANK}; {pd.unique(lab).size} orders")

# ---- qualifying pairs
depth = otu.sum(axis=0)
pairs = []
for site, g in md.groupby("plot_name"):
    w, s = g[g.water], g[~g.water]
    if len(w) and len(s):
        a, b = w.index[0], s.index[0]
        if depth[a] >= MIN_READS and depth[b] >= MIN_READS:
            pairs.append((site, a, b, g.continent.iat[0], g.biome.iat[0]))
print(f"qualifying pairs (both samples >= {MIN_READS} reads): {len(pairs)}")

# ---- within-pair rarefaction to the pair minimum, then presence per order
rng = np.random.default_rng(SEED)
orders = pd.unique(lab)
oidx = {o: i for i, o in enumerate(orders)}
code = np.array([oidx[o] for o in lab])
W = np.zeros((len(pairs), len(orders)), bool)
S = np.zeros((len(pairs), len(orders)), bool)

def rarefy(counts, n):
    tot = counts.sum()
    if tot <= n:
        return counts
    return rng.multivariate_hypergeometric(counts.astype(np.int64), int(n))

for k, (site, a, b, cont, biome) in enumerate(pairs):
    ca, cb = otu[a].to_numpy(), otu[b].to_numpy()
    n = int(min(ca.sum(), cb.sum()))
    ra, rb = rarefy(ca, n), rarefy(cb, n)
    W[k] = np.bincount(code, weights=ra, minlength=len(orders)) > 0
    S[k] = np.bincount(code, weights=rb, minlength=len(orders)) > 0
print(f"rarefied each pair to its own minimum (median {int(np.median([min(otu[a].sum(), otu[b].sum()) for _,a,b,_,_ in pairs]))} reads)")

# ---- primary test
rows = []
for j, o in enumerate(orders):
    w_only = int((W[:, j] & ~S[:, j]).sum())
    s_only = int((~W[:, j] & S[:, j]).sum())
    both = int((W[:, j] & S[:, j]).sum())
    disc = w_only + s_only
    if disc < MIN_DISCORDANT:
        continue
    p = stats.binomtest(w_only, disc, 0.5).pvalue
    frac = w_only / disc
    lo, hi = stats.binomtest(w_only, disc).proportion_ci(0.95, method="wilson")
    rows.append(dict(order=o, water_only=w_only, sediment_only=s_only, both=both,
                     discordant=disc, water_fraction=frac, ci_lo=lo, ci_hi=hi, p=p))
r = pd.DataFrame(rows)
from statsmodels.stats.multitest import multipletests
r["q"] = multipletests(r.p, method="fdr_bh")[1]
r["prefers"] = np.where(r.q >= 0.05, "—", np.where(r.water_fraction > .5, "water", "sediment"))
r = r.sort_values("q")
r.to_csv(P / "results/primary_order.tsv", sep="\t", index=False)

sig = r[r.q < 0.05]
print(f"\norders tested (>= {MIN_DISCORDANT} discordant pairs): {len(r)}")
print(f"orders with a significant preference (BH q<0.05): {len(sig)}")
print(f"  prefer water   : {int((sig.prefers=='water').sum())}")
print(f"  prefer sediment: {int((sig.prefers=='sediment').sum())}")
print("\ntop 18 by q:")
print(sig.head(18)[["order","water_only","sediment_only","both","water_fraction","q","prefers"]]
      .to_string(index=False, float_format=lambda v: f"{v:.3g}"))
json.dump({"pairs": len(pairs), "orders_tested": int(len(r)), "significant": int(len(sig)),
           "water": int((sig.prefers=='water').sum()), "sediment": int((sig.prefers=='sediment').sum())},
          open(P / "results/primary_summary.json", "w"), indent=1)
np.savez_compressed(P / "results/presence_matrices.npz", W=W, S=S, orders=orders,
                    sites=np.array([p[0] for p in pairs]),
                    continents=np.array([p[3] for p in pairs]))
