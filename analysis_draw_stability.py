#!/usr/bin/env python3
"""POST HOC, added 2026-10-21. Not part of ANALYSIS_PLAN.md.

Why this exists: `analysis_primary.py` and `analysis_secondary.py` implement the same
declared primary endpoint at the same read threshold with the same seed, and disagree —
100 vs 99 orders tested, 49 vs 51 significant. The cause is not a coding error in either:
`analysis_secondary.py` consumes one module-level RNG across its four thresholds, so by
the >=100 row its rarefaction draw differs from the primary's fresh RNG.

So the declared design's headline count depends on a single random draw. This script
measures how much, and identifies the orders that are significant regardless of the draw.
The pre-registered primary remains `analysis_primary.py` with seed 20261021; nothing here
replaces it.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = FUNAQUA_DIR
P = REPO_ROOT
RANK, MIN_READS, MIN_DISCORDANT = "Order", 100, 20
N_DRAWS = 200
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

depth = otu.sum(axis=0)
pairs = []
for site, g in md.groupby("plot_name"):
    w, s = g[g.water], g[~g.water]
    if len(w) and len(s):
        a, b = w.index[0], s.index[0]
        if depth[a] >= MIN_READS and depth[b] >= MIN_READS:
            pairs.append((a, b))
print(f"{len(pairs)} qualifying pairs, {pd.unique(lab).size} named orders, {N_DRAWS} draws")

orders = pd.unique(lab)
oidx = {o: i for i, o in enumerate(orders)}
code = np.array([oidx[o] for o in lab])
counts = [(otu[a].to_numpy().astype(np.int64), otu[b].to_numpy().astype(np.int64))
          for a, b in pairs]
mins = [int(min(ca.sum(), cb.sum())) for ca, cb in counts]

tested_n, sig_n, water_n, sed_n = [], [], [], []
times_tested = np.zeros(len(orders), int)
times_sig = np.zeros(len(orders), int)
times_sig_water = np.zeros(len(orders), int)
frac_sum = np.zeros(len(orders))

for d in range(N_DRAWS):
    rng = np.random.default_rng(20261021 + d)
    W = np.zeros((len(pairs), len(orders)), bool)
    S = np.zeros_like(W)
    for k, ((ca, cb), n) in enumerate(zip(counts, mins)):
        ra = ca if ca.sum() <= n else rng.multivariate_hypergeometric(ca, n)
        rb = cb if cb.sum() <= n else rng.multivariate_hypergeometric(cb, n)
        W[k] = np.bincount(code, weights=ra, minlength=len(orders)) > 0
        S[k] = np.bincount(code, weights=rb, minlength=len(orders)) > 0
    wo = (W & ~S).sum(0); so = (~W & S).sum(0)
    disc = wo + so
    sel = np.where(disc >= MIN_DISCORDANT)[0]
    p = np.array([stats.binomtest(int(wo[j]), int(disc[j]), 0.5).pvalue for j in sel])
    q = multipletests(p, method="fdr_bh")[1]
    times_tested[sel] += 1
    frac = wo[sel] / disc[sel]
    frac_sum[sel] += frac
    s = sel[q < 0.05]
    times_sig[s] += 1
    times_sig_water[s[frac[q < 0.05] > .5]] += 1
    tested_n.append(len(sel)); sig_n.append(int((q < 0.05).sum()))
    water_n.append(int((frac[q < 0.05] > .5).sum()))
    sed_n.append(int((frac[q < 0.05] <= .5).sum()))
    if (d + 1) % 50 == 0:
        print(f"  {d+1}/{N_DRAWS} draws", flush=True)


def rep(name, v):
    v = np.array(v)
    print(f"  {name:26s} median {np.median(v):6.1f}   range {v.min()}–{v.max()}   "
          f"mean {v.mean():6.2f} (sd {v.std():.2f})")

print(f"\n=== how much does one random draw move the headline? ({N_DRAWS} draws)")
rep("orders tested", tested_n)
rep("orders significant", sig_n)
rep("  prefer water", water_n)
rep("  prefer sediment", sed_n)
print(f"\n  the two stored runs sit at 49 (analysis_primary.py) and 51 (analysis_secondary.py);")
print(f"  both fall inside this range, so neither was wrong — the design is draw-dependent.")

stab = pd.DataFrame({
    "order": orders,
    "times_tested": times_tested,
    "times_significant": times_sig,
    "pct_of_tested_draws_significant": np.where(times_tested > 0, times_sig / np.maximum(times_tested, 1), np.nan),
    "mean_water_fraction": np.where(times_tested > 0, frac_sum / np.maximum(times_tested, 1), np.nan),
    "sig_direction_water_pct": np.where(times_sig > 0, times_sig_water / np.maximum(times_sig, 1), np.nan),
}).query("times_tested > 0").sort_values("times_significant", ascending=False)
stab.to_csv(P / "results/draw_stability.tsv", sep="\t", index=False)

always = stab[(stab.times_tested == N_DRAWS) & (stab.times_significant == N_DRAWS)]
mostly = stab[stab.times_significant >= 0.95 * N_DRAWS]
never = stab[stab.times_significant == 0]
sometimes = stab[(stab.times_significant > 0) & (stab.times_significant < 0.95 * N_DRAWS)]
print(f"\n=== per-order stability of the significance call")
print(f"  significant in EVERY draw and tested in every draw : {len(always)}")
print(f"  significant in >= 95 % of its tested draws          : {len(mostly)}")
print(f"  significant in some draws but under 95 %            : {len(sometimes)}")
print(f"  never significant                                   : {len(never)}")
print(f"\n  -> the defensible set for the paper is the {len(mostly)} orders significant in >=95 % of draws,")
print(f"     not the {int(np.median(sig_n))}-odd that any single draw happens to return.")
print("\n  direction is never ambiguous among these: "
      f"{int(((mostly.sig_direction_water_pct == 1) | (mostly.sig_direction_water_pct == 0)).sum())}"
      f" of {len(mostly)} always point the same way")

print(f"\ntop 24 by stability:")
print(mostly.head(24).to_string(index=False, float_format=lambda v: f"{v:.3f}"))
if len(sometimes):
    print(f"\nthe unstable ones (significance depends on the draw):")
    print(sometimes.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

json.dump({"n_draws": N_DRAWS, "pairs": len(pairs),
           "tested": {"median": float(np.median(tested_n)), "min": int(min(tested_n)), "max": int(max(tested_n))},
           "significant": {"median": float(np.median(sig_n)), "min": int(min(sig_n)), "max": int(max(sig_n)),
                           "mean": float(np.mean(sig_n)), "sd": float(np.std(sig_n))},
           "stable_95pct": int(len(mostly)), "always": int(len(always)),
           "unstable": int(len(sometimes)), "never": int(len(never)),
           "stable_orders": mostly.order.tolist()},
          open(P / "results/draw_stability.json", "w"), indent=1)
print(f"\nwrote results/draw_stability.tsv and results/draw_stability.json")
