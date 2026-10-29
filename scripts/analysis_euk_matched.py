#!/usr/bin/env python3
"""S4 corrected: fungi vs other eukaryotes on the SAME pairs at the SAME rarefaction depth.
The first attempt compared 239 fungal pairs against 580 eukaryote pairs, so the difference it found
was power, not biology."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
pd.set_option("display.width", 200)
from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = FUNAQUA_DIR
P = REPO_ROOT
SEED, UNASSIGNED = 20261021, {".", "", "unidentified", "unassigned", "NA"}
rng = np.random.default_rng(SEED)

md = pd.read_csv(P / "metadata_clean.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("sample_code")
md["water"] = md.sample_source.str.lower().str.contains("water")

def load(otu_csv, tax_csv):
    with open(D / otu_csv) as fh:
        cols = fh.readline().rstrip("\n").split(",")
    dt = {c: np.int32 for c in cols[1:]}; dt[cols[0]] = str
    o = pd.read_csv(D / otu_csv, index_col=0, dtype=dt)
    t = pd.read_csv(D / tax_csv, index_col=0, dtype=str, keep_default_na=False).reindex(o.index)
    return o, t

fo, ft = load("FunAqua_fungi_otu.csv", "FunAqua_fungi_tax.csv")
eo, et = load("FunAqua_euk_otu.csv", "FunAqua_euk_tax.csv")
eo, et = eo[(et.Kingdom.str.strip() != "Fungi").values], et[(et.Kingdom.str.strip() != "Fungi").values]

def named(t):
    k = ~t.Order.str.strip().isin(UNASSIGNED)
    return k.values, t.loc[k, "Order"].str.strip().values

fk, flab = named(ft); fo = fo[fk]
ek, elab = named(et); eo = eo[ek]

# pairs must qualify for BOTH groups, so the two analyses see identical sites
fd, ed = fo.sum(axis=0), eo.sum(axis=0)
pairs = []
for site, g in md.loc[[c for c in fo.columns if c in md.index]].groupby("plot_name"):
    w, s = g[g.water], g[~g.water]
    if len(w) and len(s):
        a, b = w.index[0], s.index[0]
        if min(fd[a], fd[b], ed[a], ed[b]) >= 100:
            pairs.append((site, a, b))
print(f"pairs qualifying for BOTH fungi and other eukaryotes: {len(pairs)}")

# and the same rarefaction depth per pair, taken across both groups
def run(otu, lab, depths):
    orders = pd.unique(lab); oi = {o: i for i, o in enumerate(orders)}
    code = np.array([oi[o] for o in lab])
    W = np.zeros((len(pairs), len(orders)), bool); S = np.zeros_like(W)
    for k, (_, a, b) in enumerate(pairs):
        n = depths[k]
        for col, M in ((a, W), (b, S)):
            c = otu[col].to_numpy()
            rr = c if c.sum() <= n else rng.multivariate_hypergeometric(c.astype(np.int64), int(n))
            M[k] = np.bincount(code, weights=rr, minlength=len(orders)) > 0
    rows = []
    for j, o in enumerate(orders):
        x = int((W[:, j] & ~S[:, j]).sum()); y = int((~W[:, j] & S[:, j]).sum())
        if x + y >= 20:
            rows.append(dict(order=o, w=x, s=y, disc=x + y, frac=x / (x + y),
                             p=stats.binomtest(x, x + y, .5).pvalue))
    r = pd.DataFrame(rows)
    if len(r): r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r

depths = [int(min(fd[a], fd[b], ed[a], ed[b])) for _, a, b in pairs]
print(f"common rarefaction depth per pair: median {int(np.median(depths))} reads")

fr = run(fo, flab, depths)
er = run(eo, elab, depths)
fs, es = fr[fr.q < .05], er[er.q < .05]
print(f"\n{'group':22s} {'orders tested':>13s} {'significant':>12s} {'%':>6s}  {'water':>6s} {'sediment':>9s}")
for nm, r_, s_ in (("fungi", fr, fs), ("other eukaryotes", er, es)):
    print(f"{nm:22s} {len(r_):13d} {len(s_):12d} {100*len(s_)/max(len(r_),1):5.0f}% "
          f"{int((s_.frac>.5).sum()):6d} {int((s_.frac<.5).sum()):9d}")
tbl = [[len(fs), len(fr) - len(fs)], [len(es), len(er) - len(es)]]
odds, p = stats.fisher_exact(tbl)
print(f"\nfungi vs other eukaryotes (matched pairs and depth): OR = {odds:.2f}, Fisher p = {p:.4g}")

# effect magnitude, not just the count of significant orders
print(f"\nmedian |preference - 0.5| among tested orders")
print(f"  fungi            : {(fr.frac - .5).abs().median():.3f}")
print(f"  other eukaryotes : {(er.frac - .5).abs().median():.3f}")
print(f"  Mann-Whitney p   : {stats.mannwhitneyu((fr.frac-.5).abs(), (er.frac-.5).abs()).pvalue:.4g}")
fr.to_csv(P / "results/matched_fungi.tsv", sep="\t", index=False)
er.to_csv(P / "results/matched_euk.tsv", sep="\t", index=False)
json.dump({"pairs": len(pairs), "fungi_tested": int(len(fr)), "fungi_sig": int(len(fs)),
           "euk_tested": int(len(er)), "euk_sig": int(len(es)), "fisher_p": float(p),
           "median_effect_fungi": float((fr.frac-.5).abs().median()),
           "median_effect_euk": float((er.frac-.5).abs().median())},
          open(P / "results/matched_summary.json", "w"), indent=1)
