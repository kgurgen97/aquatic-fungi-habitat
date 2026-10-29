#!/usr/bin/env python3
"""Two remaining FunAqua questions:
  B  does water/sediment partitioning depend on the type of water body?
  K  which eukaryote groups drive the stronger partitioning seen in the matched benchmark?
Same design throughout: paired sites, within-pair rarefaction to the pair minimum, order level."""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
pd.set_option("display.width", 215)
D = Path.home() / "research/datasets/funaqua"
P = Path.home() / "research/projects/funaqua-habitat"
SEED, UN = 20261021, {".", "", "unidentified", "unassigned", "NA"}
rng = np.random.default_rng(SEED)

md = pd.read_csv(P / "results/metadata_clean.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("sample_code")
md["water"] = md.sample_source.str.lower().str.contains("water")

def load(o_csv, t_csv):
    with open(D / o_csv) as fh:
        cols = fh.readline().rstrip("\n").split(",")
    dt = {c: np.int32 for c in cols[1:]}; dt[cols[0]] = str
    o = pd.read_csv(D / o_csv, index_col=0, dtype=dt)
    t = pd.read_csv(D / t_csv, index_col=0, dtype=str, keep_default_na=False).reindex(o.index)
    return o, t

def qualify(otu, mr=100):
    d = otu.sum(axis=0); out = []
    for site, g in md.loc[[c for c in otu.columns if c in md.index]].groupby("plot_name"):
        w, s = g[g.water], g[~g.water]
        if len(w) and len(s) and d[w.index[0]] >= mr and d[s.index[0]] >= mr:
            out.append((site, w.index[0], s.index[0], g.continent.iat[0], g.biome.iat[0]))
    return out

def presence(otu, lab, prs):
    orders = pd.unique(lab); oi = {o: i for i, o in enumerate(orders)}
    code = np.array([oi[o] for o in lab])
    W = np.zeros((len(prs), len(orders)), bool); S = np.zeros_like(W)
    for k, p in enumerate(prs):
        a, b = p[1], p[2]
        ca, cb = otu[a].to_numpy(), otu[b].to_numpy()
        n = int(min(ca.sum(), cb.sum()))
        ra = ca if ca.sum() <= n else rng.multivariate_hypergeometric(ca.astype(np.int64), n)
        rb = cb if cb.sum() <= n else rng.multivariate_hypergeometric(cb.astype(np.int64), n)
        W[k] = np.bincount(code, weights=ra, minlength=len(orders)) > 0
        S[k] = np.bincount(code, weights=rb, minlength=len(orders)) > 0
    return W, S, orders

def tested(W, S, orders, md_=20):
    rows = []
    for j, o in enumerate(orders):
        a = int((W[:, j] & ~S[:, j]).sum()); b = int((~W[:, j] & S[:, j]).sum())
        if a + b >= md_:
            rows.append(dict(order=o, frac=a / (a + b), disc=a + b,
                             p=stats.binomtest(a, a + b, .5).pvalue))
    r = pd.DataFrame(rows)
    if len(r): r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r

fo, ft = load("FunAqua_fungi_otu.csv", "FunAqua_fungi_tax.csv")
k = ~ft.Order.str.strip().isin(UN); fo, flab = fo[k.values], ft.loc[k, "Order"].str.strip().values
prs = qualify(fo)
W, S, orders = presence(fo, flab, prs)
overall = tested(W, S, orders).set_index("order")
sig = overall[overall.q < .05]
biome = np.array([p[4] for p in prs])

print("=== B  partitioning by type of water body")
short = lambda b: b.replace("aquatic: ", "").replace(" biome", "")
cols = {}
for b in pd.unique(biome):
    m = biome == b
    if m.sum() >= 25:
        rb = tested(W[m], S[m], orders, md_=10)
        if len(rb):
            cols[short(b)] = rb.set_index("order").frac
            nsig = int((rb.q < .05).sum()) if "q" in rb else 0
            print(f"  {short(b):22s} {int(m.sum()):4d} pairs, {len(rb):3d} orders testable, {nsig:3d} significant")
B = pd.DataFrame(cols)
common = B.dropna().index.intersection(sig.index)
if len(common) and B.shape[1] >= 2:
    sub = B.loc[common]
    agree = (sub > .5).all(axis=1) | (sub < .5).all(axis=1)
    print(f"  orders significant overall and testable in all {B.shape[1]} biomes: {len(common)}")
    print(f"  same direction in every biome: {int(agree.sum())} ({100*agree.mean():.0f}%)")
    print(f"  disagreeing: {list(sub.index[~agree])[:8]}")
    print(f"\n  correlation of preference between biomes:")
    print(sub.corr(method='spearman').round(2).to_string())
B.to_csv(P / "results/biome_fractions.tsv", sep="\t")

print("\n=== K  which eukaryote groups partition most (same pairs, same depth)")
eo, et = load("FunAqua_euk_otu.csv", "FunAqua_euk_tax.csv")
rows = []
for kingdom in ["Fungi"] + [k for k in et.Kingdom.str.strip().value_counts().index[:7] if k != "Fungi"]:
    sel = (et.Kingdom.str.strip() == kingdom) & ~et.Order.str.strip().isin(UN)
    if sel.sum() < 500:
        continue
    sub_o, sub_l = eo[sel.values], et.loc[sel, "Order"].str.strip().values
    pk = qualify(sub_o)
    if len(pk) < 50:
        rows.append(dict(kingdom=kingdom, otus=int(sel.sum()), pairs=len(pk), tested=0, sig=0,
                         pct=np.nan, median_effect=np.nan)); continue
    Wk, Sk, ok_ = presence(sub_o, sub_l, pk)
    rk = tested(Wk, Sk, ok_)
    s = rk[rk.q < .05] if len(rk) else rk
    rows.append(dict(kingdom=kingdom, otus=int(sel.sum()), pairs=len(pk), tested=len(rk), sig=len(s),
                     pct=100*len(s)/max(len(rk), 1),
                     median_effect=float((rk.frac - .5).abs().median()) if len(rk) else np.nan))
kd = pd.DataFrame(rows).sort_values("median_effect", ascending=False)
print(kd.round(3).to_string(index=False))
print("\nNOTE: pair counts differ between kingdoms (read depth differs), so the % column is not a fair")
print("cross-kingdom comparison on its own; median_effect is less sensitive to N.")
kd.to_csv(P / "results/kingdom_partitioning.tsv", sep="\t", index=False)
