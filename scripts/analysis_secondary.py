#!/usr/bin/env python3
"""FunAqua secondary analyses, as declared in ANALYSIS_PLAN.md:
  S0 sensitivity to the read threshold
  S1 taxonomic conservatism (permutation over order -> class)
  S2 regional consistency
  S3 functional reading via FungalTraits
  S4 eukaryote benchmark (same pairs, non-fungal eukaryotes)
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
pd.set_option("display.width", 215)

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = FUNAQUA_DIR
T = TRAITS_DIR
P = REPO_ROOT
SEED, UNASSIGNED = 20261021, {".", "", "unidentified", "unassigned", "NA"}
rng = np.random.default_rng(SEED)

md = pd.read_csv(P / "results/metadata_clean.tsv", sep="\t", dtype=str, keep_default_na=False).set_index("sample_code")
md["water"] = md.sample_source.str.lower().str.contains("water")

def load(otu_csv, tax_csv, rank="Order"):
    with open(D / otu_csv) as fh:
        cols = fh.readline().rstrip("\n").split(",")
    dt = {c: np.int32 for c in cols[1:]}; dt[cols[0]] = str
    o = pd.read_csv(D / otu_csv, index_col=0, dtype=dt)
    t = pd.read_csv(D / tax_csv, index_col=0, dtype=str, keep_default_na=False).reindex(o.index)
    k = ~t[rank].str.strip().isin(UNASSIGNED)
    return o[k.values], t[k.values], t.loc[k, rank].str.strip().values

def pairs_at(otu, min_reads):
    depth = otu.sum(axis=0)
    out = []
    for site, g in md.loc[[c for c in otu.columns if c in md.index]].groupby("plot_name"):
        w, s = g[g.water], g[~g.water]
        if len(w) and len(s):
            a, b = w.index[0], s.index[0]
            if depth[a] >= min_reads and depth[b] >= min_reads:
                out.append((site, a, b, g.continent.iat[0]))
    return out

def presence(otu, lab, prs):
    orders = pd.unique(lab); oi = {o: i for i, o in enumerate(orders)}
    code = np.array([oi[o] for o in lab])
    W = np.zeros((len(prs), len(orders)), bool); S = np.zeros_like(W)
    for k, (_, a, b, _) in enumerate(prs):
        ca, cb = otu[a].to_numpy(), otu[b].to_numpy()
        n = int(min(ca.sum(), cb.sum()))
        ra = ca if ca.sum() <= n else rng.multivariate_hypergeometric(ca.astype(np.int64), n)
        rb = cb if cb.sum() <= n else rng.multivariate_hypergeometric(cb.astype(np.int64), n)
        W[k] = np.bincount(code, weights=ra, minlength=len(orders)) > 0
        S[k] = np.bincount(code, weights=rb, minlength=len(orders)) > 0
    return W, S, orders

def test(W, S, orders, min_disc=20):
    rows = []
    for j, o in enumerate(orders):
        a = int((W[:, j] & ~S[:, j]).sum()); b = int((~W[:, j] & S[:, j]).sum())
        if a + b >= min_disc:
            rows.append(dict(order=o, w=a, s=b, disc=a + b, frac=a / (a + b),
                             p=stats.binomtest(a, a + b, .5).pvalue))
    r = pd.DataFrame(rows)
    if len(r):
        r["q"] = multipletests(r.p, method="fdr_bh")[1]
    return r

otu, tax, lab = load("FunAqua_fungi_otu.csv", "FunAqua_fungi_tax.csv")

print("=== S0  sensitivity to the read threshold")
for mr in (50, 100, 200, 400):
    prs = pairs_at(otu, mr)
    W, S, orders = presence(otu, lab, prs)
    r = test(W, S, orders)
    sig = r[r.q < .05] if len(r) else r
    print(f"  >= {mr:4d} reads : {len(prs):4d} pairs, {len(r):3d} orders tested, "
          f"{len(sig):3d} significant ({int((sig.frac>.5).sum()) if len(sig) else 0} water / "
          f"{int((sig.frac<.5).sum()) if len(sig) else 0} sediment)")

prs = pairs_at(otu, 100)
W, S, orders = presence(otu, lab, prs)
r = test(W, S, orders).set_index("order")
sig = r[r.q < .05]

print("\n=== S1  taxonomic conservatism: do orders in the same class agree?")
o2c = tax.drop_duplicates("Order").set_index(tax.drop_duplicates("Order")["Order"].str.strip())["Class"].str.strip()
cls = o2c.reindex(r.index)
dat = pd.DataFrame({"frac": r.frac, "cls": cls}).dropna()
dat = dat[~dat.cls.isin(UNASSIGNED)]
grp = dat.groupby("cls").filter(lambda g: len(g) >= 3)
obs = grp.groupby("cls").frac.var(ddof=0).mean()
null = []
for _ in range(10000):
    sh = grp.assign(frac=rng.permutation(grp.frac.values))
    null.append(sh.groupby("cls").frac.var(ddof=0).mean())
null = np.array(null)
pval = (null <= obs).mean()
print(f"  orders with a class and a preference: {len(dat)}; classes with >=3 orders: {grp.cls.nunique()} "
      f"covering {len(grp)} orders")
print(f"  mean within-class variance of preference: observed {obs:.4f}, "
      f"null {null.mean():.4f} (sd {null.std():.4f})")
print(f"  permutation p = {pval:.4f}  -> {'conserved' if pval < .05 else 'no evidence of conservatism'}")

print("\n=== S2  regional consistency")
cont = np.array([p[3] for p in prs])
reg = {}
for c in pd.unique(cont):
    m = cont == c
    if m.sum() >= 30:
        rc = test(W[m], S[m], orders, min_disc=10)
        if len(rc): reg[c] = rc.set_index("order").frac
        print(f"  {c:15s} {int(m.sum()):4d} pairs, {len(rc):3d} orders testable")
# the >=100 filter leaves only two regions testable; repeat at >=50 to see the wider picture
prs50 = pairs_at(otu, 50)
W50, S50, _ = presence(otu, lab, prs50)
c50 = np.array([p[3] for p in prs50])
print("  at the >=50 read threshold:")
for c in pd.unique(c50):
    m = c50 == c
    if m.sum() >= 25:
        rc = test(W50[m], S50[m], orders, min_disc=10)
        print(f"    {c:15s} {int(m.sum()):4d} pairs, {len(rc):3d} orders testable")

R = pd.DataFrame(reg)
common = R.dropna().index.intersection(sig.index)
if len(common) and R.shape[1] >= 2:
    sub = R.loc[common]
    agree = ((sub > .5).all(axis=1) | (sub < .5).all(axis=1))
    print(f"  orders significant overall AND testable in all {R.shape[1]} regions: {len(common)}")
    print(f"  same direction in every region: {int(agree.sum())} ({100*agree.mean():.0f}%)")
    print("  disagreeing:", list(sub.index[~agree])[:8])

print("\n=== S3  functional reading (FungalTraits primary_lifestyle, via genus -> order)")
ft = pd.read_csv(T / "polme2020-s1-fungal-traits-genera.csv", encoding="utf-8-sig")
ft["g"] = ft.GENUS.astype(str).str.strip()
ft = ft.dropna(subset=["primary_lifestyle"]).drop_duplicates("g")   # the table has repeated genus rows
g2l = ft.set_index("g")["primary_lifestyle"]
gt = tax.assign(g=tax.Genus.str.strip(), o=tax.Order.str.strip())
gt["life"] = gt.g.map(g2l)
dom = (gt.dropna(subset=["life"]).groupby(["o", "life"]).size().reset_index(name="n")
         .sort_values("n", ascending=False).drop_duplicates("o").set_index("o")["life"])
dl = pd.DataFrame({"frac": r.frac, "life": dom.reindex(r.index)}).dropna()
summ = dl.groupby("life").agg(orders=("frac", "size"), median_water_fraction=("frac", "median"))
summ = summ[summ.orders >= 3].sort_values("median_water_fraction", ascending=False)
print(summ.round(3).to_string())
if len(summ) >= 2:
    groups = [dl[dl.life == l].frac.values for l in summ.index]
    print(f"  Kruskal-Wallis across lifestyles: p = {stats.kruskal(*groups).pvalue:.4g}")

print("\n=== S4  eukaryote benchmark (non-fungal eukaryotes, same pairs)")
eotu, etax, elab = load("FunAqua_euk_otu.csv", "FunAqua_euk_tax.csv")
nonfungal = etax.Kingdom.str.strip() != "Fungi"
eotu2, elab2 = eotu[nonfungal.values], elab[nonfungal.values]
eprs = [p for p in pairs_at(eotu2, 100)]
EW, ES, eorders = presence(eotu2, elab2, eprs)
er = test(EW, ES, eorders)
esig = er[er.q < .05] if len(er) else er
print(f"  non-fungal eukaryote OTUs: {int(nonfungal.sum())}, orders {len(eorders)}")
print(f"  pairs {len(eprs)}, orders tested {len(er)}, significant {len(esig)} "
      f"({100*len(esig)/max(len(er),1):.0f}% of tested)")
print(f"  fungi for comparison      : orders tested {len(r)}, significant {len(sig)} "
      f"({100*len(sig)/max(len(r),1):.0f}% of tested)")
if len(er) and len(r):
    tbl = [[len(sig), len(r) - len(sig)], [len(esig), len(er) - len(esig)]]
    print(f"  fungi vs other eukaryotes, Fisher p = {stats.fisher_exact(tbl)[1]:.4g}")

r.to_csv(P / "results/primary_order_rerun.tsv", sep="\t")
er.to_csv(P / "results/euk_order.tsv", sep="\t", index=False)
R.to_csv(P / "results/regional_fractions.tsv", sep="\t")
json.dump({"conservatism_p": float(pval), "obs_within_class_var": float(obs),
           "fungi_sig_frac": float(len(sig)/max(len(r),1)),
           "euk_sig_frac": float(len(esig)/max(len(er),1))},
          open(P / "results/secondary_summary.json", "w"), indent=1)
