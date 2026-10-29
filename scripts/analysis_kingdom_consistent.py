#!/usr/bin/env python3
"""Per-kingdom comparison on ONE consistent footing.

Needed because the headline has to change. "Fungi partition less than other eukaryotes"
pools protists with plants and animals, and the Levins run shows that is misleading:
per OTU across all samples, Alveolata (2.057) is broader than fungi (1.888), Straminipila
ties, and Viridiplantae (1.264) and Metazoa (1.142) are far narrower. Against protists -
which is what the published niche-breadth work compares fungi to - fungi are equal or
narrower, so there is no real conflict with that literature.

But those numbers are per OTU across 408 samples, while our preference index is per Order
across 2 compartments. Mixing scopes to make a kingdom claim would be sloppy. This
computes, for every kingdom including Fungi, on the same matched pairs and depth:
  - Levins' B across the 2 compartments, per Order and per OTU
  - the preference index |water_fraction - 0.5|, per Order
so the ordering of kingdoms rests on one metric at a time.
"""
import collections, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = FUNAQUA_DIR
P = REPO_ROOT
UNASSIGNED = {".", "", "unidentified", "unassigned", "NA"}
MIN_READS, MIN_DISC, SEED = 100, 20, 20261021


def load(o_, t_):
    with open(D / o_) as fh:
        cols = fh.readline().rstrip("\n").split(",")
    dt = {c: np.int32 for c in cols[1:]}; dt[cols[0]] = str
    o = pd.read_csv(D / o_, index_col=0, dtype=dt)
    t = pd.read_csv(D / t_, index_col=0, dtype=str, keep_default_na=False).reindex(o.index)
    return o, t


md = pd.read_csv(P / "results/metadata_clean.tsv", sep="\t", dtype=str,
                 keep_default_na=False).set_index("sample_code")
md["water"] = md.sample_source.str.lower().str.contains("water")

fo, ft = load("FunAqua_fungi_otu.csv", "FunAqua_fungi_tax.csv")
eo, et = load("FunAqua_euk_otu.csv", "FunAqua_euk_tax.csv")
nf = (et.Kingdom.str.strip() != "Fungi").values
eo, et = eo[nf], et[nf]

fsel = ~ft.Order.str.strip().isin(UNASSIGNED)
esel = ~et.Order.str.strip().isin(UNASSIGNED)
fo = fo[fsel.values]; flab = ft.loc[fsel, "Order"].str.strip().values
eo = eo[esel.values]; elab = et.loc[esel, "Order"].str.strip().values
ekg = et.loc[esel, "Kingdom"].str.strip().values
fkg = np.array(["Fungi"] * len(flab))

# one OTU table with a kingdom label per row
otu = pd.concat([fo, eo])
lab = np.r_[flab, elab]
kng = np.r_[fkg, ekg]
print(f"combined: {otu.shape[0]} OTUs, {len(set(kng))} kingdoms")

fd, ed = fo.sum(axis=0), eo.sum(axis=0)
pairs = []
for site, g in md.loc[[c for c in fo.columns if c in md.index]].groupby("plot_name"):
    w, s = g[g.water], g[~g.water]
    if len(w) and len(s):
        a, b = w.index[0], s.index[0]
        if min(fd[a], fd[b], ed[a], ed[b]) >= MIN_READS:
            pairs.append((a, b, int(min(fd[a], fd[b], ed[a], ed[b]))))
print(f"pairs: {len(pairs)}  median depth {int(np.median([p[2] for p in pairs]))}\n")

rng = np.random.default_rng(SEED)
W = np.zeros((otu.shape[0], len(pairs)), np.int64)
S = np.zeros_like(W)
for i, (a, b, n) in enumerate(pairs):
    for col, M in ((a, W), (b, S)):
        c = otu[col].to_numpy().astype(np.int64)
        M[:, i] = c if c.sum() <= n else rng.multivariate_hypergeometric(c, n)

# NOTE: rarefaction here is per sample to the pair minimum computed over the FULL
# fungal and eukaryote tables, exactly as analysis_euk_matched.py does, so the depth
# every kingdom sees is identical.

def two_comp_B(w, s):
    tot = (w + s).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        pw = np.where(tot > 0, w / tot, np.nan)
        ps = np.where(tot > 0, s / tot, np.nan)
        return 1.0 / (pw ** 2 + ps ** 2)


rows = []
for king in sorted(set(kng)):
    m = kng == king
    o_ = lab[m]; Wk = W[m]; Sk = S[m]
    orders = pd.unique(o_)
    oi = {x: i for i, x in enumerate(orders)}
    code = np.array([oi[x] for x in o_])
    AW = np.zeros((len(orders), len(pairs)), np.int64)
    AS = np.zeros_like(AW)
    for j in range(len(pairs)):
        AW[:, j] = np.bincount(code, weights=Wk[:, j], minlength=len(orders))
        AS[:, j] = np.bincount(code, weights=Sk[:, j], minlength=len(orders))

    # per-OTU 2-compartment B
    tot_otu = Wk.sum(1) + Sk.sum(1)
    bo = two_comp_B(Wk.sum(1), Sk.sum(1))[tot_otu >= 20]
    bo = bo[np.isfinite(bo)]

    # per-Order 2-compartment B
    tot_ord = AW.sum(1) + AS.sum(1)
    bord = two_comp_B(AW.sum(1), AS.sum(1))[tot_ord >= 20]
    bord = bord[np.isfinite(bord)]

    # per-Order preference effect, same >=20 discordant rule as everywhere else
    pres_w = AW > 0; pres_s = AS > 0
    wo = (pres_w & ~pres_s).sum(1); so = (~pres_w & pres_s).sum(1)
    disc = wo + so
    sel = disc >= MIN_DISC
    eff = np.abs(wo[sel] / disc[sel] - .5)

    rows.append(dict(kingdom=king, n_otu=int((tot_otu >= 20).sum()),
                     n_order_B=len(bord), n_order_pref=int(sel.sum()),
                     B_otu=float(np.median(bo)) if len(bo) else np.nan,
                     B_order=float(np.median(bord)) if len(bord) else np.nan,
                     effect=float(np.median(eff)) if len(eff) else np.nan))

t = pd.DataFrame(rows)
t = t[(t.n_order_pref >= 5) | (t.kingdom == "Fungi")].copy()

print("=" * 86)
print("ALL ON THE SAME MATCHED PAIRS AND DEPTH, 2 COMPARTMENTS")
print("B: 1 = one compartment only (specialist), 2 = even split (generalist)")
print("effect: 0 = no preference, 0.5 = absolute preference\n")
print(f"{'kingdom':18s} {'n OTU':>7s} {'B per OTU':>10s} {'n ord':>6s} {'B per order':>12s} "
      f"{'n pref':>7s} {'median effect':>14s}")
for _, r in t.sort_values("B_otu", ascending=False).iterrows():
    print(f"{r.kingdom:18s} {r.n_otu:7d} {r.B_otu:10.3f} {r.n_order_B:6d} {r.B_order:12.3f} "
          f"{r.n_order_pref:7d} {r.effect:14.3f}")

print("\n=== do the two metrics rank the kingdoms the same way?")
tt = t.dropna(subset=["B_otu", "effect"])
rho, pv = stats.spearmanr(tt.B_otu, tt.effect)
print(f"  Spearman(B per OTU, median effect) across {len(tt)} kingdoms = {rho:+.3f}  p = {pv:.3g}")
rho2, pv2 = stats.spearmanr(tt.B_order, tt.effect)
print(f"  Spearman(B per order, median effect) = {rho2:+.3f}  p = {pv2:.3g}")

print("\n=== where do fungi sit?")
for col, nm in (("B_otu", "B per OTU"), ("B_order", "B per order"), ("effect", "median effect")):
    s = t.dropna(subset=[col]).sort_values(col, ascending=(col == "effect"))
    rank = list(s.kingdom).index("Fungi") + 1
    broader = [k for k in s.kingdom[:rank - 1]]
    print(f"  {nm:14s} fungi rank {rank} of {len(s)} "
          f"(more generalist than fungi: {broader if broader else 'none'})")

print("\n=== fungi vs protists only (Alveolata + Straminipila + Rhizaria), per OTU B")
prot = {"Alveolata", "Straminipila", "Rhizaria"}
mp = np.isin(kng, list(prot)); mf = kng == "Fungi"
bp = two_comp_B(W[mp].sum(1), S[mp].sum(1)); bp = bp[np.isfinite(bp) & ((W[mp].sum(1) + S[mp].sum(1)) >= 20)]
bf = two_comp_B(W[mf].sum(1), S[mf].sum(1)); bf = bf[np.isfinite(bf) & ((W[mf].sum(1) + S[mf].sum(1)) >= 20)]
u = stats.mannwhitneyu(bf, bp, alternative="two-sided")
print(f"  fungi    n={len(bf):5d}  median B {np.median(bf):.3f}")
print(f"  protists n={len(bp):5d}  median B {np.median(bp):.3f}")
print(f"  Mann-Whitney p = {u.pvalue:.4g}  ->  "
      f"{'fungi broader' if np.median(bf) > np.median(bp) else 'protists broader'}")
print("  This is the comparison the published niche-breadth literature actually makes.")

t.to_csv(P / "results/kingdom_consistent.tsv", sep="\t", index=False)
json.dump({"table": t.to_dict("records"),
           "spearman_B_otu_vs_effect": [float(rho), float(pv)],
           "fungi_vs_protists": {"fungi_B": float(np.median(bf)), "protist_B": float(np.median(bp)),
                                 "p": float(u.pvalue), "n_f": int(len(bf)), "n_p": int(len(bp))}},
          open(P / "results/kingdom_consistent.json", "w"), indent=1)
print("\nwrote results/kingdom_consistent.{tsv,json}")
