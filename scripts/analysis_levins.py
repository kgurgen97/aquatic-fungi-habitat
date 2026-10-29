#!/usr/bin/env python3
"""ADDED AFTER EUK_ROBUSTNESS_PLAN.md was frozen, and before its results were seen.

Why: the published niche-breadth literature says the opposite of our result. In
freshwater sediment, Levins' niche breadth runs bacteria < fungi < protists, i.e. fungi
NARROWER than protists; in alpine soils fungi are narrowest of all groups. Our paired
compartment-preference metric says fungi partition water vs sediment LESS than other
eukaryotes do.

Those may simply be different quantities: Levins' B measures how evenly a taxon is
spread across many samples, while our index measures whether a taxon favours one of two
compartments at the same site. This computes BOTH on the same data, the same matched
pairs and the same rarefaction depth, so the tension can be settled empirically instead
of argued.

Levins' B = 1 / sum(p_i^2), p_i = taxon's share of its own total in unit i.
  - across the 2 compartments: B in [1, 2]; 1 = found in one only, 2 = even split
  - across all 408 samples:    B in [1, 408]; the literature's formulation
Computed per Order (to match our analysis) and per OTU (to match how the literature
usually does it).
"""
import collections, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = FUNAQUA_DIR
P = REPO_ROOT
UNASSIGNED = {".", "", "unidentified", "unassigned", "NA"}
MIN_READS, SEED = 100, 20261021


def load(otu_csv, tax_csv):
    with open(D / otu_csv) as fh:
        cols = fh.readline().rstrip("\n").split(",")
    dt = {c: np.int32 for c in cols[1:]}; dt[cols[0]] = str
    o = pd.read_csv(D / otu_csv, index_col=0, dtype=dt)
    t = pd.read_csv(D / tax_csv, index_col=0, dtype=str,
                    keep_default_na=False).reindex(o.index)
    return o, t


md = pd.read_csv(P / "results/metadata_clean.tsv", sep="\t", dtype=str,
                 keep_default_na=False).set_index("sample_code")
md["water"] = md.sample_source.str.lower().str.contains("water")

fo, ft = load("FunAqua_fungi_otu.csv", "FunAqua_fungi_tax.csv")
eo, et = load("FunAqua_euk_otu.csv", "FunAqua_euk_tax.csv")
k = (et.Kingdom.str.strip() != "Fungi").values
eo, et = eo[k], et[k]

fk = ~ft.Order.str.strip().isin(UNASSIGNED)
ek = ~et.Order.str.strip().isin(UNASSIGNED)
fo, flab = fo[fk.values], ft.loc[fk, "Order"].str.strip().values
eo, elab = eo[ek.values], et.loc[ek, "Order"].str.strip().values
ekg = et.loc[ek, "Kingdom"].str.strip().values

fd, ed = fo.sum(axis=0), eo.sum(axis=0)
pairs = []
for site, g in md.loc[[c for c in fo.columns if c in md.index]].groupby("plot_name"):
    w, s = g[g.water], g[~g.water]
    if len(w) and len(s):
        a, b = w.index[0], s.index[0]
        if min(fd[a], fd[b], ed[a], ed[b]) >= MIN_READS:
            pairs.append((a, b, int(min(fd[a], fd[b], ed[a], ed[b]))))
print(f"pairs: {len(pairs)}   median common depth "
      f"{int(np.median([p[2] for p in pairs]))} reads")

rng = np.random.default_rng(SEED)


def rarefied_matrix(otu):
    """OTU x (2*npairs) matrix, every sample rarefied to its pair's common depth."""
    M = np.zeros((otu.shape[0], 2 * len(pairs)), dtype=np.int64)
    for i, (a, b, n) in enumerate(pairs):
        for j, col in ((0, a), (1, b)):
            c = otu[col].to_numpy().astype(np.int64)
            M[:, 2 * i + j] = c if c.sum() <= n else rng.multivariate_hypergeometric(c, n)
    return M


def levins(M):
    """B per row, across the columns given. Rows with no reads -> nan."""
    tot = M.sum(1, keepdims=True).astype(float)
    with np.errstate(invalid="ignore", divide="ignore"):
        p = np.where(tot > 0, M / tot, np.nan)
        return 1.0 / np.nansum(p ** 2, axis=1)


print("rarefying ...", flush=True)
Mf, Me = rarefied_matrix(fo), rarefied_matrix(eo)
water_cols = np.arange(0, 2 * len(pairs), 2)
sed_cols = np.arange(1, 2 * len(pairs), 2)


def two_compartment(M):
    """B across exactly 2 units: water total vs sediment total."""
    A = np.column_stack([M[:, water_cols].sum(1), M[:, sed_cols].sum(1)])
    return levins(A)


def agg_to_order(M, lab):
    orders = pd.unique(lab)
    oi = {o: i for i, o in enumerate(orders)}
    code = np.array([oi[o] for o in lab])
    A = np.zeros((len(orders), M.shape[1]), dtype=np.int64)
    for j in range(M.shape[1]):
        A[:, j] = np.bincount(code, weights=M[:, j], minlength=len(orders)).astype(np.int64)
    return A, orders


res = {}
print("\n" + "=" * 74)
print("LEVINS' B ACROSS THE 2 COMPARTMENTS   (1 = one compartment only, 2 = even)")
rows = []
for nm, M, lab in (("fungi", Mf, flab), ("other eukaryotes", Me, elab)):
    A, orders = agg_to_order(M, lab)
    keep = A.sum(1) >= 20                     # same spirit as the >=20 discordant rule
    b = two_compartment(A)[keep]
    b = b[np.isfinite(b)]
    rows.append((nm, b, orders[keep]))
    print(f"  {nm:18s} orders n={len(b):4d}   median B = {np.median(b):.4f}   "
          f"mean {b.mean():.4f}")
u = stats.mannwhitneyu(rows[0][1], rows[1][1], alternative="two-sided")
wider = "fungi" if np.median(rows[0][1]) > np.median(rows[1][1]) else "other eukaryotes"
print(f"  Mann-Whitney p = {u.pvalue:.4g}   BROADER (more generalist): {wider}")
res["two_compartment_order"] = {"fungi": float(np.median(rows[0][1])),
                                "euk": float(np.median(rows[1][1])),
                                "p": float(u.pvalue), "broader": wider,
                                "n_f": int(len(rows[0][1])), "n_e": int(len(rows[1][1]))}

print("\n" + "=" * 74)
print(f"LEVINS' B ACROSS ALL {2*len(pairs)} SAMPLES   (the literature's formulation)")
for unit, (Af, Ae) in (("per Order", (agg_to_order(Mf, flab)[0], agg_to_order(Me, elab)[0])),
                       ("per OTU", (Mf, Me))):
    out_ = []
    for nm, A in (("fungi", Af), ("other eukaryotes", Ae)):
        keep = A.sum(1) >= 20
        b = levins(A[keep]); b = b[np.isfinite(b)]
        out_.append((nm, b))
        print(f"  {unit:10s} {nm:18s} n={len(b):6d}   median B = {np.median(b):7.3f}   "
              f"mean {b.mean():7.3f}")
    uu = stats.mannwhitneyu(out_[0][1], out_[1][1], alternative="two-sided")
    wid = "fungi" if np.median(out_[0][1]) > np.median(out_[1][1]) else "other eukaryotes"
    print(f"  {unit:10s} Mann-Whitney p = {uu.pvalue:.4g}   BROADER: {wid}\n")
    res[f"all_samples_{unit.split()[-1].lower()}"] = {
        "fungi": float(np.median(out_[0][1])), "euk": float(np.median(out_[1][1])),
        "p": float(uu.pvalue), "broader": wid,
        "n_f": int(len(out_[0][1])), "n_e": int(len(out_[1][1]))}

print("=" * 74)
print("PER KINGDOM, Levins' B across all samples, per OTU")
ko = collections.defaultdict(list)
keep = Me.sum(1) >= 20
bo = levins(Me[keep]); kk = ekg[keep]
for b, kname in zip(bo, kk):
    if np.isfinite(b):
        ko[kname].append(b)
keepf = Mf.sum(1) >= 20
bf = levins(Mf[keepf]); bf = bf[np.isfinite(bf)]
ko["Fungi"] = list(bf)
print(f"  {'group':20s} {'n OTUs':>8s} {'median B':>10s}")
for kname, v in sorted(ko.items(), key=lambda x: -np.median(x[1])):
    if len(v) >= 20:
        print(f"  {kname:20s} {len(v):8d} {np.median(v):10.3f}")
res["per_kingdom"] = {kname: {"n": len(v), "median_B": float(np.median(v))}
                      for kname, v in ko.items() if len(v) >= 20}

print("\n" + "=" * 74)
print("DO THE TWO METRICS AGREE ON THE SAME DATA?")
pref = pd.read_csv(P / "results/matched_fungi.tsv", sep="\t")
pe = pd.read_csv(P / "results/matched_euk.tsv", sep="\t")
Af, of_ = agg_to_order(Mf, flab); Ae, oe_ = agg_to_order(Me, elab)
bf2 = dict(zip(of_, two_compartment(Af)))
be2 = dict(zip(oe_, two_compartment(Ae)))
per_order = []
for nm, d, bmap in (("fungi", pref, bf2), ("other eukaryotes", pe, be2)):
    d = d.copy()
    d["B"] = [bmap.get(o, np.nan) for o in d.order]
    d["effect"] = (d.frac - .5).abs()
    ok = d.dropna(subset=["B", "effect"])
    rho, pv = stats.spearmanr(ok.B, ok.effect)
    print(f"  {nm:18s} Spearman(Levins B, |preference-0.5|) = {rho:+.3f}  p = {pv:.3g}  n={len(ok)}")
    res[f"agreement_{nm.split()[0]}"] = {"rho": float(rho), "p": float(pv), "n": int(len(ok))}
    per_order.append(ok.assign(group=nm)[["order", "group", "frac", "effect", "B", "disc"]])
# Saved per order so the supplementary figure can show this agreement rather than
# only quote the two correlation coefficients.
pd.concat(per_order, ignore_index=True).to_csv(
    P / "results/levins_per_order.tsv", sep="\t", index=False)
print("  wrote results/levins_per_order.tsv")
print("  (a strong negative correlation would mean the two metrics measure the same thing)")

json.dump(res, open(P / "results/levins.json", "w"), indent=1)
print("\nwrote results/levins.json")
