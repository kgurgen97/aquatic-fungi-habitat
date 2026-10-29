#!/usr/bin/env python3
"""The five checks fixed in EUK_ROBUSTNESS_PLAN.md
(sha256 a12794f91b0df75b433df5fe37afb346ab3a603b5ed53152ca1a3dee85f879f8).

Run once. A check that fails is reported as failing.
"""
import collections, json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = FUNAQUA_DIR
P = REPO_ROOT
UNASSIGNED = {".", "", "unidentified", "unassigned", "NA"}
MIN_READS, MIN_DISC, N_DRAWS = 100, 20, 200
out = {}


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
keep = (et.Kingdom.str.strip() != "Fungi").values
eo, et = eo[keep], et[keep]
print(f"fungal OTUs {fo.shape[0]}, non-fungal eukaryote OTUs {eo.shape[0]}")


def prep(otu, tax, rank):
    k = ~tax[rank].str.strip().isin(UNASSIGNED)
    return otu[k.values], tax.loc[k, rank].str.strip().values, tax.loc[k]


def pairs_at(fo_, eo_):
    fd, ed = fo_.sum(axis=0), eo_.sum(axis=0)
    out_ = []
    for site, g in md.loc[[c for c in fo_.columns if c in md.index]].groupby("plot_name"):
        w, s = g[g.water], g[~g.water]
        if len(w) and len(s):
            a, b = w.index[0], s.index[0]
            if min(fd[a], fd[b], ed[a], ed[b]) >= MIN_READS:
                out_.append((site, a, b, int(min(fd[a], fd[b], ed[a], ed[b]))))
    return out_


def one_draw(otu, lab, prs, rng):
    orders = pd.unique(lab); oi = {o: i for i, o in enumerate(orders)}
    code = np.array([oi[o] for o in lab])
    cols = {c: otu[c].to_numpy().astype(np.int64) for _, a, b, _ in prs for c in (a, b)}
    W = np.zeros((len(prs), len(orders)), bool); S = np.zeros_like(W)
    for k, (_, a, b, n) in enumerate(prs):
        for col, M in ((a, W), (b, S)):
            c = cols[col]
            rr = c if c.sum() <= n else rng.multivariate_hypergeometric(c, n)
            M[k] = np.bincount(code, weights=rr, minlength=len(orders)) > 0
    wo = (W & ~S).sum(0); so = (~W & S).sum(0); disc = wo + so
    sel = np.where(disc >= MIN_DISC)[0]
    frac = wo[sel] / disc[sel]
    p = np.array([stats.binomtest(int(wo[j]), int(disc[j]), .5).pvalue for j in sel])
    q = multipletests(p, method="fdr_bh")[1] if len(p) else np.array([])
    return pd.DataFrame({"order": [orders[j] for j in sel], "disc": disc[sel],
                         "frac": frac, "q": q, "effect": np.abs(frac - .5)})


# ---------------------------------------------------------------- check 1
fo1, flab, ftx = prep(fo, ft, "Order")
eo1, elab, etx = prep(eo, et, "Order")
prs = pairs_at(fo1, eo1)
print(f"pairs qualifying for both groups: {len(prs)}  "
      f"median common depth {int(np.median([p[3] for p in prs]))} reads\n")

print("=" * 74)
print(f"CHECK 1  draw stability, {N_DRAWS} draws")
rev = 0; recs = []
for d in range(N_DRAWS):
    fr = one_draw(fo1, flab, prs, np.random.default_rng(20261021 + d))
    er = one_draw(eo1, elab, prs, np.random.default_rng(70261021 + d))
    ef, ee = fr.effect.median(), er.effect.median()
    sf = (fr.q < .05).mean(); se = (er.q < .05).mean()
    if ef >= ee:
        rev += 1
    recs.append((ef, ee, sf, se, len(fr), len(er)))
    if (d + 1) % 50 == 0:
        print(f"  {d+1}/{N_DRAWS}", flush=True)
r = pd.DataFrame(recs, columns=["eff_f", "eff_e", "sh_f", "sh_e", "n_f", "n_e"])
print(f"  median effect, fungi            : {r.eff_f.median():.3f}  "
      f"range {r.eff_f.min():.3f}-{r.eff_f.max():.3f}")
print(f"  median effect, other eukaryotes : {r.eff_e.median():.3f}  "
      f"range {r.eff_e.min():.3f}-{r.eff_e.max():.3f}")
print(f"  share significant, fungi        : {r.sh_f.median():.1%}  "
      f"range {r.sh_f.min():.1%}-{r.sh_f.max():.1%}")
print(f"  share significant, eukaryotes   : {r.sh_e.median():.1%}  "
      f"range {r.sh_e.min():.1%}-{r.sh_e.max():.1%}")
print(f"  orders tested, fungi {r.n_f.min()}-{r.n_f.max()}   eukaryotes {r.n_e.min()}-{r.n_e.max()}")
print(f"  ORDERING REVERSED in {rev}/{N_DRAWS} draws = {rev/N_DRAWS:.1%}")
c1 = rev / N_DRAWS <= 0.05
print(f"  -> CHECK 1 {'PASS' if c1 else 'FAIL'} (threshold: reversal in <= 5 % of draws)")
out["check1"] = {"reversals": rev, "n_draws": N_DRAWS, "pass": bool(c1),
                 "eff_f_median": float(r.eff_f.median()), "eff_e_median": float(r.eff_e.median()),
                 "eff_f_range": [float(r.eff_f.min()), float(r.eff_f.max())],
                 "eff_e_range": [float(r.eff_e.min()), float(r.eff_e.max())]}

# one reference draw for the remaining checks
fr = one_draw(fo1, flab, prs, np.random.default_rng(20261021))
er = one_draw(eo1, elab, prs, np.random.default_rng(70261021))
fr["grp"], er["grp"] = "fungi", "euk"
both = pd.concat([fr, er], ignore_index=True)

# ---------------------------------------------------------------- check 2
print("\n" + "=" * 74)
print("CHECK 2  matched on power (number of discordant pairs)")
qs = both.disc.quantile([0, .25, .5, .75, 1.0]).to_numpy()
both["bin"] = pd.cut(both.disc, bins=np.unique(qs), include_lowest=True)
print(f"  {'discordant-pair bin':24s} {'n fungi':>8s} {'n euk':>7s} "
      f"{'eff fungi':>10s} {'eff euk':>9s} {'p':>9s}")
ps, diffs = [], []
for b, g in both.groupby("bin", observed=True):
    f_, e_ = g[g.grp == "fungi"].effect, g[g.grp == "euk"].effect
    if len(f_) >= 4 and len(e_) >= 4:
        u = stats.mannwhitneyu(f_, e_, alternative="two-sided")
        ps.append(u.pvalue); diffs.append(f_.median() - e_.median())
        print(f"  {str(b):24s} {len(f_):8d} {len(e_):7d} {f_.median():10.3f} "
              f"{e_.median():9.3f} {u.pvalue:9.4g}")
    else:
        print(f"  {str(b):24s} {len(f_):8d} {len(e_):7d}   too few")
# stratified test: van Elteren via combining, plus a single rank test with bin as stratum
both["rank_in_bin"] = both.groupby("bin", observed=True).effect.rank(pct=True)
u_all = stats.mannwhitneyu(both[both.grp == "fungi"].rank_in_bin,
                           both[both.grp == "euk"].rank_in_bin, alternative="two-sided")
print(f"  stratified rank test across bins: p = {u_all.pvalue:.4g}")
print(f"  bins where fungi lower: {sum(d < 0 for d in diffs)}/{len(diffs)}")
c2 = u_all.pvalue < 0.05 and sum(d < 0 for d in diffs) >= max(1, len(diffs) - 1)
print(f"  -> CHECK 2 {'PASS' if c2 else 'FAIL'} (difference must survive holding power fixed)")
out["check2"] = {"stratified_p": float(u_all.pvalue), "bins": len(diffs),
                 "bins_fungi_lower": int(sum(d < 0 for d in diffs)), "pass": bool(c2)}

# ---------------------------------------------------------------- check 3
print("\n" + "=" * 74)
print("CHECK 3  fungi against each kingdom separately")
kof = etx.Kingdom.str.strip()
o2k = {}
for o, k in zip(elab, kof.values):
    o2k.setdefault(o, collections.Counter())[k] += 1
er["kingdom"] = [o2k[o].most_common(1)[0][0] if o in o2k else "?" for o in er.order]
print(f"  {'kingdom':20s} {'orders':>7s} {'median effect':>14s} {'vs fungi p':>11s} {'lower?':>8s}")
print(f"  {'FUNGI':20s} {len(fr):7d} {fr.effect.median():14.3f}")
ks = []
for k, g in er.groupby("kingdom"):
    if len(g) < 5:
        continue
    u = stats.mannwhitneyu(fr.effect, g.effect, alternative="two-sided")
    lower = fr.effect.median() < g.effect.median()
    ks.append((k, len(g), g.effect.median(), u.pvalue, lower))
    print(f"  {k:20s} {len(g):7d} {g.effect.median():14.3f} {u.pvalue:11.4g} "
          f"{'yes' if lower else 'NO':>8s}")
nlow = sum(1 for *_, l in ks if l)
nsig = sum(1 for *_, p_, l in ks if l and p_ < 0.05)
c3 = nlow == len(ks) and nsig >= 2
print(f"  fungi lower than {nlow}/{len(ks)} kingdoms, significantly so against {nsig}")
print(f"  -> CHECK 3 {'PASS' if c3 else 'FAIL'} (must be lower than every kingdom, "
      f"significantly against at least two)")
out["check3"] = {"kingdoms": [{"kingdom": k, "n": n, "effect": float(e), "p": float(p_),
                               "fungi_lower": bool(l)} for k, n, e, p_, l in ks],
                 "n_lower": nlow, "n_sig": nsig, "pass": bool(c3)}

# ---------------------------------------------------------------- check 4
print("\n" + "=" * 74)
print("CHECK 4  the same comparison at class and family level")
print(f"  {'rank':8s} {'n fungi':>8s} {'n euk':>7s} {'eff fungi':>10s} {'eff euk':>9s} {'p':>10s}")
rank_rows = []
for rank in ("Class", "Order", "Family"):
    fa, fl, _ = prep(fo, ft, rank)
    ea, el, _ = prep(eo, et, rank)
    pr = pairs_at(fa, ea)
    a = one_draw(fa, fl, pr, np.random.default_rng(20261021))
    b = one_draw(ea, el, pr, np.random.default_rng(70261021))
    u = stats.mannwhitneyu(a.effect, b.effect, alternative="two-sided")
    print(f"  {rank:8s} {len(a):8d} {len(b):7d} {a.effect.median():10.3f} "
          f"{b.effect.median():9.3f} {u.pvalue:10.4g}")
    rank_rows.append({"rank": rank, "n_f": len(a), "n_e": len(b),
                      "eff_f": float(a.effect.median()), "eff_e": float(b.effect.median()),
                      "p": float(u.pvalue), "fungi_lower": bool(a.effect.median() < b.effect.median())})
c4 = sum(r_["fungi_lower"] for r_ in rank_rows) >= 2
print(f"  -> CHECK 4 {'PASS' if c4 else 'FAIL'} (fungi lower at a majority of ranks)")
out["check4"] = {"ranks": rank_rows, "pass": bool(c4)}

# ---------------------------------------------------------------- check 5
print("\n" + "=" * 74)
print("CHECK 5  adjust for OTUs per order")
fper = pd.Series(collections.Counter(flab))
eper = pd.Series(collections.Counter(elab))
fr["n_otu"] = [fper.get(o, np.nan) for o in fr.order]
er["n_otu"] = [eper.get(o, np.nan) for o in er.order]
bb = pd.concat([fr.assign(grp="fungi"), er.assign(grp="euk")], ignore_index=True).dropna(subset=["n_otu"])
rk = lambda x: stats.rankdata(x)
Z = np.column_stack([np.ones(len(bb)), rk(bb.n_otu), rk(bb.disc)])
res_g = rk((bb.grp == "fungi").astype(float)) - Z @ np.linalg.lstsq(Z, rk((bb.grp == "fungi").astype(float)), rcond=None)[0]
res_e = rk(bb.effect) - Z @ np.linalg.lstsq(Z, rk(bb.effect), rcond=None)[0]
pr_, pp_ = stats.pearsonr(res_g, res_e)
print(f"  partial Spearman(being fungal, effect | OTUs per order, discordant pairs) = "
      f"{pr_:+.3f}  p = {pp_:.4g}")
lo = max(fr.n_otu.min(), er.n_otu.min()); hi = min(fr.n_otu.max(), er.n_otu.max())
f2 = fr[(fr.n_otu >= lo) & (fr.n_otu <= hi)]; e2 = er[(er.n_otu >= lo) & (er.n_otu <= hi)]
u2 = stats.mannwhitneyu(f2.effect, e2.effect, alternative="two-sided")
print(f"  restricted to the overlapping OTU-count range [{lo:.0f}, {hi:.0f}]: "
      f"fungi n={len(f2)} {f2.effect.median():.3f}  euk n={len(e2)} {e2.effect.median():.3f}  "
      f"p={u2.pvalue:.4g}")
c5 = pp_ < 0.05 and pr_ < 0
print(f"  -> CHECK 5 {'PASS' if c5 else 'FAIL'} (gap must survive the adjustment)")
out["check5"] = {"partial_r": float(pr_), "partial_p": float(pp_),
                 "overlap_p": float(u2.pvalue), "pass": bool(c5)}

# ---------------------------------------------------------------- verdict
print("\n" + "=" * 74)
for k in ("check1", "check2", "check3", "check4", "check5"):
    print(f"  {k}: {'PASS' if out[k]['pass'] else 'FAIL'}")
headline_ok = out["check1"]["pass"] and out["check2"]["pass"] and out["check3"]["pass"]
print(f"\n  CHECKS 1-3 (the gate for the title): {'PASS' if headline_ok else 'FAIL'}")
out["headline_allowed"] = bool(headline_ok)
both.to_csv(P / "results/euk_robustness.tsv", sep="\t", index=False)
json.dump(out, open(P / "results/euk_robustness.json", "w"), indent=1)
print("wrote results/euk_robustness.{tsv,json}")
