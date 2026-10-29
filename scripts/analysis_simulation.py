#!/usr/bin/env python3
"""When is a single-draw rarefaction unstable, and why does depth change the answer?

Added 2026-10-22. This is deductive, not an analysis of the data: it asks what the
estimator does under known truth, so there is no outcome to pre-register. Its job is to
turn two observations from the data into general statements a reader can apply to their
own dataset.

The model. In a paired design, both samples of a pair are rarefied to a common depth D.
A taxon present at relative abundance a in a sample is detected with probability

    P(detect) = 1 - (1 - a)^D  ~  1 - exp(-lambda),   lambda = a*D expected reads

Writing P_w and P_s for the two compartments, the discordant pairs - the only ones the
test uses - split as

    water-only     P_w (1 - P_s)
    sediment-only  P_s (1 - P_w)

so the expected preference index is

    E[index] = P_w(1 - P_s) / [P_w(1 - P_s) + P_s(1 - P_w)]

E[index] depends on D even at fixed true abundances, which suggested a prediction:
deeper data should show more extreme indices. **That prediction is wrong on the real
data** - NEON's indices are 1.02x FunAqua's (Mann-Whitney p = 0.69) despite NEON being
~150x deeper - and the reason it is wrong turned out to be the useful result.

A taxon contributes discordant pairs only when its detection probability is away from
0 and 1. The ">= 20 discordant pairs" rule is therefore a filter that selects taxa whose
lambda sits near 1, AT WHATEVER DEPTH. Raise the depth and the taxa meeting that
condition move down the abundance distribution: section 3 below shows the median lambda
of the testable set holding at 0.22-0.50 across a 1000-fold change in depth while the
median abundance of that set falls 447-fold.

So the marginal distribution of indices barely moves with depth - which is why aggregate
contrasts survive the comparison - while the identity of the taxa occupying that
position changes completely, which is why per-order values do not transfer. In the real
data the two testable sets overlap at Jaccard 0.51: 31 orders are testable only in
FunAqua and 34 only in NEON.

The failed prediction is kept in section 1 rather than deleted, because the contrast
between it and section 3 is the argument.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
P = REPO_ROOT
R = P / "results"
RNG = np.random.default_rng(20261022)
out = {}

# ---------------------------------------------------------------- 1. analytic
print("=" * 74)
print("1  The expected preference index moves with depth, at fixed true abundances")
print("=" * 74)


def expected_index(a_w, a_s, D):
    pw = 1 - np.exp(-a_w * D)
    ps = 1 - np.exp(-a_s * D)
    num = pw * (1 - ps)
    den = num + ps * (1 - pw)
    return np.where(den > 0, num / den, np.nan)


RATIO = 4.0           # the taxon is 4x more abundant in water than in sediment
A_S = 1e-4            # its relative abundance in sediment
depths = np.array([50, 131, 211, 500, 1_000, 5_000, 31_131, 100_000])
rows = []
for D in depths:
    e = float(expected_index(A_S * RATIO, A_S, D))
    rows.append(dict(depth=int(D), expected_index=e,
                     lambda_water=A_S * RATIO * D, lambda_sed=A_S * D))
    print(f"  D = {D:7,d}   lambda_w = {A_S*RATIO*D:8.3f}   "
          f"E[index] = {e:.3f}")
print(f"\n  A taxon whose true abundance ratio is fixed at {RATIO:.0f}x appears to have a")
print(f"  preference of {rows[1]['expected_index']:.2f} at FunAqua's depth and "
      f"{rows[6]['expected_index']:.2f} at NEON's.")
print("  Same organism, same truth, different number. No amount of extra sampling")
print("  removes this: it is what presence/absence means at that depth.")
out["depth_dependence"] = rows

# ---------------------------------------------------------------- 1b. the test
print("\n" + "=" * 74)
print("1b  That prediction, tested against the real data")
print("=" * 74)
fa_e = (pd.read_csv(R / "primary_order.tsv", sep="\t").water_fraction - .5).abs()
ne_e = (pd.read_csv(R / "neon_order.tsv", sep="\t").water_fraction - .5).abs()
u = stats.mannwhitneyu(ne_e, fa_e, alternative="two-sided")
print(f"  FunAqua (depth ~211)    n={len(fa_e):3d}  median |index-0.5| {fa_e.median():.3f}")
print(f"  NEON    (depth 31,131)  n={len(ne_e):3d}  median |index-0.5| {ne_e.median():.3f}")
print(f"  ratio {ne_e.median()/fa_e.median():.2f}x,  Mann-Whitney p = {u.pvalue:.3g}")
print("\n  PREDICTION NOT SUPPORTED. The estimator's algebra says indices should inflate")
print("  with depth; they do not. Section 3 shows why, and that is the real finding.")
out["depth_prediction"] = dict(funaqua_median=float(fa_e.median()),
                               neon_median=float(ne_e.median()),
                               ratio=float(ne_e.median() / fa_e.median()),
                               p=float(u.pvalue), supported=False)

# ---------------------------------------------------------------- 2. instability band
print("\n" + "=" * 74)
print("2  Where the significance call is unstable")
print("=" * 74)
N_PAIRS, N_DRAWS, MIN_DISC = 239, 200, 20
lams = np.geomspace(0.02, 50, 26)
ratios = [2.0, 4.0, 8.0]
band = []
for ratio in ratios:
    for lam_s in lams:
        lam_w = lam_s * ratio
        pw, ps = 1 - np.exp(-lam_w), 1 - np.exp(-lam_s)
        # simulate N_DRAWS independent realisations of the whole design
        w = RNG.random((N_DRAWS, N_PAIRS)) < pw
        s = RNG.random((N_DRAWS, N_PAIRS)) < ps
        wo = (w & ~s).sum(1)
        so = (~w & s).sum(1)
        disc = wo + so
        testable = disc >= MIN_DISC
        sig = np.zeros(N_DRAWS, bool)
        for i in np.where(testable)[0]:
            sig[i] = stats.binomtest(int(wo[i]), int(disc[i]), 0.5).pvalue < 0.05
        band.append(dict(ratio=ratio, lambda_sed=float(lam_s),
                         mean_discordant=float(disc.mean()),
                         frac_testable=float(testable.mean()),
                         frac_significant=float(sig.mean())))
bd = pd.DataFrame(band)
for ratio in ratios:
    g = bd[bd.ratio == ratio]
    unstable = g[(g.frac_significant > 0.05) & (g.frac_significant < 0.95)]
    if len(unstable):
        print(f"  ratio {ratio:.0f}x : unstable for lambda_sed in "
              f"{unstable.lambda_sed.min():.2f} – {unstable.lambda_sed.max():.2f} "
              f"({len(unstable)} of {len(g)} points)")
    else:
        print(f"  ratio {ratio:.0f}x : no unstable region at this resolution")
bd.to_csv(R / "sim_instability_band.tsv", sep="\t", index=False)

# ---------------------------------------------------------------- 3. self-selection
print("\n" + "=" * 74)
print("3  Why the prediction failed: the threshold pins the testable set")
print("=" * 74)
abund = np.geomspace(1e-6, 1e-1, 400)
RATIO3 = 3.0
print(f"  {'depth':>8s} {'tested':>7s} {'med lambda':>11s} {'med abundance':>14s} "
      f"{'med |index-0.5|':>16s}")
sel_rows = []
for D in (100, 211, 1_000, 5_000, 31_131, 100_000):
    lam_s = abund * D
    pw, ps = 1 - np.exp(-lam_s * RATIO3), 1 - np.exp(-lam_s)
    w = RNG.random((N_PAIRS, len(abund))) < pw
    s = RNG.random((N_PAIRS, len(abund))) < ps
    wo, so = (w & ~s).sum(0), (~w & s).sum(0)
    disc = wo + so
    sel = disc >= MIN_DISC
    if not sel.any():
        continue
    eff = np.abs(wo[sel] / disc[sel] - .5)
    sel_rows.append(dict(depth=int(D), n_tested=int(sel.sum()),
                         med_lambda=float(np.median(lam_s[sel])),
                         med_abundance=float(np.median(abund[sel])),
                         med_effect=float(np.median(eff))))
    r = sel_rows[-1]
    print(f"  {D:8,d} {r['n_tested']:7d} {r['med_lambda']:11.3f} "
          f"{r['med_abundance']:14.2e} {r['med_effect']:16.3f}")
sr = pd.DataFrame(sel_rows)
print(f"\n  median lambda of the testable set: {sr.med_lambda.min():.2f}–"
      f"{sr.med_lambda.max():.2f} across a 1000-fold depth range")
print(f"  median abundance of that set falls {sr.med_abundance.max()/sr.med_abundance.min():.0f}-fold")
print(f"  median |index-0.5| holds at {sr.med_effect.min():.3f}–{sr.med_effect.max():.3f}")
fa_set = set(pd.read_csv(R / "primary_order.tsv", sep="\t").order)
ne_set = set(pd.read_csv(R / "neon_order.tsv", sep="\t").order)
jac = len(fa_set & ne_set) / len(fa_set | ne_set)
print(f"\n  real data: FunAqua tests {len(fa_set)}, NEON {len(ne_set)}, "
      f"both {len(fa_set & ne_set)}, Jaccard {jac:.2f}")
print(f"  {len(fa_set - ne_set)} orders testable only in FunAqua, "
      f"{len(ne_set - fa_set)} only in NEON")
sr.to_csv(R / "sim_self_selection.tsv", sep="\t", index=False)
out["self_selection"] = dict(rows=sel_rows, jaccard=float(jac),
                             only_funaqua=len(fa_set - ne_set),
                             only_neon=len(ne_set - fa_set))

# ---------------------------------------------------------------- 4. community
print("\n" + "=" * 74)
print("4  A FunAqua-like community reproduces the observed spread")
print("=" * 74)
# Abundance distribution taken from the real data rather than assumed: the per-order
# share of reads in the fungal table.
agg = R.parent / "agg_order.csv"
if agg.exists():
    a = pd.read_csv(agg, index_col=0)
    tot = a.sum(axis=1)
    share = (tot / tot.sum()).sort_values(ascending=False)
    share = share[share > 0].to_numpy()[:120]
    print(f"  order abundance shares taken from agg_order.csv: {len(share)} orders, "
          f"median {np.median(share):.2e}")
else:
    share = np.sort(RNG.lognormal(-7, 2.0, 120))[::-1]
    share /= share.sum()
    print("  agg_order.csv absent; using a lognormal stand-in")

D_MED = 211
TRUE_RATIO = RNG.choice([1.0, 1.0, 1.0, 2.0, 3.0, 5.0], size=len(share))
sig_counts, tested_counts = [], []
for d in range(N_DRAWS):
    lam_s = share * D_MED
    lam_w = lam_s * TRUE_RATIO
    pw, ps = 1 - np.exp(-lam_w), 1 - np.exp(-lam_s)
    w = RNG.random((N_PAIRS, len(share))) < pw
    s = RNG.random((N_PAIRS, len(share))) < ps
    wo, so = (w & ~s).sum(0), (~w & s).sum(0)
    disc = wo + so
    sel = np.where(disc >= MIN_DISC)[0]
    if not len(sel):
        continue
    p = np.array([stats.binomtest(int(wo[j]), int(disc[j]), .5).pvalue for j in sel])
    q = multipletests(p, method="fdr_bh")[1]
    sig_counts.append(int((q < 0.05).sum()))
    tested_counts.append(len(sel))
sc, tc = np.array(sig_counts), np.array(tested_counts)
print(f"  simulated: {tc.min()}–{tc.max()} orders tested, "
      f"{sc.min()}–{sc.max()} significant (median {np.median(sc):.0f}, sd {sc.std():.1f})")
obs = json.load((R / "draw_stability.json").open())["significant"]
print(f"  observed : 96–103 tested, {obs['min']}–{obs['max']} significant "
      f"(median {obs['median']:.0f}, sd {obs['sd']:.1f})")
print("\n  The spread in the real data is not a quirk of FunAqua. A community with this")
print("  abundance distribution, at this depth and this many pairs, produces a spread of")
print("  the same order from the rarefaction draw alone.")
out["community"] = dict(sim_sig_min=int(sc.min()), sim_sig_max=int(sc.max()),
                        sim_sig_median=float(np.median(sc)), sim_sig_sd=float(sc.std()),
                        sim_tested_min=int(tc.min()), sim_tested_max=int(tc.max()),
                        obs_sig_min=obs["min"], obs_sig_max=obs["max"],
                        obs_sig_sd=obs["sd"], depth=D_MED, n_pairs=N_PAIRS)
pd.DataFrame({"draw": np.arange(len(sc)), "n_tested": tc,
              "n_significant": sc}).to_csv(R / "sim_community.tsv", sep="\t", index=False)

json.dump(out, open(R / "simulation.json", "w"), indent=1)
print("\nwrote results/sim_instability_band.tsv, results/sim_community.tsv, "
      "results/simulation.json")
