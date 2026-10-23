#!/usr/bin/env python3
"""Which estimator should be used? As BENCHMARK_PLAN.md specifies.

The plan's sha256 is in PLAN_HASHES.md and is asserted here before anything runs.

The generator is the occupancy-plus-mean-abundance model validated in
PREVALENCE_RESULTS.md - the only one that reproduced the measured re-selection
exponent. Five estimators see identical simulated counts and must decide, per
order, whether a water/sediment preference exists and in which direction.
"""
import hashlib
import itertools
import json
import pathlib
import re
import warnings
import zlib

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore", category=RuntimeWarning)

P = pathlib.Path.home() / "research/projects/funaqua-habitat"
R = P / "results"

PLAN = P / "BENCHMARK_PLAN.md"
plan_sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
rec = re.search(r"`BENCHMARK_PLAN\.md` \| `([0-9a-f]{64})`", (P / "PLAN_HASHES.md").read_text())
assert rec and plan_sha == rec.group(1), "BENCHMARK_PLAN.md changed since freezing"
print(f"plan verified: {plan_sha[:16]}…")

N_PAIRS, N_REP = 102, 40
DEPTHS = [100, 1_000, 10_000]
RATIOS = [2, 4, 8]
FRAC_PREFER = 0.30
MIN_DISCORDANT = 20
N_DRAWS_REPEATED = 25
STABLE_FRAC = 0.95
ALPHA = 0.05

emp = pd.read_csv(R / "sweep_empirical_abundance.tsv", sep="\t")
A_POOLED = emp.a_hat.to_numpy()
T = len(A_POOLED)

mat = pd.read_csv(R / "order_abundance_by_pair.tsv.gz", sep="\t", index_col=0)
assert list(mat.columns) == emp.order.tolist(), "order lists disagree"
OCC = np.clip((mat.to_numpy() > 0).mean(0), 1e-6, 1.0)
A_IF_PRESENT = A_POOLED / OCC
print(f"orders {T}, occupancy median {np.median(OCC):.3f}, pairs {N_PAIRS}")


def simulate(ratio, depth, rng):
    """-> counts (water, sediment), each (N_PAIRS, T), and the true index per order."""
    prefers = rng.random(T) < FRAC_PREFER
    r = np.where(prefers, ratio, 1.0)
    true_index = r / (1.0 + r)          # the estimand: the abundance-ratio preference
    pw = rng.random((N_PAIRS, T)) < OCC
    ps = rng.random((N_PAIRS, T)) < OCC
    lam_w = np.where(pw, A_IF_PRESENT * r, 0.0) * depth
    lam_s = np.where(ps, A_IF_PRESENT, 0.0) * depth
    return rng.poisson(lam_w), rng.poisson(lam_s), true_index, prefers


def _bh(p):
    out = np.full(len(p), np.nan)
    ok = np.isfinite(p)
    if ok.any():
        out[ok] = multipletests(p[ok], alpha=ALPHA, method="fdr_bh")[1]
    return out


def _pa_once(cw, cs, rng):
    """One rarefaction draw -> (index, discordant count) per order."""
    nw, ns = cw.sum(1), cs.sum(1)
    n = np.minimum(nw, ns)
    W = np.zeros_like(cw, dtype=bool)
    S = np.zeros_like(cs, dtype=bool)
    for j in range(cw.shape[0]):
        if n[j] <= 0:
            continue
        W[j] = rng.multivariate_hypergeometric(cw[j], int(n[j])) > 0
        S[j] = rng.multivariate_hypergeometric(cs[j], int(n[j])) > 0
    wo = (W & ~S).sum(0)
    so = (~W & S).sum(0)
    disc = wo + so
    idx = np.full(T, np.nan)
    sel = disc >= MIN_DISCORDANT
    np.divide(wo, disc, out=idx, where=sel)
    return idx, disc, wo, sel


def est_pa_single(cw, cs, rng):
    idx, disc, wo, sel = _pa_once(cw, cs, rng)
    p = np.full(T, np.nan)
    for i in np.flatnonzero(sel):
        p[i] = stats.binomtest(int(wo[i]), int(disc[i]), 0.5).pvalue
    return idx, p, sel


def est_pa_repeated(cw, cs, rng):
    sig = np.zeros(T)
    idxs = np.full((N_DRAWS_REPEATED, T), np.nan)
    testable = np.zeros(T)
    for d in range(N_DRAWS_REPEATED):
        idx, disc, wo, sel = _pa_once(cw, cs, rng)
        idxs[d] = idx
        testable += sel
        p = np.full(T, np.nan)
        for i in np.flatnonzero(sel):
            p[i] = stats.binomtest(int(wo[i]), int(disc[i]), 0.5).pvalue
        q = _bh(p)
        sig += (q < ALPHA) & np.isfinite(q)
    sel = testable >= STABLE_FRAC * N_DRAWS_REPEATED
    # q is already applied inside; express the stable rule as a pseudo-p so the
    # common BH step downstream cannot double-correct: 0 if stable-significant.
    pseudo = np.where(sig >= STABLE_FRAC * N_DRAWS_REPEATED, 0.0, 1.0)
    pseudo[~sel] = np.nan
    return np.nanmedian(idxs, axis=0), pseudo, sel


def est_prop_wilcoxon(cw, cs, rng):
    fw = cw / np.maximum(cw.sum(1, keepdims=True), 1)
    fs = cs / np.maximum(cs.sum(1, keepdims=True), 1)
    d = fw - fs
    p = np.full(T, np.nan)
    sel = (d != 0).sum(0) >= 5          # Wilcoxon needs some non-zero differences
    for i in np.flatnonzero(sel):
        p[i] = stats.wilcoxon(d[:, i], zero_method="wilcox",
                              alternative="two-sided").pvalue
    # direction only; this estimator returns no index on the 0-1 scale
    direction = np.where(np.nanmedian(d, axis=0) > 0, 1.0, 0.0)
    return np.where(sel, direction, np.nan), p, sel


def _clr(c):
    x = c + 0.5
    return np.log(x) - np.log(x).mean(1, keepdims=True)


def est_clr_ttest(cw, cs, rng):
    d = _clr(cw) - _clr(cs)
    p = np.full(T, np.nan)
    sel = np.isfinite(d).all(0) & (d.std(0) > 0)
    for i in np.flatnonzero(sel):
        p[i] = stats.ttest_rel(_clr(cw)[:, i], _clr(cs)[:, i]).pvalue
    direction = np.where(d.mean(0) > 0, 1.0, 0.0)
    return np.where(sel, direction, np.nan), p, sel


def est_betabinom(cw, cs, rng):
    """Intercept-only quasi-binomial, solved analytically.

    p_hat = sum(k)/sum(n); Pearson dispersion phi inflates the SE. Not a
    hierarchical model - the plan says so - but it is an honest overdispersed
    binomial and it returns an estimate on the same 0-1 scale as the estimand.
    """
    k = cw.astype(float)
    n = (cw + cs).astype(float)
    tot_n = n.sum(0)
    tot_k = k.sum(0)
    sel = tot_n >= MIN_DISCORDANT      # same evidence floor, on reads not pairs
    phat = np.full(T, np.nan)
    np.divide(tot_k, tot_n, out=phat, where=sel)
    p = np.full(T, np.nan)
    for i in np.flatnonzero(sel):
        ph = phat[i]
        if not (0 < ph < 1):
            p[i] = 1.0
            continue
        nj = n[:, i]
        use = nj > 0
        m = int(use.sum())
        if m < 2:
            continue
        resid = (k[use, i] - nj[use] * ph) / np.sqrt(nj[use] * ph * (1 - ph))
        phi = max(float((resid ** 2).sum() / (m - 1)), 1.0)
        se = np.sqrt(phi * ph * (1 - ph) / nj.sum())
        if se <= 0:
            continue
        p[i] = 2 * stats.norm.sf(abs(ph - 0.5) / se)
    return phat, p, sel


ESTIMATORS = {
    "pa_single": est_pa_single,
    "pa_repeated": est_pa_repeated,
    "prop_wilcoxon": est_prop_wilcoxon,
    "clr_ttest": est_clr_ttest,
    "betabinom": est_betabinom,
}
ON_INDEX_SCALE = {"pa_single", "pa_repeated", "betabinom"}

rows = []
evaluable_sets = {}
cells = list(itertools.product(DEPTHS, RATIOS))
print(f"\n=== {len(cells)} cells x {N_REP} replicates x {len(ESTIMATORS)} estimators")
for depth, ratio in cells:
    acc = {k: [] for k in ESTIMATORS}
    for rep in range(N_REP):
        seed = zlib.crc32(f"{depth}|{ratio}|{rep}".encode())
        cw, cs, true_index, prefers = simulate(ratio, depth, np.random.default_rng(seed))
        for name, fn in ESTIMATORS.items():
            idx, p, sel = fn(cw, cs, np.random.default_rng(seed + 1))
            q = _bh(p)
            called = np.isfinite(q) & (q < ALPHA)
            # direction: index > 0.5 means water, or the 1/0 flag for rank methods
            est_water = (idx > 0.5) if name in ON_INDEX_SCALE else (idx > 0.5)
            correct_dir = est_water == (true_index > 0.5)
            tp = called & prefers & correct_dir
            fp = called & (~prefers | ~correct_dir)
            power = tp.sum() / max(int((prefers & sel).sum()), 1)
            fdr = fp.sum() / max(int(called.sum()), 1)
            diracc = (called & correct_dir).sum() / max(int(called.sum()), 1)
            if name in ON_INDEX_SCALE:
                m = sel & prefers & np.isfinite(idx)
                bias = float(np.mean(idx[m] - true_index[m])) if m.any() else np.nan
                rmse = float(np.sqrt(np.mean((idx[m] - true_index[m]) ** 2))) if m.any() else np.nan
            else:
                bias = rmse = np.nan
            acc[name].append(sel)
            rows.append(dict(depth=depth, ratio=ratio, rep=rep, estimator=name,
                             n_evaluable=int(sel.sum()), power=float(power),
                             fdr=float(fdr), dir_acc=float(diracc),
                             bias=bias, rmse=rmse))
    for name in ESTIMATORS:
        evaluable_sets[(depth, ratio, name)] = np.array(acc[name]).mean(0) >= 0.5
    print(f"  depth {depth:6,d}  r={ratio}: " + "  ".join(
        f"{n}={int(np.median([r['n_evaluable'] for r in rows if r['depth']==depth and r['ratio']==ratio and r['estimator']==n]))}"
        for n in ESTIMATORS))

bm = pd.DataFrame(rows)
bm.to_csv(R / "benchmark.tsv", sep="\t", index=False)

# depth_stability: Jaccard of the evaluable set between the extreme depths
stab = {}
for name in ESTIMATORS:
    js = []
    for ratio in RATIOS:
        lo = evaluable_sets[(DEPTHS[0], ratio, name)]
        hi = evaluable_sets[(DEPTHS[-1], ratio, name)]
        u = int((lo | hi).sum())
        js.append(int((lo & hi).sum()) / u if u else np.nan)
    stab[name] = float(np.median(js))

med = bm.groupby(["estimator", "depth", "ratio"]).median(numeric_only=True)
print(f"\n{'estimator':15s} {'evaluable':>10s} {'power':>7s} {'FDR':>7s} "
      f"{'dir':>6s} {'bias':>8s} {'rmse':>7s} {'stability':>10s}")
for name in ESTIMATORS:
    s = bm[bm.estimator == name]
    print(f"{name:15s} {s.n_evaluable.median():10.0f} {s.power.median():7.3f} "
          f"{s.fdr.median():7.3f} {s.dir_acc.median():6.3f} "
          f"{s.bias.median():8.3f} {s.rmse.median():7.3f} {stab[name]:10.3f}")

verdicts = {}


def verdict(name, text, confirmed, failed, value):
    v = "CONFIRMED" if confirmed else ("FAILED" if failed else "INCONCLUSIVE")
    verdicts[name] = dict(prediction=text, value=value, verdict=v)
    print(f"  {name}  {v:13s} {text}\n        observed: {value}")


COUNT_BASED = ["prop_wilcoxon", "clr_ttest", "betabinom"]
lo_d = DEPTHS[0]
ev_lo = {n: bm[(bm.estimator == n) & (bm.depth == lo_d)].n_evaluable.median()
         for n in ESTIMATORS}
best_cb = max(COUNT_BASED, key=lambda n: ev_lo[n])
ratio_ev = ev_lo[best_cb] / max(ev_lo["pa_single"], 1)

print("\n=== pre-committed predictions")
verdict("B1", "count-based estimators evaluate far more orders",
        ratio_ev >= 1.5, ratio_ev < 1.1,
        f"at depth {lo_d}: {best_cb} evaluates {ev_lo[best_cb]:.0f} vs "
        f"pa_single {ev_lo['pa_single']:.0f} = {ratio_ev:.2f}x")

best_stab = max(COUNT_BASED, key=lambda n: stab[n])
verdict("B2", "and their evaluable set is far less depth-dependent",
        stab[best_stab] >= 0.80 and stab["pa_single"] <= 0.50,
        stab[best_stab] < stab["pa_single"],
        f"{best_stab} stability {stab[best_stab]:.3f} vs "
        f"pa_single {stab['pa_single']:.3f}")

fdr_by_depth = {n: {d: bm[(bm.estimator == n) & (bm.depth == d)].fdr.median()
                    for d in DEPTHS} for n in ESTIMATORS}
diffs = {d: fdr_by_depth["pa_repeated"][d] - fdr_by_depth["pa_single"][d] for d in DEPTHS}
verdict("B3", "repeated rarefaction controls FDR better than a single draw",
        all(v <= -0.02 for v in diffs.values()), any(v > 0 for v in diffs.values()),
        "repeated - single FDR by depth: " +
        ", ".join(f"{d}: {v:+.3f}" for d, v in diffs.items()))

cell = bm[(bm.depth == 1000) & (bm.ratio == 4)]
base = cell[cell.estimator == "pa_single"]
beats = [n for n in ESTIMATORS if n != "pa_single"
         and cell[cell.estimator == n].power.median() > base.power.median()
         and cell[cell.estimator == n].fdr.median() <= base.fdr.median()]
verdict("B4", "some estimator beats pa_single on power at equal or lower FDR",
        bool(beats), not beats,
        f"at depth 1000, r=4: pa_single power {base.power.median():.3f} "
        f"FDR {base.fdr.median():.3f}; beaten by {beats or 'none'}")

worst = bm.groupby(["estimator", "depth", "ratio"]).fdr.median()
all_ok = {n: bool((worst.loc[n] <= 0.05).all()) for n in ESTIMATORS}
any_cell_all_fail = False
for d, r in cells:
    if all(bm[(bm.estimator == n) & (bm.depth == d) & (bm.ratio == r)].fdr.median() > 0.05
           for n in ESTIMATORS):
        any_cell_all_fail = True
verdict("B5", "no estimator achieves nominal FDR everywhere",
        any_cell_all_fail, all(all_ok.values()),
        "holds FDR<=0.05 in every cell: " +
        ", ".join(f"{n}={'yes' if v else 'no'}" for n, v in all_ok.items()))

out = dict(plan_sha256=plan_sha, n_pairs=N_PAIRS, n_orders=T, n_rep=N_REP,
           depths=DEPTHS, ratios=RATIOS,
           per_estimator={n: dict(
               n_evaluable=float(bm[bm.estimator == n].n_evaluable.median()),
               power=float(bm[bm.estimator == n].power.median()),
               fdr=float(bm[bm.estimator == n].fdr.median()),
               dir_acc=float(bm[bm.estimator == n].dir_acc.median()),
               bias=(None if not np.isfinite(bm[bm.estimator == n].bias.median())
                     else float(bm[bm.estimator == n].bias.median())),
               rmse=(None if not np.isfinite(bm[bm.estimator == n].rmse.median())
                     else float(bm[bm.estimator == n].rmse.median())),
               depth_stability=stab[n]) for n in ESTIMATORS},
           verdicts=verdicts)
(R / "benchmark.json").write_text(json.dumps(out, indent=1))
n_fail = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
print(f"\nwrote results/benchmark.{{tsv,json}} — {n_fail} of {len(verdicts)} FAILED")
