#!/usr/bin/env python3
"""Where the threshold mechanism operates, as SWEEP_PLAN.md specifies.

The plan's sha256 is recorded in PLAN_HASHES.md, not copied here; this script
asserts at startup that the plan on disk still matches that record.

Two questions. Descriptively, over what effect sizes, abundance shapes, thresholds
and pair counts does a minimum-evidence threshold re-select taxa with depth.
Predictively, and this is the point: fed the EMPIRICAL abundance vector, does the
detection model reproduce the exponent the NEON depth ladder measured (b = 0.485)
without being told it.
"""
import hashlib
import itertools
import json
import pathlib
import re
import sys
import zlib

import numpy as np
import pandas as pd

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
P = REPO_ROOT
R = P / "results"

PLAN = P / "plans/SWEEP_PLAN.md"
plan_sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
rec = re.search(r"`SWEEP_PLAN\.md` \| `([0-9a-f]{64})`", (P / "plans/PLAN_HASHES.md").read_text())
assert rec, "the sweep plan's hash is not recorded in PLAN_HASHES.md"
assert plan_sha == rec.group(1), (
    f"SWEEP_PLAN.md changed since freezing:\n  now      {plan_sha}\n"
    f"  recorded {rec.group(1)}")
print(f"plan verified: {plan_sha[:16]}…")

# ------------------------------------------------- plan section 3: the grid, fixed
RUNGS = [100, 300, 1_000, 3_000, 10_000]
RATIOS = [1, 2, 4, 8, 16]
DISTS = ["loguniform", "lognormal", "powerlaw", "empirical"]
THRESHOLDS = [10, 20, 40]
NPAIRS = [102, 239]
SIGMAS = [0.0, 1.0]
N_DRAWS = 30
N_TAXA = 254
MEASURED_B = 0.485          # from results/depth_ladder.json, section 3.4a

emp = pd.read_csv(R / "sweep_empirical_abundance.tsv", sep="\t")
EMPIRICAL = emp.a_hat.to_numpy()
assert len(EMPIRICAL) == N_TAXA, f"expected {N_TAXA} empirical abundances, got {len(EMPIRICAL)}"
print(f"empirical abundances: {len(EMPIRICAL)}  "
      f"range {EMPIRICAL.min():.2e}-{EMPIRICAL.max():.2e}  "
      f"median {np.median(EMPIRICAL):.2e}")


def abundances(kind, rng):
    """One community of N_TAXA relative abundances."""
    if kind == "empirical":
        return EMPIRICAL.copy()
    if kind == "loguniform":                      # the section 3.4 community
        return np.geomspace(1e-6, 1e-1, N_TAXA)
    if kind == "lognormal":                       # sigma = 2 in log10
        a = 10.0 ** rng.normal(-3.0, 2.0, N_TAXA)
        return np.clip(a, 1e-9, 1.0)
    if kind == "powerlaw":                        # Pareto, alpha = 1.5, heavy-tailed
        a = (rng.pareto(1.5, N_TAXA) + 1) * 1e-5
        return np.clip(a, 1e-9, 1.0)
    raise ValueError(kind)


def one_draw(a, ratio, depth, n_pairs, sigma, threshold, rng):
    """Simulate n_pairs water/sediment pairs; return the per-taxon selection."""
    t = len(a)
    if sigma > 0:
        # site heterogeneity: a lognormal multiplier per pair per taxon
        mult = np.exp(rng.normal(0.0, sigma, (n_pairs, t)))
    else:
        mult = 1.0
    lam_s = a * depth * mult
    lam_w = lam_s * ratio
    p_s = 1.0 - np.exp(-lam_s)
    p_w = 1.0 - np.exp(-lam_w)
    p_s = np.broadcast_to(p_s, (n_pairs, t))
    p_w = np.broadcast_to(p_w, (n_pairs, t))
    w = rng.random((n_pairs, t)) < p_w
    s = rng.random((n_pairs, t)) < p_s
    wo = (w & ~s).sum(0)
    so = (~w & s).sum(0)
    disc = wo + so
    sel = disc >= threshold
    idx = np.full(t, np.nan)
    np.divide(wo, disc, out=idx, where=sel)
    return sel, idx


def fit_exponent(depths, values):
    """b in  value ~ D^-b, plus R^2, in log-log. nan where it cannot be fitted."""
    d = np.asarray(depths, float)
    v = np.asarray(values, float)
    ok = np.isfinite(v) & (v > 0)
    if ok.sum() < 3:
        return np.nan, np.nan
    x, y = np.log10(d[ok]), np.log10(v[ok])
    m, c = np.polyfit(x, y, 1)
    ss = np.sum((y - y.mean()) ** 2)
    r2 = 1 - np.sum((y - (m * x + c)) ** 2) / ss if ss > 0 else np.nan
    return -m, r2


cells = list(itertools.product(DISTS, RATIOS, THRESHOLDS, NPAIRS, SIGMAS))
print(f"\n=== {len(cells)} cells x {len(RUNGS)} rungs x {N_DRAWS} draws")

rows = []
for n, (dist, ratio, thr, npair, sigma) in enumerate(cells, 1):
    # Seeds shared across rungs within a cell, so depth differences are not draw luck.
    # crc32, not hash(): Python randomises string hashing per process, so hash()
    # would make this script irreproducible between runs - the one property the
    # depth ladder was verified to have.
    tag = f"{dist}|{ratio}|{thr}|{npair}|{sigma}".encode()
    seeds = [zlib.crc32(tag + b"|" + str(d).encode()) for d in range(N_DRAWS)]
    a = abundances(dist, np.random.default_rng(12345))
    per_rung = {}
    sel_at = {}
    for depth in RUNGS:
        med_a, med_lam, med_eff, n_test, sels = [], [], [], [], []
        for sd in seeds:
            rng = np.random.default_rng(sd)
            sel, idx = one_draw(a, ratio, depth, npair, sigma, thr, rng)
            sels.append(sel)
            n_test.append(int(sel.sum()))
            if sel.any():
                med_a.append(float(np.median(a[sel])))
                med_lam.append(float(np.median(a[sel] * depth)))
                med_eff.append(float(np.median(np.abs(idx[sel] - 0.5))))
        per_rung[depth] = dict(
            n_testable=float(np.median(n_test)),
            median_a=float(np.median(med_a)) if med_a else np.nan,
            median_lambda=float(np.median(med_lam)) if med_lam else np.nan,
            median_effect=float(np.median(med_eff)) if med_eff else np.nan)
        sel_at[depth] = np.array(sels).mean(0) >= 0.5      # majority of draws

    b, r2 = fit_exponent(RUNGS, [per_rung[d]["median_a"] for d in RUNGS])
    lam = np.array([per_rung[d]["median_lambda"] for d in RUNGS], float)
    eff = np.array([per_rung[d]["median_effect"] for d in RUNGS], float)
    lo, hi = sel_at[RUNGS[0]], sel_at[RUNGS[-1]]
    union = int((lo | hi).sum())
    nt0, nt1 = per_rung[RUNGS[0]]["n_testable"], per_rung[RUNGS[-1]]["n_testable"]

    rows.append(dict(
        dist=dist, ratio=ratio, threshold=thr, n_pairs=npair, sigma=sigma,
        b=b, r2=r2,
        lambda_span=(np.nanmax(lam) / np.nanmin(lam)) if np.isfinite(lam).all()
                    and np.nanmin(lam) > 0 else np.nan,
        jaccard=(int((lo & hi).sum()) / union) if union else np.nan,
        effect_range=(np.nanmax(eff) - np.nanmin(eff)) if np.isfinite(eff).any() else np.nan,
        growth=(nt1 / nt0) if nt0 > 0 else np.nan,
        n_testable_lo=nt0, n_testable_hi=nt1))
    if n % 40 == 0:
        print(f"  {n}/{len(cells)} cells")

sw = pd.DataFrame(rows)
sw.to_csv(R / "sweep.tsv", sep="\t", index=False)
print(f"wrote results/sweep.tsv ({len(sw)} cells)")


def pick(**kw):
    q = sw
    for k, v in kw.items():
        q = q[q[k] == v]
    return q


print("\n=== exponent b by abundance distribution (threshold 20, 102 pairs, sigma 0)")
print(f"  {'dist':12s} " + "  ".join(f"r={r:<6d}" for r in RATIOS))
for dist in DISTS:
    cellrow = [pick(dist=dist, threshold=20, n_pairs=102, sigma=0.0, ratio=r)
               for r in RATIOS]
    bs = "  ".join(f"{c.b.iloc[0]:+.3f}  " if len(c) and np.isfinite(c.b.iloc[0])
                   else "   n/a   " for c in cellrow)
    print(f"  {dist:12s} {bs}")

print("\n=== R^2 of those fits")
for dist in DISTS:
    r2s = "  ".join(
        f"{pick(dist=dist, threshold=20, n_pairs=102, sigma=0.0, ratio=r).r2.iloc[0]:.3f}  "
        for r in RATIOS)
    print(f"  {dist:12s} {r2s}")

print("\n=== jaccard(100, 10000) vs effect size (empirical, threshold 20, 102 pairs, sigma 0)")
for r in RATIOS:
    c = pick(dist="empirical", threshold=20, n_pairs=102, sigma=0.0, ratio=r)
    print(f"  r={r:<3d} jaccard {c.jaccard.iloc[0]:.3f}  growth {c.growth.iloc[0]:.2f}  "
          f"testable {c.n_testable_lo.iloc[0]:.0f} -> {c.n_testable_hi.iloc[0]:.0f}")

# ------------------------------------------------- plan section 5: the verdicts
verdicts = {}


def verdict(name, text, confirmed, failed, value):
    v = "CONFIRMED" if confirmed else ("FAILED" if failed else "INCONCLUSIVE")
    verdicts[name] = dict(prediction=text, value=value, verdict=v)
    print(f"  {name}  {v:13s} {text}\n        observed: {value}")


print("\n=== pre-committed predictions")

# S1 - shape sets the exponent
diffs = []
for r in (2, 4):
    bu = pick(dist="loguniform", ratio=r, threshold=20, n_pairs=102, sigma=0.0).b.iloc[0]
    be = pick(dist="empirical", ratio=r, threshold=20, n_pairs=102, sigma=0.0).b.iloc[0]
    diffs.append(bu - be)
dmin = float(np.nanmin(diffs))
verdict("S1", "the exponent is set by the abundance distribution's shape",
        dmin >= 0.25, (dmin < 0.10) or (dmin < 0),
        f"b(loguniform) - b(empirical) = {diffs[0]:+.3f} at r=2, {diffs[1]:+.3f} at r=4")

# S2 - the predictive test
emp_b = {r: pick(dist="empirical", ratio=r, threshold=20, n_pairs=102,
                 sigma=0.0).b.iloc[0] for r in RATIOS}
in_band = [r for r in (2, 4, 8) if abs(emp_b[r] - MEASURED_B) <= 0.15]
nearest = min(RATIOS, key=lambda r: abs(emp_b[r] - MEASURED_B))
verdict("S2", f"the empirical vector reproduces the measured b = {MEASURED_B}",
        bool(in_band),
        all(abs(emp_b[r] - MEASURED_B) > 0.25 for r in RATIOS),
        f"b(empirical) = " + ", ".join(f"{emp_b[r]:+.3f} (r={r})" for r in RATIOS)
        + f"; within 0.15 for r in {in_band or 'none'}; nearest r={nearest} "
          f"at {emp_b[nearest]:+.3f}, off by {abs(emp_b[nearest]-MEASURED_B):.3f}")

# S3 - large effects escape the selection
jac = [pick(dist="empirical", ratio=r, threshold=20, n_pairs=102,
            sigma=0.0).jaccard.iloc[0] for r in RATIOS]
rising = all(b >= a - 1e-9 for a, b in zip(jac, jac[1:]))
verdict("S3", "large effects escape the lambda~1 selection",
        rising and (jac[-1] - jac[0] >= 0.15),
        (jac[-1] <= jac[0]),
        f"jaccard across r=1..16: " + ", ".join(f"{j:.3f}" for j in jac)
        + f"; monotone rising: {rising}; r16-r1 = {jac[-1]-jac[0]:+.3f}")

# S4 - threshold sets the level, not the exponent
bs = [pick(dist="empirical", ratio=4, threshold=t, n_pairs=102, sigma=0.0).b.iloc[0]
      for t in THRESHOLDS]
js = [pick(dist="empirical", ratio=4, threshold=t, n_pairs=102, sigma=0.0).jaccard.iloc[0]
      for t in THRESHOLDS]
b_var = float(np.nanmax(bs) - np.nanmin(bs))
j_var = float(np.nanmax(js) - np.nanmin(js))
verdict("S4", "the threshold sets the level, not the exponent",
        (b_var <= 0.15) and (j_var > b_var), b_var > 0.30,
        f"b across thresholds 10/20/40: " + ", ".join(f"{x:+.3f}" for x in bs)
        + f" (spread {b_var:.3f}); jaccard spread {j_var:.3f}")

# S5 - depth-invariance of apparent effect is general
good = sw[sw.r2 >= 0.5]
frac = float((good.effect_range <= 0.10).mean()) if len(good) else np.nan
verdict("S5", "depth-invariance of apparent effect is general",
        frac >= 0.80, frac <= 0.50,
        f"{frac:.0%} of {len(good)} cells with R^2>=0.5 have effect_range <= 0.10")

out = dict(plan_sha256=plan_sha, measured_b=MEASURED_B, n_cells=len(sw),
           n_draws=N_DRAWS, rungs=RUNGS,
           empirical_b_by_ratio={int(r): float(emp_b[r]) for r in RATIOS},
           verdicts=verdicts)
(R / "sweep.json").write_text(json.dumps(out, indent=1))
n_fail = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
print(f"\nwrote results/sweep.json — {n_fail} of {len(verdicts)} predictions FAILED")
if n_fail:
    print("A failure is a result. Report it; do not adjust the plan.")
