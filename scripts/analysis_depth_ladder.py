#!/usr/bin/env python3
"""Controlled depth experiment inside NEON, exactly as DEPTH_LADDER_PLAN.md specifies.

The plan this implements is DEPTH_LADDER_PLAN.md. Its sha256 is recorded in
PLAN_HASHES.md rather than copied here, and this script asserts at startup that the
plan on disk still matches that record, so a plan edited after freezing stops the
run instead of silently changing what was promised.

Ecology, samples, reference database, pairing rule and test are held fixed; only the
rarefaction depth changes. This is the test the FunAqua-vs-NEON Jaccard cannot be:
that comparison confounds depth with geography, water-body type and what "sediment"
means, so it cannot say whether depth re-selects the taxa a survey can measure.
"""
import ast
import collections
import hashlib
import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
P = REPO_ROOT
R = P / "results"

# ----------------------------------------------------------------- frozen plan
PLAN = P / "plans/DEPTH_LADDER_PLAN.md"
plan_sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
recorded = re.search(r"`DEPTH_LADDER_PLAN\.md` \| `([0-9a-f]{64})`",
                     (P / "plans/PLAN_HASHES.md").read_text())
assert recorded, "the plan's hash is not recorded in PLAN_HASHES.md"
assert plan_sha == recorded.group(1), (
    f"DEPTH_LADDER_PLAN.md has changed since it was frozen:\n"
    f"  now      {plan_sha}\n  recorded {recorded.group(1)}")
print(f"plan verified against its frozen hash: {plan_sha[:16]}…")

# --------------------------------------------- plan section 3: the ladder, fixed
RUNGS = [100, 300, 1_000, 3_000, 10_000]
SUBSET_MIN_DEPTH = 10_000      # plan 2
MIN_DISCORDANT = 20            # plan 3
N_DRAWS = 50                   # plan 3
SEEDS = [20261023 + i for i in range(N_DRAWS)]   # same seeds at every rung

# ------------------------------------------------- borrow the frozen NEON loader
# analysis_neon_validation.py runs the whole validation at import, and it is not
# ours to restructure: its outputs are published and a sibling script once gave a
# wrong answer because of a module-level RNG consumed across calls. So take only
# its imports, constants and function definitions, and execute nothing else.
src = (P / "scripts/analysis_neon_validation.py").read_text()
tree = ast.parse(src)
# Constants sit before the first function; every top-level assignment AFTER the
# functions is the validation actually running (water = load(...), and so on), so
# those must be dropped. Keeping all Assign nodes ran the loaders and then died on
# an undefined name - the filter has to be positional, not just by node type.
first_def = next(i for i, n in enumerate(tree.body)
                 if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)))
keep = [n for i, n in enumerate(tree.body)
        if isinstance(n, (ast.Import, ast.ImportFrom))
        or isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
        or (isinstance(n, ast.Assign) and i < first_def)]
ns: dict = {"__name__": "neon_io_borrowed"}
exec(compile(ast.Module(body=keep, type_ignores=[]), "<neon_io>", "exec"), ns)
load, build_pairs = ns["load"], ns["build_pairs"]
MIN_CONF, MIN_READS, WINDOW_DAYS = ns["MIN_CONF"], ns["MIN_READS"], ns["WINDOW_DAYS"]
print(f"borrowed loader: MIN_CONF={MIN_CONF} MIN_READS={MIN_READS} "
      f"WINDOW_DAYS={WINDOW_DAYS}")

print("\n=== loading NEON ITS tables")
water = load("water", MIN_CONF)
benthic = load("benthic", MIN_CONF)
pairs_all = build_pairs(water, benthic, WINDOW_DAYS)
print(f"candidate pairs under the primary rule: {len(pairs_all)}")

# ------------------------------------------- plan section 2: the fixed 102 pairs
# Take the pair set from the frozen validation's own output, so this experiment
# cannot quietly re-pair or re-filter, then keep those deep enough for the top rung.
frozen = pd.read_csv(R / "neon_pairs.tsv", sep="\t")
frozen_keys = {(r.site, r.water_sample, r.benthic_sample) for r in frozen.itertuples()}
depth_of = dict(zip(zip(frozen.site, frozen.water_sample, frozen.benthic_sample),
                    frozen.rarefied_depth))

pairs = [p for p in pairs_all if (p[0], p[2], p[3]) in frozen_keys
         and depth_of[(p[0], p[2], p[3])] >= SUBSET_MIN_DEPTH]
print(f"pairs in the frozen validation set: "
      f"{sum(1 for p in pairs_all if (p[0], p[2], p[3]) in frozen_keys)}")
print(f"of those, within-pair depth >= {SUBSET_MIN_DEPTH:,}: {len(pairs)}")

# The plan surveyed 102; a different number means the pair set moved underneath us.
assert len(pairs) == 102, (
    f"the plan fixed a 102-pair subset; this run finds {len(pairs)}. "
    f"Do not proceed - establish why before computing anything.")

# Every rung must be reachable for every pair, or a pair would be upsampled.
avail = np.array([min(sum(water[p[2]]["otu"].values()),
                      sum(benthic[p[3]]["otu"].values())) for p in pairs])
assert avail.min() >= max(RUNGS), (
    f"a pair carries only {avail.min()} reads, below the top rung {max(RUNGS)}")
print(f"available within-pair depth: min {avail.min():,} "
      f"median {int(np.median(avail)):,} max {avail.max():,}")

pd.DataFrame([{"site": p[0], "water_date": p[1], "water_sample": p[2],
               "benthic_sample": p[3], "day_offset": p[4],
               "available_depth": int(a)}
              for p, a in zip(pairs, avail)]).to_csv(
    R / "depth_ladder_pairs.tsv", sep="\t", index=False)
print(f"wrote the frozen subset to results/depth_ladder_pairs.tsv")

# ------------------------------------------------------------ order bookkeeping
orders = sorted({o for s in list(water.values()) + list(benthic.values())
                 for _, o in s["otu"]})
oidx = {o: i for i, o in enumerate(orders)}
print(f"orders present across NEON: {len(orders)}")


def vec(sample):
    ks = list(sample["otu"])
    return (np.array([oidx[o] for _, o in ks]),
            np.array([sample["otu"][k] for k in ks], dtype=np.int64))


VECS = [(vec(water[p[2]]), vec(benthic[p[3]])) for p in pairs]

# ---- plan section 4: a-hat, estimated ONCE at full depth and held fixed.
# For a fixed order lambda-hat = a-hat * D is proportional to D by construction.
# Every prediction below is about the MEDIAN over the set that passes the
# threshold, whose membership is free to change. The plan says so explicitly.
tot = np.zeros(len(orders))
denom = 0
for (kw, cw), (kb, cb) in VECS:
    tot += np.bincount(kw, weights=cw, minlength=len(orders))
    tot += np.bincount(kb, weights=cb, minlength=len(orders))
    denom += cw.sum() + cb.sum()
a_hat = tot / denom
print(f"a-hat estimated at full depth over {len(pairs)} pairs, "
      f"{denom:,} reads; nonzero orders {int((a_hat > 0).sum())}")

# Exported for SWEEP_PLAN.md, which asks whether this vector, fed to the detection
# model, reproduces the exponent this script measures. Adding an output changes no
# computation here; the run is deterministic on fixed seeds, and the commit adding
# it records that results/depth_ladder.json came back byte-identical.
pd.DataFrame({"order": orders, "a_hat": a_hat})[a_hat > 0].sort_values(
    "a_hat", ascending=False).to_csv(
    R / "sweep_empirical_abundance.tsv", sep="\t", index=False)


def one_draw(rung, seed):
    """Rarefy both members of every pair to `rung`; return the per-order table."""
    rng = np.random.default_rng(seed)
    W = np.zeros((len(VECS), len(orders)), dtype=bool)
    S = np.zeros((len(VECS), len(orders)), dtype=bool)
    for i, ((kw, cw), (kb, cb)) in enumerate(VECS):
        rw = rng.multivariate_hypergeometric(cw, rung)
        rb = rng.multivariate_hypergeometric(cb, rung)
        W[i] = np.bincount(kw, weights=rw, minlength=len(orders)) > 0
        S[i] = np.bincount(kb, weights=rb, minlength=len(orders)) > 0
    wo = (W & ~S).sum(0)
    bo = (~W & S).sum(0)
    disc = wo + bo
    sel = disc >= MIN_DISCORDANT
    idx = np.full(len(orders), np.nan)
    idx[sel] = wo[sel] / disc[sel]
    return sel, idx


print(f"\n=== the ladder: {RUNGS}, {N_DRAWS} draws each, same seeds at every rung")
rows = []
# sel_sets[rung][draw] = boolean mask of testable orders; idx likewise
sel_by_rung, idx_by_rung = {}, {}
for rung in RUNGS:
    sels, idxs = [], []
    for seed in SEEDS:
        sel, idx = one_draw(rung, seed)
        sels.append(sel)
        idxs.append(idx)
        eff = np.abs(idx[sel] - 0.5)
        rows.append(dict(rung=rung, seed=seed, n_testable=int(sel.sum()),
                         median_abundance=float(np.median(a_hat[sel])),
                         median_lambda=float(np.median(a_hat[sel] * rung)),
                         median_effect=float(np.median(eff))))
    sel_by_rung[rung] = np.array(sels)
    idx_by_rung[rung] = np.array(idxs)
    t = pd.DataFrame([r for r in rows if r["rung"] == rung])
    print(f"  rung {rung:6,d}: testable {t.n_testable.median():5.0f}   "
          f"median a-hat {t.median_abundance.median():.3e}   "
          f"median lambda {t.median_lambda.median():7.3f}   "
          f"median |index-0.5| {t.median_effect.median():.3f}")

lad = pd.DataFrame(rows)
lad.to_csv(R / "depth_ladder.tsv", sep="\t", index=False)

# ------------------------------------------------- turnover and persistent taxa
top = max(RUNGS)
# A order counts as testable at a rung if it is testable in the majority of draws,
# so turnover is not just draw noise.
maj = {r: sel_by_rung[r].mean(0) >= 0.5 for r in RUNGS}
jac = {}
for r in RUNGS:
    a, b = maj[r], maj[top]
    inter = int((a & b).sum())
    union = int((a | b).sum())
    jac[r] = inter / union if union else np.nan

persistent = np.ones(len(orders), dtype=bool)
for r in RUNGS:
    persistent &= maj[r]
print(f"\norders testable (majority of draws) at every rung: {int(persistent.sum())}")

# median |index(rung) - index(top)| over persistent orders, per draw then median
delta = {}
for r in RUNGS:
    per_draw = []
    for d in range(N_DRAWS):
        a = idx_by_rung[r][d][persistent]
        b = idx_by_rung[top][d][persistent]
        ok = np.isfinite(a) & np.isfinite(b)
        if ok.sum():
            per_draw.append(float(np.median(np.abs(a[ok] - b[ok]))))
    delta[r] = float(np.median(per_draw)) if per_draw else np.nan

print(f"\n{'rung':>8s} {'jaccard vs top':>15s} {'median |d index|':>18s}")
for r in RUNGS:
    print(f"{r:8,d} {jac[r]:15.3f} {delta[r]:18.3f}")

# ------------------------------------------------- plan section 5: the verdicts
med = lad.groupby("rung").median(numeric_only=True)
lo, hi = min(RUNGS), max(RUNGS)

ab_fall = med.loc[lo, "median_abundance"] / med.loc[hi, "median_abundance"]
lam_span = med["median_lambda"].max() / med["median_lambda"].min()
eff_range = med["median_effect"].max() - med["median_effect"].min()

verdicts = {}


def verdict(name, text, confirmed, failed, value):
    if confirmed:
        v = "CONFIRMED"
    elif failed:
        v = "FAILED"
    else:
        v = "INCONCLUSIVE"
    verdicts[name] = {"prediction": text, "value": value, "verdict": v}
    print(f"  {name}  {v:13s} {text}")
    print(f"        observed: {value}")


print("\n=== pre-committed predictions")
verdict("P1", "testable set moves down the abundance distribution",
        ab_fall >= 10, (ab_fall < 2) or (ab_fall < 1),
        f"median a-hat falls {ab_fall:.1f}-fold from rung {lo} to {hi}")
verdict("P2", "testable set holds its position on the detection curve",
        lam_span <= 3, lam_span > 10,
        f"median lambda spans {lam_span:.2f}-fold across rungs "
        f"({med['median_lambda'].min():.3f}-{med['median_lambda'].max():.3f})")
verdict("P3", "apparent effect size is roughly depth-invariant",
        eff_range <= 0.10, eff_range > 0.10,
        f"median |index-0.5| range {eff_range:.3f} "
        f"({med['median_effect'].min():.3f}-{med['median_effect'].max():.3f})")
verdict("P4", "the testable set turns over",
        jac[lo] < 0.6, jac[lo] > 0.8,
        f"jaccard(rung {lo}, rung {hi}) = {jac[lo]:.3f}")
verdict("P5", "depth changes membership rather than a given taxon's estimate",
        delta[lo] <= 0.10, delta[lo] > 0.20,
        f"median |index({lo}) - index({hi})| over "
        f"{int(persistent.sum())} persistent orders = {delta[lo]:.3f}")

out = {
    "plan_sha256": plan_sha,
    "rungs": RUNGS,
    "n_pairs": len(pairs),
    "n_draws": N_DRAWS,
    "min_discordant": MIN_DISCORDANT,
    "per_rung": {int(r): {
        "n_testable": float(med.loc[r, "n_testable"]),
        "median_abundance": float(med.loc[r, "median_abundance"]),
        "median_lambda": float(med.loc[r, "median_lambda"]),
        "median_effect": float(med.loc[r, "median_effect"]),
        "jaccard_vs_top": float(jac[r]),
        "delta_index_persistent": float(delta[r]),
    } for r in RUNGS},
    "n_persistent_orders": int(persistent.sum()),
    "verdicts": verdicts,
}
(R / "depth_ladder.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/depth_ladder.{tsv,json} and depth_ladder_pairs.tsv")

n_fail = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
print(f"\n{len(verdicts) - n_fail} of {len(verdicts)} predictions not failed; "
      f"{n_fail} FAILED")
if n_fail:
    print("A failure is a result. Report it; do not adjust the plan.")
