#!/usr/bin/env python3
"""Does the estimator ranking depend on how preference is encoded?

As BENCHMARK2_PLAN.md specifies; its sha256 is in PLAN_HASHES.md and is asserted
before anything runs.

The first benchmark could not answer "which estimator is best" because its
generator encoded preference as a read ratio, which one estimator measures and the
others cannot see. This asks the question that has an answer: whether the ranking
changes with the encoding.

NOTE ON THE PLAN'S CELL COUNT. Section 3 lists three arms x three depths x three
ratios and then says "81 cells". That product is 27. The axes are unambiguous, so
this runs the 27 cells they define; the 81 is an arithmetic slip in the frozen
plan's prose and is recorded in BENCHMARK2_RESULTS.md rather than corrected there.
"""
import ast
import hashlib
import itertools
import json
import pathlib
import re
import warnings
import zlib

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=RuntimeWarning)

P = pathlib.Path.home() / "research/projects/funaqua-habitat"
R = P / "results"

PLAN = P / "BENCHMARK2_PLAN.md"
plan_sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
rec = re.search(r"`BENCHMARK2_PLAN\.md` \| `([0-9a-f]{64})`", (P / "PLAN_HASHES.md").read_text())
assert rec and plan_sha == rec.group(1), "BENCHMARK2_PLAN.md changed since freezing"
print(f"plan verified: {plan_sha[:16]}…")

# Borrow the five estimators and the loaded abundance/occupancy from the first
# benchmark, without running it: imports, pre-definition assignments and function
# definitions only.
src = (P / "analysis_benchmark.py").read_text()
tree = ast.parse(src)
first_def = next(i for i, n in enumerate(tree.body)
                 if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)))
keep = [n for i, n in enumerate(tree.body)
        if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.AsyncFunctionDef))
        or (isinstance(n, ast.Assign) and i < first_def)]
ns: dict = {"__name__": "benchmark_borrowed", "__file__": str(P / "analysis_benchmark.py")}
exec(compile(ast.Module(body=keep, type_ignores=[]), "<bench_io>", "exec"), ns)

OCC, A_IF_PRESENT, T = ns["OCC"], ns["A_IF_PRESENT"], ns["T"]
N_PAIRS, ALPHA, FRAC_PREFER = ns["N_PAIRS"], ns["ALPHA"], ns["FRAC_PREFER"]
_bh = ns["_bh"]
ESTIMATORS = {k: ns[f"est_{k}"] for k in
              ("pa_single", "pa_repeated", "prop_wilcoxon", "clr_ttest", "betabinom")}
print(f"borrowed {len(ESTIMATORS)} estimators; orders {T}, pairs {N_PAIRS}, "
      f"occupancy median {np.median(OCC):.3f}")

ARMS = ["abundance", "occupancy", "both"]
DEPTHS = [100, 1_000, 10_000]
RATIOS = [2, 4, 8]
N_REP = 30


def odds_scale(p, r):
    """Scale the odds of presence by r, keeping the result a probability."""
    return p * r / (1.0 + p * (r - 1.0))


def simulate(arm, ratio, depth, rng):
    prefers = rng.random(T) < FRAC_PREFER
    r = np.where(prefers, float(ratio), 1.0)

    occ_w = odds_scale(OCC, r) if arm in ("occupancy", "both") else OCC
    ab_w = A_IF_PRESENT * (r if arm in ("abundance", "both") else 1.0)

    pw = rng.random((N_PAIRS, T)) < occ_w
    ps = rng.random((N_PAIRS, T)) < OCC
    cw = rng.poisson(np.where(pw, ab_w, 0.0) * depth)
    cs = rng.poisson(np.where(ps, A_IF_PRESENT, 0.0) * depth)
    return cw, cs, prefers


rows = []
cells = list(itertools.product(ARMS, DEPTHS, RATIOS))
print(f"\n=== {len(cells)} cells x {N_REP} replicates x {len(ESTIMATORS)} estimators")
for arm, depth, ratio in cells:
    for rep in range(N_REP):
        seed = zlib.crc32(f"{arm}|{depth}|{ratio}|{rep}".encode())
        cw, cs, prefers = simulate(arm, ratio, depth, np.random.default_rng(seed))
        for name, fn in ESTIMATORS.items():
            idx, p, sel = fn(cw, cs, np.random.default_rng(seed + 1))
            q = _bh(p)
            called = np.isfinite(q) & (q < ALPHA)
            est_water = idx > 0.5
            correct = est_water == prefers        # every preferring taxon prefers WATER
            tp = called & prefers & correct
            fp = called & (~prefers | ~correct)
            denom = max(int((prefers & sel).sum()), 1)
            rows.append(dict(
                arm=arm, depth=depth, ratio=ratio, rep=rep, estimator=name,
                n_evaluable=int(sel.sum()),
                power=float(tp.sum() / denom),
                fdr=float(fp.sum() / max(int(called.sum()), 1)),
                dir_acc=float((called & correct).sum() / max(int(called.sum()), 1))))
    done = len([r for r in rows if r["arm"] == arm]) // (N_REP * len(ESTIMATORS))
    if done == len(DEPTHS) * len(RATIOS):
        print(f"  arm {arm} done")

bm = pd.DataFrame(rows)
bm.to_csv(R / "benchmark2.tsv", sep="\t", index=False)

print(f"\n{'arm':11s} {'estimator':15s} {'evaluable':>10s} {'power':>7s} "
      f"{'FDR':>7s} {'dir':>7s}")
for arm in ARMS:
    for name in ESTIMATORS:
        s = bm[(bm.arm == arm) & (bm.estimator == name)]
        print(f"{arm:11s} {name:15s} {s.n_evaluable.median():10.0f} "
              f"{s.power.median():7.3f} {s.fdr.median():7.3f} {s.dir_acc.median():7.3f}")

best = {arm: max(ESTIMATORS, key=lambda n: bm[(bm.arm == arm) &
                                              (bm.estimator == n)].power.median())
        for arm in ARMS}
print("\nbest estimator by median power: " +
      ", ".join(f"{a}={best[a]}" for a in ARMS))

verdicts = {}


def verdict(name, text, confirmed, failed, value):
    v = "CONFIRMED" if confirmed else ("FAILED" if failed else "INCONCLUSIVE")
    verdicts[name] = dict(prediction=text, value=value, verdict=v)
    print(f"  {name}  {v:13s} {text}\n        observed: {value}")


def pw(arm, name, depth=None, ratio=None, col="power"):
    q = bm[(bm.arm == arm) & (bm.estimator == name)]
    if depth is not None:
        q = q[q.depth == depth]
    if ratio is not None:
        q = q[q.ratio == ratio]
    return float(q[col].median())


print("\n=== pre-committed predictions")
verdict("C1", "the ranking depends on the encoding",
        best["abundance"] != best["occupancy"],
        len(set(best.values())) == 1,
        f"best by power: " + ", ".join(f"{a}={best[a]}" for a in ARMS))

p_occ = pw("occupancy", "pa_single", 1000, 4)
verdict("C2", "presence/absence recovers under occupancy-encoded preference",
        p_occ >= 0.30, p_occ < 0.10,
        f"power(pa_single, occupancy arm, D=1000, r=4) = {p_occ:.3f}")

adv_ab = pw("abundance", "betabinom") - pw("abundance", "pa_single")
adv_oc = pw("occupancy", "betabinom") - pw("occupancy", "pa_single")
verdict("C3", "betabinom's dominance was an artefact of the first generator",
        (adv_ab - adv_oc) >= 0.20, adv_oc >= adv_ab,
        f"betabinom - pa_single power: abundance arm {adv_ab:+.3f}, "
        f"occupancy arm {adv_oc:+.3f}, shrinkage {adv_ab - adv_oc:+.3f}")

fdrs = {a: pw(a, "clr_ttest", 1000, col="fdr") for a in ARMS}
verdict("C4", "no evidence floor means poor FDR whatever the encoding",
        all(v > 0.20 for v in fdrs.values()), any(v <= 0.10 for v in fdrs.values()),
        "clr_ttest FDR at D=1000: " + ", ".join(f"{a} {v:.3f}" for a, v in fdrs.items()))

base_p = pw("both", "pa_single", 1000, 4)
base_f = pw("both", "pa_single", 1000, 4, col="fdr")
beats = [n for n in ESTIMATORS if n != "pa_single"
         and pw("both", n, 1000, 4) > base_p
         and pw("both", n, 1000, 4, col="fdr") <= base_f]
verdict("C5", "under the realistic arm something beats the paper's own method",
        bool(beats), not beats,
        f"both arm, D=1000, r=4: pa_single power {base_p:.3f} FDR {base_f:.3f}; "
        f"beaten by {beats or 'none'}")

out = dict(plan_sha256=plan_sha, n_cells=len(cells), n_rep=N_REP,
           arms=ARMS, depths=DEPTHS, ratios=RATIOS,
           best_by_arm=best,
           per_arm_estimator={
               f"{a}|{n}": dict(
                   n_evaluable=float(bm[(bm.arm == a) & (bm.estimator == n)].n_evaluable.median()),
                   power=pw(a, n), fdr=pw(a, n, col="fdr"), dir_acc=pw(a, n, col="dir_acc"))
               for a in ARMS for n in ESTIMATORS},
           verdicts=verdicts)
(R / "benchmark2.json").write_text(json.dumps(out, indent=1))
n_fail = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
print(f"\nwrote results/benchmark2.{{tsv,json}} — {n_fail} of {len(verdicts)} FAILED")
