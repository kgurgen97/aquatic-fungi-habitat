#!/usr/bin/env python3
"""A detection-explicit occupancy model, as OCCUPANCY_PLAN.md specifies.

Its sha256 is in PLAN_HASHES.md and is asserted before anything runs.

Occupancy models need repeated visits to separate occupancy from detection, and a
metabarcoding survey gives one observation per sample. K independent subsamples of
one library are K detection trials, which supplies them at no cost in data.

What psi means here is presence in the LIBRARY, not at the site: the model corrects
subsampling, not sequencing. The plan says so and the write-up repeats it.
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
from scipy import stats

warnings.filterwarnings("ignore", category=RuntimeWarning)

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
P = REPO_ROOT
R = P / "results"

PLAN = P / "OCCUPANCY_PLAN.md"
plan_sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
rec = re.search(r"`OCCUPANCY_PLAN\.md` \| `([0-9a-f]{64})`", (P / "PLAN_HASHES.md").read_text())
assert rec and plan_sha == rec.group(1), "OCCUPANCY_PLAN.md changed since freezing"
print(f"plan verified: {plan_sha[:16]}…")

# Borrow the generator and constants from benchmark2, which borrows them from
# benchmark: imports, pre-definition assignments and function definitions only, so
# nothing executes. Same seeds, same arms, same cells - the point is comparability.
def borrow(path, names):
    """Run a sibling script's setup and definitions, but not its analysis.

    Everything before the first `def` is kept whatever its node type, not just
    assignments: analysis_benchmark2.py does its own borrowing through a bare
    `exec(...)`, which is an Expr, and dropping it would leave its namespace empty
    and the lookups below raising KeyError. Keeping the whole prelude also runs
    that script's own frozen-plan hash check, which is wanted.
    """
    src = (P / path).read_text()
    tree = ast.parse(src)
    fd = next(i for i, n in enumerate(tree.body)
              if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)))
    keep = [n for i, n in enumerate(tree.body)
            if i < fd or isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))]
    ns = {"__name__": "borrowed", "__file__": str(P / path)}
    exec(compile(ast.Module(body=keep, type_ignores=[]), f"<{path}>", "exec"), ns)
    missing = [k for k in names if k not in ns]
    assert not missing, f"{path} did not define {missing}"
    return {k: ns[k] for k in names}


b2 = borrow("analysis_benchmark2.py", ["simulate", "odds_scale", "ARMS", "DEPTHS",
                                       "RATIOS", "N_REP", "T", "N_PAIRS", "ALPHA", "_bh"])
simulate, _bh = b2["simulate"], b2["_bh"]
ARMS, DEPTHS, RATIOS, N_REP = b2["ARMS"], b2["DEPTHS"], b2["RATIOS"], b2["N_REP"]
T, N_PAIRS, ALPHA = b2["T"], b2["N_PAIRS"], b2["ALPHA"]
K = 8                       # subsample replicates per sample, plan section 3
MIN_DETECTED = 5            # evaluable floor, plan section 3
print(f"borrowed generator: arms {ARMS}, {T} orders, {N_PAIRS} pairs, K={K}")

# E[s | s>0] = K p / (1 - (1-p)^K). Invert once on a grid, then look up - this runs
# 254 times per sample group per dataset and must not call a root finder each time.
_PG = np.linspace(1e-4, 1 - 1e-6, 20000)
_MG = K * _PG / (1.0 - (1.0 - _PG) ** K)


def p_from_truncated_mean(mbar):
    """Invert the zero-truncated binomial mean. mbar in (1, K]."""
    out = np.full_like(mbar, np.nan, dtype=float)
    ok = np.isfinite(mbar) & (mbar > 1.0 + 1e-9) & (mbar <= K)
    out[ok] = np.interp(mbar[ok], _MG, _PG)
    # mbar at or below 1 means every detected sample was detected exactly once:
    # detection is very low, so take the smallest p the grid carries.
    out[np.isfinite(mbar) & (mbar <= 1.0 + 1e-9)] = _PG[0]
    return out


def detection_histories(counts, depth, rng):
    """-> s[sample, taxon], detections out of K subsamples of size `depth`."""
    n, t = counts.shape
    s = np.zeros((n, t), dtype=np.int16)
    tot = counts.sum(1)
    for j in range(n):
        if tot[j] <= 0:
            continue
        pvec = counts[j] / tot[j]
        # With replacement, so the K trials are exchangeable given the library.
        draws = rng.multinomial(depth, pvec, size=K)
        s[j] = (draws > 0).sum(0)
    return s


def occupancy_test(sw, ss):
    """Detection-explicit occupancy model per taxon.

    psi_g = (d_g/n_g) / pi_g with pi_g = 1-(1-p_g)^K, p_g estimated per compartment
    from the zero-truncated mean. Likelihood ratio of free (psi_w, psi_s) against a
    common psi, with pi_g held at its estimate, referred to chi2 with 1 df.
    """
    out_psi_w = np.full(T, np.nan)
    out_psi_s = np.full(T, np.nan)
    pval = np.full(T, np.nan)

    res = {}
    for tag, s in (("w", sw), ("s", ss)):
        det = s > 0
        d = det.sum(0).astype(float)                 # samples with >=1 detection
        n = float(s.shape[0])
        with np.errstate(invalid="ignore", divide="ignore"):
            mbar = np.where(d > 0, s.sum(0) / np.maximum(d, 1), np.nan)
        p = p_from_truncated_mean(mbar)
        pi = 1.0 - (1.0 - p) ** K
        res[tag] = (d, n, pi)

    dw, nw, piw = res["w"]
    ds, ns, pis = res["s"]
    sel = (dw + ds >= MIN_DETECTED) & np.isfinite(piw) & np.isfinite(pis) \
        & (piw > 0) & (pis > 0)

    psi_w = np.clip(np.divide(dw / nw, piw, out=np.full(T, np.nan), where=sel), 0, 1)
    psi_s = np.clip(np.divide(ds / ns, pis, out=np.full(T, np.nan), where=sel), 0, 1)
    psi0 = np.divide(dw + ds, nw * piw + ns * pis,
                     out=np.full(T, np.nan), where=sel)
    psi0 = np.clip(psi0, 0, 1)

    def loglik(d, n, pi, psi):
        q = np.clip(psi * pi, 1e-12, 1 - 1e-12)
        return d * np.log(q) + (n - d) * np.log1p(-q)

    with np.errstate(invalid="ignore", divide="ignore"):
        ll1 = loglik(dw, nw, piw, psi_w) + loglik(ds, ns, pis, psi_s)
        ll0 = loglik(dw, nw, piw, psi0) + loglik(ds, ns, pis, psi0)
        lr = 2.0 * (ll1 - ll0)
    ok = sel & np.isfinite(lr) & (lr >= 0)
    pval[ok] = stats.chi2.sf(lr[ok], 1)
    out_psi_w[sel] = psi_w[sel]
    out_psi_s[sel] = psi_s[sel]
    return out_psi_w, out_psi_s, pval, sel


rows = []
cells = list(itertools.product(ARMS, DEPTHS, RATIOS))
print(f"\n=== {len(cells)} cells x {N_REP} replicates")
for arm, depth, ratio in cells:
    for rep in range(N_REP):
        seed = zlib.crc32(f"{arm}|{depth}|{ratio}|{rep}".encode())
        cw, cs, prefers = simulate(arm, ratio, depth, np.random.default_rng(seed))
        rng = np.random.default_rng(seed + 2)
        sw = detection_histories(cw, depth, rng)
        ss = detection_histories(cs, depth, rng)
        psi_w, psi_s, p, sel = occupancy_test(sw, ss)
        q = _bh(p)
        called = np.isfinite(q) & (q < ALPHA)
        correct = (psi_w > psi_s) == prefers
        tp = called & prefers & correct
        fp = called & (~prefers | ~correct)
        rows.append(dict(
            arm=arm, depth=depth, ratio=ratio, rep=rep, estimator="occupancy",
            n_evaluable=int(sel.sum()),
            power=float(tp.sum() / max(int((prefers & sel).sum()), 1)),
            fdr=float(fp.sum() / max(int(called.sum()), 1)),
            dir_acc=float((called & correct).sum() / max(int(called.sum()), 1))))
    print(f"  {arm:10s} D={depth:6,d} r={ratio}: "
          f"evaluable {rows[-1]['n_evaluable']:3d}  "
          f"power {np.median([r['power'] for r in rows[-N_REP:]]):.3f}  "
          f"fdr {np.median([r['fdr'] for r in rows[-N_REP:]]):.3f}")

oc = pd.DataFrame(rows)
oc.to_csv(R / "occupancy_model.tsv", sep="\t", index=False)

b2r = pd.read_csv(R / "benchmark2.tsv", sep="\t")
allr = pd.concat([b2r, oc], ignore_index=True)
ESTS = list(dict.fromkeys(list(b2r.estimator.unique()) + ["occupancy"]))

print(f"\n{'arm':11s} {'estimator':15s} {'evaluable':>10s} {'power':>7s} "
      f"{'FDR':>7s} {'dir':>7s}")
for arm in ARMS:
    for name in ESTS:
        s = allr[(allr.arm == arm) & (allr.estimator == name)]
        print(f"{arm:11s} {name:15s} {s.n_evaluable.median():10.0f} "
              f"{s.power.median():7.3f} {s.fdr.median():7.3f} {s.dir_acc.median():7.3f}")


def med(arm, name, col="power", depth=None):
    q = allr[(allr.arm == arm) & (allr.estimator == name)]
    if depth is not None:
        q = q[q.depth == depth]
    return float(q[col].median())


verdicts = {}


def verdict(name, text, confirmed, failed, value):
    v = "CONFIRMED" if confirmed else ("FAILED" if failed else "INCONCLUSIVE")
    verdicts[name] = dict(prediction=text, value=value, verdict=v)
    print(f"  {name}  {v:13s} {text}\n        observed: {value}")


print("\n=== pre-committed predictions")
pw_occ_arm = {n: med("occupancy", n) for n in ESTS}
order = sorted(ESTS, key=lambda n: -pw_occ_arm[n])
rank_of_occ = order.index("occupancy") + 1
verdict("O1", "highest power of all six in the occupancy arm",
        rank_of_occ == 1, rank_of_occ > 3,
        "power in the occupancy arm: " +
        ", ".join(f"{n} {pw_occ_arm[n]:.3f}" for n in order) +
        f"  -> occupancy ranks {rank_of_occ} of {len(ESTS)}")

p_ab, f_ab = med("abundance", "occupancy"), med("abundance", "occupancy", "fdr")
verdict("O2", "near-blind in the abundance arm, by design",
        (p_ab <= 0.15) and (f_ab <= 0.10), f_ab > 0.25,
        f"abundance arm: power {p_ab:.3f} (bar <=0.15), FDR {f_ab:.3f} (bar <=0.10)")

f_both = med("both", "occupancy", "fdr")
f_pa = med("both", "pa_single", "fdr")
verdict("O3", "controls FDR better than a single presence/absence draw",
        f_both <= 0.08, f_both > f_pa,
        f"both arm: occupancy FDR {f_both:.3f} vs pa_single {f_pa:.3f}")

ev = {a: (med(a, "occupancy", "n_evaluable"), med(a, "pa_single", "n_evaluable"))
      for a in ARMS}
n_lower = sum(1 for a in ARMS if ev[a][0] < ev[a][1])
verdict("O4", "evaluates at least as many orders as pa_single",
        n_lower == 0, n_lower >= 2,
        "evaluable, occupancy vs pa_single: " +
        ", ".join(f"{a} {ev[a][0]:.0f}/{ev[a][1]:.0f}" for a in ARMS))

out = dict(plan_sha256=plan_sha, K=K, n_cells=len(cells), n_rep=N_REP,
           per_arm={a: dict(n_evaluable=med(a, "occupancy", "n_evaluable"),
                            power=med(a, "occupancy"),
                            fdr=med(a, "occupancy", "fdr"),
                            dir_acc=med(a, "occupancy", "dir_acc")) for a in ARMS},
           occupancy_arm_ranking=order,
           verdicts=verdicts)
(R / "occupancy_model.json").write_text(json.dumps(out, indent=1))
n_fail = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
print(f"\nwrote results/occupancy_model.{{tsv,json}} — {n_fail} of {len(verdicts)} FAILED")
