#!/usr/bin/env python3
"""Does per-site abundance explain the exponent pooled abundance could not?

As PREVALENCE_PLAN.md specifies; the plan's sha256 is recorded in PLAN_HASHES.md and
asserted here before anything is computed.

The sweep found that a pooled abundance vector makes independent-Poisson detection
predict b ~ 0.96 against the ladder's measured 0.485. The untested diagnosis was
that pooling destroys local structure. Four variants separate the parts.
"""
import ast
import hashlib
import json
import pathlib
import re
import zlib

import numpy as np
import pandas as pd

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
P = REPO_ROOT
R = P / "results"

PLAN = P / "PREVALENCE_PLAN.md"
plan_sha = hashlib.sha256(PLAN.read_bytes()).hexdigest()
rec = re.search(r"`PREVALENCE_PLAN\.md` \| `([0-9a-f]{64})`", (P / "PLAN_HASHES.md").read_text())
assert rec and plan_sha == rec.group(1), "PREVALENCE_PLAN.md changed since freezing"
print(f"plan verified: {plan_sha[:16]}…")

RUNGS = [100, 300, 1_000, 3_000, 10_000]
MIN_DISCORDANT = 20
N_DRAWS = 30
MEASURED_B, MEASURED_J = 0.485, 0.281
RATIO = 1.0   # no true habitat effect is imposed: the question is which taxa are
              # TESTABLE at each depth, and that is what b and Jaccard measure.

# ------------------------------------------------- borrow the frozen NEON loader
# Same device as analysis_depth_ladder.py: take imports, constants and function
# definitions only. Keeping every top-level Assign would run the validation itself.
src = (P / "analysis_neon_validation.py").read_text()
tree = ast.parse(src)
first_def = next(i for i, n in enumerate(tree.body)
                 if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef)))
keep = [n for i, n in enumerate(tree.body)
        if isinstance(n, (ast.Import, ast.ImportFrom, ast.FunctionDef, ast.AsyncFunctionDef))
        or (isinstance(n, ast.Assign) and i < first_def)]
ns: dict = {"__name__": "neon_io_borrowed"}
exec(compile(ast.Module(body=keep, type_ignores=[]), "<neon_io>", "exec"), ns)
load, build_pairs = ns["load"], ns["build_pairs"]

print("\n=== loading NEON ITS tables")
water = load("water", ns["MIN_CONF"])
benthic = load("benthic", ns["MIN_CONF"])
pairs_all = build_pairs(water, benthic, ns["WINDOW_DAYS"])

# The exact 102 pairs the ladder used, read from its own output so the two
# analyses cannot drift apart.
frozen = pd.read_csv(R / "depth_ladder_pairs.tsv", sep="\t")
want = set(zip(frozen.site, frozen.water_sample, frozen.benthic_sample))
pairs = [p for p in pairs_all if (p[0], p[2], p[3]) in want]
assert len(pairs) == len(frozen) == 102, f"{len(pairs)} pairs, expected 102"
print(f"pairs: {len(pairs)} (the depth ladder's own subset)")

# Same 254 orders the sweep used, in the same order.
emp = pd.read_csv(R / "sweep_empirical_abundance.tsv", sep="\t")
ORDERS = emp.order.tolist()
A_POOLED = emp.a_hat.to_numpy()
oidx = {o: i for i, o in enumerate(ORDERS)}
T = len(ORDERS)
print(f"orders: {T}  pooled a-hat median {np.median(A_POOLED):.2e}")


def sample_abundance(sample):
    """Relative abundance per order in one sample, over the fixed order list."""
    v = np.zeros(T)
    for (_, order), c in sample["otu"].items():
        j = oidx.get(order)
        if j is not None:
            v[j] += c
    tot = v.sum()
    return v / tot if tot > 0 else v


A_w = np.array([sample_abundance(water[p[2]]) for p in pairs])     # (102, T)
A_s = np.array([sample_abundance(benthic[p[3]]) for p in pairs])
print(f"abundance matrices: {A_w.shape}, nonzero "
      f"{100*(A_w > 0).mean():.1f}% water / {100*(A_s > 0).mean():.1f}% benthic")

pd.DataFrame(np.vstack([A_w, A_s]),
             columns=ORDERS,
             index=[f"w:{p[2]}" for p in pairs] + [f"b:{p[3]}" for p in pairs]
             ).to_csv(R / "order_abundance_by_pair.tsv.gz", sep="\t",
                      compression="gzip")
print("wrote results/order_abundance_by_pair.tsv.gz")

# Occupancy across the 204 samples of these 102 pairs. Detection is simulated per
# sample, so sample-level occupancy is the analogue of the plan's "across the 102
# pairs"; stated here and in the write-up so the choice is visible.
OCC = ((A_w > 0).sum(0) + (A_s > 0).sum(0)) / (2 * len(pairs))
OCC = np.clip(OCC, 1e-6, 1.0)
print(f"occupancy: median {np.median(OCC):.3f}, "
      f"{int((OCC >= 0.95).sum())} orders in >=95% of samples, "
      f"{int((OCC <= 0.05).sum())} in <=5%")

# A-hat where present, preserving pooled abundance: a/p * p = a.
A_IF_PRESENT = A_POOLED / OCC

# M2b: zeros replaced by that order's smallest non-zero observed value.
def no_zeros(A):
    out = A.copy()
    for j in range(T):
        col = A[:, j]
        nz = col[col > 0]
        if nz.size:
            out[col == 0, j] = nz.min()
    return out


A_w_nz, A_s_nz = no_zeros(A_w), no_zeros(A_s)

VARIANTS = ["M1_pooled", "M2_per_pair", "M2b_per_pair_no_zeros", "M3_prevalence"]


def draw(variant, depth, rng):
    """-> (testable mask, preference index). Returns per-order selection."""
    n = len(pairs)
    if variant == "M1_pooled":
        aw = np.broadcast_to(A_POOLED * RATIO, (n, T))
        as_ = np.broadcast_to(A_POOLED, (n, T))
    elif variant == "M2_per_pair":
        aw, as_ = A_w * RATIO, A_s
    elif variant == "M2b_per_pair_no_zeros":
        aw, as_ = A_w_nz * RATIO, A_s_nz
    elif variant == "M3_prevalence":
        pw = rng.random((n, T)) < OCC
        ps = rng.random((n, T)) < OCC
        aw = np.where(pw, A_IF_PRESENT * RATIO, 0.0)
        as_ = np.where(ps, A_IF_PRESENT, 0.0)
    else:
        raise ValueError(variant)
    w = rng.random((n, T)) < (1.0 - np.exp(-aw * depth))
    s = rng.random((n, T)) < (1.0 - np.exp(-as_ * depth))
    wo = (w & ~s).sum(0)
    so = (~w & s).sum(0)
    disc = wo + so
    sel = disc >= MIN_DISCORDANT
    idx = np.full(T, np.nan)
    np.divide(wo, disc, out=idx, where=sel)
    return sel, idx


def fit_b(depths, vals):
    d, v = np.asarray(depths, float), np.asarray(vals, float)
    ok = np.isfinite(v) & (v > 0)
    if ok.sum() < 3:
        return np.nan, np.nan
    x, y = np.log10(d[ok]), np.log10(v[ok])
    m, c = np.polyfit(x, y, 1)
    ss = np.sum((y - y.mean()) ** 2)
    return -m, (1 - np.sum((y - (m * x + c)) ** 2) / ss if ss > 0 else np.nan)


rows, summary = [], {}
print(f"\n=== {len(VARIANTS)} variants x {len(RUNGS)} rungs x {N_DRAWS} draws")
for variant in VARIANTS:
    seeds = [zlib.crc32(f"{variant}|{i}".encode()) for i in range(N_DRAWS)]
    per, sel_at = {}, {}
    for depth in RUNGS:
        ma, me, nt, sels = [], [], [], []
        for sd in seeds:
            sel, idx = draw(variant, depth, np.random.default_rng(sd))
            sels.append(sel)
            nt.append(int(sel.sum()))
            if sel.any():
                ma.append(float(np.median(A_POOLED[sel])))
                me.append(float(np.median(np.abs(idx[sel] - 0.5))))
        per[depth] = dict(n_testable=float(np.median(nt)),
                          median_a=float(np.median(ma)) if ma else np.nan,
                          median_effect=float(np.median(me)) if me else np.nan)
        sel_at[depth] = np.array(sels).mean(0) >= 0.5
        rows.append(dict(variant=variant, rung=depth, **per[depth]))

    b, r2 = fit_b(RUNGS, [per[d]["median_a"] for d in RUNGS])
    lo, hi = sel_at[RUNGS[0]], sel_at[RUNGS[-1]]
    union = int((lo | hi).sum())
    jac = int((lo & hi).sum()) / union if union else np.nan
    summary[variant] = dict(b=b, r2=r2, jaccard=jac,
                            n_lo=per[RUNGS[0]]["n_testable"],
                            n_hi=per[RUNGS[-1]]["n_testable"])
    print(f"  {variant:24s} b {b:+.3f} (R2 {r2:.3f})  jaccard {jac:.3f}  "
          f"testable {per[RUNGS[0]]['n_testable']:.0f} -> {per[RUNGS[-1]]['n_testable']:.0f}")

pd.DataFrame(rows).to_csv(R / "prevalence.tsv", sep="\t", index=False)
print(f"\nmeasured on real data: b {MEASURED_B}, jaccard {MEASURED_J}")

verdicts = {}


def verdict(name, text, confirmed, failed, value):
    v = "CONFIRMED" if confirmed else ("FAILED" if failed else "INCONCLUSIVE")
    verdicts[name] = dict(prediction=text, value=value, verdict=v)
    print(f"  {name}  {v:13s} {text}\n        observed: {value}")


b2, j2 = summary["M2_per_pair"]["b"], summary["M2_per_pair"]["jaccard"]
b1 = summary["M1_pooled"]["b"]
b3 = summary["M3_prevalence"]["b"]
b2b = summary["M2b_per_pair_no_zeros"]["b"]

print("\n=== pre-committed predictions")
verdict("D1", "per-pair abundance reproduces the measured exponent",
        abs(b2 - MEASURED_B) <= 0.15, abs(b2 - MEASURED_B) > 0.25,
        f"b(M2) = {b2:+.3f} vs measured {MEASURED_B}, off by {abs(b2-MEASURED_B):.3f}")
verdict("D2", "and the measured turnover",
        abs(j2 - MEASURED_J) <= 0.10, abs(j2 - MEASURED_J) > 0.20,
        f"jaccard(M2) = {j2:.3f} vs measured {MEASURED_J}, off by {abs(j2-MEASURED_J):.3f}")
verdict("D3", "CONTROL: the pooled model still fails",
        b1 > 0.80, b1 < 0.65, f"b(M1) = {b1:+.3f} (sweep gave ~0.96)")
verdict("D4", "structure alone suffices, without the real matrix",
        b3 < 0.75, b3 > 0.85, f"b(M3) = {b3:+.3f}")
verdict("D5", "absences, not graded variation, are what matter",
        (b2b - b2) >= 0.15, (b2b - b2) <= 0.05,
        f"b(M2b) - b(M2) = {b2b:+.3f} - {b2:+.3f} = {b2b-b2:+.3f}")

out = dict(plan_sha256=plan_sha, measured_b=MEASURED_B, measured_jaccard=MEASURED_J,
           n_pairs=len(pairs), n_orders=T, n_draws=N_DRAWS, rungs=RUNGS,
           summary={k: {kk: (None if not np.isfinite(vv) else float(vv))
                        for kk, vv in v.items()} for k, v in summary.items()},
           verdicts=verdicts)
(R / "prevalence.json").write_text(json.dumps(out, indent=1))
n_fail = sum(1 for v in verdicts.values() if v["verdict"] == "FAILED")
print(f"\nwrote results/prevalence.{{tsv,json}} — {n_fail} of {len(verdicts)} FAILED")
if verdicts["D3"]["verdict"] != "CONFIRMED":
    print("!! the control did not pass: this script disagrees with the sweep.")
    print("   Nothing else here may be believed until that is explained.")
