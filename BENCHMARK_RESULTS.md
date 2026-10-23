# The benchmark as frozen cannot compare these estimators — my generator decided the winner

2026-10-23. Analysis as fixed in `BENCHMARK_PLAN.md` (sha256 `ada8cb3a0403dcce…`, full value in
`PLAN_HASHES.md`, commit `828f8da`, frozen before `analysis_benchmark.py` existed).

**The design is flawed, the flaw is in the plan I wrote, and the numbers below should not be used to
choose an estimator.** They are reported because the plan says outcomes are reported, and because the
flaw is itself worth recording.

## What came out

9 cells (depths 100/1 000/10 000 × ratios 2/4/8), 40 replicates, 254 orders, 102 pairs. Medians:

| estimator | evaluable | power | FDR | direction | bias | RMSE | depth stability |
|---|---|---|---|---|---|---|---|
| `pa_single` | 79 | **0.000** | 0.000 | 0.000 | −0.218 | 0.242 | 0.345 |
| `pa_repeated` | 75 | **0.000** | 0.000 | 0.000 | −0.226 | 0.245 | 0.340 |
| `prop_wilcoxon` | 152 | 0.209 | 0.826 | 0.808 | — | — | 0.444 |
| `clr_ttest` | 254 | 0.184 | 0.883 | 0.813 | — | — | 1.000 |
| `betabinom` | 115 | **0.850** | **0.053** | 0.976 | −0.000 | 0.058 | 0.304 |

| # | prediction | verdict | observed |
|---|---|---|---|
| B1 | count-based estimators evaluate far more orders | CONFIRMED | `clr_ttest` 254 vs `pa_single` 44 at depth 100 = **5.77×** |
| B2 | and their evaluable set is less depth-dependent | CONFIRMED | `clr_ttest` stability 1.000 vs `pa_single` 0.345 |
| B3 | repeated rarefaction controls FDR better | INCONCLUSIVE | −0.092 at depth 100, 0.000 at 1 000 and 10 000 |
| B4 | some estimator beats `pa_single` on power at equal or lower FDR | **FAILED** | `pa_single` power was 0.000, so nothing could beat it "at equal or lower FDR" |
| B5 | no estimator achieves nominal FDR everywhere | INCONCLUSIVE | none holds ≤ 0.05 in every cell, but no single cell defeats all five |

## Why `pa_single` scoring zero is a red flag, not a finding

The manuscript's own primary analysis is `pa_single`, and it finds 49 significant orders in FunAqua. An
estimator that detects nothing at all under simulation is not a weak estimator; it is a sign the
simulation is wrong. It is.

**`BENCHMARK_PLAN.md` § 2 gives both compartments the same occupancy *p*ᵢ and makes preference an
abundance ratio applied only where the taxon is present.** Presence/absence estimators read occupancy.
Probing the generator directly, for a taxon whose true index is 0.800:

| depth | observed index, preferring taxa | observed index, null taxa | realised occupancy, water vs sediment |
|---|---|---|---|
| 100 | 0.636 | 0.381 | 0.434 vs 0.310 |
| 1 000 | 0.576 | 0.424 | 0.407 vs 0.366 |
| 10 000 | **0.500** | 0.485 | 0.378 vs 0.368 |

The signal does not merely weaken with depth — it **disappears**. Once detection saturates, two
compartments with identical occupancy are indistinguishable, so the index collapses to 0.5. That is the
opposite of how presence/absence behaves on real data, where deeper sequencing makes more orders testable
(§ 3.4a). The null taxa are also off-centre at low depth (0.381, not 0.5), which inflates FDR for every
method that uses this index.

And the mirror image: `betabinom` estimates the water share of reads, which is **exactly** the quantity
the generator varies. Its bias of −0.000 and RMSE of 0.058 are not evidence of a good estimator; they are
evidence that it was handed its own estimand.

So the generator pre-decided the comparison. One family was given the truth in the form it measures, the
other in a form it cannot see. B1 and B2 "confirmed" are hollow for the same reason — `clr_ttest`
evaluates all 254 orders and has perfect depth stability because it tests every order regardless of
evidence, which is a property of applying no evidence floor rather than a virtue.

## What is not wrong

The code implements the plan. The failure is in § 2 of the plan, written by me, which specified a
preference acting only on abundance-given-presence. Ecologically that is not what "prefers water" means:
a water-preferring taxon is both **present at more water sites** and **more abundant there**.

Under this project's rule, a frozen plan is not edited because of a result. The fix is a new,
separately labelled analysis with a corrected generator, pre-registered afresh: `BENCHMARK2_PLAN.md`.

## What the corrected version must do differently

1. **Preference must act on occupancy as well as abundance.** Otherwise presence/absence methods are
   being asked to detect something invisible to them.
2. **Encode preference three ways and compare across them** — abundance-only (this generator, kept as a
   known-biased arm), occupancy-only, and both — because the honest question is not "which estimator
   wins" but "which estimator wins under which kind of preference". A single generator always flatters
   whichever estimator shares its parameterisation, which is the lesson of this attempt.
3. **Make direction-based metrics primary.** With preference acting on two quantities there is no single
   scalar estimand that all five estimators target, so bias on an index scale is only meaningful for
   named estimators against a named estimand.
4. **Keep an evidence floor comparison honest.** `clr_ttest` testing all 254 orders is not coverage in
   the sense § 3.4a cares about; it is an absence of filtering, and it shows up as an FDR of 0.883.

## Files

`BENCHMARK_PLAN.md` (frozen, flawed) · `analysis_benchmark.py` → `results/benchmark.tsv`,
`results/benchmark.json`. Superseded by `BENCHMARK2_PLAN.md` and `BENCHMARK2_RESULTS.md`.
