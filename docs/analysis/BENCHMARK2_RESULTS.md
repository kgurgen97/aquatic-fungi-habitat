# There is no best estimator, the ranking flips with how preference is encoded

2026-10-23. Analysis as fixed in `BENCHMARK2_PLAN.md` (sha256 `9b101499c7779afa…`, full value in
`PLAN_HASHES.md`, commit `6a94646`, frozen before `analysis_benchmark2.py` existed).

**All five pre-registered predictions confirmed.** 27 cells (3 encodings × 3 depths × 3 effect sizes),
30 datasets each, 254 orders, 102 pairs, five estimators on identical simulated counts.

> **A slip in the frozen plan, not corrected there.** § 3 lists three arms × three depths × three ratios
> and then says "81 cells". The product is 27. The axes are unambiguous, so the run covers the 27 cells
> they define; the frozen text keeps its arithmetic error, as the rule requires.

## The result

Medians over each arm:

| encoding | estimator | evaluable | power | FDR | direction |
|---|---|---|---|---|---|
| **abundance** | `pa_single` | 79 | 0.000 | 0.000 | 0.000 |
| | `pa_repeated` | 76 | 0.000 | 0.000 | 0.000 |
| | `prop_wilcoxon` | 153 | 0.212 | 0.819 | 0.813 |
| | `clr_ttest` | 254 | 0.180 | 0.880 | 0.812 |
| | **`betabinom`** | 117 | **0.852** | 0.056 | 0.973 |
| **occupancy** | **`pa_single`** | 86 | **0.385** | 0.000 | 1.000 |
| | `pa_repeated` | 83 | 0.317 | 0.000 | 1.000 |
| | `prop_wilcoxon` | 154 | 0.062 | 0.615 | 0.486 |
| | `clr_ttest` | 254 | 0.129 | 0.853 | 0.836 |
| | `betabinom` | 114 | 0.383 | 0.052 | 0.981 |
| **both** | **`pa_single`** | 87 | **0.920** | 0.125 | 1.000 |
| | `pa_repeated` | 82 | 0.875 | **0.037** | 1.000 |
| | `prop_wilcoxon` | 161 | 0.376 | 0.755 | 0.742 |
| | `clr_ttest` | 254 | 0.441 | 0.813 | 0.884 |
| | `betabinom` | 126 | 0.908 | 0.050 | 0.980 |

| # | prediction | observed | verdict |
|---|---|---|---|
| C1 | the ranking depends on the encoding | best by power: abundance → `betabinom`, occupancy → `pa_single`, both → `pa_single` | **CONFIRMED** |
| C2 | presence/absence recovers under occupancy-encoded preference | `pa_single` power 0.754 at *D* = 1 000, *r* = 4 (bar ≥ 0.30) | **CONFIRMED** |
| C3 | `betabinom`'s dominance was an artefact | advantage over `pa_single` +0.852 in the abundance arm, **−0.002** in the occupancy arm; shrinkage 0.854 | **CONFIRMED** |
| C4 | no evidence floor means poor FDR | `clr_ttest` FDR 0.865 / 0.835 / 0.813 across arms at *D* = 1 000 | **CONFIRMED** |
| C5 | something beats the paper's method under the realistic arm | at *D* = 1 000, *r* = 4: `pa_single` power 0.947 FDR 0.186, beaten by `betabinom` | **CONFIRMED** |

## What this settles

**C3 vindicates the decision to throw the first benchmark away.** `betabinom` beat `pa_single` by 0.852
in power under the first generator and by **−0.002**, nothing at all, once preference acts on occupancy
instead of only on abundance. Its apparent dominance was the generator handing it its own estimand. Had
the first benchmark been published as it stood, it would have recommended an estimator on the strength of
a modelling choice.

**There is no best estimator, and that is the finding.** Which method wins is decided by how the
preference is encoded, which is a property of the biology and is not known in advance. Any benchmark
using one generator will recommend whichever estimator shares its parameterisation.

## What a practitioner should take from it

Under `both`, the only arm where a water-preferring taxon behaves as one ecologically, being both more
often present and more abundant in water, three things hold:

1. **A single rarefaction draw inflates FDR.** `pa_single` reaches power 0.920 but FDR **0.125**, against
 a nominal 0.05.
2. **The draw-stable rule fixes that, and costs little.** `pa_repeated`, significant in ≥ 95 % of draws
gives FDR **0.037** for power 0.875. Trading 0.045 of power for a threefold reduction in false
 discoveries. **This is a second, independent reason for § 3.2's recommendation**, which was argued
 there only from reproducibility; it also controls error.
3. **An overdispersed binomial on read counts is the other defensible choice**: power 0.908 at FDR 0.050,
 and it evaluates 126 orders against 87. It is not uniformly better; it loses badly in the occupancy
 arm's low-depth cells, but it is competitive under realistic encoding and sees more taxa.

And a clear negative: **applying no evidence floor is not coverage.** `clr_ttest` evaluates all 254
orders at every depth, which looked like a virtue in the first benchmark, and runs an FDR of 0.81–0.88.
`prop_wilcoxon` is poor throughout, with direction accuracy of 0.486 in the occupancy arm, no better
than a coin.

## Limits

- Three encodings are not all encodings. Preference could act on dispersion or seasonality; neither is
 modelled.
- Occupancy is odds-scaled; a logit-additive effect would give different numbers.
- `betabinom` is intercept-only, not hierarchical. **A detection-explicit occupancy model remains the
 clearest gap in this project**; it is the natural estimator for exactly the structure
 `PREVALENCE_RESULTS.md` found, and it is not tested here.
- Poisson sampling, independent taxa, independent compartments within a pair; real libraries are
 compositional with a fixed total, which is the setting `clr_ttest` is designed for and this generator
 does not produce, so its showing here is not a fair test of CLR methods in general.
- Nothing was run on the real FunAqua or NEON data. These are recommendations to test, not a re-analysis.

## Files

`BENCHMARK2_PLAN.md` (frozen) · `analysis_benchmark2.py` → `results/benchmark2.tsv`,
`results/benchmark2.json`. Supersedes `BENCHMARK_RESULTS.md`, whose generator decided its own outcome.
