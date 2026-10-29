# Which estimator wins depends on how preference is encoded — testing that directly

**Written 2026-10-23, before `analysis_benchmark2.py` existed and before any outcome of it was computed.**
Nothing below may be changed because of a result. The hash is recorded in `PLAN_HASHES.md` with the commit
that introduced it, and that commit precedes any result commit.

## 1. Why there is a second plan

`BENCHMARK_RESULTS.md` records that the first benchmark could not compare these estimators. Its generator
(`BENCHMARK_PLAN.md` § 2) gave both compartments identical occupancy and made preference an abundance
ratio applied only where a taxon was present. Presence/absence estimators read occupancy, so they were
asked to detect something invisible to them: for a taxon whose true index was 0.800, the observed index
fell from 0.636 at depth 100 to exactly 0.500 at depth 10 000. Meanwhile `betabinom` estimates the water
share of reads, which is precisely what that generator varied, and scored a bias of −0.000.

The generator decided the winner. The frozen plan is not edited; this is the separately labelled
replacement, and the correction is to stop asking "which estimator is best" and ask the question that has
an answer.

## 2. The question

**Does the ranking of estimators depend on how a habitat preference is encoded?** If it does — and the
first attempt is reason to think so — then no benchmark with one generator can recommend an estimator,
and the practical advice must be conditional on what kind of preference is expected.

## 3. The generator, fixed

Base structure as validated in `PREVALENCE_RESULTS.md`: 254 orders, 102 pairs, each order with occupancy
*p*ᵢ and pooled abundance *â*ᵢ from `results/sweep_empirical_abundance.tsv` and
`results/order_abundance_by_pair.tsv.gz`; where present, expected relative abundance *â*ᵢ/*p*ᵢ.

A fixed 30 % of orders, drawn once per dataset by seed, prefer water with strength *r*. **Three encodings,
each a separate arm:**

| arm | occupancy in water | abundance in water where present |
|---|---|---|
| `abundance` | *p*ᵢ, unchanged | *â*ᵢ/*p*ᵢ × *r* |
| `occupancy` | odds-scaled: *p*ᵢ*r* / (1 + *p*ᵢ(*r* − 1)) | *â*ᵢ/*p*ᵢ, unchanged |
| `both` | odds-scaled as above | *â*ᵢ/*p*ᵢ × *r* |

Sediment keeps *p*ᵢ and *â*ᵢ/*p*ᵢ in every arm. Reads are Poisson with mean (relative abundance × depth).
`abundance` is deliberately the first benchmark's generator, carried over as a known-biased control so the
three arms are comparable on one axis.

Grid: arm ∈ {abundance, occupancy, both}; *D* ∈ {100, 1 000, 10 000}; *r* ∈ {2, 4, 8}; **30 datasets per
cell**. 81 cells.

## 4. Estimators

The same five as `BENCHMARK_PLAN.md` § 3, unchanged and re-used from `analysis_benchmark.py`:
`pa_single`, `pa_repeated`, `prop_wilcoxon`, `clr_ttest`, `betabinom`. Significance is BH *q* < 0.05
across the orders each estimator can evaluate.

## 5. Metrics

Direction-based metrics are primary, because with preference acting on two quantities there is no single
scalar estimand all five target. Per cell per estimator:

1. `n_evaluable` — orders the estimator can test at all.
2. `power` — truly preferring orders called significant in the correct direction, over truly preferring
   orders the estimator could evaluate.
3. `fdr` — significant calls that are null or wrong-direction.
4. `dir_acc` — significant calls with the correct direction.
5. `rank` — the estimator's rank by `power` within its cell, 1 = best.

Bias and RMSE are **not** reported. The first attempt showed they reward whichever estimator shares the
generator's parameterisation, which is the artefact this plan exists to avoid.

## 6. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| C1 | the ranking depends on the encoding | the estimator with the best median `power` differs between the `abundance` and `occupancy` arms | the same estimator is best in all three arms |
| C2 | presence/absence recovers under occupancy-encoded preference | `power`(`pa_single`) in the `occupancy` arm ≥ 0.30 at *D* = 1 000, *r* = 4 | < 0.10, i.e. still near-blind |
| C3 | `betabinom`'s first-attempt dominance was an artefact of that generator | its `power` advantage over `pa_single` in the `occupancy` arm is smaller than in the `abundance` arm by ≥ 0.20 | the advantage is as large or larger |
| C4 | no evidence floor means poor FDR, whatever the encoding | `fdr`(`clr_ttest`) > 0.20 in all three arms at *D* = 1 000 | ≤ 0.10 in any arm |
| C5 | under the most realistic arm, something beats the paper's own method | in the `both` arm at *D* = 1 000, *r* = 4, some estimator has higher `power` than `pa_single` at `fdr` no higher | none does |

**C1 is the point.** If it fails — one estimator best everywhere — that is a clean recommendation and a
better outcome than expected, and will be reported as such.

**C5 decides what the manuscript may advise.** If it fails, the paper's own presence/absence method is
not improved on under a realistic encoding and the recommendation stays "report the draw-stable set and
state the evidence floor", with no change of estimator proposed.

## 7. What this cannot do

- Three encodings are not all encodings; real preference may act on dispersion, seasonality or
  co-occurrence, none of which is modelled.
- Poisson sampling, independent taxa, independent compartments within a pair — as throughout this paper.
  Real libraries are compositional with a fixed total.
- `betabinom` remains intercept-only, not hierarchical; a detection-explicit occupancy model is still
  absent and is the clearest remaining gap in this project.
- The odds-scaling of occupancy is one choice among several; a logit-additive effect would differ.
- Nothing here touches the real FunAqua or NEON data. Any recommendation this supports is a
  recommendation to test on real data, not a re-analysis of it.

## 8. Deliverables

`analysis_benchmark2.py` → `results/benchmark2.tsv`, `results/benchmark2.json`.
Outcome written to `BENCHMARK2_RESULTS.md`, reporting every verdict including failures.
