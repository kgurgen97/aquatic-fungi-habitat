# Which estimator should be used, given that the threshold decides what is measurable?

**Written 2026-10-23, before `analysis_benchmark.py` existed and before any outcome of it was computed.**
Nothing below may be changed because of a result. If a decision here turns out to be wrong, the fix is a
new, separately labelled analysis — not an edit to this file. The hash is recorded in `PLAN_HASHES.md`
with the commit that introduced it, and that commit precedes any result commit.

## 1. Why, and why only now

§§ 3.4a–3.4c establish a problem: a minimum-evidence threshold on presence/absence makes *which* taxa a
survey can assess depend on its depth, while leaving the apparent strength of what it does measure almost
unchanged. That is a warning. It becomes useful only if some estimator avoids it, so this benchmark asks
which one to use.

It is deliberately last. A benchmark is only as good as the data generator behind it, and until
`PREVALENCE_RESULTS.md` there was no validated generator: pooled abundance reproduced neither the measured
re-selection exponent (0.958 against 0.485) nor the turnover. The occupancy-plus-mean-abundance model does
(*b* = 0.432, Jaccard 0.366), from two numbers per taxon, so it is the generator used here.

## 2. The generator, fixed

Per simulated dataset: **254 orders, 102 pairs**, matching the ladder.

- Order *i* has occupancy *p*ᵢ and pooled abundance *â*ᵢ taken from
  `results/sweep_empirical_abundance.tsv` and `results/order_abundance_by_pair.tsv.gz`, the same values
  the prevalence analysis used.
- In each pair and compartment the order is present with probability *p*ᵢ; where present its expected
  relative abundance is *â*ᵢ/*p*ᵢ, so pooled abundance is preserved.
- **True preference.** A fixed 30 % of orders, chosen by a seeded draw once per dataset, are given a true
  water:sediment abundance ratio *r*; the remaining 70 % have *r* = 1 and are true nulls. The true index
  for a preferring taxon is *r*/(1 + *r*); for a null it is exactly 0.5.
- Reads are drawn as Poisson with mean (relative abundance × depth) per order per sample, then the sample
  is treated as that count vector. Depth *D* is the per-sample read budget.

Grid: *D* ∈ {100, 1 000, 10 000}; *r* ∈ {2, 4, 8}; **40 simulated datasets per cell**.

## 3. The five estimators

All are given the same simulated count table and must decide, per order, whether there is a water/sediment
preference and in which direction. Significance is Benjamini–Hochberg *q* < 0.05 across the orders that
estimator can evaluate.

| key | estimator | what it does |
|---|---|---|
| `pa_single` | presence/absence, one rarefaction draw | the paper's own method: rarefy both members to the pair minimum, require ≥ 20 discordant pairs, exact binomial on water-only vs sediment-only |
| `pa_repeated` | presence/absence, 25 draws | as above, but an order counts as significant only if it is significant in ≥ 95 % of draws, and its index is the median across draws — the draw-stable rule of § 3.2 |
| `prop_wilcoxon` | paired Wilcoxon on proportions | per order, the within-sample relative abundance in water against sediment across pairs, Wilcoxon signed-rank; no thresholding on discordance |
| `clr_ttest` | paired *t*-test on CLR | centred log-ratio per sample with a 0.5 pseudocount, paired *t*-test across pairs |
| `betabinom` | overdispersed binomial | per order, water reads out of water + sediment reads per pair, intercept-only binomial GLM with a Pearson overdispersion scale, testing the intercept against 0.5 |

`pa_single` and `pa_repeated` are what the manuscript currently does and recommends. The other three use
counts and never threshold on discordant pairs, so they are the candidates for escaping the problem.

## 4. Metrics, per cell per estimator

1. `n_evaluable` — orders the estimator could test at all. **This is the metric the rest of the paper
   implies and that the usual benchmark omits.**
2. `power` — of the truly preferring orders, the fraction called significant in the correct direction.
3. `fdr` — of the orders called significant, the fraction that are true nulls or wrong-direction.
4. `dir_acc` — of significant calls, the fraction whose direction is correct.
5. `bias` — mean (estimated index − true index) over evaluable preferring orders, on the index scale;
   estimators not returning an index on that scale contribute `nan` and are excluded from this metric
   only, which is stated rather than worked around.
6. `rmse` — root mean square of the same difference.
7. `depth_stability` — Jaccard of the evaluable set between *D* = 100 and *D* = 10 000 at matched *r*.

## 5. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| B1 | count-based estimators evaluate far more orders than presence/absence | median `n_evaluable` for the best count-based estimator ≥ 1.5× that of `pa_single` at *D* = 100 | < 1.1× |
| B2 | and their evaluable set is far less depth-dependent | best count-based `depth_stability` ≥ 0.80 while `pa_single` is ≤ 0.50 | the count-based value is below `pa_single`'s |
| B3 | repeated rarefaction controls FDR better than a single draw | `fdr`(`pa_repeated`) ≤ `fdr`(`pa_single`) − 0.02 at every depth | higher at any depth |
| B4 | at least one estimator beats `pa_single` on power at equal or lower FDR | such an estimator exists at *D* = 1 000, *r* = 4 | none does at any cell |
| B5 | no estimator achieves nominal FDR everywhere | some cell has `fdr` > 0.05 for every estimator | every estimator holds ≤ 0.05 in every cell |

**B5 is deliberately a prediction of failure.** Compositional data, overdispersion and a sparse occupancy
structure are hard, and claiming any of these methods is calibrated would be the easy, wrong conclusion.
If B5 fails — if something does hold nominal FDR throughout — that is a stronger and more useful result
than expected and will be reported as such.

**B1 and B2 are the ones that matter for the paper's recommendation.** If they fail, presence/absence
thresholding is not costing coverage after all, and §§ 3.4a–3.4c remain a caution with no remedy to offer.

## 6. What this cannot do

- The generator is validated on one statistic of one dataset (the re-selection exponent and turnover of
  NEON fungal orders). Agreement there does not make it right in every respect, and an estimator could
  win here for reasons specific to it.
- Poisson sampling, independent taxa, and independent compartments within a pair — the same
  simplifications as every simulation in this paper. Real libraries are compositional with a fixed total,
  which `clr_ttest` is designed for and the generator does not produce.
- `betabinom` here is an intercept-only overdispersed binomial, not a full hierarchical model; a
  hierarchical occupancy model with explicit detection is out of scope and its absence is a real gap.
- One pairing design, one threshold value for the presence/absence methods, one occupancy structure.
- Nothing here is run on the real FunAqua or NEON data. This benchmark is about estimator behaviour under
  known truth, and any recommendation it supports is a recommendation to test, not a re-analysis.

## 7. Deliverables

`analysis_benchmark.py` → `results/benchmark.tsv` (per cell per estimator per replicate),
`results/benchmark.json` (medians and each prediction's verdict).
Outcome written to `BENCHMARK_RESULTS.md`, reporting every verdict including failures.
