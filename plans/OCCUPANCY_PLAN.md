# A detection-explicit occupancy model: the estimator this project kept saying was missing

**Written 2026-10-23, before `analysis_occupancy.py` existed and before any outcome of it was computed.**
Nothing below may be changed because of a result. The hash is recorded in `PLAN_HASHES.md` with the
commit that introduced it, and that commit precedes any result commit.

## 1. Why, and why it was possible to omit until now

`BENCHMARK2_RESULTS.md`, `PREVALENCE_RESULTS.md` and § 3.4d all name the same gap: a model with an
explicit detection process is the natural estimator for the structure § 3.4c found — taxa that are
locally present and widely absent, detected or not depending on depth — and none of the five estimators
benchmarked has one.

Occupancy models need **repeated visits** to separate occupancy ψ from detection *p*, and a metabarcoding
survey gives one observation per sample. That is why this was left out.

**Rarefaction supplies the repeated visits.** Drawing *K* independent subsamples of depth *D* from one
library produces *K* detection trials for that sample: a taxon present in the library is detected in some
fraction of them, a taxon absent in none. That is exactly the detection history the model needs, and it
costs no new data.

## 2. What ψ means here, stated before any result

The occupancy this recovers is **presence in the sequenced library**, not ecological presence at the
site. The library is itself a sample of the community at whatever depth it was sequenced, so ψ here
corrects for the subsampling detection process and not for field detection. Any claim beyond that would
be unsupported, and the write-up must not make one.

## 3. The model

Per order, per compartment group *g* ∈ {water, sediment}, over the *n_g* = 102 samples of that group:

- *K* = 8 subsamples of depth *D* are drawn per sample, with replacement from the sample's observed
  composition, so the trials are exchangeable given the library. Detection is ≥ 1 read in a subsample.
- *s_j* ∈ {0…K} is the detection count for sample *j*.
- Detection given presence, *p_g*, is estimated from the samples with *s_j* > 0 by matching the mean of
  the zero-truncated binomial: E[*s* | *s* > 0] = *Kp*/(1 − (1 − *p*)^*K*).
- π_g = 1 − (1 − *p_g*)^*K* is the probability a present taxon is detected at least once.
- *d_g*, the number of samples with *s_j* > 0, is Binomial(*n_g*, ψ_g π_g), giving
  ψ̂_g = (*d_g*/*n_g*)/π_g, clipped at 1.

**Why *p* must be estimated per compartment.** With a shared *p* the correction cancels and the test
reduces to comparing raw detection frequencies, which adds nothing. It is precisely because a taxon more
abundant in water is *also detected more readily there* that *p_w* ≠ *p_s*, and separating that from a
genuine occupancy difference is the whole point of the model.

Test: likelihood ratio of free (ψ_w, ψ_s) against a common ψ, with π_g held at its estimate; under the
null ψ̂ = (*d_w* + *d_s*)/(*n_w*π_w + *n_s*π_s). The statistic is referred to χ²₁, then BH across orders.
Direction is ψ̂_w > ψ̂_s. An order is evaluable if *d_w* + *d_s* ≥ 5 and both π_g are estimable.

## 4. Grid and comparison

Identical to `BENCHMARK2_PLAN.md` so the numbers are directly comparable: arms {abundance, occupancy,
both} × depths {100, 1 000, 10 000} × ratios {2, 4, 8}, 30 datasets per cell, 27 cells, the same
generator and the same seeds. Metrics as there: `n_evaluable`, `power`, `fdr`, `dir_acc`.

## 5. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| O1 | in the `occupancy` arm it has the highest power of all six estimators | its median power exceeds every estimator in `BENCHMARK2_RESULTS.md` for that arm | it is not in the top three |
| O2 | **by design it should be near-blind in the `abundance` arm**, where true occupancy does not differ | its power there is ≤ 0.15 **and** its FDR there is ≤ 0.10 | FDR there > 0.25 |
| O3 | it controls FDR better than a single presence/absence draw under the realistic arm | `fdr` in the `both` arm ≤ 0.08, against `pa_single`'s 0.125 | > 0.125 |
| O4 | it evaluates at least as many orders as `pa_single` | median `n_evaluable` ≥ that of `pa_single` in every arm | lower in two or more arms |

**O2 is the unusual one and is deliberately a prediction of low power.** In the `abundance` arm the
taxon's occupancy genuinely does not differ between compartments — only its abundance where present does.
An occupancy model *should* decline to call those taxa. Under this benchmark's definition of truth that
counts as missed power, and it is the right behaviour for the estimand, not a defect. Scoring it as a
failure would be scoring the model for answering the question it was asked. If instead it fires in that
arm, its detection correction is not working and O2 fails on FDR.

**O1 is the claim that justifies the extra machinery.** If the model is not best in the arm built for it,
there is no case for recommending it over the simpler draw-stable rule, and the manuscript will say so.

## 6. What this cannot do

- ψ is occupancy in the library, as § 2 states, not in the environment.
- Subsamples drawn with replacement from one library are exchangeable but not independent ecological
  visits; they cannot reveal a taxon the library missed entirely. The model therefore corrects
  subsampling, not sequencing.
- *p* is constant across samples within a compartment, which is wrong: *p* depends on a taxon's abundance
  in each particular sample. A per-sample detection covariate is the obvious extension and is not done.
- Single-season, two-group, no covariates, no random effects.
- The same generator limits as `BENCHMARK2_PLAN.md` § 7, which this inherits unchanged.
- Nothing is run on real FunAqua or NEON data.

## 7. Deliverables

`analysis_occupancy.py` → `results/occupancy_model.tsv`, `results/occupancy_model.json`.
Outcome in `OCCUPANCY_RESULTS.md`, reporting every verdict including failures, with a combined
six-estimator table against `BENCHMARK2_RESULTS.md`.
