# Where the threshold mechanism operates, and whether it predicts the exponent we measured

**Written 2026-10-23, before any sweep outcome was computed and before `analysis_sweep.py` existed.**
Nothing below may be changed because of a result. If a decision here turns out to be wrong, the fix is a
new, separately labelled analysis — not an edit to this file. The hash of this file is recorded in
`PLAN_HASHES.md` with the commit that introduced it, and that commit precedes any result commit.

## 1. Why, and what is new since the last two analyses

§ 3.4 derives that a minimum-evidence threshold admits taxa whose expected read count is near one, so
depth re-selects which taxa are measurable. Two things are now known and they do not agree:

- **Simulation** (§ 3.4, a community of 400 abundances spread log-uniformly over six orders of magnitude):
  the testable set's median λ holds at 0.22–0.49 across a 1000-fold depth change. Complete re-selection.
- **Real data** (§ 3.4a, the NEON depth ladder, ecology held constant): median *â* ∝ *D*^−0.485 and median
  *λ̂* ∝ *D*^+0.515, R² 0.92. Re-selection about half complete, and the testable set *grows* 2.15-fold
  rather than exchanging members.

So the idealised form is a property of the simulated community. The obvious candidate explanation is the
**shape of the abundance distribution**: log-uniform over six decades is flat in log space, so as depth
rises the set slides down it and stays the same size, whereas a steeper or truncated distribution runs out
of rare taxa to recruit and the set grows instead.

This sweep asks two questions, one descriptive and one predictive:

1. **Where does the mechanism operate?** The reviewer's point: the λ ≈ 1 statement was derived at one
   effect size and one abundance distribution, and a taxon reliably detected in one compartment and
   essentially never in the other accumulates discordant pairs at λ ≫ 1, so it can pass the threshold
   without sitting mid-curve. Map the domain instead of conceding it.
2. **Does the mechanism predict the exponent we measured?** If the detection model is right, feeding it
   the *empirical* abundance vector should reproduce *b* ≈ 0.485 without being told the answer.

## 2. The generative model, fixed

For a taxon with sediment abundance *a* and water:sediment ratio *r* ≥ 1, at depth *D*:

    lambda_water = a * r * D        lambda_sediment = a * D
    P(detect) = 1 - exp(-lambda)

Per pair, water and sediment detection are drawn independently as Bernoulli at those probabilities. With
site heterogeneity enabled, *a* is multiplied per pair by a lognormal draw with σ = `sigma_site`
(σ = 0 reproduces the homogeneous case). Presence/absence only, as in the real analysis.

Discordant pairs = water-only + sediment-only. A taxon is **testable** at a depth if discordant pairs
≥ `threshold`. The preference index is water-only / discordant, and effect is |index − 0.5|.

This mirrors § 3.4a's measurement exactly: select by simulated discordant counts at each depth, then take
the median *a* of the selected set, and fit *b* in median *a* ∝ *D*^−*b* by least squares in log–log over
the same five rungs. **The same estimator is applied to simulated and real data**, so the comparison in
question 2 is a comparison of one number computed one way.

## 3. The grid, fixed now

| axis | values |
|---|---|
| effect size *r* | 1, 2, 4, 8, 16 |
| abundance distribution | `loguniform` (1e−6…1e−1, the § 3.4 community), `lognormal` (σ = 2 in log10), `powerlaw` (Pareto, α = 1.5, heavy-tailed), `empirical` (the 254 nonzero *â* from `results/depth_ladder.json`'s run, written to `results/sweep_empirical_abundance.tsv`) |
| threshold | 10, 20, 40 discordant pairs |
| pairs | 102 (the ladder), 239 (FunAqua's primary) |
| site heterogeneity σ | 0, 1 |
| depth rungs | 100, 300, 1 000, 3 000, 10 000 — the ladder's rungs |

30 draws per cell, seeds shared across depths within a cell so that depth differences are not draw luck.
Taxa per simulated community: 254, matching the empirical vector's length, so cell sizes are comparable.

## 4. Outcomes per cell

1. `b` — the exponent in median *a* ∝ *D*^−*b*. *b* = 1 is complete re-selection, *b* = 0 is none.
2. `r2` — R² of that log–log fit. A cell with R² < 0.5 has no power law and its *b* is not interpreted.
3. `lambda_span` — max/min of median λ across rungs.
4. `jaccard` — Jaccard between the testable sets at 100 and at 10 000.
5. `effect_range` — max − min of median |index − 0.5| across rungs.
6. `growth` — n_testable at 10 000 divided by n_testable at 100.

## 5. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| S1 | the exponent is set by the abundance distribution's shape | `b` for `loguniform` exceeds `b` for `empirical` by ≥ 0.25, at *r* = 2–4, threshold 20, 102 pairs, σ = 0 | the difference is < 0.10, or reversed |
| S2 | **the predictive test:** the empirical abundance vector reproduces the measured exponent | `b` for `empirical` lies within **0.485 ± 0.15** for at least one *r* in 2–8 at threshold 20, 102 pairs | no *r* in 1–16 brings `b` within 0.25 of 0.485 |
| S3 | large effects escape the λ ≈ 1 selection | `jaccard` rises monotonically with *r* over 1 → 16 for `empirical`, threshold 20, and is ≥ 0.15 higher at *r* = 16 than at *r* = 1 | no increase, or a decrease |
| S4 | the threshold sets the level, not the exponent | across thresholds 10/20/40, `b` varies by ≤ 0.15 while `jaccard` varies by more | `b` varies by > 0.30 |
| S5 | depth-invariance of apparent effect is general | `effect_range` ≤ 0.10 in ≥ 80 % of cells with R² ≥ 0.5 | ≤ 50 % of such cells |

**S2 is the point of this analysis and it can fail cleanly.** What is shared with the real measurement is
the abundance vector; what is tested is whether the independent-Poisson detection model, given that
vector, selects the same taxa the real rarefaction selected. Failure would mean the detection model is
wrong — site heterogeneity, non-independence between compartments, or compositional effects from
rarefying a fixed total — and § 3.4a's power law would then stand as an empirical description with no
mechanism behind it. That is a publishable outcome and will be reported as readily as success.

**S3 is the reviewer's counterexample, made testable.** If it fails, the "λ near 1" framing survives
at large effect sizes too and the domain caveat now in § 3.4 is too cautious and must be relaxed.

## 6. What this cannot do

- It is a simulation. Confirming S2 shows the detection model suffices to reproduce one exponent on one
  dataset; it does not prove that model is what operates in NEON.
- The empirical abundance vector is NEON's, pooled over compartments and pairs, and is itself estimated
  from finite data at full depth.
- Presence/absence only. Abundance-based estimators are a separate question, deliberately out of scope
  here; the comparison of estimators is a different analysis and is not pre-registered by this plan.
- Independent taxa throughout. Real communities have correlated abundances, and rarefying to a fixed total
  makes presence calls weakly dependent across taxa within a sample; neither is modelled.
- *r* is one number per taxon per cell. Real effect sizes vary across taxa; a mixture is not swept.

## 7. Deliverables

`analysis_sweep.py` → `results/sweep.tsv` (one row per cell),
`results/sweep.json` (each prediction's verdict), `results/sweep_empirical_abundance.tsv`.
Outcome written to `SWEEP_RESULTS.md`, reporting every verdict including failures.
