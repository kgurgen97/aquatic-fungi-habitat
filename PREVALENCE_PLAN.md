# Does per-site abundance explain the exponent that pooled abundance could not?

**Written 2026-10-23, before `analysis_prevalence.py` existed and before any outcome of it was computed.**
Nothing below may be changed because of a result. If a decision here turns out to be wrong, the fix is a
new, separately labelled analysis — not an edit to this file. The hash of this file is recorded in
`PLAN_HASHES.md` with the commit that introduced it, and that commit precedes any result commit.

## 1. The specific failure this tests

`SWEEP_RESULTS.md` recorded two failed predictions. Fed the **pooled** empirical abundance vector — one
*â* per order, averaged over 102 pairs and both compartments — an independent-Poisson detection model
predicts a re-selection exponent *b* ≈ 0.96 and near-total turnover (Jaccard ≈ 0.000). The NEON depth
ladder measured *b* = 0.485 and Jaccard 0.281. Neither the shape of the abundance distribution (S1) nor
per-pair lognormal site heterogeneity (σ = 1, which moves *b* by ≤ 0.08) closes a gap of about 0.48.

`SWEEP_RESULTS.md` offered a diagnosis and explicitly labelled it as untested: real low-depth testable
sets contain rarer taxa than independent detection from a pooled abundance allows, because **a taxon
locally abundant at a subset of sites is detected there at every depth while its pooled abundance stays
low**. Pooling destroys exactly that structure. The σ axis could not supply it, because it redraws a
multiplier per pair per taxon, which averages out over 102 pairs instead of producing persistent local
specialists.

This analysis tests that diagnosis. It is a test of the *measurement model*, not of the biology.

## 2. The four model variants, fixed

All use the same 102 pairs, the same five rungs (100, 300, 1 000, 3 000, 10 000), 30 draws on shared
seeds, the same ≥ 20-discordant-pair threshold, and the same estimator as § 3.4a — select by simulated
discordant count at each rung, take the median abundance of the selected set, fit *b* in
median *a* ∝ *D*^−*b* by least squares in log–log across the five rungs. Detection in a given pair and
compartment is Bernoulli at 1 − exp(−*a* · *D*).

| variant | abundance used for order *i* in pair *j*, compartment *c* |
|---|---|
| **M1 pooled** | *â*ᵢ, one value per order, as in the sweep. The known-failing control. |
| **M2 per-pair** | *â*ᵢⱼᶜ, the observed relative abundance of order *i* in that exact sample at full depth. Zeros stay zero. |
| **M2b per-pair, no zeros** | as M2, but a zero is replaced by that order's smallest non-zero *â*ᵢⱼᶜ across the 102 pairs. |
| **M3 parametric prevalence** | order *i* is present in a pair with probability *p*ᵢ = its observed occupancy across the 102 pairs; where present its abundance is *â*ᵢ / *p*ᵢ, so the expected pooled abundance is preserved. Drawn fresh per draw. |

**Why M2b and M3 are both here.** M2 uses the real matrix and so carries both graded variation *and* true
absences. M2b removes the absences while keeping the variation. M3 keeps the absences in expectation,
discards the real matrix, and preserves pooled abundance by construction. Comparing the three separates
"abundance varies across sites" from "taxa are absent from many sites" from "the exact observed matrix".

**The abundance matrix is observed, not fitted.** *â*ᵢⱼᶜ is a read fraction from the full-depth data for
that sample. It is shared with the real measurement, exactly as the pooled vector was in the sweep. What
is tested is whether detection simulated from that matrix selects the taxa real rarefaction selected — a
test of the detection model, not of the abundance estimate.

## 3. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| D1 | per-pair abundance reproduces the measured exponent | *b*(M2) within **0.485 ± 0.15** | \|*b*(M2) − 0.485\| > 0.25 |
| D2 | and the measured turnover | Jaccard(M2) within **0.281 ± 0.10** | off by > 0.20 |
| D3 | **control:** the pooled model still fails | *b*(M1) > 0.80 | *b*(M1) < 0.65 |
| D4 | the structure alone suffices, without the real matrix | *b*(M3) < 0.75 | *b*(M3) > 0.85 |
| D5 | absences, not graded variation, are what matter | *b*(M2b) − *b*(M2) ≥ 0.15 | ≤ 0.05 |

**D3 is a harness control and must pass.** M1 here is the sweep's own configuration recomputed through
new code. If M1 comes out near 0.485, this script disagrees with the sweep and nothing else in it may be
believed; that would be reported as a defect to be found, not as a result.

**D1 is the diagnosis.** If it fails, pooling was not the problem, and the ladder's *D*^−0.485 keeps
standing as an empirical law with no validated mechanism — in which case the honest conclusion is that the
independent-Poisson detection account is wrong in some way this analysis has not identified, and § 3.4b
should say so without a replacement story.

**D5 can fail while D1 succeeds**, which would mean graded variation alone is enough and the zeros are
incidental. Either way it constrains what the manuscript may claim about *why*.

## 4. What this cannot do

- Confirming D1 shows that detection simulated from the observed per-sample abundances reproduces one
  exponent and one turnover statistic. It does not prove that independent Bernoulli detection is the true
  data-generating process, only that it suffices once abundance is not pooled.
- *â*ᵢⱼᶜ is itself a finite-sample estimate from the full-depth library, so rare taxa have noisy values
  and exact zeros may be sampling zeros rather than true absences. M2b is the sensitivity check for that
  and is the reason it is in the plan.
- Compartments within a pair are simulated independently, as in every earlier analysis here. Real water
  and sediment samples from one site may share taxa for reasons this does not model.
- Rarefying a real library to a fixed total makes presence calls weakly dependent across taxa within a
  sample; none of the four variants models that.
- One dataset, one amplicon, one threshold, one habitat contrast.

## 5. Deliverables

`analysis_prevalence.py` → `results/prevalence.tsv` (per variant per rung),
`results/prevalence.json` (fitted *b*, Jaccard and each prediction's verdict),
`results/order_abundance_by_pair.tsv.gz` (the observed *â*ᵢⱼᶜ matrix).
Outcome written to `PREVALENCE_RESULTS.md`, reporting every verdict including failures.
