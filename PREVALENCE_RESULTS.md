# Pooling abundance was the whole discrepancy, and the informative number is not the exact match

2026-10-23. Analysis as fixed in `PREVALENCE_PLAN.md` (sha256 `35e1f6be2c337ffd…`, full value in
`PLAN_HASHES.md`, commit `4ac0440`, frozen before `analysis_prevalence.py` existed). The script asserts
the plan still matches that hash before computing anything.

**All five pre-registered predictions confirmed, including the control.** The diagnosis `SWEEP_RESULTS.md`
left untested is correct: the sweep's failure came from pooling abundance across sites, not from the
detection account.

## The four variants

Same 102 pairs as the depth ladder, same five rungs, same ≥ 20-discordant threshold, same estimator.
Measured on real data: *b* = **0.485**, Jaccard = **0.281**, testable set **36 → 77** orders.

| variant | *b* | R² | Jaccard | testable, 100 → 10 000 |
|---|---|---|---|---|
| **M1** pooled *â* (the control) | **+0.958** | 1.000 | 0.000 | 40 → 96 |
| **M2** observed per-sample *â*ᵢⱼᶜ | **+0.484** | 0.918 | **0.300** | **36 → 78** |
| **M2b** as M2, zeros replaced by each order's smallest non-zero value | **+1.183** | 0.989 | 0.121 | 39 → 234 |
| **M3** parametric prevalence: present with probability = observed occupancy, abundance *â*/*p* where present | **+0.432** | 0.961 | 0.366 | 45 → 109 |

| # | prediction | observed | verdict |
|---|---|---|---|
| D1 | per-pair abundance reproduces the exponent | *b* = 0.484 vs 0.485, off by **0.001** | **CONFIRMED** |
| D2 | and the turnover | Jaccard 0.300 vs 0.281, off by 0.019 | **CONFIRMED** |
| D3 | control: pooled model still fails | *b* = 0.958 (sweep gave 0.96) | **CONFIRMED** |
| D4 | structure alone suffices without the real matrix | *b*(M3) = 0.432, bar was < 0.75 | **CONFIRMED** |
| D5 | absences, not graded variation, matter | *b*(M2b) − *b*(M2) = **+0.699**, bar was ≥ 0.15 | **CONFIRMED** |

## The exact match in D1 is close to a tautology, and should not be quoted as the result

*â*ᵢⱼᶜ is the read fraction of order *i* in that exact sample at full depth, and M2 simulates detection as
Bernoulli at 1 − exp(−*â*ᵢⱼᶜ·*D*). That is a Poisson approximation to what real within-pair rarefaction
does: subsample a library of that composition down to *D* reads. **M2 therefore reproduces the real
procedure almost by construction**, and agreeing with it to 0.001 in *b* is near-arithmetic rather than
evidence about biology. A reviewer would spot this, so it is said here first.

What is *not* tautological is that the same harness, fed the same orders and pairs, produces 0.958, 0.484,
1.183 and 0.432 depending only on how abundance is represented. The content of this analysis is in those
contrasts:

1. **Pooling is the entire discrepancy.** M1 and M2 differ in nothing but whether abundance varies across
 samples, and that single change moves *b* from 0.958 to 0.484, across the whole gap the sweep could
 not close. Neither distribution shape (S1) nor lognormal per-pair noise (σ = 1) came close, because
 neither creates persistent local presence.
2. **Sparse occupancy is the mechanism, and absences dominate.** The median order is present in only
 **9.8 %** of these 204 samples; 91 of 254 orders are in ≤ 5 % and only 6 in ≥ 95 %. Filling the zeros
 in (M2b) while keeping all the graded variation sends *b* to 1.183, past complete re-selection, and
 inflates the testable set from 78 to 234 orders at the top rung. So it is **where a taxon is absent**,
 not how much its abundance varies where present, that governs which taxa are measurable at low depth.
3. **Two numbers per taxon are enough.** M3 discards the observed matrix entirely and keeps only each
 order's occupancy *p*ᵢ and its pooled abundance, rescaled to *â*ᵢ/*p*ᵢ where present so that pooled
 abundance is preserved by construction. It recovers *b* = 0.432 against the measured 0.485 and
 Jaccard 0.366 against 0.281. **This is the result worth reporting**: the depth-dependence of which taxa
 a survey can measure is predictable from occupancy and mean abundance alone, with no appeal to the
 observed per-sample table.

## What this changes

`SWEEP_RESULTS.md` concluded that *D*^−0.485 was "an empirical measurement without a validated mechanism".
**That is now too pessimistic and is corrected there.** The mechanism is validated, with the amendment
that abundance must be represented per site: an independent-detection account reproduces the measured
exponent and turnover once taxa are allowed to be locally present and widely absent, and a two-parameter
occupancy model does most of the work.

The earlier failure was a modelling artefact. Pooling a sparse community into one abundance per taxon
makes every taxon look uniformly rare instead of locally common, which makes low-depth testability look
impossible and so predicts near-total re-selection with depth.

## Limits

- The D1 match is near-tautological, as stated above. Quote M1-vs-M3 instead.
- *â*ᵢⱼᶜ is a finite-sample estimate; an exact zero may be a sampling zero rather than a true absence.
 M2b is the sensitivity check and it moves the answer a great deal, so the distinction matters. A taxon
 with one read in a 69 440-read sample is "present" here, and the boundary is arbitrary at that scale.
- Occupancy is computed over the 204 samples of these 102 pairs, which are NEON's deepest; occupancy
 measured on shallower samples would itself be depth-dependent, and the circularity in that is not
 resolved here.
- Water and sediment are simulated independently within a pair, and no variant models the weak dependence
 rarefying to a fixed total induces across taxa within a sample.
- No true habitat effect is imposed (ratio 1.0): the question was which taxa become *testable*, which is
 what *b* and Jaccard measure. Effect sizes were swept separately in `SWEEP_RESULTS.md`.
- One dataset, one amplicon, one threshold, one habitat contrast.

## Files

`PREVALENCE_PLAN.md` (frozen) · `analysis_prevalence.py` → `results/prevalence.tsv`,
`results/prevalence.json`, `results/order_abundance_by_pair.tsv.gz` (the observed 204 × 254 matrix).
