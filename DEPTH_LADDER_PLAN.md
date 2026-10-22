# Controlled depth experiment inside NEON: does the threshold re-select taxa as depth changes?

**Written 2026-10-23, before any outcome of this experiment was computed.** Nothing below may be changed
because of a result. If a decision here turns out to be wrong, the fix is a new, separately labelled
analysis — not an edit to this file. The hash of this file is recorded in `PLAN_HASHES.md` together with
the commit that introduced it, and that commit precedes the commit carrying any result.

## 1. Why this experiment exists

§ 3.4 of the manuscript derives a mechanism: a minimum-evidence threshold admits taxa whose expected read
count is near one, so raising depth moves the set that qualifies *down* the abundance distribution rather
than sharpening the estimates for taxa already in it. The evidence offered for it in real data is that
FunAqua and NEON share only 69 of 134 testable orders (Jaccard 0.51).

**That evidence is confounded and we have said so.** The two datasets differ in depth, but also in
geography, water-body type and in what "sediment" means — lake and river sediment against stream benthic
biofilm. An order testable in one and not the other may be ecologically absent rather than filtered out by
depth, and the observed turnover cannot separate the two.

This experiment removes the confound by staying inside NEON. The ecosystem, the samples, the reference
database, the pairing rule and the statistical test are all held fixed; **only the rarefaction depth
changes**. If the mechanism operates, it must show itself here. If it does not show itself here, the
§ 3.4 claim about real data is not supported and the section must be reduced to a statement about
simulation.

## 2. The fixed pair subset — and why it is fixed

`exploratory_neon_depth.py` (the post hoc probe behind Figure 4) capped depth but let every pair
participate at every cap, so a pair shallower than the cap simply used its own minimum. The pair set
therefore changed between caps and pair composition was confounded with depth. That is not repeated here.

- The subset is **every qualifying NEON pair whose within-pair minimum depth is ≥ 10 000 reads**, taken
  from `results/neon_pairs.tsv` (`rarefied_depth`), which the frozen NEON validation produced. A survey of
  that column gives **102 of 149 pairs**; the exact membership is written to
  `results/depth_ladder_pairs.tsv` before any statistic is computed.
- **The same 102 pairs are used at every rung.** No pair is ever upsampled, because every rung is at or
  below the subset's minimum available depth.
- Pairs are not re-selected, re-paired or re-filtered between rungs for any reason.

## 3. The ladder

Rungs: **100, 300, 1 000, 3 000, 10 000** reads per sample. A 100-fold range, with the bottom rung below
FunAqua's own within-pair depth (median 131–179 reads per sample) and the top rung at NEON's own scale.

- Within-pair rarefaction is at OTU level by `rng.multivariate_hypergeometric`, presence then aggregated
  to order — identical to `analysis_neon_validation.py`, not reimplemented.
- At each rung, both members of a pair are drawn to exactly that rung.
- **50 draws per rung**, with the same 50 seeds reused at every rung so that draw noise is paired across
  rungs and a difference between rungs cannot be a difference in draw luck.
- An order is **testable** at a rung if it reaches **≥ 20 discordant pairs**, the same threshold used
  everywhere else in this paper. Significance, where reported, is the same two-sided exact binomial
  against 0.5 with Benjamini–Hochberg correction across testable orders.

## 4. Quantities measured at each rung

Each is the median across the 50 draws of the per-draw value.

1. `n_testable` — orders reaching the threshold.
2. `median_abundance` — median, over the testable set, of each order's estimated relative abundance *â*.
3. `median_lambda` — median, over the testable set, of *λ̂ = â · D*.
4. `median_effect` — median, over the testable set, of |preference index − 0.5|.
5. `jaccard_vs_top` — Jaccard similarity between the testable set at this rung and at 10 000.
6. `delta_index_persistent` — among orders testable at **every** rung, the median
   |index at this rung − index at 10 000|.

**How *â* is estimated, and what that does and does not license.** *â* for an order is its share of reads
in the uncapped data, averaged over the 102 pairs and both compartments. It is therefore estimated once,
at full depth, and held fixed across rungs. Consequently *λ̂* for a *fixed order* is proportional to *D* by
construction — that is arithmetic, not a finding, and no prediction below is about it. Every prediction
about *λ̂* is about the **median over the set that passes the threshold**, whose membership is free to
change. Reporting *λ̂* this way is the only way to ask whether the selected set keeps its position on the
detection curve; it is not an independent measurement of abundance.

## 5. Pre-committed predictions, with failure criteria

| # | prediction | confirmed if | **failed if** |
|---|---|---|---|
| P1 | the testable set moves down the abundance distribution | `median_abundance` falls ≥ 10-fold from rung 100 to rung 10 000 | falls < 2-fold, or rises |
| P2 | the testable set holds its position on the detection curve | `median_lambda` stays inside a 3-fold band across all five rungs | spans > 10-fold |
| P3 | apparent effect size is roughly depth-invariant | `median_effect` range across rungs ≤ 0.10 | a monotone trend larger than 0.10 |
| P4 | the testable set turns over | `jaccard_vs_top` at rung 100 is < 0.6 | > 0.8 |
| P5 | **the discriminator:** depth changes membership rather than the estimate for a given taxon | `delta_index_persistent` at rung 100 is ≤ 0.10 | > 0.20 |

**P5 is the test that can sink the section's wording.** § 3.4 currently says depth changes which taxa are
measurable "not how specialised they appear". That is a claim about taxa held fixed. If orders testable at
every rung nevertheless have substantially different indices at 100 and at 10 000, then depth changes both
things, and the "not" is wrong. It will be reported either way, and if it fails the manuscript sentence
changes rather than the plan.

**P2 and P4 can also fail in a way that matters.** If the set barely turns over (P4 fails) then the
mechanism is not operating at these depths in this dataset, and § 3.4's real-data claim must be withdrawn
and replaced by the simulation-only statement.

**Interpretation is fixed now:** confirmation of P1, P2 and P4 together establishes that the threshold
re-selects taxa with depth *in real data, with ecology held constant*. It does not establish that this is
what caused the FunAqua/NEON disagreement, which remains confounded and will continue to be reported as
such.

## 6. What this experiment cannot do

- It is one dataset, one amplicon (ITS), one habitat contrast and one threshold value (20 discordant
  pairs). The threshold sweep and the effect-size sweep that § 3.4's domain caveat calls for are separate
  work and are not in this plan.
- The 102 deepest-sequenced pairs are not a random sample of NEON pairs; deeply sequenced samples may
  differ systematically. Restricting to them is what makes the design clean, and the cost is that the
  result is about those pairs.
- *â* is estimated from the same data, as § 4 states.
- Subsampling a finished library is not the same as having sequenced it shallower: it reproduces
  depth-driven sampling but not any depth-dependent step earlier in the protocol.
- A null result here would constrain § 3.4's real-data claim, not the algebra in § 3.4's first paragraph,
  which stands or falls on its own derivation.

## 7. Deliverables

`analysis_depth_ladder.py` → `results/depth_ladder.tsv` (per rung per draw),
`results/depth_ladder.json` (the medians and each prediction's verdict),
`results/depth_ladder_pairs.tsv` (the frozen 102-pair subset).
Outcome written to `DEPTH_LADDER_RESULTS.md`, reporting every prediction's verdict including failures.
