# Depth ladder inside NEON: the consequence holds, the idealised form does not

2026-10-23. Analysis as fixed in `DEPTH_LADDER_PLAN.md`
(sha256 `b53dd38de0be59e304aa0844545d54ed…`, full value in `PLAN_HASHES.md`, commit `464e481`, frozen
before the analysis script existed). `analysis_depth_ladder.py` asserts the plan still matches that hash
before computing anything.

**Four of five pre-registered predictions confirmed; the fifth came out inconclusive and against the
idealised form of the mechanism.** That fifth result is the most informative thing here.

## Design, as frozen

102 of the 149 qualifying NEON pairs, namely every pair whose within-pair minimum depth reaches 10 000
reads. **The same 102 pairs at every rung**, available depth min 10 641, median 69 440, max 425 195, so
no pair is ever upsampled. Rungs 100 / 300 / 1 000 / 3 000 / 10 000 reads per sample, a 100-fold range
whose bottom sits below FunAqua's own within-pair depth. 50 draws per rung on the same 50 seeds, so draw
noise is paired across rungs. Loading, within-pair rarefaction and the ≥ 20-discordant-pair threshold are
taken from `analysis_neon_validation.py` unchanged.

This is what the FunAqua-versus-NEON comparison cannot be. Ecology, samples, reference database, pairing
rule and test are all fixed; only depth moves.

## What happens along the ladder

| rung | testable orders | median *â* | median *λ̂* | median \|index − 0.5\| | Jaccard vs 10 000 | median \|Δindex\| |
|---|---|---|---|---|---|---|
| 100 | 36 | 1.11 × 10⁻² | 1.109 | 0.129 | **0.281** | 0.094 |
| 300 | 51 | 3.52 × 10⁻³ | 1.057 | 0.145 | 0.488 | 0.077 |
| 1 000 | 66 | 2.10 × 10⁻³ | 2.104 | 0.133 | 0.671 | 0.062 |
| 3 000 | 72 | 1.51 × 10⁻³ | 4.523 | 0.122 | 0.783 | 0.044 |
| 10 000 | 77 | 1.03 × 10⁻³ | 10.296 | 0.125 | 1.000 | 0.000 |

## The verdicts, as pre-committed

| # | prediction | observed | verdict |
|---|---|---|---|
| P1 | the testable set moves down the abundance distribution | median *â* falls **10.8-fold** (bar: ≥ 10-fold) | **CONFIRMED**, with little margin |
| P2 | the set holds its position on the detection curve | median *λ̂* spans **9.74-fold**, 1.06 → 10.30 (confirm ≤ 3-fold, fail > 10-fold) | **INCONCLUSIVE**, and against |
| P3 | apparent effect size is roughly depth-invariant | range **0.024** (0.122–0.145); bar was ≤ 0.10 | **CONFIRMED**, comfortably |
| P4 | the testable set turns over | Jaccard(100, 10 000) = **0.281**; bar was < 0.6 | **CONFIRMED**, comfortably |
| P5 | depth changes membership rather than a given taxon's estimate | median \|Δindex\| = **0.094** over 25 persistent orders; bar was ≤ 0.10 | **CONFIRMED**, with little margin |

P2 landed at 9.74-fold against a failure threshold of 10-fold. It is reported as inconclusive because that
is what the frozen rule says, but nobody should read 9.74 as support: the prediction was a 3-fold band and
the observation is three times wider than that. Had the ladder been one rung longer it would very likely
have failed outright.

## What P2 actually shows, post hoc, and labelled as such

The following log–log fits are **not** in the frozen plan. They describe a prediction that came out
inconclusive and are post hoc.

Median *â* of the testable set falls as **D^−0.485** (R² 0.92), so median *λ̂* rises as **D^+0.515**
(R² 0.93). Complete re-selection, the simulation's behaviour, would be *â* ∝ D^−1 and *λ̂* flat. So in
real data **the set re-selects about half the way.**

Part of the reason is visible in the first column: the testable set **grows 2.15-fold**, from 36 orders to
77, as depth rises. It is not a pure exchange of members. Deeper sequencing keeps the orders that were
marginal at 100 reads, which by 10 000 reads sit far above λ = 1, and *adds* newly marginal ones. The
median of a growing, upward-drifting set rises.

> **Correction, 2026-10-23, after `SWEEP_PLAN.md` was executed.** This section originally continued: "In
> the simulation of § 3.4 the community was 400 abundances spread log-uniformly over six orders of
> magnitude, so the testable count stayed roughly stable and membership was exchanged rather than
> extended; the real fungal abundance distribution is steeper, and the exponent changes with it."
>
> **The sweep refutes that explanation.** Fed the empirical abundance vector, the same detection model
> gives *b* ≈ 0.96, statistically indistinguishable from the log-uniform community's 1.00, the shapes
> behave alike (S1 FAILED), and neither reproduces the measured 0.485 (S2 FAILED, nearest 0.852, off by
> 0.367). Site heterogeneity does not close the gap either.
>
> What survives is narrower and is kept above: distribution shape *does* explain the **growth**.
> Log-uniform grows 0.96–1.01-fold while the empirical vector grows 2.15-fold at *r* = 4, the same
> 2.15-fold seen here. Shape accounts for the set getting bigger, not for it re-selecting only half way.
>
> The exponent *D*^−0.485 therefore stands as an **empirical measurement with no validated model behind
> it**. `SWEEP_RESULTS.md` records what the failure points at: real low-depth testable sets contain rarer
> taxa than independent detection from a pooled abundance allows, which is what locally abundant taxa
> would do. That is a diagnosis awaiting a test, not a result.

The turnover itself is the reviewer's point about the mechanism's domain, now measured rather than
conceded.

## What this does and does not establish

**Establishes, with ecology held constant:**

- Depth strongly changes *which* taxa a survey can measure. The testable sets at 100 and 10 000 reads
 share 28 % of their union.
- Depth barely changes *how specialised the measured taxa appear*. Median |index − 0.5| moves by 0.024
 across a 100-fold depth change, a sixth of the spread the same quantity shows between FunAqua and NEON.
- For taxa measurable at every depth, the estimate itself is stable: median |Δindex| 0.094 between the
 extreme rungs.

Those three together are the consequence § 3.4 claims, and they now rest on real data with the ecological
confound removed, not on simulation plus a confounded cross-dataset Jaccard.

**Does not establish:**

- That the selected set holds a fixed position on the detection curve. It does not; *λ̂* rises as the
 square root of depth. § 3.4's "median λ holds at 0.22–0.49" is a property of the simulated community and
 must be presented as such.
- That depth caused the FunAqua/NEON disagreement. That comparison remains confounded and is still
 reported as confounded. This experiment shows the mechanism operates in real data; it does not attribute
 a particular disagreement to it.
- Anything about other thresholds, other amplicons, other habitat contrasts or other abundance
 distributions. One dataset, ITS, water-versus-benthic, threshold fixed at 20 discordant pairs.

**Weak spots, stated plainly:**

- P1 (10.8-fold against a 10-fold bar) and P5 (0.094 against a 0.10 bar) both cleared their thresholds
 narrowly. Neither should be quoted as a wide margin.
- P5 rests on only **25 orders** testable at every rung.
- The 102 deepest-sequenced pairs are not a random sample of NEON pairs.
- *â* is estimated once at full depth from these same data, so *λ̂* for a fixed order is proportional to
 depth by construction. Every statement above is about the median over the *set that passes the
 threshold*, whose membership changes. The plan said this before the run.

## Files

`DEPTH_LADDER_PLAN.md` (frozen) · `analysis_depth_ladder.py` →
`results/depth_ladder.tsv` (per rung per draw), `results/depth_ladder.json` (medians and verdicts),
`results/depth_ladder_pairs.tsv` (the frozen 102-pair subset).
