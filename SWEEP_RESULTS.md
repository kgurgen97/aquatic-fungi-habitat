# The sweep: one claim is universal, and my explanation for the exponent was wrong

2026-10-23. Analysis as fixed in `SWEEP_PLAN.md` (sha256 `2af9a5a9e4002fc0…`, full value in
`PLAN_HASHES.md`, commit `7a0d398`, frozen before `analysis_sweep.py` existed). The script asserts the
plan still matches that hash before computing anything.

240 cells — 4 abundance distributions × effect sizes 1–16 × thresholds 10/20/40 × 102 and 239 pairs ×
site heterogeneity σ = 0 and 1 — each at five depth rungs with 30 draws on shared seeds.

**Two of five pre-registered predictions FAILED, and they are the two that carried my explanation of the
depth-ladder result. One prediction was confirmed so completely that it is now the paper's firmest
statement.**

## S5 CONFIRMED, without exception: apparent effect size does not depend on depth

In **239 of 239** interpretable cells (R² ≥ 0.5) the median |index − 0.5| moves by ≤ 0.10 across a
100-fold depth change. The largest range anywhere in the grid is **0.058**; the median is **0.0109**. The
real NEON ladder gave **0.024**, sitting inside that distribution.

This holds across every effect size from 1× to 16×, all four abundance distributions including a
heavy-tailed Pareto, all three thresholds, both pair counts, and with or without site heterogeneity. It is
the claim the manuscript relies on — *changing depth leaves the apparent strength of what you measure
almost untouched* — and it is now supported over the whole parameter space rather than at one point in it.

## S1 FAILED: the abundance distribution does not set the exponent

Exponent *b* in median *â* ∝ *D*^−*b*, at threshold 20, 102 pairs, σ = 0 (*b* = 1 is complete
re-selection):

| distribution | r = 1 | r = 2 | r = 4 | r = 8 | r = 16 |
|---|---|---|---|---|---|
| log-uniform | +1.000 | +1.001 | +1.002 | +0.996 | +0.997 |
| lognormal | +0.978 | +1.008 | +1.024 | +1.060 | +1.015 |
| power law | +0.937 | +0.961 | +0.914 | +0.838 | +0.574 |
| **empirical** | +0.961 | +0.989 | +0.963 | +0.931 | +0.852 |

The prediction required *b*(log-uniform) − *b*(empirical) ≥ 0.25. Observed: **+0.012** at *r* = 2 and
**+0.039** at *r* = 4. The two distributions behave the same.

## S2 FAILED: the detection model does not reproduce the measured exponent

This was the point of the analysis. Fed the real abundance vector, the model predicts *b* ≈ **0.96** —
essentially complete re-selection — against the **0.485** the NEON ladder measured. No effect size in
1–16 brings it within 0.25; the nearest, *r* = 16, gives 0.852, off by **0.367**. Simulated turnover is
also near-total (Jaccard ≈ 0.000) where the real ladder gave **0.281**.

So **reality re-selects about half as much as independent-Poisson detection says it should.** Real taxa
persist in the testable set across depths more than the model allows.

### What I claimed in `DEPTH_LADDER_RESULTS.md`, and why it is withdrawn

That write-up attributed the half-re-selection to the shape of the abundance distribution: log-uniform is
flat in log space so the set slides down it at constant size, while a steeper real distribution runs out
of rare taxa to recruit. **S1 refutes the exponent half of that and S2 refutes it as an explanation of
0.485.** The sentence has been corrected there and in the manuscript.

One part of it survives and is worth keeping separate. The **growth** of the testable set *is* explained
by distribution shape: log-uniform grows 0.96–1.01-fold (constant size, pure exchange of members) while
the empirical vector grows **2.15-fold at r = 4**, which is exactly the 2.15-fold growth the real ladder
showed. So the shape accounts for the set getting bigger; it does not account for the set re-selecting
only half way.

### What the failure points at

Site heterogeneity was the plan's first named candidate and is in the grid. It is **not** the answer:
σ = 1 moves *b* from 0.961 to 0.970 at *r* = 1 and from 0.852 to 0.772 at *r* = 16 — too little, in a
grid where the gap to close is 0.48.

Scanning all 240 cells, only **3** reproduce both *b* and Jaccard within 0.15 of the measured values, and
all three use the empirical vector with **239 pairs** — more than the 102 the real experiment used — with
a large effect and a low threshold (closest: empirical, *r* = 16, threshold 10, 239 pairs, σ = 0, giving
*b* = 0.565 and Jaccard 0.315). Raising the pair count at a fixed threshold lowers the effective bar, so
rarer taxa qualify at low depth and the set is less depth-determined.

Read together with the σ result, that says the real low-depth testable sets contain **rarer taxa than
independent detection from a pooled abundance allows**. The most likely reason is structure the model does
not have: a taxon locally abundant at a subset of sites is detected there at every depth and keeps
supplying discordant pairs, while its pooled *â* is low. The σ axis does not create that structure,
because it draws a fresh multiplier per pair per taxon, which averages out over 102 pairs instead of
producing persistent site specialists.

**That is a diagnosis, not a demonstration.** Testing it needs abundance estimated per site rather than
pooled, and a generative model with site-level taxon prevalence. Not in this plan and not claimed here.

## S3 INCONCLUSIVE: large effects barely escape the selection

The reviewer's counterexample — a taxon reliably detected in one compartment and essentially never in the
other passes the threshold at λ ≫ 1 — predicts turnover should fall as effect size rises. Direction is
right, magnitude is small: Jaccard goes 0.000 → 0.062 from *r* = 1 to *r* = 16 at σ = 0, and 0.000 →
0.128 at σ = 1. The confirmation bar was a rise of ≥ 0.15.

The exponent does respond: *b* falls from 0.961 at *r* = 1 to 0.852 at *r* = 16 for the empirical vector,
and the power-law community drops to 0.574. So large effects partially escape, but at these parameters the
λ ≈ 1 selection is more robust to effect size than § 3.4's caveat implies. The caveat is kept — it is
correct in principle, and the power-law column shows it can matter — but it is no longer presented as the
main limit on the mechanism.

## S4 INCONCLUSIVE: the threshold mostly sets the level

*b* across thresholds 10/20/40 at *r* = 4, empirical, 102 pairs: **+0.866, +0.963, +0.993**, a spread of
0.127, inside the ≤ 0.15 confirmation band. But the prediction also required Jaccard to vary *more* than
*b* did, and Jaccard's spread was only 0.024. So the threshold moves the exponent about as much as it
moves the turnover, and the "level not exponent" framing is not supported as stated.

## Verdict table

| # | prediction | verdict | observed |
|---|---|---|---|
| S1 | abundance shape sets the exponent | **FAILED** | difference +0.012 to +0.039, bar was ≥ 0.25 |
| S2 | empirical vector reproduces *b* = 0.485 | **FAILED** | *b* ≈ 0.96; nearest 0.852, off by 0.367 |
| S3 | large effects escape the selection | INCONCLUSIVE | Jaccard +0.062, bar was ≥ 0.15 |
| S4 | threshold sets level, not exponent | INCONCLUSIVE | *b* spread 0.127, Jaccard spread 0.024 |
| S5 | depth-invariance of apparent effect is general | **CONFIRMED** | 239/239 cells, max range 0.058 |

## What the paper can now say, and what it cannot

**Can say, over the whole swept space and in real data:** a minimum-evidence threshold makes the set of
measurable taxa depend on depth, while the apparent strength of the preferences you do measure is
essentially depth-invariant. S5 plus the ladder's 0.024 support this without qualification.

**Cannot say:** that the independent-Poisson account is quantitatively right for this dataset. It
over-predicts re-selection by about a factor of two and over-predicts turnover badly. The mechanism is
*sufficient* to produce depth-dependent membership — that was § 3.4's simulation — but it is not
calibrated, and the ladder's *D*^−0.485 remains an empirical measurement without a validated model behind
it.

**Must not say:** that the abundance distribution's shape explains the measured exponent. It does not.

## Files

`SWEEP_PLAN.md` (frozen) · `analysis_sweep.py` → `results/sweep.tsv` (240 cells),
`results/sweep.json` (verdicts), `results/sweep_empirical_abundance.tsv` (the 254 *â* values,
exported from the depth-ladder run, which reproduced `depth_ladder.json` byte-identically).
