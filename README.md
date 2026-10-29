# The threshold decides which taxa you can measure

Every metabarcoding analysis that tests taxa for habitat preference has to decide first which taxa
carry enough evidence to test at all. That rule is usually a minimum count of informative samples, it is
usually treated as a neutral piece of bookkeeping, and it is not neutral. It admits taxa whose expected
read count sits near one at whatever depth the study happens to use. Change the depth and you change
*which* taxa the survey can measure, not how specialised they look.

That is what this repository is about. The paper built from it is in `SUBMISSION.md`.

## How I got here, including the part that failed

I started somewhere else. FunAqua is a global paired survey, the same sites sampled in both water and
sediment, and nobody had asked whether the fungi in it specialise between the two. They do: at order
level, 100 orders have enough discordant pairs to test and about 49 show a significant preference, split
coherently between zoosporic lineages going to the water column and soil and dung organisms going to the
sediment. Nothing in the method knows any fungal biology, so that split is a real check on it.

Then I tried to replicate the per-order table in NEON and it did not replicate. Agreement in sign came
out at 14 of 27, 51.9 %, exact binomial p = 1.0, against a rule I had written down in advance. A coin
flip.

NEON is not a null dataset, which is what makes the failure interesting. On its own it finds 37 of 103
orders with a significant preference. So habitat partitioning is real in both datasets, measured by
different people with different primers. What does not transfer is which order prefers which habitat.

The two datasets differ about 150-fold in sequencing depth. Working out why that should matter is what
turned an empirical failure into the paper.

## The mechanism

A taxon at relative abundance *a* in a sample rarefied to depth *D* is detected with probability
1 − (1 − *a*)^*D*, which is about 1 − e^−λ with λ = *aD*. The preference index is built from discordant
pairs only, so its expectation depends on *D* even when the true abundances are fixed. A taxon whose
true abundance ratio is a constant 4× has an expected index of 0.80 at FunAqua's depth and 1.00 at
NEON's. Presence at 211 reads and presence at 31 131 reads are not noisy measurements of the same thing.

That algebra predicts deeper data should show more extreme indices. **It does not happen.** NEON's
median |index − 0.5| is 0.151 against FunAqua's 0.149, a ratio of 1.02 at p = 0.69, despite the depth
gap. The prediction failed and the reason it failed is the result: the minimum-evidence rule re-selects
the testable set at every depth, so the marginal distribution of indices stays put while the membership
of the set turns over underneath it.

## What I got wrong along the way, all of it still in here

This is the part I would want to read in someone else's repository, so it stays.

**An explanation I gave was refuted by my own next analysis.** I had an account of why the median
abundance of the testable set falls with a particular exponent. The parameter sweep in
`SWEEP_RESULTS.md` tested it and both predictions, S1 and S2, failed. The explanation is withdrawn in
two files with the original sentence quoted, rather than quietly edited out.

**A whole benchmark was discarded.** The first estimator comparison encoded preference only as
abundance given presence, which handed one estimator its own estimand and scored a competitor at power
0.000. The generator was deciding the winner. It was thrown away and replaced by the three-arm design in
`BENCHMARK2_PLAN.md`, where no estimator dominates because the ranking flips with how preference is
encoded. That is the honest answer to "which estimator should I use" and it is not a satisfying one.

**The headline count was unstable and I nearly published one draw of it.** Two scripts implementing the
same endpoint at the same threshold with the same seed disagreed, 49 significant orders against 51.
Neither was miscoded: one consumed a module-level RNG across four thresholds. Running the primary over
200 draws gives a range of 43 to 57, median 51, sd 2.5. What is stable is the set, not the count: 35
orders are significant in at least 95 % of draws. Single-draw rarefaction is not a reporting unit.

**A result I had written down evaporated.** Three orders I reported as disagreeing between Europe and
North America are significant in only 27, 136 and 185 of 200 draws. That was draw noise, not
biogeography.

**A claim as first written was too broad.** "Fungi partition less than other eukaryotes" pooled protists
with plants and animals. Against Alveolata specifically, fungi are indistinguishable. The group-level
contrast in the paper is stated at the resolution the data actually support.

`results/quarantine/` holds a third results file whose provenance I cannot establish. No script in here
regenerates it and it contradicts the deterministic primary. Kept so the discrepancy stays on record.

## Pre-registration

Eleven analysis plans were committed before the corresponding outcomes were computed. Their SHA-256
hashes live in `PLAN_HASHES.md`, outside the plan files, together with the commit that introduced each
one, and the analysis scripts assert the hash at run time before computing anything.
`verify_plan_hashes.py` checks all eleven.

The rule I held to: a frozen plan is never edited because of a result. Where a frozen decision turned
out to be wrong, the fix is a new analysis under a new label, not a revision of the old one.

## Layout

Plans and their outcomes pair up by name: `ANALYSIS_PLAN.md`, `NEON_VALIDATION_PLAN.md`,
`EUK_ROBUSTNESS_PLAN.md`, `ZOOSPORIC_PLAN.md`, `METAPR2_REPLICATION_PLAN.md`, `DEPTH_LADDER_PLAN.md`,
`SWEEP_PLAN.md`, `PREVALENCE_PLAN.md`, `BENCHMARK_PLAN.md`, `BENCHMARK2_PLAN.md`, `OCCUPANCY_PLAN.md`,
each with a matching `*_RESULTS.md`.

`SUBMISSION.md` is the paper. `SUPPLEMENTARY.md` is its supplement. `MANUSCRIPT.md` is the same text
with my working notes attached, and `SUBMISSION_NOTES.md` holds what did not belong in a paper.
`RESULTS.md` and `PAPER.md` are earlier write-ups from before the reframing, kept for the record.

Analyses are plain Python needing numpy, pandas, scipy and statsmodels. Each script is standalone and
writes into `results/`.

    analysis_primary.py          the pre-registered FunAqua primary, seed 20261021, deterministic
    analysis_secondary.py        thresholds, conservatism, regions, traits, eukaryote benchmark
    analysis_draw_stability.py   200 draws; post hoc, added after the RNG discrepancy
    analysis_neon_validation.py  the pre-registered NEON test
    analysis_depth_ladder.py     the controlled within-NEON depth experiment
    analysis_simulation.py       detection, re-selection and the index under known truth
    analysis_sweep.py            the parameter sweep that refuted my exponent explanation
    analysis_prevalence.py       local prevalence against global abundance
    analysis_benchmark2.py       the three-arm estimator comparison
    analysis_occupancy.py        the detection-explicit occupancy model
    analysis_euk_matched.py      the depth- and pair-matched eukaryote comparison
    analysis_euk_robustness.py   five robustness checks, written down before they were run
    analysis_levins.py           niche breadth as an alternative metric
    analysis_zoosporic.py        the pre-specified motile-stage prediction
    analysis_biome_kingdom.py    rivers against lakes, and the per-kingdom breakdown
    make_figures.py              the five main figures and the supplementary ones
    make_supplementary.py        tables S1 to S14
    verify_plan_hashes.py        asserts all eleven frozen plans

In the supplementary tables and the per-order results, an em dash in a `prefers` column means the taxon
showed no significant preference. It is a placeholder, not a missing value.

## Data

FunAqua's OTU tables go in `~/research/datasets/funaqua/` and the trait table in
`~/research/datasets/fungaltraits/`. The NEON tables come from the API and need a personal token in
`~/.neon_token`, read from there and sent as a header, never written into any script here. Use
DP1.20141.002 for surface water and DP1.20086.002 for benthic. Do not use the `.001` versions: they are
superseded and assign only about 5 % of reads to an order, which I found out the slow way.

Nothing larger than the metadata is committed.

## Licence

Code under the MIT licence, see `LICENSE`. Text, figures and derived tables under CC BY 4.0.
