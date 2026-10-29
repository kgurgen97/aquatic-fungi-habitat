# Stress-testing the fungi-vs-other-eukaryotes result before it becomes a headline

**Written 2026-10-21, before any of the checks below were run.** The primary result already
exists (`analysis_euk_matched.py`, `results/matched_summary.json`) and is **not** being re-chosen here;
this fixes in advance which checks it must survive and what would count as failing them.

## The result under test

On 204 pairs that qualify for both groups, rarefied to a common per-pair depth (median 211 reads):

| | orders tested | significant (q < 0.05) | share | median &#124;preference − 0.5&#124; |
|---|---|---|---|---|
| fungi | 98 | 42 | 43 % | 0.149 |
| other eukaryotes | 122 | 77 | 63 % | 0.246 |

Fisher OR 0.44, p = 0.0029; Mann–Whitney on the effect magnitudes p = 0.0009.
Reading: fungi are comparatively **generalist** across the water–sediment boundary.

## Two things already settled, recorded so they are not re-litigated

**The two tables are one assay.** FunAqua's own README: eukaryotic long-read metabarcoding with the
**ITS** region targeted; the pipeline's variables are `fungiotu` and `eukotu` — "non-fungal eukaryotes
OTU table". So both groups come from the same primers, the same reaction and the same DNA extract on the
same samples. That is the strongest available design for this contrast, and it is why the comparison is
worth making at all.

**The obvious saturation artefact runs the wrong way to explain the result.** Within fungi, orders with
more OTUs do show weaker apparent preference (Spearman −0.246, p = 0.0144). But non-fungal eukaryote
orders carry *more* OTUs per order (median 256 vs 142) while showing *stronger* preference (0.246 vs
0.149). If aggregation-driven saturation were producing the gap, the direction would be reversed.
Adjusting for it should widen the gap, not close it.

## Primary metric is changed, and why — decided before the checks

The count of significant orders is power-dependent: it answers "how many orders cleared a threshold",
which depends on how many discordant pairs each order has. The median effect magnitude does not. So for
any write-up the primary becomes **median &#124;preference − 0.5&#124;**, with the count reported as
secondary. This is a change of emphasis between two already-computed numbers that agree, not a swap to a
better-looking one.

## Checks, fixed now

1. **Draw stability.** 200 independent rarefaction draws. Report the distribution of both metrics and,
   crucially, **how often the fungi < eukaryote ordering reverses.** `analysis_secondary.py` and
   `analysis_primary.py` already disagreed on the fungal primary because of a single draw, so a result
   that flips across draws cannot be a headline.
   *Fails if* the ordering reverses in more than 5 % of draws.
2. **Matched on power.** Compare the two groups' effect magnitudes among orders **matched on the number
   of discordant pairs**, which is the direct determinant of power, by stratifying into quantile bins of
   discordant count and testing within bins (Cochran–Mantel–Haenszel style for the count metric, and a
   stratified Mann–Whitney / rank test for the effect metric).
   *Fails if* the difference disappears once discordant count is held fixed.
3. **Per kingdom, not just "other eukaryotes".** Fungi against each of Alveolata, Straminipila,
   Metazoa, Viridiplantae, Rhizaria separately.
   *Fails as a general claim if* fungi are only lower than one kingdom; then the honest statement names
   that kingdom instead.
4. **Other ranks.** The same comparison at class and family level.
   *Fails if* the result is specific to the order level with no trace at neighbouring ranks.
5. **Adjust for OTUs per order.** Partial Spearman of group membership against effect magnitude given
   OTUs per order, plus the comparison restricted to overlapping OTU-count ranges.
   *Fails if* the gap closes.

A check that fails is reported as failing. If checks 1–3 fail the claim does not go in the title.

## Limitation that stays regardless of the outcome

**ITS is a fungal barcode.** Non-fungal eukaryotes here are off-target amplification: their recovery
efficiency is taxon-dependent and their ITS reference coverage is far thinner than UNITE's is for fungi.
The design controls sample, day, assay and depth, but it cannot make the two groups equally well served
by the primers or the reference database. Any write-up states this in the abstract's limitations, not
only in a methods footnote. The claim is about what one ITS assay sees in paired water and sediment, and
the correct generalisation of it is a hypothesis for a multi-marker test, not a settled fact.

Published niche-breadth work also points the other way — fungi *narrower* than protists in freshwater
sediment (Levins' B) and narrowest of all groups in alpine soils. Those measure occupancy across samples
and gradients, not paired compartment preference, so they are not a direct contradiction, but the tension
is real and must be discussed rather than omitted.

## Deliverable

`analysis_euk_robustness.py`, `results/euk_robustness.{tsv,json}`, and a rewritten headline section in
`RESULTS.md` only if checks 1–3 pass.
