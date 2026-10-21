# FunAqua analysis plan — written before computing any outcome

Short by design: this is an ecological re-analysis, not a clinical trial. Its purpose is to fix the
primary endpoint and the multiple-testing rule in advance, because with 346 orders × several ranks ×
6 continents × 10 biomes something "significant" will appear by chance if the analysis is left open.

## Question

Within a single sampling site, do fungal clades differ in occurrence between the water column and the
sediment; is that preference shared among related clades; and is it consistent across regions?

## Unit and sample

- **Analytic unit:** taxonomic **Order** (346 named orders, 82.8 % of reads). Family and Genus are
  pre-specified sensitivity levels.
- **Sites:** the 920 sites with both a water and a sediment sample. A site enters the primary analysis
  when **both** of its samples have ≥ 100 reads.
- `'.'` is FunAqua's unassigned marker and is excluded at every rank.

## The confounder, handled by design

Sediment samples carry about twice the reads of water samples (median 196 vs 94, p 1.5e-13) and richness
tracks depth (Spearman 0.84). Therefore **within each pair both samples are rarefied without replacement
to the minimum depth of that pair** before anything is counted. Detection power is then equal within a
pair by construction, and a difference cannot be attributed to sequencing effort.

## Primary endpoint

For each order, across qualifying pairs, classify each pair as water-only, sediment-only, both, or neither.
Test the discordant pairs (water-only vs sediment-only) with an **exact binomial test against 0.5**
(the paired, McNemar-style test), two-sided.

- Orders tested: those with **≥ 20 discordant pairs**.
- Multiple testing: **Benjamini–Hochberg FDR at 0.05** across the orders tested.
- Effect size: the water-only fraction of discordant pairs, with a 95 % Wilson interval.

## Secondary, declared now

1. **Taxonomic conservatism.** Do orders within the same class share a preference more than expected?
   Tested by permuting order→class labels (10 000 permutations) and comparing the within-class variance of
   the preference index. Reported as *taxonomic*, never phylogenetic: the reads carry only ~61 bp of LSU,
   so no defensible phylogeny exists (ITSx, this dataset).
2. **Regional consistency.** Repeat the primary test within each continent with ≥ 50 qualifying pairs and
   report the sign agreement across continents for the orders significant overall.
3. **Functional reading.** Map orders to FungalTraits `primary_lifestyle` via genus (91.5 % of genera
   covered) and ask whether preference differs between lifestyles.
4. **Eukaryote benchmark.** Repeat the primary endpoint on the non-fungal eukaryote table (87 657 OTUs,
   same samples, same pairs) to ask whether fungi are more or less habitat-partitioned than other
   eukaryotes in the same water bodies.

## Negative-result policy

If no order passes FDR, that is the result. Thresholds, ranks and the rarefaction rule are not revised
afterwards; any change is reported as a labelled sensitivity analysis.

Written 2026-10-21, before any water-vs-sediment outcome was computed.
