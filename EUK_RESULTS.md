# Where fungi sit among eukaryotes across the water–sediment boundary

2026-10-21. Checks fixed in `EUK_ROBUSTNESS_PLAN.md`
(sha256 `a12794f91b0df75b433df5fe37afb346ab3a603b5ed53152ca1a3dee85f879f8`, commit `a5c9fe6`)
before any of them was run. Levins' niche breadth was added after that plan was frozen but before its
results were seen (commit `9e01ec4`), to confront the published literature on its own metric.

## The claim, stated precisely

> Across 204 paired water/sediment sites at matched sequencing depth, **fungi sit at the generalist end of
> the eukaryotic range**: their habitat preference is weaker than that of every other eukaryotic kingdom
> measured here, significantly so against Metazoa and Straminipila, and **indistinguishable from
> Alveolata**, the most generalist protists.

The earlier phrasing — "fungi partition less than other eukaryotes" — is withdrawn as imprecise. It pools
protists with plants and animals, and the per-kingdom analysis below shows that is what produced the gap.
Against protists specifically there is no significant difference.

## Design

One assay, one reaction, one DNA extract. FunAqua's README: eukaryotic long-read metabarcoding with the
**ITS** region targeted; the pipeline's two tables are `fungiotu` and `eukotu` ("non-fungal eukaryotes").
So fungi and other eukaryotes here are the same samples, the same day, the same primers and the same
sequencing run, split only by taxonomy. 204 sites qualify for both groups; both samples of a pair are
rarefied to that pair's common minimum across all four of fungi/eukaryotes × water/sediment
(median 211 reads), so detection power is equal by construction.

**Primary metric is the median effect magnitude**, |preference − 0.5|, not the count of significant
orders. The count answers "how many orders cleared a threshold" and so depends on each order's number of
discordant pairs; the effect magnitude does not. Both were already computed and agree; this is a change of
emphasis, decided in the plan before the checks were run, not a swap to a better-looking number.

## Primary result, over 200 rarefaction draws

| | median effect | range over 200 draws | share of orders significant | orders tested |
|---|---|---|---|---|
| fungi | **0.156** | 0.140 – 0.173 | 47.4 % (39.6–53.7) | 91 – 99 |
| other eukaryotes | **0.228** | 0.202 – 0.257 | 61.4 % (54.8–67.5) | 115 – 127 |

**The two ranges do not overlap, and the ordering reversed in 0 of 200 draws.** This matters because the
fungal primary endpoint elsewhere in this project *is* draw-dependent (43–57 significant orders); this
comparison is not.

## The five pre-fixed checks: all pass

| check | what it tests | result | verdict |
|---|---|---|---|
| 1 | draw stability, 200 draws | ordering reversed 0/200 times | **PASS** |
| 2 | difference survives holding power fixed | stratified rank test p = 0.0023; fungi lower in 4 of 4 discordant-count bins | **PASS** |
| 3 | fungi lower than each kingdom separately | lower than 5 of 5; significantly against Metazoa (p = 4.4e−05) and Straminipila (p = 0.015) | **PASS** |
| 4 | holds at neighbouring ranks | Class p = 0.018, Order p = 0.00037, Family p = 0.0011 | **PASS** |
| 5 | survives adjusting for OTUs per order | partial Spearman −0.264, p = 7.3e−05; on the overlapping OTU-count range p = 1.4e−04 | **PASS** |

Check 2 is the one that mattered most. Fungal and non-fungal orders differ in how many discordant pairs
they have, and that alone sets statistical power. Stratifying into quantile bins of discordant count and
testing within bins leaves the difference intact in every bin.

Check 5 disposes of the aggregation worry. Within fungi, orders built from more OTUs do show weaker
apparent preference (Spearman −0.246, p = 0.0144) — real saturation. But non-fungal eukaryote orders carry
*more* OTUs (median 256 vs 142) while showing *stronger* preference, so saturation cannot be producing the
gap; adjusting for it widens the gap rather than closing it.

## Per kingdom, on one metric at a time

Same 204 pairs, same depth, two compartments. Levins' B: 1 = found in one compartment only,
2 = even split. Effect: 0 = no preference, 0.5 = absolute preference.

| kingdom | orders (preference) | median effect | orders (B) | Levins B per order |
|---|---|---|---|---|
| **Fungi** | 62 | **0.150** | 111 | 1.528 |
| Alveolata | 39 | 0.157 | 54 | **1.665** |
| Rhizaria | 7 | 0.208 | 12 | 1.530 |
| Viridiplantae | 17 | 0.227 | 45 | 1.208 |
| Straminipila | 19 | 0.271 | 35 | 1.478 |
| Metazoa | 13 | **0.294** | 37 | **1.106** |

On the preference metric fungi rank first of six, most generalist. On Levins' B per order they rank third,
behind Alveolata and Rhizaria — but Rhizaria rests on 12 orders and 5 OTUs, so it carries little weight.
The two metrics rank the kingdoms similarly (Spearman −0.78) though with six kingdoms that is p ≈ 0.07.

The biologically coherent part: **animals and land-plant-related lineages are the most
compartment-bound** — structurally, an animal or a plant is benthic or planktonic — while fungi and the
protist groups, which disperse as spores or zoospores and use dissolved and particulate substrate in both
compartments, are the most spread out. Nothing in the method knows this.

## Levins' niche breadth: the literature's metric, on the same data

The published niche-breadth work reports the opposite ordering — in freshwater sediment,
bacteria < fungi < protists, i.e. fungi *narrower* than protists; in alpine soils fungi narrowest of all
groups. Those studies use Levins' B across samples, not paired compartment preference, so the first
question is whether the two metrics are even measuring the same thing here.

**They are.** On this data, Levins' B and the preference index are strongly correlated:
Spearman **−0.590** (p = 1.7e−10) within fungi and **−0.782** (p = 2.2e−26) within other eukaryotes.
So "these are different quantities" is not available as an explanation.

Levins' B across all 408 samples, the literature's formulation:

| unit | fungi | other eukaryotes | p | broader |
|---|---|---|---|---|
| per OTU | **1.888** | 1.652 | 5.6e−04 | fungi |
| per Order | 7.328 | 6.736 | 0.082 | fungi (n.s.) |

Per kingdom, per OTU across all samples:

| group | n OTUs | median B |
|---|---|---|
| Alveolata | 587 | **2.057** |
| Straminipila | 199 | 1.888 |
| **Fungi** | 1315 | **1.888** |
| Cryptista | 24 | 1.639 |
| Viridiplantae | 244 | 1.264 |
| Metazoa | 190 | 1.142 |

**This is where the apparent conflict with the literature dissolves.** Fungi tie with Straminipila and sit
below Alveolata. Pooled against protists (Alveolata + Straminipila + Rhizaria) the difference is not
significant (p = 0.163). So this dataset does not reproduce the reported fungi-narrower-than-protists
ordering, but neither does it significantly contradict it: fungi and protists are indistinguishable here.
What the data does show clearly is fungi broader than **plants and animals**, which those studies did not
include.

Methodological note: the *two-compartment* per-OTU Levins B is floored — most OTUs with ≥ 20 reads occur
in only one compartment, so B = 1.000 for four of six kingdoms and the statistic cannot discriminate. Only
the per-order two-compartment version and the across-samples version are informative. The per-OTU
two-compartment column is reported above for completeness and should not be read as a ranking.

## Limitations

**ITS is a fungal barcode.** This is the central caveat and belongs in any abstract, not a methods
footnote. Non-fungal eukaryotes here are off-target amplification: their recovery is taxon-dependent and
their ITS reference coverage is far thinner than UNITE's is for fungi. The design controls sample, day,
assay and depth; it cannot make the two groups equally well served by the primers or the reference
database. The claim is about what one ITS assay sees in paired water and sediment. Generalising it is a
hypothesis for a multi-marker test, not a settled fact.

**No external validation is possible in NEON.** Checked: NEON's paired water/benthic taxonomy products
carry only 16S and ITS, and the ITS tables are 100 % Fungi (37 493 and 16 919 rows, no other kingdom).
Macroinvertebrate metabarcoding is benthic-only and zooplankton metabarcoding covers 7 sites and water
only — different organisms, different sampling, different primers. They cannot replicate this contrast,
and presenting them as if they could would be dishonest. This result therefore rests on one dataset.

**Orders are not independent.** The checks address power and aggregation, not phylogenetic
non-independence among the orders being compared.

**But the two groups ARE treated at comparable taxonomic resolution** — established while checking
metaPR2 (`METAPR2_RESULTS.md`). At the rank this analysis uses, a fungal order aggregates 0.96× as much
diversity as a non-fungal one (145 vs 151 OTUs per label; Class 0.95, Family 0.76). So the comparison is
not between a coarse group and a fine one. This matters because coarser aggregation demonstrably weakens
apparent preference (Spearman −0.246 within fungi), and it is the first thing a reviewer should ask about.

**Attempted external replication in metaPR2 was abandoned, deliberately.** Its 18S has the opposite
primer bias, which would have made it the ideal test, but PR2 assigns fungi systematically coarser ranks
than protists (1.5–2.2× more diversity per label at every usable rank), which would have manufactured
this result rather than tested it. Full reasoning in `METAPR2_RESULTS.md`.

**Rhizaria and Cryptista are thin** (5–24 OTUs) and their positions should not be quoted.

**Shallow data.** Median common depth 211 reads. Presence/absence at order level is defensible;
abundance-based inference is not.

## Files

`EUK_ROBUSTNESS_PLAN.md` (frozen) · `analysis_euk_matched.py` (the primary) ·
`analysis_euk_robustness.py` → `results/euk_robustness.{tsv,json}` ·
`analysis_levins.py` → `results/levins.json` ·
`analysis_kingdom_consistent.py` → `results/kingdom_consistent.{tsv,json}` ·
`results/matched_fungi.tsv`, `results/matched_euk.tsv`, `results/matched_summary.json`
