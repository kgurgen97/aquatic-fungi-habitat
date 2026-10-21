# FunAqua — first results

2026-10-21. Analysis as declared in `ANALYSIS_PLAN.md`, written before any water-vs-sediment outcome was
computed. Unit: taxonomic Order. Within each site both samples are rarefied to that pair's minimum depth,
so detection power is equal by construction and a difference cannot be attributed to sequencing effort.

## Primary: fungal orders do partition between water and sediment

239 qualifying pairs (both samples ≥ 100 reads), 100 orders testable (≥ 20 discordant pairs).
**49 orders show a significant preference** (Benjamini–Hochberg q < 0.05): 30 water, 19 sediment.

> **Corrected 2026-10-21.** This read "99 orders testable, 51 significant, 30 water / 21 sediment".
> Those came from the ≥ 100 row of `analysis_secondary.py`, not from `analysis_primary.py`. The two
> implement the same declared endpoint with the same seed but differ in rarefaction draw, because
> `analysis_secondary.py` consumes one module-level RNG across its four thresholds. The pre-registered
> primary is `analysis_primary.py` (seed 20261021), which reproduces `results/primary_order.tsv`
> bit-identically on rerun: **100 tested, 49 significant, 30 / 19**.
>
> **The count is draw-dependent.** Over 200 draws (`analysis_draw_stability.py`, post hoc) the number of
> significant orders ranges **43–57** (median 51, sd 2.5), orders tested 96–103. Both stored runs fall
> inside that range, so neither was miscoded — the single-draw design is unstable at about ±5 orders.
> **The defensible claim is the 35 orders significant in ≥ 95 % of draws** (21 water, 14 sediment;
> 28 in every draw), in `results/draw_stability.tsv`. Use that set, not one draw's count, in anything
> written up.

The pattern is biologically coherent, which is the strongest evidence that the pipeline works — nothing in
the method knows about fungal biology:

| prefers the water column | prefers the sediment |
|---|---|
| Zygophlyctidales (94 vs 14 pairs), Chytridiales, Zygorhizidiales, Rozellales — **zoosporic lineages with a motile flagellate stage** | Sordariales, Microascales, Zoopagales — soil and dung saprotrophs, amoeba parasites, i.e. lineages of terrestrial origin |
| Taphrinales, Pucciniales, Dothideales, Mycosphaerellales | Lulworthiales (marine, wood-degrading) |

**Sensitivity to the read threshold**

| threshold | pairs | orders tested | significant |
|---|---|---|---|
| ≥ 50 | 353 | 105 | 54 |
| ≥ 100 | 239 | 99 | 51 |
| ≥ 200 | 127 | 84 | 26 |
| ≥ 400 | 69 | 56 | 14 |

These are `analysis_secondary.py`'s rows, each computed from a different point in one sequential RNG, so
every count here carries the ±5-order draw uncertainty measured above; the ≥ 100 row is **not** the
pre-registered primary (which is 100 tested / 49 significant). The count falls with the threshold because
N falls; the direction and identity of the strongest orders do not.

## Habitat preference is taxonomically conserved

Orders within the same class share a preference more than chance allows: mean within-class variance of the
preference index is **0.0243 observed vs 0.0350 expected** (10 000 permutations of order→class),
**p = 0.0041**. Based on 9 classes with ≥ 3 tested orders, covering 63 orders.

Reported as **taxonomic**, not phylogenetic. Three routes to a real phylogeny were tried and all are closed:

| route | outcome |
|---|---|
| ask the authors for representative sequences | would not help: ITSx shows the amplicon carries only ~61 bp of LSU |
| derive sequences ourselves from the raw reads (downloaded, 5.7 GB) | same amplicon limit |
| map orders onto a published tree (Open Tree of Life) | 70 orders retrieved, 47 with a preference value, but the synthetic fungal tree is almost entirely polytomous — only 4 resolved terminal sister pairs |

On that Open Tree topology, phylogenetic distance does not predict difference in preference
(Spearman +0.025, p = 0.41). This is **not** evidence against conservatism: the tree has no branch lengths,
patristic distance reduces to counting nodes, and the resolution is too poor for the test to be
informative. The class-level permutation above remains the best available evidence, and the claim stays
taxonomic.

## Consistency: Europe and North America, not yet global

Of 33 orders significant overall and testable in both regions, **30 (91 %) point the same way**.
Disagreeing: Auriculariales, Pezizales, Trechisporales.

> **Caveat added 2026-10-21.** All three "disagreeing" orders are in the draw-unstable set: across 200
> draws they are significant in only 27, 136 and 185 of 200 respectively. Their disagreement is in
> substantial part rarefaction noise rather than biogeography, and the 91 % figure should be read as an
> agreement among orders of mixed stability. Restricting to the draw-stable 35 would be the cleaner test
> and has not been run.

**This is not a global claim.** Only Europe (152 pairs) and North America (37) clear the threshold; at
≥ 50 reads Asia joins with 30 pairs. Africa, South America, Oceania and Antarctica have too few
sufficiently deep pairs. The dataset is global in sampling but not in usable paired depth.

## The pattern is identical in rivers and lakes

| water body | pairs | orders testable | significant |
|---|---|---|---|
| freshwater river | 95 | 103 | 20 |
| freshwater lake | 116 | 96 | 32 |

Of the 45 orders significant overall and testable in both, **45 of 45 (100 %) point the same way**, and the
preference indices correlate at Spearman **0.79** between the two biome types.

This is stronger evidence of reproducibility than the continental comparison (91 % over 33 orders) and it
is the better framing: the water–sediment preference is a property of the compartment, not of the kind of
water body. Marine and brackish sites have too few deep pairs to join the test.

## Fungi partition less than other eukaryotes

On the **same 204 pairs at the same rarefaction depth** (median 211 reads):

| | orders tested | significant | share |
|---|---|---|---|
| fungi | 98 | 42 | 43 % |
| other eukaryotes | 122 | 77 | 63 % |

Fisher OR 0.44, **p = 0.0029**. Effect magnitude agrees: median |preference − 0.5| is 0.149 for fungi
versus 0.246 for other eukaryotes (Mann–Whitney **p = 0.0009**).

A first run compared 239 fungal pairs against 580 eukaryote pairs and was discarded: eukaryotes carry more
reads per sample, so more sites passed the filter and the apparent difference was power. The matched
comparison above replaces it.

So aquatic fungi are comparatively **generalist** across the water–sediment boundary, against a benchmark
measured in the same water bodies on the same day by the same assay. Plausible reading: spore dispersal
places many fungi in both compartments, whereas many protists and small metazoans are obligately planktonic
or benthic. This is the opposite of a "fungi are special" story, and it is the more informative result.

Broken down by group, using median |preference − 0.5|, which is less sensitive to sample size than the
count of significant orders:

| group | OTUs | pairs | orders tested | median effect |
|---|---|---|---|---|
| Metazoa | 10 928 | 145 | 19 | **0.340** |
| Straminipila | 13 849 | 137 | 28 | 0.230 |
| Viridiplantae | 11 090 | 161 | 29 | 0.196 |
| Alveolata | 17 705 | 389 | 55 | 0.192 |
| **Fungi** | 33 370 | 204 | 98 | **0.149** |

An interpretable gradient: animals separate the two compartments most sharply — they are typically either
benthic or planktonic — and fungi least. Caveat: read depth differs by group, so the number of qualifying
pairs differs (Rhizaria, Amoebozoa and Euglenozoa yielded none); the percentage of significant orders is
therefore not comparable across groups and only the median effect is reported here.

## Not supported

Functional lifestyle (FungalTraits `primary_lifestyle`) does **not** significantly predict preference
(Kruskal–Wallis p = 0.32). The ordering is suggestive — algal parasites 0.73 toward water, litter
saprotrophs 0.38 toward sediment — but the data do not support it. Reported as an observation only.

## Limitations

- Shallow data: median rarefaction depth 211 reads. Presence/absence at the order level is defensible;
  abundance-based inference is not.
- 239 of 920 paired sites survive the depth filter. The paper must state the usable N, not the headline one.
- Some "orders" are environmental clade labels (`Branch03`, `GS05`, `*_incertae_sedis`), not formal taxa.
- Only the primary model's multiple-testing correction is applied; the regional and functional analyses are
  secondary and uncorrected across their own comparisons.
- **The single-draw design is draw-dependent** (43–57 significant orders over 200 draws). Only the 35
  orders significant in ≥ 95 % of draws should be treated as established.
- **This result did not replicate in NEON.** See `NEON_RESULTS.md`: the phenomenon reproduces (37 of 103
  orders significant in NEON's own data) but the per-order direction does not (sign agreement
  14/27 = 51.9 %, p = 1.0). The zoosporic lineages do replicate, pre-specified and strongly. Any write-up
  must carry this; a per-order preference table cannot be presented as a general result.

## Files

`results/primary_order.tsv` · `results/regional_fractions.tsv` · `results/matched_fungi.tsv` ·
`results/matched_euk.tsv` · `results/secondary_summary.json` · `results/matched_summary.json` ·
`results/draw_stability.tsv` · external validation in `NEON_RESULTS.md`
