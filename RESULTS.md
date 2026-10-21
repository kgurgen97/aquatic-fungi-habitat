# Fungi are the generalists among aquatic eukaryotes

2026-10-21. Unit: taxonomic Order. Within each site both samples are rarefied to that pair's minimum
depth, so detection power is equal by construction and a difference cannot be attributed to sequencing
effort. The water-vs-sediment analysis follows `ANALYSIS_PLAN.md`, written before any outcome was
computed; the eukaryote comparison's robustness checks follow `EUK_ROBUSTNESS_PLAN.md`, likewise frozen
before they were run.

## Headline

> Across 204 paired water/sediment sites at matched sequencing depth and in one ITS assay, **fungi sit at
> the generalist end of the eukaryotic range**. Their habitat preference is weaker than that of every
> other eukaryotic kingdom measured, significantly so against Metazoa (p = 4.4e−05) and Straminipila
> (p = 0.015), and indistinguishable from Alveolata, the most generalist protists.

Median effect |preference − 0.5| is **0.156** for fungi against **0.228** for other eukaryotes. Over 200
rarefaction draws the two ranges do not overlap (0.140–0.173 vs 0.202–0.257) and **the ordering never
reverses**. Five robustness checks, fixed before they were run, all pass — including the decisive one,
that the difference survives holding statistical power fixed.

Animals and plants are the most compartment-bound, which is biologically coherent: an animal or a plant is
structurally benthic or planktonic, while fungi and protists disperse as spores or zoospores and use
substrate in both compartments. Nothing in the method knows any of this.

Full account, per-kingdom table, the Levins' niche-breadth cross-check and the limitations:
**`EUK_RESULTS.md`**. The central caveat travels with the claim: ITS is a fungal barcode, so non-fungal
eukaryotes are off-target amplification, and no external validation is possible in NEON because its paired
products carry no eukaryote marker.

## Supporting: fungal orders do partition between water and sediment

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

## The headline result in full

Moved to **`EUK_RESULTS.md`**, which carries the primary metric, the five frozen robustness checks, the
per-kingdom table on a consistent metric, the Levins' niche-breadth cross-check against the published
literature, and the limitations. Summary of what changed since this section was first written:

- **The metric changed.** Primary is now the median effect magnitude, not the count of significant
  orders, because the count depends on each order's discordant-pair count and so measures power as much
  as biology. Decided in `EUK_ROBUSTNESS_PLAN.md` before the checks ran. Both numbers agree.
- **The claim was narrowed and is now more precise.** "Fungi partition less than other eukaryotes" pooled
  protists with plants and animals. On a consistent per-kingdom metric fungi are indistinguishable from
  **Alveolata** (effect 0.150 vs 0.157; Levins B per order 1.528 vs 1.665) and the gap against "other
  eukaryotes" is driven by **Metazoa** (0.294) and **Straminipila** (0.271). The defensible statement is
  that fungi sit at the generalist end, not that they beat every eukaryote.
- **The apparent conflict with the literature dissolved.** Published work reports fungi *narrower* than
  protists. On this data, pooled against protists the difference is not significant (p = 0.163) — fungi
  and protists are indistinguishable — while fungi are clearly broader than plants and animals, which
  those studies did not include. And the two metrics are not measuring different things here: Levins' B
  and the preference index correlate at Spearman −0.59 and −0.78.
- **It survived every check**, including matching on power (stratified p = 0.0023, fungi lower in 4 of 4
  bins) and adjusting for OTUs per order (partial Spearman −0.264, p = 7.3e−05).

An earlier run compared 239 fungal pairs against 580 eukaryote pairs and was discarded: eukaryotes carry
more reads per sample, so more sites passed the filter and the apparent difference was power. The matched
comparison replaced it.

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
