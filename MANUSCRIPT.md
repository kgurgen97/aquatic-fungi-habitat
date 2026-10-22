# Draft manuscript

**Status: first full draft, 2026-10-22.** Every number below is quoted from a file in `results/`; none is
from memory. The reference list is incomplete and flagged as such in § References — do not submit before
it is completed and every citation checked against the actual paper.

Target: *Molecular Ecology Resources* (methods-leaning, Q1). Alternatives: *Environmental Microbiology*,
*ISME Communications*.

---

## Title

**What transfers between amplicon datasets and what does not: habitat specialisation of aquatic fungi
measured twice**

Alternative, biology-led: *Fungi are generalists among aquatic eukaryotes across the water–sediment
boundary*

## Running title

Transferability of habitat-specialisation estimates

## Authors

Gurgen Karapetyan¹

¹ *affiliation to be completed*

Correspondence: karagurgen23@gmail.com

---

## Abstract

*(248 words)*

Metabarcoding surveys are increasingly used to ask which taxa specialise between adjacent habitats, but
estimates of specialisation are rarely tested for transferability between datasets. We measured habitat
preference across the water–sediment boundary in FunAqua, a global paired survey of aquatic fungi
(50 152 OTUs, 2 375 samples), using a design in which both members of every site pair are rarefied to that
pair's own minimum depth so that detection power is equal by construction. At order level, 239 qualifying
pairs yielded 49 orders with a significant habitat preference, and the pattern was biologically coherent,
with zoosporic lineages in the water column and soil-derived saprotrophs in the sediment.

We then tested how far these estimates travel. Three findings follow. First, the count of significant
orders is unstable to the rarefaction draw alone: across 200 draws it ranges from 43 to 57, and only 35
orders are significant in at least 95 % of draws. Second, the per-order preferences do not transfer. In an
independent paired dataset (NEON surface water and benthos, 149 pairs, 23 US streams) habitat partitioning
is reproduced — 37 of 103 orders show a significant preference — but the direction of preference agrees
with FunAqua for only 14 of 27 comparable orders (51.9 %, p = 1.0). Third, what does transfer is a
group-level pattern: fungi sit at the generalist end of the eukaryotic range, with a median preference
magnitude of 0.156 against 0.228 for non-fungal eukaryotes in the same assay and samples, a difference
that survives five pre-registered robustness checks and never reverses across 200 draws.

We conclude that per-order specialisation estimates from single amplicon surveys should not be treated as
transferable, and we identify matched read depth and matched taxonomic rank as prerequisites for any
cross-dataset comparison.

**Keywords:** metabarcoding; habitat specialisation; aquatic fungi; reproducibility; niche breadth;
rarefaction

---

## 1 Introduction

*(draft — tighten and add citations)*

Amplicon surveys now routinely report which taxa prefer one habitat over an adjacent one, and such
statements are used to infer dispersal strategy, substrate use and ecological strategy more broadly. The
inference is attractive because it is cheap: a presence/absence table and a habitat label are enough. It
is also fragile in ways that are not usually examined. Sequencing depth differs between habitats, between
samples and between datasets; taxonomic resolution differs between reference databases and between the
groups being compared; and presence/absence at any fixed depth is a threshold applied to a continuous
abundance distribution.

Aquatic fungi are a good case in which to examine this. They are abundant in both the water column and the
sediment, they have well-described dispersal differences — notably the flagellated zoospore of the
chytrids and their relatives — and a global paired survey now exists. FunAqua samples the same sites in
both compartments across 87 countries, which removes the site-level confounding that bedevils comparisons
assembled from separate studies. Crucially, it also sequences the non-fungal eukaryotes that the same ITS
primers capture, which provides a within-assay comparison group: the same samples, the same day, the same
reaction.

The closest existing work on habitat conservatism in aquatic microbial communities explicitly discarded
all eukaryotic reads [ref: Malfertheiner et al. 2026], leaving the fungal question open. Published
niche-breadth comparisons that do include fungi report them as *narrower* than protists, in freshwater
sediment [ref: snail/sediment study] and in soils [ref: alpine soils study], which makes any opposite
result worth examining carefully rather than announcing.

We set out to answer one question — do fungal clades specialise between water and sediment — and found
that the more informative question was how far the answer travels. We report three things: the extent to
which a single rarefaction draw determines the headline count; the failure of per-order preferences to
replicate in an independent paired dataset; and the one result that does replicate and survives
stress-testing, namely where fungi sit relative to other eukaryotes.

**Hypotheses, and when they were fixed.** All analysis plans were written and committed to version control
before the corresponding outcomes were computed; commit identifiers are given in § 2.7 so that the order of
events is checkable rather than asserted.

---

## 2 Materials and methods

### 2.1 FunAqua data

FunAqua [ref: data descriptor] is a global survey of aquatic eukaryotes in which the internal transcribed
spacer (ITS) region was targeted by long-read metabarcoding. The release used here contains 50 152 fungal
OTUs and 87 656 non-fungal eukaryote OTUs across 2 375 samples (1 076 sediment, 1 299 filtered water) from
1 452 sites in 87 countries, with water and sediment chemistry. The table is 99.87 % zeros and the median
sample depth is 125 reads.

The fungal and non-fungal tables derive from **one amplicon, one reaction and one DNA extract**, split only
by taxonomy: the project's own pipeline names them `fungiotu` and `eukotu` ("non-fungal eukaryotes OTU
table"). This is the basis of the within-assay comparison in § 3.4 and also its principal limitation,
since ITS is a fungal barcode and non-fungal eukaryotes are off-target amplification (§ 4.4).

Samples join to metadata on `sample_code` (2 375 of 2 375). A full stop (`.`) is FunAqua's marker for an
unassigned rank and is excluded; treating it as a taxon name inflates the usable pair count from 577 to 625
and falsely implies complete naming at every rank.

### 2.2 Pair construction and the depth confound

Sediment samples are systematically deeper than water samples (median 179 against 131 reads,
p = 1.5 × 10⁻¹³) and richness tracks depth closely (Spearman 0.84). An unadjusted comparison would
therefore measure sequencing effort. We handle this in the design rather than by covariate adjustment:
**within each site pair, both samples are rarefied to that pair's own minimum depth** (multivariate
hypergeometric sampling, seed 20261021) before anything is counted. Detection power is then equal inside
every pair by construction, and a difference cannot be attributed to depth.

920 sites have both habitats. Requiring both samples of a pair to carry at least 100 reads after
rank filtering leaves **239 qualifying pairs** (median rarefaction depth 211 reads for the matched
eukaryote analysis).

### 2.3 Preference index and the primary test

Rarefied counts are reduced to presence/absence per taxonomic Order. For each order we count the pairs in
which it is present in water only (*w*) and in sediment only (*s*); concordant pairs carry no information
about preference. Orders with at least 20 discordant pairs are tested by a two-sided exact binomial
against 0.5, with Benjamini–Hochberg correction across tested orders. The **preference index** is
*w*/(*w* + *s*), and the **effect magnitude** is |index − 0.5|.

### 2.4 Draw stability

Because the primary endpoint depends on one random rarefaction, we repeated the primary analysis over
**200 independent draws** (`analysis_draw_stability.py`) and recorded, per order, the proportion of draws
in which it was significant.

### 2.5 The eukaryote comparison and its robustness checks

Fungi and non-fungal eukaryotes were compared on **the same 204 site pairs at the same per-pair depth**,
the depth being the minimum across all four of fungi/eukaryotes × water/sediment. An earlier version
compared 239 fungal pairs against 580 eukaryote pairs and was discarded: eukaryotes carry more reads, so
more sites passed the filter and the apparent difference was power.

Five checks were fixed in `EUK_ROBUSTNESS_PLAN.md` before any of them was run: (i) stability over 200
draws, with a declared failure threshold of reversal in more than 5 % of draws; (ii) the difference must
survive stratification on the number of discordant pairs, which is the direct determinant of power;
(iii) fungi must be lower than each eukaryotic kingdom separately; (iv) the result must hold at
neighbouring taxonomic ranks; (v) it must survive adjustment for the number of OTUs per order. The plan
also moved the primary metric from the count of significant orders to the median effect magnitude, on the
grounds that the count measures power as much as biology.

Levins' niche breadth *B* = 1/Σ*p*ᵢ² was computed on the same rarefied data, across the two compartments
and across all 408 samples, per order and per OTU, to allow comparison with the published niche-breadth
literature on its own metric.

### 2.6 External datasets

**NEON.** Surface water (DP1.20141.002) and benthic (DP1.20086.002) microbe community taxonomy, RELEASE
2025 or later. The superseded `.001` products were measured and rejected: they assign only 5.2 % of reads
to an order (104 837 of 2 032 229 reads over 123 files). The `.002` products emit a rank name for every
sequence regardless of support, so a per-rank confidence cutoff is mandatory; we used `orderConfidence`
≥ 0.80, the RDP/DADA2 convention. The frozen universe is the **298 site-months where both products have
data** (23 wadeable-stream sites, 2018-03 to 2024-10). Pairing was by site and calendar date; because
exact-date pairing yielded only 20 qualifying pairs, below the trigger of 60 fixed in advance, the
pre-declared ±7-day rule applied, giving **149 pairs** across 23 sites and 128 dates. All other
measurement was identical to § 2.2–2.3.

**metaPR2 — examined and deliberately not used.** metaPR2 v4.0 (187 datasets, 18S rRNA, PR2 5.1.0) was
obtained and assessed as a second external set, attractive because 18S has the opposite primer bias to
ITS. It was rejected on a measurable ground: PR2 assigns fungi systematically coarser ranks than protists.
Of 43 distinct fungal labels in its `order` column only 19 are genuinely order-rank; twelve are subphyla,
and Ascomycota and Basidiomycota together collapse into six subphylum labels. ASVs per label run 1.5–2.2×
higher for fungi at every rank coarse enough for the available 44 paired coastal samples. Since coarser
aggregation demonstrably weakens apparent preference (§ 3.4), a test on that data would have produced a
spurious confirmation. For contrast, the same audit on FunAqua gives near parity — a fungal order
aggregates 0.96× as much diversity as a non-fungal one (Class 0.95, Family 0.76).

### 2.7 Pre-registration and availability

Analysis plans were committed before the corresponding results: `bbc9892` (NEON validation plan) precedes
`3d130c3` (its result); `a12b677` (flagellum-trait plan) precedes `9dec095` (its result); `a5c9fe6`
(eukaryote robustness checks) precedes `ab5a24e` (their outcome); Levins' *B* was added at `9e01ec4`,
before those outcomes were seen. Code and results: `github.com/kgurgen97/aquatic-fungi-habitat`.

---

## 3 Results

### 3.1 Fungal orders partition between water and sediment

Of 100 orders with enough discordant pairs to test, **49 show a significant habitat preference**
(BH *q* < 0.05): 30 toward water, 19 toward sediment. The pattern is biologically coherent, which matters
because nothing in the method has access to fungal biology. Zoosporic lineages with a motile flagellate
stage go to the water column — Zygophlyctidales (94 water-only against 14 sediment-only pairs),
Chytridiales, Zygorhizidiales, Rozellales — while Sordariales, Microascales and Zoopagales, all soil and
dung organisms, go to the sediment. Lulworthiales, marine and wood-degrading, is strongly sedimentary
(index 0.074).

Preference is conserved within taxonomic classes more than chance allows: mean within-class variance of
the preference index is 0.0243 observed against 0.0350 expected under 10 000 permutations of order to
class, **p = 0.0041**, over 9 classes containing 63 tested orders. We report this as *taxonomic* and not
phylogenetic conservatism; three independent routes to a resolved fungal phylogeny were attempted and all
are closed (§ 4.5).

### 3.2 The headline count is unstable to the rarefaction draw

Repeating the primary analysis over 200 draws gives a count of significant orders ranging from **43 to 57**
(median 51, sd 2.5), with 96 to 103 orders tested. Two independently written implementations of the same
declared endpoint, with the same seed, returned 49 and 51 significant orders; neither was miscoded, and
both fall inside this range.

Per order, the instability is concentrated:

| | orders |
|---|---|
| significant in every draw | 28 |
| significant in ≥ 95 % of draws | **35** |
| significant in some draws, under 95 % | 37 |
| never significant | 35 |

**The defensible unit of reporting is therefore the set of 35 draw-stable orders** (21 water, 14 sediment),
not any single draw's count. This has a direct consequence for our own earlier reading: the three orders
we had identified as disagreeing in direction between Europe and North America — Auriculariales, Pezizales
and Trechisporales — are significant in only 27, 136 and 185 of 200 draws respectively, so that
disagreement was substantially draw noise rather than biogeography.

### 3.3 Per-order preferences do not transfer to an independent dataset

NEON is not a null dataset. On its own it reproduces the phenomenon: of 103 orders testable across 149
paired water/benthic samples, **37 show a significant preference** (18 water, 19 benthic), against
FunAqua's 49 of 100.

What fails to transfer is which order prefers which habitat. Of the 49 FunAqua-significant orders, 27 were
testable in NEON, and the direction agrees for **14 of 27 — 51.9 %** (95 % CI 34.0–69.3 %, two-sided exact
binomial **p = 1.0**). This was the pre-registered primary endpoint of the validation, with replication
declared only at agreement above 50 % with p < 0.05; the agreement is at chance. The preference index
correlates weakly and non-significantly between datasets (Spearman +0.281, p = 0.156 over the 27 compared
orders; +0.186, p = 0.126 over all 69 testable in both). Agreement is flat at 44–52 % across every read
threshold and confidence cutoff examined, so this is not a knife-edge result.

Two orders reverse significantly: Filobasidiales (0.646 → 0.167) and Amphisphaeriales (0.676 → 0.292),
water-preferring in FunAqua and benthic in NEON; and GS11 (0.277 → 0.721) in the other direction.

**A candidate explanation, tested and found insufficient.** Within-pair rarefaction equalises detection
power inside a pair but cannot make "present" mean the same thing in two datasets. FunAqua's pairs carry a
median of 131–179 reads per sample; NEON's within-pair rarefaction depth is a median of **31 131** reads, a
roughly 200-fold difference. At 130 reads presence means "among the most abundant sequences"; at 31 000 it
means "detectable at all". Capping NEON's pairs at FunAqua-like depths raises agreement monotonically from
51.9 % to about 69 %, the direction the hypothesis predicts, but it never reaches significance (best
p = 0.0755) and capping shrinks the comparison to 16–26 orders. Depth mismatch is therefore a credible but
unproven contributor, and a design requirement for future comparisons rather than an established
explanation.

**A trait-level hypothesis also failed to replicate.** Because the per-order comparison was underpowered at
the level of individual named orders, we pre-registered a trait test: orders in flagellated (zoosporic)
fungal phyla should show a higher water preference than those in aflagellate phyla, assigned by phylum
membership from fungal biology before any group value was computed. NEON supports it emphatically
(flagellated median 0.858 against aflagellate 0.379, Mann–Whitney p = 5.1 × 10⁻⁸, rank-biserial +0.823;
16 of 18 flagellated orders prefer water against 18 of 83 aflagellate, Fisher OR 28.9, p = 1.3 × 10⁻⁷).
FunAqua does not, and what direction it has runs the other way (0.439 against 0.548, p = 0.179). Under the
pre-declared rule — predicted direction with p < 0.05 in both datasets independently — the hypothesis is
**not replicated**.

Post hoc, the two flagellated phyla behave oppositely, and this is itself consistent across datasets:
Chytridiomycota prefer water more than Rozellomycota in both (p = 0.012 and p = 0.0013). Narrowing the
trait to chytrids does not rescue it (NEON p = 6.0 × 10⁻⁷ with rank-biserial 0.971; FunAqua p = 0.269 with
a 0.04 difference in medians). We note that 89–100 % of the Rozellomycota "orders" in these datasets are
undescribed environmental clade labels (GS05, GS07, GS10, GS11, Branch02, Branch03) whose dispersal stage
has never been observed, so assigning them a flagellum on the strength of *Rozella* alone was weaker
evidence than the Chytridiomycota assignment and should not have been treated as equivalent.

### 3.4 What does transfer: fungi sit at the generalist end of the eukaryotic range

On the same 204 pairs at the same depth, the median effect magnitude is **0.156 for fungi against 0.228 for
non-fungal eukaryotes**. Over 200 rarefaction draws the two ranges are 0.140–0.173 and 0.202–0.257 — they
do not overlap — and **the ordering reverses in 0 of 200 draws**. By the count metric, 42 of 98 fungal
orders are significant (43 %) against 77 of 122 non-fungal (63 %), Fisher p = 0.0029.

All five pre-registered checks pass:

| check | result |
|---|---|
| draw stability (200 draws) | 0 reversals; threshold was ≤ 5 % |
| matched on power (discordant-count bins) | stratified rank test p = 0.0023; fungi lower in 4 of 4 bins |
| each kingdom separately | fungi lower than 5 of 5; significantly against Metazoa (p = 4.4 × 10⁻⁵) and Straminipila (p = 0.015) |
| neighbouring ranks | Class p = 0.018, Order p = 3.7 × 10⁻⁴, Family p = 0.0011 |
| adjusted for OTUs per order | partial Spearman −0.264, p = 7.3 × 10⁻⁵; on the overlapping OTU-count range p = 1.4 × 10⁻⁴ |

The power check is the one that mattered. Fungal and non-fungal orders differ in how many discordant pairs
they have, and that alone sets statistical power; stratifying into quantile bins of discordant count leaves
the difference intact in every bin.

The aggregation check disposes of the obvious artefact. Within fungi, orders built from more OTUs do show
weaker apparent preference (Spearman −0.246, p = 0.0144) — saturation is real. But non-fungal eukaryote
orders carry *more* OTUs per order (median 256 against 142) while showing *stronger* preference, so
saturation cannot produce the gap; adjusting for it widens the gap.

**The claim must be stated narrowly.** On a consistent per-kingdom metric fungi are indistinguishable from
Alveolata, and the gap against "non-fungal eukaryotes" as a pooled group is driven by Metazoa and
Straminipila:

| kingdom | orders | median effect | Levins *B* per order |
|---|---|---|---|
| **Fungi** | 62 | **0.150** | 1.528 |
| Alveolata | 39 | 0.157 | **1.665** |
| Rhizaria | 7 | 0.208 | 1.530 |
| Viridiplantae | 17 | 0.227 | 1.208 |
| Straminipila | 19 | 0.271 | 1.478 |
| Metazoa | 13 | **0.294** | **1.106** |

Animals and plants are the most compartment-bound, which is biologically coherent: an animal or a plant is
structurally benthic or planktonic, whereas fungi and the protist groups disperse as spores or zoospores
and use substrate in both compartments.

### 3.5 The result is not an artefact of our choice of metric

Published niche-breadth work reports the opposite ordering — fungi *narrower* than protists in freshwater
sediment and narrowest of all groups in alpine soils — using Levins' *B* across samples rather than paired
compartment preference. The first question is therefore whether the two metrics measure the same thing
here. **They do**: Spearman between Levins' *B* and the preference index is **−0.590** (p = 1.7 × 10⁻¹⁰)
within fungi and **−0.782** (p = 2.2 × 10⁻²⁶) within non-fungal eukaryotes. "These are different
quantities" is not available as an explanation.

On Levins' own formulation across all 408 samples, fungi are broader per OTU (1.888 against 1.652,
p = 5.6 × 10⁻⁴) and non-significantly broader per order (7.328 against 6.736, p = 0.082). Per kingdom,
per OTU: Alveolata 2.057, Straminipila 1.888, **Fungi 1.888**, Cryptista 1.639, Viridiplantae 1.264,
Metazoa 1.142. **Pooled against protists the difference is not significant (p = 0.163)** — fungi and
protists are indistinguishable here — while fungi are clearly broader than plants and animals, which those
studies did not include. The apparent conflict with the literature therefore dissolves rather than
requiring explanation.

A methodological note: the *two-compartment* per-OTU Levins *B* is floored, since most OTUs with ≥ 20 reads
occur in only one compartment, giving *B* = 1.000 for four of six kingdoms. Only the per-order
two-compartment version and the across-samples version are informative.

---

## 4 Discussion

### 4.1 Specialisation estimates are less transferable than they appear

The usual way to present an analysis like § 3.1 is a table of taxa with their preferred habitat. Our
results argue against that presentation. The identity and direction of per-order preferences agreed with
an independent paired dataset at chance (51.9 %, p = 1.0), even though both datasets independently detect
habitat partitioning at similar rates (49/100 and 37/103). A per-order preference table from a single
survey describes that survey.

What did transfer was coarser: a group-level contrast with a large effect, measured within one assay on
one set of samples. The practical reading is that cross-dataset claims should be pitched at the level of
aggregate contrasts, and that per-taxon claims need either replication or an explicit statement that they
have not been replicated.

### 4.2 Two prerequisites for cross-dataset comparison

Our two failed attempts at external validation identify the same two requirements, from opposite
directions.

**Matched depth.** FunAqua and NEON differ ~200-fold in within-pair depth. Within-pair rarefaction, which
correctly removes the depth confound *inside* a dataset, does nothing about this. Presence/absence is a
threshold on an abundance distribution, and moving the threshold by two orders of magnitude changes which
slice of the community carries the signal. Capping NEON toward FunAqua's depth moved agreement from 51.9 %
to ~69 %, consistent with this being a contributor, though underpowered to confirm.

**Matched rank.** The metaPR2 attempt failed for a different reason that is just as general: a reference
database may assign systematically coarser ranks to one of the groups being compared. PR2 collapses
Ascomycota and Basidiomycota into six subphylum labels while resolving chytrids to real orders. Because
coarser aggregation weakens apparent preference (§ 3.4), such a test would have confirmed our hypothesis
for the wrong reason. We report the rank-parity audit (fungi 0.96× non-fungal in FunAqua at order level)
as a check that should accompany any between-group comparison.

### 4.3 Single-draw rarefaction is not a reporting unit

The count of significant orders moved between 43 and 57 on the same data under nothing but a different
random draw, and our own two implementations landed at 49 and 51. Where rarefaction is used — and for
paired designs with unequal depth it is the cleanest available fix — the result should be reported over
many draws, and the reporting unit should be the set of taxa that are stable across them. We found 35 of
100 orders stable at the 95 % level.

### 4.4 Limitations

**ITS is a fungal barcode.** This is the central limitation of § 3.4 and belongs with the claim wherever it
is quoted. Non-fungal eukaryotes here are off-target amplification: their recovery efficiency is
taxon-dependent and their ITS reference coverage is far thinner than UNITE's is for fungi. The design
controls sample, day, assay and depth; it cannot make the two groups equally well served by the primers or
the reference database. The claim is about what one ITS assay sees in paired water and sediment.
Generalising it is a hypothesis for a multi-marker test.

**No external validation of § 3.4 was possible.** NEON's paired products carry only 16S and ITS, and its
ITS tables are 100 % Fungi (37 493 and 16 919 rows, no other kingdom); metaPR2 fails the rank-comparability
requirement above. Macroinvertebrate metabarcoding is benthic-only and zooplankton metabarcoding covers 7
sites and water only — different organisms, sampling and primers. The result therefore rests on one
dataset, and the obvious next step is a paired survey carrying ITS and 18S on the same DNA.

**Shallow data.** Median common depth 211 reads. Presence/absence at order level is defensible;
abundance-based inference is not.

**Non-independence.** Orders are not independent observations. The checks address power and aggregation,
not phylogenetic non-independence.

**Thin groups.** Rhizaria (7 orders) and Cryptista (24 OTUs) positions should not be quoted.

### 4.5 Why the conservatism result is taxonomic and not phylogenetic

Three routes to a resolved fungal phylogeny were attempted and all are closed. ITSx shows the amplicon
carries a median of only ~61 bp of LSU, so neither the authors' representative sequences nor the raw reads
(downloaded, 5.7 GB) can resolve deep fungal relationships; deep fungal phylogeny needs D1/D2, ~600 bp.
Mapping orders onto the Open Tree of Life synthetic tree retrieved 70 orders, 47 with a preference value,
but that tree is almost entirely polytomous — 4 resolved terminal sister pairs — and carries no branch
lengths, so patristic distance reduces to counting nodes. On that topology phylogenetic distance does not
predict difference in preference (Spearman +0.025, p = 0.41), which is **not** evidence against
conservatism but a statement that the test is uninformative. The class-level permutation (p = 0.0041)
remains the best available evidence and the claim stays taxonomic.

---

## 5 Data availability

FunAqua is publicly available under CC-BY [ref: data descriptor and Zenodo DOI to be inserted]. NEON data
products DP1.20141.002 and DP1.20086.002 were obtained through the NEON API. metaPR2 v4.0 is archived at
Zenodo (10.5281/zenodo.21106507, MIT). Analysis code, all frozen analysis plans with their SHA-256 hashes,
and all result tables: `github.com/kgurgen97/aquatic-fungi-habitat`.

---

## 6 References — INCOMPLETE, DO NOT SUBMIT AS IS

These are the works the analysis actually leaned on. **Each entry still needs its full citation verified
against the paper itself; several are recorded here only by topic and first author.** Nothing in this list
should be cited without that check.

1. FunAqua data descriptor — *needs full citation.* University of Tartu / Estonian University of Life
   Sciences; dataset release used here distributed with `pipeline.zip` and a README describing ITS
   long-read metabarcoding of 1 076 sediment and 1 299 filtered water samples.
2. Malfertheiner et al. (2026) *Nature Ecology & Evolution* — community conservatism in aquatic microbial
   communities; eukaryotic reads excluded. *Verify authors, title, volume, pages.*
3. Reñé et al. (2023) *Marine Ecology* — vertical and temporal distribution of chytrids infecting diatoms,
   Gulf of Naples; reports higher chytrid diversity and relative abundance in sediment than water column.
   *Verify.*
4. Vaulot, Geisen, Mahé, Bass — metaPR² database of eukaryotic 18S rRNA metabarcodes. *Verify journal and
   year; the version used here is v4.0, Zenodo 10.5281/zenodo.21106507.*
5. PR² reference database v5.1.0 — *needs citation.*
6. Madin, J.S. et al. (2020) A synthesis of bacterial and archaeal phenotypic trait data.
   *Scientific Data* **7**, 170. doi:10.1038/s41597-020-0497-4. *(Used in the companion pseudogene study,
   not here — remove unless cited.)*
7. Freshwater sediment niche-breadth comparison across bacteria, fungi and protists (snail *Bellamya
   purificata* sediment study, MDPI *Microorganisms*, 2024) — reports bacteria < fungi < protists.
   *Verify.*
8. Alpine soil comparison of bacterial, archaeal, fungal and protist niche breadth
   (*Soil Biology & Biochemistry*, 2022) — reports fungi narrowest. *Verify.*
9. Benjamini & Hochberg (1995) — FDR control. *Standard; add full citation.*
10. Levins, R. (1968) *Evolution in Changing Environments* — niche breadth *B*. *Add edition and pages.*
11. ITSx — Bengtsson-Palme et al. *Add citation.*
12. Open Tree of Life synthetic tree — *add citation and the version/date retrieved.*
13. NEON data product documentation for DP1.20141.002 and DP1.20086.002, including the note that the
    `.001` products were superseded from RELEASE-2025. *Add accession dates.*

---

## Author's notes on what still has to happen before submission

1. **Complete and verify every reference.** Nothing above is safe to cite as written.
2. **Decide the framing.** The draft leads with transferability. The alternative leads with the fungal
   generalist result and relegates the failed replications to a later section; that reads better but is
   more exposed, because the headline would then rest on one dataset with an off-target amplicon.
3. **Figures — drawn.** `make_figures.py` produces all four as vector PDF and 600 dpi PNG in `figures/`,
   sized to 89 mm (single column) and 183 mm (double). Palette is categorical slots 1–2 of a validated
   default (blue `#2a78d6`, orange `#eb6834`); the validator passes all checks on all pairs, and the
   aqua slot was dropped because it warns below 3:1 against the light surface. Two layout bugs and one
   numerical inconsistency were found by looking at the rendered output and are fixed:
   Fig 3's mirrored layout clipped "Straminipila" to "Straminip", Fig 1 printed a blue annotation on top
   of blue bars, and Fig 1 counted stability against the draws in which an order was *testable* rather
   than against all 200, which made it say 36 where the text says 35. An assertion in the script now
   fails if the figure and `draw_stability.json` ever disagree again.

   | figure | panel | file |
   |---|---|---|
   | 1 | per-order significance frequency across 200 draws, 35-order stable set marked | `fig1_draw_stability` |
   | 2 | FunAqua against NEON preference index, 27 orders, agreement quadrants | `fig2_transfer_failure` |
   | 3 | median effect magnitude and Levins *B* by kingdom, two panels | `fig3_kingdoms` |
   | 4 | agreement against NEON depth cap, both comparison sets | `fig4_depth_probe` |

   Still to decide: whether Fig 4 belongs in the main text or the supplement, given that it is post hoc
   and reaches significance at no cap.
4. **Abstract word count** is 248; check the target journal's limit.
5. **Statement on the corrected numbers.** § 3.2 reports that our own two implementations disagreed. Decide
   whether that belongs in the Results as written, or in a methods note. It is honest either way and a
   reviewer will respect it, but it does draw attention.
