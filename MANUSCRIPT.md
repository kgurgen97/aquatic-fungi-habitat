# Sequencing depth changes which taxa meet metabarcoding evidence thresholds: implications for habitat-preference analysis

**Authors:** Gurgen Karapetyan¹; Astghik Z. Pepoyan¹

¹ Department of Food Safety and Biotechnology, Armenian National Agrarian University, 74 Teryan Street, Yerevan 0009, Armenia.

**Corresponding author:** Gurgen Karapetyan (karagurgen23@gmail.com)

**ORCID:** Gurgen Karapetyan — https://orcid.org/0009-0005-5509-2576; Astghik Z. Pepoyan — https://orcid.org/0000-0002-1935-5341

## Abstract

Amplicon metabarcoding studies often infer habitat preferences from the presence or absence of taxa. Before a taxon can be tested, it must pass a minimum-evidence threshold. We examined how that threshold interacts with sequencing depth, rarefaction and taxonomic assignment in two aquatic microbial datasets. In FunAqua, 49 of 100 testable fungal orders showed a significant water–sediment preference in the primary draw, but this count ranged from 43 to 57 over 200 rarefaction draws; 35 orders remained significant in at least 95% of draws. In NEON, a controlled depth ladder used the same 102 paired samples at 100–10,000 reads per sample. The sets of testable orders at the shallowest and deepest depths shared only 28% of their union, although the median absolute preference effect changed by 0.024. Simulations reproduced this contrast across 239 interpretable parameter combinations. Between FunAqua and NEON, a prespecified sign-agreement test was underpowered, while post hoc count-based tests detected differences for 26 of 69 shared testable orders after multiple-testing correction. These differences cannot be attributed to sequencing depth alone because the datasets also differ in habitat definitions and taxonomy. An estimator benchmark showed that performance depends on whether preference affects occurrence, abundance or both. Minimum-evidence thresholds therefore affect which taxa can be assessed, even when aggregate estimates appear stable. Comparisons should report taxon eligibility, uncertainty, sequencing depth and taxonomic resolution.

**Keywords:** metabarcoding; fungal ecology; sequencing depth; rarefaction; reproducibility; niche breadth

## 1 Introduction

Amplicon surveys now routinely report which taxa prefer one habitat over an adjacent one, and such
statements are used to infer dispersal strategy, substrate use and ecological strategy more broadly. The
inference is attractive because it is cheap: a presence/absence table and a habitat label are enough. It
is also fragile in ways that are not usually examined. Sequencing depth differs between habitats, between
samples and between datasets; taxonomic resolution differs between reference databases and between the
groups being compared; and presence/absence at any fixed depth is a threshold applied to a continuous
abundance distribution.

One consequence of that last point turns out to govern everything else, and it is the main result of this
paper. Every such analysis needs a minimum-evidence rule — a taxon must appear in enough samples, or
enough discordant pairs, before it can be tested at all. This filter also affects which taxa enter the analysis. A taxon
contributes evidence only while its detection probability sits away from both 0 and 1, which happens when
its expected read count is near one; so the rule admits taxa whose abundance multiplied by depth is near
one, at whatever depth is used. Raise the depth and the set that passes does not become more extreme — it
moves down the abundance distribution to rarer taxa that now sit in the same place on the detection curve.
Depth therefore changes *which* taxa a survey can measure rather than how specialised they appear.

This predicts a specific and testable pattern: the distribution of specialisation estimates should
reproduce between datasets of very different depth, while the estimates for individually named taxa should
not, because the two datasets are measuring different taxa. It also says which claims from such a survey
are safe to carry forward — aggregate contrasts between large groups, which average over whichever taxa
happened to pass — and which are not, namely any statement about a named taxon. We derive the mechanism
under known truth, confirm the pattern it predicts in two independent paired datasets, and show that a
prediction we first drew from the same algebra fails, for a reason that is itself the result.

That depth and rarefaction affect diversity metrics is not itself news: Shirazi et al. (2021) show that
read sampling depth moves alpha diversity, beta diversity and ordination position in eDNA metabarcoding.
The closest result to ours is García-Machado et al. (2023), who varied sequencing depth in lake fish eDNA
and found that depth strongly affected species occupancy *before* a minimum-read-proportion threshold was
applied, and had negligible effect on community assessment after it. They read that as the threshold
curing depth sensitivity. We read the same phenomenon the other way: the threshold conceals depth
sensitivity by changing *which* taxa are measured, so the stability it buys is stability of the marginal
distribution, not of any named taxon. That distinction is this paper's novelty claim, and §§ 3.4a–3.4c
are the test of it.

Aquatic fungi are a good case in which to examine this. They are abundant in both the water column and the
sediment, they have well-described dispersal differences — notably the flagellated zoospore of the
chytrids and their relatives — and a global paired survey now exists. FunAqua samples the same sites in
both compartments across 87 countries, which removes the site-level confounding that bedevils comparisons
assembled from separate studies. Crucially, it also sequences the non-fungal eukaryotes that the same ITS
primers capture, which provides a within-assay comparison group: the same samples, the same day, the same
reaction.

The closest existing work on habitat conservatism in microbial communities is prokaryotic by design:
Malfertheiner et al. (2026) "removed all eukaryotic reads to solely focus on the prokaryotic diversity",
and their phylum-level analysis covers 3 archaeal and 16 bacterial phyla. They do not discuss fungi, so
whether the pattern they describe extends to fungi is untested rather than disputed. Published
niche-breadth comparisons that do include fungi report them as *narrower* than protists, in pond sediment
(Hou et al., 2024) and in alpine soils, where fungi have the smallest mean niche breadth of four microbial
groups (Malard et al., 2022). That makes any opposite result worth examining carefully rather than
announcing.

We set out to answer one question — do fungal clades specialise between water and sediment — and found
that the measurement was more interesting than the answer. We report four things, in the order in which
they constrain each other. The mechanism above is the first: a derivation under known truth, with the
depth range over which the testable set holds its position on the detection curve. The second and third
are what it accounts for in real data: a single rarefaction draw moves the headline count of specialised
orders across a range of 14, and per-order preferences fail to replicate in an independent paired dataset
even though that dataset reproduces habitat partitioning in aggregate. The fourth is what does carry
across — a contrast between large groups, in which animals and plants are compartment-bound while fungi
and the protist groups are not — together with the limits on how narrowly that contrast may be stated.

**Hypotheses, and when they were fixed.** All analysis plans were written and committed to version control
before the corresponding outcomes were computed; commit identifiers are given in § 2.7 so that the order of
events is checkable rather than asserted.

---

## 2 Materials and methods

### 2.1 FunAqua data

FunAqua (Prins et al., 2026) is a global survey of aquatic eukaryotes in which the internal transcribed
spacer (ITS) region was targeted by long-read metabarcoding. The release used here contains 50 152 fungal
OTUs and 87 656 non-fungal eukaryote OTUs across 2 375 samples (1 076 sediment, 1 299 filtered water) from
1 452 sites in 87 countries, with water and sediment chemistry. The table is 99.87 % zeros and the median
sample depth is 125 reads across all 2 375 samples. (The water and sediment depths quoted in § 2.2 are
computed on the paired subset, which is deeper; the two figures describe different sample sets.)

The fungal and non-fungal tables derive from **one amplicon, one reaction and one DNA extract**, split only
by taxonomy: the project's own pipeline names them `fungiotu` and `eukotu` ("non-fungal eukaryotes OTU
table"). This is the basis of the within-assay comparison in § 3.5 and also its principal limitation,
since ITS is a fungal barcode and non-fungal eukaryotes are off-target amplification (§ 4.4).

Samples join to metadata on `sample_code` (2 375 of 2 375). A full stop (`.`) is FunAqua's marker for an
unassigned rank and is excluded; treating it as a taxon name inflates the usable pair count from 577 to 625
and falsely implies complete naming at every rank.

### 2.2 Pair construction and the depth confound

Across the 920 sites carrying both habitats, sediment samples are systematically deeper than water
samples (median 179 against 131 reads, p = 1.5 × 10⁻¹³) and richness tracks depth closely (Spearman 0.84). An unadjusted comparison would
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
against 0.5, with Benjamini–Hochberg correction (Benjamini and Hochberg, 1995) across tested orders. The **preference index** is
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

Levins' niche breadth *B* = 1/Σ*p*ᵢ² (Levins, 1968) was computed on the same rarefied data, across the two compartments
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

**metaPR2 — examined and deliberately not used.** metaPR2 v4.0 (Vaulot et al., 2022; 187 datasets, 18S
rRNA, annotated against PR² 5.1.0, Guillou et al., 2013) was
obtained and assessed as a second external set, attractive because 18S has the opposite primer bias to
ITS. It was rejected on a measurable ground: PR2 assigns fungi systematically coarser ranks than protists.
Of 43 distinct fungal labels in its `order` column only 19 are genuinely order-rank; twelve are subphyla,
and Ascomycota and Basidiomycota together collapse into six subphylum labels. ASVs per label run 1.5–2.2×
higher for fungi at every rank coarse enough for the available 44 paired coastal samples. Since coarser
aggregation demonstrably weakens apparent preference (§ 3.5), a test on that data would have produced a
spurious confirmation. For contrast, the same audit on FunAqua gives near parity — a fungal order
aggregates 0.96× as much diversity as a non-fungal one (Class 0.95, Family 0.76).

### 2.7 Pre-registration and availability

Analysis plans were committed before the corresponding results, and each plan's SHA-256 is recorded in
`PLAN_HASHES.md` rather than inside the plan file, so that a frozen plan is never edited. Eleven plans
were frozen; the later analysis scripts assert that the plan still matches its recorded hash before
computing anything.

| plan | frozen at | result at | section |
|---|---|---|---|
| `ANALYSIS_PLAN.md` | `bbc9892` | — | § 3.1 |
| `NEON_VALIDATION_PLAN.md` | `bbc9892` | `3d130c3` | § 3.3 |
| `ZOOSPORIC_PLAN.md` | `a12b677` | `9dec095` | § 3.3 |
| `EUK_ROBUSTNESS_PLAN.md` | `a5c9fe6` | `ab5a24e` | § 3.5 |
| `METAPR2_REPLICATION_PLAN.md` | `eba1edf` | `b298db6` (not executed) | § 2.6 |
| `DEPTH_LADDER_PLAN.md` | `464e481` | `5058008` | § 3.4a |
| `SWEEP_PLAN.md` | `7a0d398` | `fa995db` | § 3.4b |
| `PREVALENCE_PLAN.md` | `4ac0440` | `96be558` | § 3.4c |
| `BENCHMARK_PLAN.md` | `828f8da` | `6a94646` | § 3.4d |
| `BENCHMARK2_PLAN.md` | `6a94646` | `8c6df32` | § 3.4d |
| `OCCUPANCY_PLAN.md` | `be9490c` | `6fd306a` | § 3.4d |

Levins' *B* was added at `9e01ec4`, before those outcomes were seen. Analyses labelled **post hoc** in the
text are not covered by any plan: `analysis_draw_stability.py` (§ 3.2), `exploratory_neon_depth.py`
(Figure 4), `exploratory_zoosporic_byphylum.py` (§ 3.3), the log–log fits in § 3.4a, and
`analysis_agreement_power.py` (§ 3.3a). Code and results:
`github.com/kgurgen97/aquatic-fungi-habitat`.

### 2.8 The mechanism analyses

Five analyses support §§ 3.4–3.4d; each has its own frozen plan (§ 2.7) and its own results file.

**Simulation (§ 3.4, Figure 5).** `analysis_simulation.py` draws a community of 400 taxa log-uniformly
over 10⁻⁶ to 10⁻¹, applies a fixed true abundance ratio between compartments, detects each taxon in each
sample with probability 1 − e^(−*aD*), and runs the § 2.3 estimator on the result at each depth.

**Depth ladder (§ 3.4a, Figure 6a–b).** The 102 of 149 NEON pairs whose within-pair minimum depth reaches
10 000 reads, subsampled to 100 / 300 / 1 000 / 3 000 / 10 000 reads per sample, the same 102 pairs at
every rung, 50 draws per rung on shared seeds. Loading, within-pair rarefaction and the
≥ 20-discordant-pair threshold are taken unchanged from `analysis_neon_validation.py`.

**Parameter sweep (§ 3.4b, Figure 6c).** 240 cells: four abundance distributions (log-uniform, lognormal,
power law, and the empirical vector exported from the ladder run) × effect sizes 1, 2, 4, 8, 16 ×
thresholds of 10, 20 and 40 discordant pairs × 102 and 239 pairs × site heterogeneity σ = 0 and 1, each
at five depth rungs with 30 draws on shared seeds. A cell is *interpretable* when the log–log fit for the
re-selection exponent reaches R² ≥ 0.5; 239 of 240 do.

**Abundance representation (§ 3.4c, Figure 6a, 6d).** The same 102 pairs, rungs and threshold, with
abundance supplied to the detection model four ways: pooled to one value per order; the observed
per-sample read fractions; those fractions with every zero replaced by the order's smallest non-zero
value; and a two-parameter form keeping only each order's occupancy *p*ᵢ and its pooled abundance
rescaled to *â*ᵢ/*p*ᵢ where present.

**Estimator benchmark and occupancy model (§ 3.4d).** Six estimators on identical simulated counts over
three encodings of preference (abundance-only, occupancy-only, both) × three depths × three effect sizes,
30 datasets per cell. The occupancy model obtains its repeated visits from *K* = 8 rarefaction subsamples
of each library, so ψ is occupancy in the sequenced library rather than at the site.

**Post hoc agreement analyses (§ 3.3a).** `analysis_agreement_power.py` reads only the supplementary
tables and computes exact binomial power, the minimum detectable agreement, and intervals bootstrapped
over phyla rather than orders.

---

## 3 Results

### 3.1 Fungal orders partition between water and sediment

Of 100 orders with enough discordant pairs to test, **49 show a significant habitat preference**
(BH *q* < 0.05): 30 toward water, 19 toward sediment (Table S1). The pattern is biologically coherent, which matters
because nothing in the method has access to fungal biology. Zoosporic lineages with a motile flagellate
stage go to the water column — Zygophlyctidales (94 water-only against 14 sediment-only pairs),
Chytridiales, Zygorhizidiales, Rozellales — while Sordariales, Microascales and Zoopagales, all soil and
dung organisms, go to the sediment. Lulworthiales, marine and wood-degrading, is strongly sedimentary
(index 0.074).

One published vertical pattern runs the other way. Reñé et al. (2023) report that in the Gulf of Naples
the sediment chytrid community is both more diverse and more abundant than the water-column one. That
study is marine, single-site and based on 18S, and it measures standing diversity rather than paired
preference, so it is not a contradiction of a paired test across 239 site pairs; but it is a reminder that
the water-column signal we recover for zoosporic lineages is an average over those pairs, not a claim
about any one site.

Preference is conserved within taxonomic classes more than chance allows: mean within-class variance of
the preference index is 0.0243 observed against 0.0350 expected under 10 000 permutations of order to
class, **p = 0.0041**, over 9 classes containing 63 tested orders. We report this as *taxonomic* and not
phylogenetic conservatism; three independent routes to a resolved fungal phylogeny were attempted and all
are closed (§ 4.5).

### 3.2 The headline count is unstable to the rarefaction draw

Repeating the primary analysis over 200 draws gives a count of significant orders ranging from **43 to 57**
(median 51, sd 2.5), with 96 to 103 orders tested and 107 tested in at least one draw (Figure 1,
Table S2). Two independently written implementations of the same
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
not any single draw's count. That set is internally reproducible: of the 45 orders significant overall and
testable in both freshwater rivers (95 pairs) and lakes (116 pairs), **45 of 45 point the same way**, and
the preference indices correlate at Spearman 0.79 between the two water-body types — so the water–sediment
preference is a property of the compartment rather than of the kind of water body, within FunAqua. This has a direct consequence for our own earlier reading: the three orders
we had identified as disagreeing in direction between Europe and North America — Auriculariales, Pezizales
and Trechisporales — are significant in only 27, 136 and 185 of 200 draws respectively, so that
disagreement was substantially draw noise rather than biogeography.

### 3.3 Directional agreement with an independent dataset is not demonstrated

NEON is not a null dataset. On its own it reproduces the phenomenon: of 103 orders testable across 149
paired water/benthic samples, **37 show a significant preference** (18 water, 19 benthic), against
FunAqua's 49 of 100.

What fails to transfer is which order prefers which habitat. Of the 49 FunAqua-significant orders, 27 were
testable in NEON. One of those 27, Coniochaetales, has a NEON index of exactly 0.500 — no direction at all
— and is therefore excluded as uninformative rather than counted as an agreement, which is how
`results/neon_vs_funaqua.tsv` had recorded it. On the remaining 26 the direction agrees for **13 — 50.0 %**
(95 % Wilson CI 32.1–67.9 %, two-sided exact binomial **p = 1.0**). This was the pre-registered primary
endpoint of the validation, with replication declared only at agreement above 50 % with p < 0.05.
The preference index correlates weakly and non-significantly between datasets (Spearman +0.281, p = 0.156
over the compared orders; +0.186, p = 0.126 over all 69 testable in both). Agreement is flat at 44–52 %
across every read threshold and confidence cutoff examined (Table S7), so this is not a knife-edge result.
Per-order values are in Tables S7 and S8; the comparison is drawn in Figure 2.

**The pre-specified zoosporic lineages did transfer, and that is reported here rather than only in the
trait test below.** Secondary endpoint 3 of the validation plan named four zoosporic orders in advance.
Two were testable in NEON and both agree in direction with the strongest *q*-values in the whole NEON
table — Zygophlyctidales 0.870 → 0.875 (*q* = 1.8 × 10⁻⁴) and Chytridiales 0.752 → 0.933
(*q* = 2.3 × 10⁻¹¹) — as does Microascales, 0.150 → 0.323 (*q* = 0.022), while Sordariales shows no NEON
preference (0.168 → 0.491). With two zoosporic orders there is no formal test and this is a consistent
observation rather than a demonstrated effect, but it is a pre-registered endpoint that came out positive
and it is the part of the FunAqua pattern that survived.

Two orders reverse significantly: Filobasidiales (0.646 → 0.167) and Amphisphaeriales (0.676 → 0.292),
water-preferring in FunAqua and benthic in NEON; and GS11 (0.277 → 0.721) in the other direction.

**A candidate explanation, tested and found insufficient.** Within-pair rarefaction equalises detection
power inside a pair but cannot make "present" mean the same thing in two datasets. FunAqua's pairs carry a
median of 131–179 reads per sample; NEON's within-pair rarefaction depth is a median of **31 131** reads. That is
roughly 200-fold the 131–179 reads a FunAqua sample carries, and ~150-fold the 211-read matched depth at
which the FunAqua analyses are run; both denominators are used below and are named wherever quoted. At 130 reads presence means "among the most abundant sequences"; at 31 000 it
means "detectable at all". Capping NEON's pairs toward FunAqua-like depths raises agreement from 50.0 %
to a maximum of 69.2 % at a cap of 300 reads — the direction the hypothesis predicts — but the rise is
**not monotonic**: at a cap of 130 reads, the only rung inside FunAqua's own 131–179, agreement falls back
to 60.9 % (Figure 4, Table S11). It never reaches significance at any cap (best p = 0.0755), and capping
shrinks the comparison to 16–26 orders. That the trend reverses at the matched rung is awkward for the
hypothesis, which predicts its best agreement exactly there. Depth mismatch is therefore a credible but
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

### 3.3a What that endpoint could have detected, and a better-powered version of it

*(Post hoc; Supplementary Methods § S.6, Figure S4, Tables S15–S17.)* A non-replication is
only informative if the test could have detected a replication, and the frozen plan never asked what it
could detect. It cannot detect much. With 26 comparable orders, a two-sided exact binomial at α = 0.05
first reaches significance at **19 of 26 — 73.1 % agreement** — so the endpoint was declaring replication
only for near-total transfer. Its power is **0.46** against a true agreement of 70 % and 0.26 against
65 %; reaching 80 % power would have needed **49** comparable orders against a true 70 %, 90 against 65 %
and 199 against 60 % (Table S16, Figure S4b). Equivalence to chance is not established either: against a
± 15 percentage-point margin the 90 % interval for 13/26 is 34.6–65.4 %, which breaches both bounds. The
pre-registered endpoint therefore excludes neither replication nor its absence, and 13/26 is uninformative
*as a sign test*. It does not follow that the question is unanswerable: reducing each order to the sign of
its index discards the discordant-pair counts that carry the evidence, and § 3.3b recovers them.

**The same data support a better-powered endpoint.** Restricting to orders significant in FunAqua
discards 39 of the 65 orders that are testable in both datasets and conditions on significance in a design
§ 3.2 shows to be draw-dependent. Across all 65 (four further orders are excluded as exact-0.500 ties),
agreement is **35 of 65 — 53.8 %** (95 % Wilson CI 41.9–65.4 %, p = 0.62), with power 0.91 against a true
70 %. Its 90 % interval, 43.7–63.7 %, has an upper bound below 65 %: **strong per-order reproducibility is
excluded**, though chance is not. That is the one positive statement this comparison supports.

| selection set | *n* | agree | agreement | 95 % Wilson | 95 % phylum-clustered |
|---|---|---|---|---|---|
| **all orders testable in both** | **65** | 35 | **53.8 %** | 41.9–65.4 % | 38.2–66.3 % |
| FunAqua-significant (pre-registered) | 26 | 13 | 50.0 % | 32.1–67.9 % | 30.8–72.2 % |
| FunAqua draw-stable (§ 3.2) | 21 | 12 | 57.1 % | 36.6–75.5 % | 37.5–76.5 % |
| significant in both datasets | 9 | 6 | 66.7 % | 35.4–87.9 % | 0–85.7 % |

Agreement rises as the set is narrowed, which is what conditioning on significance in both datasets does
mechanically; no selection makes it demonstrable, and the last row rests on nine orders whose clustered
interval spans the scale. **Taxonomic non-independence does not overturn this.** Bootstrapping phyla
rather than orders — resampling whole phyla so that related orders travel together — widens the 65-order
interval only to 38.2–66.3 % (Figure S4a). Within the 26-order endpoint, agreement by phylum is
Ascomycota 7/12, Basidiomycota 2/8, Chytridiomycota 2/2 and Rozellomycota 2/4, so no phylum carries the
result. On the clustered interval the exclusion of strong reproducibility is marginal (90 % clustered
upper bound 65.3 %), and the claim should be stated as holding when orders are treated as independent.

### 3.3b Testing the counts rather than the signs: the disagreement is real

*(Post hoc; Supplementary Methods § S.7, Tables S18–S19.)* The sign test asks whether two datasets fall on
the same side of 0.5. A direct test asks whether they give the *same value*: for each order, a 2 × 2 Fisher
exact test of (water-only, sediment-only) in FunAqua against (water-only, benthic-only) in NEON, BH-corrected
across the 69 orders testable in both. This uses the counts the sign test throws away, and it answers the
question the sign test was too blunt for.

**They do not give the same value.** Twenty-six of 69 orders — **37.7 %** — differ significantly between
datasets. The disagreement is concentrated exactly where it matters:

| | significantly different (*q* < 0.05) | not distinguishable |
|---|---|---|
| agree in direction | 7 | 28 |
| **disagree in direction** | **16** | 14 |

Of the 30 orders that disagree in direction, **16 do so significantly**; within the pre-registered 26-order
endpoint, **11 of its 13 direction-disagreements are significant**. A further 7 orders agree in direction
while still differing significantly in magnitude. So the per-order disagreement between FunAqua and NEON is
**not** an artefact of imprecision: for a substantial minority of orders the two datasets genuinely measure
different preferences, and the sign test's failure to detect it was a property of the sign test.

This is the paper's claim stated positively rather than as a failure, and it is the form the claim should
take. §§ 3.3–3.3a establish that the pre-registered endpoint could not have demonstrated replication;
§ 3.3b establishes that, tested properly, the two datasets demonstrably disagree for at least 26 of 69
comparable orders. Heterogeneity tests have slightly more power where orders carry more discordant pairs,
and heterogeneous orders do carry marginally more (median minimum discordant count 38 against 33,
Mann–Whitney *p* = 0.081), so the 37.7 % is a lower bound on disagreement among well-measured orders
rather than an unbiased estimate across all of them.

**The result does not depend on the evidence threshold.** Recomputing both datasets' endpoints at
discordant-pair thresholds of 20, 30 and 40 gives agreement of 50.0 %, 52.6 % and 50.0 % over 26, 19 and 10
comparable orders; at 50 and 60 the NEON comparison collapses to 3 and 2 orders and is uninformative
(Table S19). The non-replication is therefore not an artefact of the ≥ 20 rule, which
§ 3.4a lists as an untested limitation.

**One asymmetry is worth recording, because it bounds § 3.4's claim.** Raising the evidence threshold at
*fixed* depth barely moves FunAqua's median |index − 0.5| (0.149 → 0.134 from threshold 20 to 60) but moves
NEON's a great deal (0.151 → 0.437). Depth-invariance of apparent effect size (§ 3.4b, S5) is therefore
specific to depth: the evidence *threshold* is a second axis, and in a deeply sequenced dataset raising it
does select more extreme taxa. The claim this paper makes — that changing depth changes which taxa are
measurable rather than how specialised they look — is unaffected, but it should not be generalised to
changes in the threshold itself.

### 3.4 The testability threshold decides which taxa are measured

Section 3.3 is an empirical failure. Simulation under known truth turns it into a
mechanism (`analysis_simulation.py`; Figure 5).

A taxon at relative abundance *a*, in a sample rarefied to depth *D*, is detected with probability
1 − (1 − *a*)^*D* ≈ 1 − e^(−λ), λ = *aD*. Writing *P*_w and *P*_s for the two compartments, the discordant
pairs — the only ones the test uses — split as *P*_w(1 − *P*_s) and *P*_s(1 − *P*_w), so

> E[index] = *P*_w(1 − *P*_s) / [ *P*_w(1 − *P*_s) + *P*_s(1 − *P*_w) ]

**which depends on *D* even when the true abundances are fixed** (Figure 5a). A taxon whose true
abundance ratio is a constant 4× has an expected index of 0.80 at FunAqua's depth and 1.00 at NEON's.
Presence/absence at 211 reads and at 31 131 reads are not noisy measurements of one quantity.

That algebra predicts that deeper data should show more extreme indices. **It does not.** NEON's median
|index − 0.5| is 0.151 against FunAqua's 0.149 — a ratio of 1.02, Mann–Whitney *p* = 0.69 — despite the
~150-fold difference in matched depth (211 against 31 131 reads). We report the failed prediction because the reason it fails is the result.

**The ≥ 20-discordant-pairs requirement is a filter.** A taxon contributes discordant pairs only when its
detection probability lies away from 0 and 1, so the rule favours taxa whose λ is near 1 *at whatever
depth is used*. In a simulated community spanning six orders of magnitude in abundance, log-uniformly
(Figure 5b): across a 1000-fold change in depth the median λ of the testable set holds at **0.22–0.49**,
the median |index − 0.5| holds at **0.304–0.350**, and the median *abundance* of that set falls
**454-fold**. Those three numbers are properties of that simulated community, not measurements; § 3.4a
reports what the same quantities do in real data, where the third holds, the second does not, and the
reason is informative.

**The domain of this result is narrower than the argument above.** The simulation fixes a single
true abundance ratio (3×) and one abundance distribution, and the selection it describes is strongest for
taxa with modest habitat effects. For a taxon with a large effect — reliably detected in one compartment
and essentially never in the other — discordant pairs accumulate over a much wider range of λ, including
λ ≫ 1, so such a taxon can pass the threshold without sitting near the middle of the detection curve. The
testable set is therefore determined jointly by both detection probabilities, the true effect, the number
of pairs and the abundance distribution. § 3.4b sweeps those axes. The escape route is real but weaker
than this argument suggests: raising the effect from 1× to 16× lowers the re-selection exponent only from
0.96 to 0.85 and raises turnover by 0.06, although a heavy-tailed community does drop to 0.57. The λ ≈ 1
statement is best read as a tendency that effect size erodes slowly, not as a law with a sharp domain.

So, within that regime, depth does not change how extreme the indices look. It changes **which taxa
occupy that position**.
The two datasets test names drawn from different parts of their own abundance distributions: they share
only **69 of 134** testable orders, that is, of the 134 orders testable in at least one dataset only 69
are testable in both (Jaccard 0.51), with 31 unique to FunAqua and 34 to NEON.

**This overlap is consistent with the mechanism but does not establish it.** FunAqua and NEON differ in
depth, but also in geography, water-body type and in what "sediment" means — river and lake sediment
against stream benthic biofilm. Orders testable in one dataset and not the other may therefore be absent
for ecological reasons rather than filtered out by depth, and the observed turnover cannot separate the
two. The clean test holds ecosystem, samples and reference database fixed and varies only depth; it is
§ 3.4a. What the simulation establishes is that the mechanism is *sufficient* to
produce such turnover under known truth, not that it is what produced this particular turnover.

This resolves the apparent paradox in § 3.3 and § 3.5 — the marginal distribution of indices is
reproducible between datasets while the per-order values are not, so aggregate contrasts survive a
comparison that per-taxon claims cannot.

Finally, the spread in § 3.2 is not peculiar to FunAqua. Simulating 239 pairs at depth 211 with the
observed per-order abundance distribution gives 48–57 significant orders across 200 draws (median 52,
sd 1.7) against the observed 43–57 (median 51, sd 2.5) — the same order of variation from the draw alone
(Figure 5c). The observed spread is somewhat wider, as expected given that the simulation assumes
independent taxa and a simplified truth.

### 3.4a The mechanism in real data, with ecology held constant

To test the mechanism without the cross-dataset confound we subsampled NEON's own libraries down a depth
ladder. The design was fixed in `DEPTH_LADDER_PLAN.md` before the analysis script existed, with five
predictions and explicit failure criteria. We used the **102 of 149 qualifying pairs whose within-pair
depth reaches 10 000 reads, the same 102 at every rung** (available depth median 69 440, minimum 10 641,
so no pair is upsampled), at rungs of 100, 300, 1 000, 3 000 and 10 000 reads per sample — a 100-fold
range whose lower end sits below FunAqua's own depth — with 50 draws per rung on shared seeds so that
draw noise is paired across rungs. Ecosystem, samples, reference database, pairing rule and test are
identical throughout; only depth changes.

| rung | testable orders | median *â* | median *λ̂* | median \|index − 0.5\| | Jaccard vs 10 000 | median \|Δindex\| |
|---|---|---|---|---|---|---|
| 100 | 36 | 1.11 × 10⁻² | 1.11 | 0.129 | **0.281** | 0.094 |
| 300 | 51 | 3.52 × 10⁻³ | 1.06 | 0.145 | 0.488 | 0.077 |
| 1 000 | 66 | 2.10 × 10⁻³ | 2.10 | 0.133 | 0.671 | 0.062 |
| 3 000 | 72 | 1.51 × 10⁻³ | 4.52 | 0.122 | 0.783 | 0.044 |
| 10 000 | 77 | 1.03 × 10⁻³ | 10.30 | 0.125 | 1.000 | 0.000 |

Figure 6a and 6b show the ladder: the abundance of the measurable set falling with depth, the set turning
over, and the apparent effect size staying put.

**Three results confirm the mechanism's consequence.** Depth strongly changes *which* orders are
measurable: the testable sets at 100 and at 10 000 reads share only 28 % of their union (P4, bar < 0.6).
The set moves down the abundance distribution, median *â* falling 10.8-fold (P1, bar ≥ 10-fold). Yet
depth barely changes *how specialised the measured taxa appear*: median |index − 0.5| moves by 0.024
across the whole 100-fold range, from 0.122 to 0.145 (P3, bar ≤ 0.10) — a sixth of the difference the
same quantity shows between FunAqua and NEON. And for the 25 orders testable at every rung the estimate
itself is stable, median |Δindex| 0.094 between the extreme rungs (P5, bar ≤ 0.10).

**One prediction did not hold, and it bounds the claim.** The selected set does *not* keep a fixed
position on the detection curve. Median *λ̂* rises from 1.06 to 10.30, a 9.74-fold span against a
pre-committed confirmation band of 3-fold; the frozen rule records this as inconclusive because the
failure threshold was 10-fold, but 9.74 is three times the band and should not be read as support.
Post hoc, the drift is a clean power law: median *â* ∝ *D*^−0.485 and median *λ̂* ∝ *D*^+0.515 (R² 0.92 and
0.93), where complete re-selection would give *â* ∝ *D*^−1 and flat *λ̂*. **The set re-selects about half
the way.** Part of the reason is in the first column: the testable set *grows* 2.15-fold, from 36 orders
to 77. Deeper sequencing is not a pure exchange of members — it retains the orders that were marginal at
100 reads, which by 10 000 reads sit far above λ = 1, and adds newly marginal ones, so the median of a
growing set drifts upward.

So the idealised form of the mechanism — a testable set pinned at λ ≈ 1 — is a property of that simulated
community. What holds in real data is the consequence that matters for interpretation: **changing depth
changes which taxa you can say anything about, while leaving the apparent strength of what you measure
almost untouched.** P1 and P5 cleared their thresholds narrowly, P5 rests on 25 orders, and the 102
deepest-sequenced pairs are not a random sample of NEON; the experiment is one dataset, one amplicon and
one threshold value. Full verdicts, including the two narrow margins, are in the depth-ladder results report in the
repository (§ 5).

### 3.4b How general is this, and does the model explain the exponent?

A sweep over 240 parameter combinations — four abundance distributions, effect sizes 1× to 16×, thresholds
of 10, 20 and 40 discordant pairs, 102 and 239 pairs, with and without site heterogeneity — was
pre-registered in `SWEEP_PLAN.md` with five predictions. It answers two questions and the answers differ
in kind.

**The depth-invariance of apparent effect size is universal.** In **239 of 239** interpretable cells the
median |index − 0.5| moves by ≤ 0.10 across the 100-fold depth range; the largest range anywhere is
**0.058** and the median is **0.0109**, with the ladder's observed 0.024 inside that distribution
(S5; Figure 6c).
This is the claim §§ 3.3–3.4a rest on, and it holds over the whole swept space rather than at one point.

**The model does not explain the measured exponent, and the explanation we first offered was wrong.** Fed
the *empirical* abundance vector, the detection model predicts *b* ≈ 0.96 — near-complete re-selection —
against the 0.485 measured, with simulated turnover near-total (Jaccard ≈ 0.000) where the ladder gave
0.281. No effect size in 1–16 brings *b* within 0.25 of 0.485; the nearest, *r* = 16, gives **0.852**,
off by 0.367 (S2 failed). Nor does the abundance
distribution's shape set the exponent: log-uniform gives 1.00 and the empirical vector 0.96, a difference
of 0.012–0.039 against a required 0.25 (S1 failed). Site heterogeneity moves *b* by at most 0.08. What
distribution shape *does* explain is the **growth** of the testable set — log-uniform grows 0.96–1.01-fold
while the empirical vector grows 2.15-fold, matching § 3.4a exactly.

Pooled, then, the model re-selects about twice as much as reality does. Of 240 cells only three reproduce
both *b* and Jaccard within 0.15, all with a lower effective bar than the experiment used, which lets
rarer taxa qualify at low depth. Two further predictions came out inconclusive and are reported in the sweep results report in the
repository (§ 5).

### 3.4c What was wrong was the pooling

The discrepancy turned out to be an artefact of representing each order by a single abundance. A
pre-registered follow-up (`PREVALENCE_PLAN.md`, frozen before its analysis existed) compared four ways of
supplying abundance to the same detection model, on the same 102 pairs, rungs and threshold.

| variant | *b* | Jaccard | testable, 100 → 10 000 |
|---|---|---|---|
| pooled *â*, one value per order (control) | +0.958 | 0.000 | 40 → 96 |
| the observed per-sample abundances | **+0.484** | **0.300** | **36 → 78** |
| per-sample, but with every zero filled in | +1.183 | 0.121 | 39 → 234 |
| occupancy + mean abundance only, two parameters | +0.432 | 0.366 | 45 → 109 |

Measured on real data: *b* = 0.485, Jaccard = 0.281, 36 → 77 orders. Figure 6a plots the four against the
measurement and Figure 6d their exponents. All five pre-registered predictions
were confirmed, including the control, which reproduced the pooled result (0.958 against the sweep's 0.96)
and so establishes that the two analyses agree where they should.

**The near-exact match in row two is close to arithmetic and is not the result.** Simulating detection
from a sample's own read fractions is a Poisson approximation to the rarefaction being modelled, so it
reproduces the real procedure by construction. The informative facts are the contrasts. Changing nothing
but whether abundance varies between samples moves the exponent from 0.958 to 0.484 — the whole gap.
Filling in the zeros while keeping all the graded variation sends it to 1.183 and inflates the testable
set three-fold, so it is **where a taxon is absent, not how much its abundance varies where present**,
that decides what is measurable at low depth. The median order here is present in only **9.8 %** of the
204 samples, 91 of 254 orders in 5 % or fewer. And two numbers per taxon suffice: occupancy plus mean
abundance, with the observed table discarded, recovers *b* = 0.432 against 0.485.

So the mechanism is quantitatively calibrated after all, with one amendment — abundance must be
represented per site. Pooling a sparse community into one value per taxon makes every taxon look uniformly
rare rather than locally common, which makes low-depth testability look impossible and predicts
near-total re-selection. The limits, including that an exact zero may be a sampling zero and that
occupancy measured on deep samples is itself depth-dependent, are set out in the repository (§ 5).

### 3.4d Which estimator to use, and why the question has no single answer

A warning is useful only if some estimator avoids the problem, so five were compared on the calibrated
generator: the single-draw presence/absence test used here, the draw-stable rule of § 3.2, a paired
Wilcoxon on proportions, a paired *t*-test on centred log-ratios, and an intercept-only overdispersed
binomial on read counts. Each saw identical simulated counts (`BENCHMARK2_PLAN.md`, frozen beforehand).

**The ranking depends on how the preference is encoded, which settles the question in the negative.**
When a preferring taxon is only more *abundant* in water, the read-count estimator wins by 0.852 in power
and presence/absence detects nothing at all. When it is instead more often *present* in water, that
advantage is **−0.002** — it vanishes. Any benchmark built on one generator recommends whichever
estimator shares its parameterisation; our own first attempt did exactly that and was discarded for it, as
the repository records (§ 5).

Under the ecologically realistic encoding, where a water-preferring taxon is both more often present and
more abundant in water, the practical ordering is clear:

| estimator | orders evaluable | power | FDR |
|---|---|---|---|
| single draw, presence/absence | 87 | 0.920 | **0.125** |
| draw-stable rule, ≥ 95 % of draws | 82 | 0.875 | **0.037** |
| overdispersed binomial on counts | 126 | 0.908 | 0.050 |
| paired *t*-test on CLR | 254 | 0.441 | 0.813 |
| paired Wilcoxon on proportions | 161 | 0.376 | 0.755 |

Three things follow. A single rarefaction draw runs at an FDR of 0.125 against a nominal 0.05. The
draw-stable rule brings that to 0.037 for 0.045 of power — **a second and independent reason for the
recommendation in § 3.2**, which was argued there only from reproducibility and turns out also to control
error. And applying no evidence floor is not the same as coverage: the CLR test evaluates every one of
the 254 orders at every depth, which looks like the remedy § 3.4a calls for, and runs an FDR above 0.8.

**A detection-explicit occupancy model was then added**, since it is the natural estimator for the
structure § 3.4c identifies. Occupancy models need repeated visits to separate occupancy ψ from detection
*p*, which one sample per site does not provide; *K* = 8 rarefaction subsamples of each library supply
them, giving a detection history per sample at no cost in data. ψ is then occupancy *in the sequenced
library*, so the model corrects subsampling and not sequencing. Its plan (`OCCUPANCY_PLAN.md`) committed
in advance to the reading of a failure, and it did fail its primary prediction: in the arm built for it,
where preference acts on occupancy alone, it ranks **fourth of six** on power (0.286 against the single
draw's 0.385). It must estimate *p* from the detection histories themselves, and where detections are
sparse that estimate is noisy — so the correction helps where data are already plentiful and hurts where
they are thin. **On power there is no case for it over the draw-stable rule.**

It earns two things that do not depend on that. Under the realistic encoding it evaluates **155 orders
against the draw-stable rule's 82**, 1.9 times as many, at power 0.833 against 0.875 and FDR 0.030
against 0.037 — and unlike the CLR test's nominal 254 that is coverage *with* an evidence floor rather
than coverage from having none. And in the abundance-only arm it is almost silent (power 0.016,
FDR 0.000) where `betabinom` fires at 0.852: alone among the six, its answer tracks its own estimand,
correctly reporting that occupancy does not differ when only abundance does.

The honest summary is that no estimator dominates and the choice depends on an encoding one does not know
in advance. Within the realistic encoding, the draw-stable rule gives the fewest false discoveries for
the least machinery; the occupancy model and the overdispersed binomial assess roughly twice and 1.5
times as many taxa at comparable error; and no method should be used without an evidence floor. The
occupancy model's case rests on coverage, not power, and it weakens below about 1 000 reads per sample.
Its obvious extension — a per-sample detection covariate, since detectability depends on a taxon's
abundance in each particular sample — is untested and would most likely repair exactly the low-depth
weakness that defeated it here.

### 3.5 A group-level contrast within FunAqua: animals and plants are compartment-bound, fungi and protists are not

On the same 204 pairs at the same depth, the median effect magnitude over 200 rarefaction draws is
**0.156 for fungi against 0.228 for non-fungal eukaryotes**; the two ranges are 0.140–0.173 and
0.202–0.257 — they do not overlap — and **the ordering reverses in 0 of 200 draws**. In the single
pre-registered draw tabulated in Tables S3 and S4 the same quantities are 0.149 and 0.246, both inside
those ranges; figures quoted in this section are the 200-draw values unless a table is named. By the count
metric, 42 of 98 fungal orders are significant (43 %) against 77 of 122 non-fungal (63 %),
Fisher p = 0.0029.

On the single pre-registered draw the difference in medians is **+0.097** (95 % bootstrap CI
+0.039 to +0.133), the Hodges–Lehmann shift is +0.062 (95 % CI +0.026 to +0.103), and the rank-biserial
correlation is **+0.261** (Mann–Whitney *p* = 8.7 × 10⁻⁴, 98 against 122 orders). The interval excludes
zero, but the rank-biserial places this as a **moderate** separation of two broadly overlapping
distributions, not a clean split: many fungal orders are more compartment-bound than many non-fungal ones.
The claim is about the central tendency of two groups, and should not be read as a statement about any
individual order.

All five pre-registered checks pass (Table S13, with per-order values in Table S12):

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

**The claim must be stated narrowly.** On a consistent per-kingdom metric (Table S5, Figure 3) fungi are
indistinguishable from Alveolata, and the gap against "non-fungal eukaryotes" as a pooled group is driven
by Metazoa and Straminipila:

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

**Why § 3.4's mechanism does not undermine this comparison.** The threshold that governs §§ 3.3–3.4a
applies identically to both groups here: the same 204 pairs, the same per-pair depth, the same
≥ 20-discordant-pair rule, so both groups are filtered to the same position on the detection curve, and
robustness check (ii) leaves the difference intact in 4 of 4 discordant-count bins. The sweep result is
the stronger argument: apparent effect magnitude — the quantity compared here — is the one thing that
proved depth-invariant, moving by ≤ 0.10 in 239 of 239 swept cells (§ 3.4b). A comparison of apparent
effect magnitudes between two groups measured at one common depth is therefore the kind of quantity the
mechanism predicts should be stable, which is why this contrast is reported while the per-order claims of
§ 3.3 are not. What the mechanism does not license is extrapolating the contrast to the taxa each group
would contribute at a different depth.

### 3.6 The result is not an artefact of our choice of metric

Published niche-breadth work reports the opposite ordering — fungi *narrower* than protists in pond
sediment (Hou et al., 2024) and smallest of four microbial groups in alpine soils (Malard et al., 2022).
Neither study measures what we measure. Malard et al. compute niche breadth along topoclimatic and edaphic
gradients, not between paired compartments; Hou et al. report a habitat niche-breadth measure and use
Levins' index only to split generalists from specialists, and their system is an indoor glass-tank
experiment with pond-sourced sediment rather than a field survey. The first question is therefore whether
a breadth measure and a paired-preference measure order taxa the same way on our own data. **They do**: Spearman between Levins' *B* and the preference index is **−0.590** (p = 1.7 × 10⁻¹⁰)
within fungi and **−0.782** (p = 2.2 × 10⁻²⁶) within non-fungal eukaryotes (Figure S1, Table S6).
"These are different quantities" is not available as an explanation.

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

### 4.1 The threshold, not the biology, decides what a survey can measure

The mechanism in § 3.4 is the result we would carry to another dataset, because it does not depend on
fungi, on this assay or on these habitats. Any analysis that tests taxa for habitat preference must first
decide which taxa carry enough evidence to test, and that decision admits taxa sitting near the middle of
the detection curve. Deeper sequencing moves the set that qualifies down the abundance distribution
instead of sharpening the estimates for the taxa already in it. Two datasets at different depths therefore
produce similar-looking distributions of specialisation over largely different taxa — here, 69 shared
names out of 134, although that turnover is confounded with the ecological differences between the two
datasets; the depth ladder inside NEON (§ 3.4a) isolates the effect, and finds the same turnover with
ecology held constant.

This reframes what looked like a replication failure. The usual way to present an analysis like § 3.1 is a
table of taxa with their preferred habitat, and our results argue against that presentation: the identity
and direction of per-order preferences showed no demonstrable agreement with an independent paired dataset
(50.0 % on the pre-registered 26 orders; 53.8 %, 90 % CI 43.7–63.7 %, across all 65 testable in both) even
though both datasets independently detect habitat partitioning at similar rates (49/100 and 37/103). The
pre-registered endpoint is uninformative by construction — it needed 73 % agreement to declare replication
and had 46 % power against a true 70 % (§ 3.3a). The defensible statements come from the two better-powered
analyses: across all 65 comparable orders agreement is 53.8 % with strong reproducibility excluded, and
testing the counts rather than the signs shows that 26 of 69 orders differ significantly between datasets,
including 11 of the 13 direction-disagreements inside the pre-registered set (§ 3.3b). Non-transfer of
per-order preferences is therefore demonstrated, not merely undemonstrated — but by estimators the
validation plan did not specify, and that is itself part of the lesson: a sign test on a filtered subset is
the wrong instrument for a replication question. Those two facts are not in tension once the threshold is taken into account — the first is a
statement about named taxa, the second about the marginal distribution, and the mechanism says the second
can hold while the first does not. A per-order preference table from a single survey describes that survey.

What did transfer was coarser: a group-level contrast with a large effect, measured within one assay on
one set of samples. The practical reading is that cross-dataset claims should be pitched at the level of
aggregate contrasts, and that per-taxon claims need either replication or an explicit statement that they
have not been replicated.

### 4.2 Two prerequisites for cross-dataset comparison

Our two failed attempts at external validation identify the same two requirements, from opposite
directions.

**Matched depth.** FunAqua and NEON differ ~200-fold in per-sample depth and ~150-fold in matched
per-pair depth. Within-pair rarefaction, which
correctly removes the depth confound *inside* a dataset, does nothing about this. Presence/absence is a
threshold on an abundance distribution, and moving the threshold by two orders of magnitude changes which
slice of the community carries the signal. Capping NEON toward FunAqua's depth moved agreement from 50.0 %
to a maximum of 69.2 % at a 300-read cap, consistent with this being a contributor, though underpowered to
confirm and non-monotonic — agreement falls back to 60.9 % at the 130-read cap that actually matches
FunAqua's depth.

**Matched rank.** The metaPR2 attempt failed for a different reason that is just as general: a reference
database may assign systematically coarser ranks to one of the groups being compared. PR2 collapses
Ascomycota and Basidiomycota into six subphylum labels while resolving chytrids to real orders. Because
coarser aggregation weakens apparent preference (§ 3.5), such a test would have confirmed our hypothesis
for the wrong reason. We report the rank-parity audit (fungi 0.96× non-fungal in FunAqua at order level)
as a check that should accompany any between-group comparison.

### 4.3 Single-draw rarefaction is not a reporting unit

The count of significant orders moved between 43 and 57 on the same data under nothing but a different
random draw, and our own two implementations landed at 49 and 51. Where rarefaction is used — and for
paired designs with unequal depth it is the cleanest available fix — the result should be reported over
many draws, and the reporting unit should be the set of taxa that are stable across them. We found 35 orders stable at the 95 % level, of the 107 testable in at least one draw.

That recommendation was made on reproducibility grounds alone, and § 3.4d supplies a second reason that
is arguably stronger: under simulation with known truth, a single draw runs at a false discovery rate of
**0.125** against a nominal 0.05, while requiring significance in ≥ 95 % of draws brings it to **0.037**
at a cost of 0.045 in power. The stable set is not merely the more reproducible thing to report; it is
the one that keeps the error rate near what is claimed for it.

### 4.4 Limitations

**ITS is a fungal barcode.** This is the central limitation of § 3.5 and belongs with the claim wherever it
is quoted. Non-fungal eukaryotes here are off-target amplification: their recovery efficiency is
taxon-dependent and their ITS reference coverage is far thinner than UNITE's is for fungi. The design
controls sample, day, assay and depth; it cannot make the two groups equally well served by the primers or
the reference database. The claim is about what one ITS assay sees in paired water and sediment.
Generalising it is a hypothesis for a multi-marker test.

**No external validation of § 3.5 was possible.** NEON's paired products carry only 16S and ITS, and its
ITS tables are 100 % Fungi (37 493 and 16 919 rows, no other kingdom); metaPR2 fails the rank-comparability
requirement above. Macroinvertebrate metabarcoding is benthic-only and zooplankton metabarcoding covers 7
sites and water only — different organisms, sampling and primers. The result therefore rests on one
dataset, and the obvious next step is a paired survey carrying ITS and 18S on the same DNA.

**Shallow data.** Median common depth 211 reads. Presence/absence at order level is defensible;
abundance-based inference is not.

**Multiplicity across analyses.** Benjamini–Hochberg correction is applied across orders within the
primary endpoint (§ 3.1) and within NEON (§ 3.3). It is *not* applied across the secondary analyses — the
regional and water-body comparisons, the functional-lifestyle test, the five robustness checks of § 3.5,
the per-kingdom contrasts and the Levins' *B* comparisons each stand uncorrected relative to one another.
Those are reported as described rather than as a family, and the per-kingdom positions in particular
should be read as descriptive.

**Non-independence.** Orders are not independent observations. For the agreement endpoint this is now
addressed — bootstrapping phyla rather than orders leaves the conclusion standing (§ 3.3a) — but the
§ 3.5 robustness checks address power and aggregation only, not phylogenetic non-independence among the
orders being compared.

**Thin groups.** Rhizaria (7 orders with a preference value, 12 with a Levins *B*) and Cryptista
(24 OTUs) positions should not be quoted.

### 4.5 Why the conservatism result is taxonomic and not phylogenetic

Three routes to a resolved fungal phylogeny were attempted and all are closed. ITSx (Bengtsson-Palme
et al., 2013) shows the amplicon carries a median of only ~61 bp of LSU, so neither the authors'
representative sequences nor the raw reads (downloaded, 5.7 GB) can resolve deep fungal relationships;
deep fungal phylogeny needs D1/D2, ~600 bp. Mapping orders onto the Open Tree of Life synthetic tree
(OpenTree, 2025: synthesis `opentree16.1`, taxonomy 3.7draft3; Hinchliff et al., 2015;
Redelings and Holder, 2017)
matched 77 of our 100 tested order names and placed 70 of them in the tree, 47 with a preference value,
but that tree is almost entirely polytomous — 4 resolved terminal sister pairs — and carries no branch
lengths, so patristic distance reduces to counting nodes. On that topology phylogenetic distance does not
predict difference in preference (Spearman +0.025, p = 0.41), which is **not** evidence against
conservatism but a statement that the test is uninformative. The class-level permutation (p = 0.0041)
remains the best available evidence and the claim stays taxonomic.

---

## 4a Figure captions

**Figure 1. Significance is partly a property of the rarefaction draw.** Each taxonomic order tested in
at least one of 200 independent within-pair rarefaction draws, ranked by the proportion of draws in which
it reached a significant habitat preference (Benjamini–Hochberg *q* < 0.05). Blue marks the 35 orders
significant in at least 95 % of draws; grey the remainder, of which 35 are never significant. Across the
200 draws the count of significant orders ranges from 43 to 57 (median 51, sd 2.5) and the number of
testable orders from 96 to 103. 239 site pairs, FunAqua.

**Figure 2. Per-order habitat preferences do not transfer between datasets.** Water preference index
(water-only discordant pairs / all discordant pairs) for the 27 orders significant in FunAqua and testable
in NEON. Blue points agree in direction between the two datasets, orange disagree; numbers give the count
in each quadrant. The dashed line is 1:1. One order plotted here, Coniochaetales, sits exactly on NEON's
0.500 and is excluded from the agreement count as having no direction (§ 3.3), giving 13 of 26 (50.0 %,
95 % Wilson CI 32.1–67.9 %, two-sided exact binomial *p* = 1.0); Spearman +0.28 (*p* = 0.16). The quadrant
counts shown therefore sum to 27 while the endpoint is computed on 26. FunAqua 239 pairs, NEON 149 pairs
across 23 US streams. Figure S4 gives the endpoint's power and the better-powered 65-order version.

**Figure 3. Animals and plants are the most compartment-bound; fungi are not distinguishable from
Alveolata.** (**a**) Median effect magnitude |preference − 0.5| and (**b**) median Levins' niche breadth
*B* across the two compartments, by kingdom, on the same 204 site pairs at the same per-pair rarefaction
depth. Lower preference and higher *B* both indicate a more generalist distribution. Fungi (blue) are
lowest on preference and third on *B*; the gap against non-fungal eukaryotes as a pooled group is driven
by Metazoa and Straminipila. Rhizaria rests on 7 orders and should not be quoted.

**Figure 4. Agreement with FunAqua under capping of NEON's sequencing depth.** Direction agreement for
the pre-registered comparison set (blue) and for the 35 draw-stable orders (orange), as NEON's within-pair
rarefaction depth is capped at successively lower values; "none" is the uncapped median of 31 131 reads.
FunAqua's own pairs carry a median of 131–179 reads per sample, the regime at the left-hand end. Agreement
rises as the cap falls to 300 reads and then **falls back at the 130-read cap**, so the rise is not
monotonic and does not peak at the matched depth; no cap reaches significance (best *p* = 0.076).
Post hoc; not part of the frozen validation plan.

**Figure 5. The testability threshold, not the depth itself, decides what is measured.**
(**a**) Expected preference index against rarefaction depth, from
E[index] = *P*_w(1 − *P*_s)/[*P*_w(1 − *P*_s) + *P*_s(1 − *P*_w)] with *P* = 1 − e^(−*aD*), for three fixed
true abundance ratios; the sediment abundance is 10⁻⁴ throughout. Dotted lines mark the two datasets'
depths. The index a study reports depends on how deeply it sequenced, at fixed biological truth.
(**b**) Median relative abundance of the orders that clear the ≥ 20-discordant-pairs requirement, against
depth, **in a simulated community** spanning 10⁻⁶ to 10⁻¹ drawn log-uniformly. In *this* community the
threshold holds the testable set at a fixed position on the detection curve (median λ 0.22–0.49) while the
abundance of the taxa occupying it falls 454-fold. Real data behave differently, for a reason that turns
out to be instructive: see § 3.4a and Figure 6. (**c**) Significant orders per draw in 200 simulated draws
of a FunAqua-like community (239 pairs, depth 211, observed abundance distribution), against the range
observed in the real data (shaded, median marked).

**Figure 6. The mechanism in real data, and what a model has to contain to reproduce it.**
(**a**) Median relative abundance of the orders clearing the ≥ 20-discordant-pairs threshold, against
rarefaction depth, for the 102 NEON pairs of § 3.4a. Heavy black rings are the measurement; coloured lines
are the same detection model supplied with abundance four different ways. The *observed per sample* line
lies under the rings — that match is close to arithmetic, as § 3.4c explains, and is not the result. The
result is the distance between the lines: pooling one abundance per order (orange squares) falls far too
steeply, and filling in the zeros (orange triangles) falls more steeply still.
(**b**) What depth changes and what it does not, on the same 102 pairs. Blue, left axis: overlap of the
testable set with the set at the deepest rung, rising from 0.281 to 1 as depth increases. Orange, right
axis: median |preference − 0.5| over the testable set, which moves only from 0.122 to 0.145 across the
same 100-fold range. **The two axes differ on purpose**; drawn on one 0–1 axis the orange line is pinned
to the floor and looks flat because Jaccard dominates the scale, which would argue the paper's claim by
squashing rather than by evidence.
(**c**) How general that flatness is. Distribution of the depth-range of median |preference − 0.5| across
the 239 interpretable cells of the § 3.4b sweep — four abundance distributions, effect sizes 1× to 16×,
three thresholds, two pair counts, site heterogeneity on and off. Solid line, the real ladder's 0.024;
dashed line, the pre-registered bar of 0.10, which no cell approaches.
(**d**) The re-selection exponent implied by each way of supplying abundance, against the measured 0.485
(solid line). Occupancy plus mean abundance — two numbers per taxon, with the observed table discarded —
lands at 0.43.

## 5 Data availability

FunAqua is described by Prins et al. (2026) and archived at Zenodo (10.5281/zenodo.17855266) under CC-BY.
The two NEON products, *Surface water microbe community taxonomy* (DP1.20141.002,
doi:10.48443/4xwg-p967) and *Benthic microbe community taxonomy* (DP1.20086.002, doi:10.48443/fe8q-qh09),
were retrieved through the NEON Data API on 21 October 2026. The request did not pin a release, so it
returned RELEASE-2026 (generated 23 January 2026) together with any provisional data then available; in
the event every site-month retrieved falls in 2018-03 to 2024-10, which RELEASE-2026 covers, so no
provisional month entered the analysis. metaPR2 v4.0 is archived at Zenodo (10.5281/zenodo.21106507, MIT).
The Open Tree query used synthesis `opentree16.1` (created 20 December 2025, taxonomy 3.7draft3), recorded
in `results/opentree_version.json`. Analysis code, all frozen analysis plans with their SHA-256 hashes
listed in `PLAN_HASHES.md`, and all result tables: `github.com/kgurgen97/aquatic-fungi-habitat`. The post
hoc analyses of § 3.3a are reproduced by `analysis_agreement_power.py`, which reads only the supplementary
tables distributed with this paper and writes Tables S15–S17 and Figure S4.

---

## Author contributions (CRediT)

Gurgen Karapetyan: Conceptualization, Methodology, Software, Formal analysis, Data curation, Visualization, Writing – original draft. Astghik Z. Pepoyan: Supervision.

## Funding

This research was supported by the Higher Education and Science Committee of the Republic of Armenia (grant numbers 21AG-4D065 and 25AA-4A009).

## Declaration of competing interests

The authors must review and complete Elsevier's Declaration of Interests form before submission.

## Declaration of generative AI and AI-assisted technologies in manuscript preparation

During manuscript preparation, the authors used OpenAI ChatGPT and Anthropic Claude to assist with language editing, manuscript organization and computational work. The authors reviewed the content and take responsibility for the analysis, interpretation and final manuscript.

## References

Benjamini, Y. and Hochberg, Y. (1995) Controlling the false discovery rate: a practical and powerful
approach to multiple testing. *Journal of the Royal Statistical Society Series B* **57**, 289–300.
doi:10.1111/j.2517-6161.1995.tb02031.x

Bengtsson-Palme, J., Ryberg, M., Hartmann, M., Branco, S., Wang, Z., Godhe, A. *et al.* (2013) Improved
software detection and extraction of ITS1 and ITS2 from ribosomal ITS sequences of fungi and other
eukaryotes for analysis of environmental sequencing data. *Methods in Ecology and Evolution* **4**,
914–919. doi:10.1111/2041-210X.12073  *(21 authors)*

Guillou, L., Bachar, D., Audic, S., Bass, D., Berney, C., Bittner, L. *et al.* (2013) The Protist
Ribosomal Reference database (PR2): a catalog of unicellular eukaryote Small Sub-Unit rRNA sequences with
curated taxonomy. *Nucleic Acids Research* **41**, D597–D604. doi:10.1093/nar/gks1160  *(32 authors;
version 5.1.0 used here)*

García-Machado, E., Normandeau, E., Côté, G., Caza-Allard, I., Babin, C. and Bernatchez, L. (2023) How eDNA data filtration, sequence
coverage, and primer selection influence assessment of fish communities in northern temperate lakes.
*Environmental DNA* **5**, 1216–1233. doi:10.1002/edn3.444

Hinchliff, C.E., Smith, S.A., Allman, J.F., Burleigh, J.G., Chaudhary, R., Coghill, L.M. *et al.* (2015)
Synthesis of phylogeny and taxonomy into a comprehensive tree of life. *Proceedings of the National
Academy of Sciences* **112**, 12764–12769. doi:10.1073/pnas.1423041112  *(22 authors)*

Hou, Y., Zhang, Y., Jia, R., Zhou, L., Li, B. and Zhu, J. (2024) Variations in the bacterial, fungal, and
protist communities and their interactions within sediment affected by the benthic organism, snail
*Bellamya purificata*. *Microorganisms* **12**, 2550. doi:10.3390/microorganisms12122550

Levins, R. (1968) *Evolution in Changing Environments: Some Theoretical Explorations*. Monographs in
Population Biology 2. Princeton, NJ: Princeton University Press, ix + 120 pp.

Malard, L.A., Mod, H.K., Guex, N., Broennimann, O., Yashiro, E., Lara, E. *et al.* (2022) Comparative
analysis of diversity and environmental niches of soil bacterial, archaeal, fungal and protist communities
reveal niche divergences along environmental gradients in the Alps. *Soil Biology and Biochemistry*
**169**, 108674. doi:10.1016/j.soilbio.2022.108674  *(9 authors)*

Malfertheiner, L., Tackmann, J., Matias Rodrigues, J.F. and von Mering, C. (2026) Community conservatism
is widespread across microbial phyla and environments. *Nature Ecology & Evolution* **10**, 232–245.
doi:10.1038/s41559-025-02957-4

National Ecological Observatory Network (2026) *Benthic microbe community taxonomy* (DP1.20086.002),
RELEASE-2026. doi:10.48443/fe8q-qh09. Retrieved through the NEON Data API, 21 October 2026.

National Ecological Observatory Network (2026) *Surface water microbe community taxonomy*
(DP1.20141.002), RELEASE-2026. doi:10.48443/4xwg-p967. Retrieved through the NEON Data API,
21 October 2026.

OpenTree (2025) *Open Tree of Life synthetic tree*, synthesis release `opentree16.1` (created
20 December 2025; Open Tree Taxonomy 3.7draft3; 2 064 source trees from 1 931 studies). Queried
22 October 2026. doi:10.5281/zenodo.3937741

Prins, V., Tedersoo, L., Mikryukov, V., Paiste, P., Sepp, M., Grossart, H.-P. *et al.* (2026) The FunAqua
dataset of global fungal biodiversity in aquatic ecosystems. *Scientific Data* **13**, 1346.
doi:10.1038/s41597-026-07729-8  *(240 authors; data archived at doi:10.5281/zenodo.17855266)*

Redelings, B.D. and Holder, M.T. (2017) A supertree pipeline for summarizing phylogenetic and taxonomic
information for millions of species. *PeerJ* **5**, e3058. doi:10.7717/peerj.3058

Reñé, A., Timoneda, N., Sarno, D., Piredda, R., Zampicinini, G., Zingone, A. *et al.* (2023) Vertical and
temporal distribution of chytrids infecting diatoms in the Gulf of Naples (Italy, Mediterranean Sea).
*Marine Ecology* **44**, e12726. doi:10.1111/maec.12726  *(8 authors; published online 13 December 2022,
issue dated 2023)*

Shirazi, S., Meyer, R.S. and Shapiro, B. (2021) Revisiting the effect of PCR replication and sequencing
depth on biodiversity metrics in environmental DNA metabarcoding. *Ecology and Evolution* **11**,
15766–15779. doi:10.1002/ece3.8239

Vaulot, D., Sim, C.W.H., Ong, D., Teo, B., Biwer, C., Jamy, M. and Lopes dos Santos, A. (2022) metaPR²: a
database of eukaryotic 18S rRNA metabarcodes with an emphasis on protists. *Molecular Ecology Resources*
**22**, 3188–3201. doi:10.1111/1755-0998.13674  *(version 4.0 used here,
doi:10.5281/zenodo.21106507)*

---

## 7 Author contributions

G.K. conceived the study, wrote all analysis code and analysis plans, performed the analyses and wrote the
manuscript.

## 8 Acknowledgements

*[to be completed]* — the FunAqua consortium for releasing the dataset under CC-BY, and the National
Ecological Observatory Network, a programme sponsored by the US National Science Foundation and operated
under cooperative agreement by Battelle, for the aquatic microbe community products used here.

## 9 Competing interests

The author declares no competing interests.

## 10 Funding

*[to be completed]*

## 11 ORCID

*[to be completed]*