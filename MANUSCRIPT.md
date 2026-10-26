# Draft manuscript

**Status: first full draft, references verified, 2026-10-22.** Every number below is quoted from a file in
`results/`; none is from memory. Every reference was checked against Crossref or the publisher's record on
2026-10-22, and each cited finding was read in the paper itself rather than inferred from its title; the
cross-check that no citation lacks an entry and no entry lacks a citation is in § 6.

Target: *Molecular Ecology Resources* (methods-leaning, Q1). Alternatives: *Environmental Microbiology*,
*ISME Communications*.

---

## Title

**Minimum-evidence thresholds decide which taxa a metabarcoding survey can measure, not how specialised
they look: why per-taxon habitat preferences do not transfer between datasets**

Shorter alternative, same claim: *Depth changes which taxa are measurable, not how specialised they
appear: a mechanism for the non-transferability of per-taxon habitat preferences*

Rejected, biology-led: *Fungi are generalists among aquatic eukaryotes across the water–sediment boundary*
— § 3.5 shows fungi are indistinguishable from Alveolata and that the pooled gap is driven by Metazoa and
Straminipila, so this title would overstate the paper's own table.

## Running title

Thresholds decide what metabarcoding can measure

## Authors

Gurgen Karapetyan¹

¹ *affiliation to be completed*

Correspondence: karagurgen23@gmail.com

---

## Abstract

*(250 words)*

Metabarcoding is used to ask which taxa specialise between adjacent habitats, but such estimates are
rarely tested for transferability. Any such analysis needs a minimum-evidence threshold, which is not
neutral: it admits taxa whose expected read count is near one at any depth, so depth should change which
taxa are measurable rather than how specialised they appear. Subsampling NEON's libraries down a 100-fold
depth ladder, with ecosystem and samples fixed, confirms the consequence: the testable sets at 100 and
10 000 reads share 28 % of their union, while median apparent effect shifts by 0.024 — the latter in all
239 swept parameter sets.

In FunAqua, a global paired survey of aquatic fungi (2 375 samples), rarefying both members of each site
pair to its own minimum depth, 49 of 100 testable orders show a significant water–sediment preference —
zoosporic lineages in water, soil saprotrophs in sediment — but that count ranges from 43 to 57 across
200 draws.

An independent dataset 150-fold deeper (NEON, 149 pairs) reproduces habitat partitioning, with
37 of 103 orders significant, yet directional agreement with FunAqua was **not demonstrated**: 14 of 27
comparable orders, 51.9 % (95 % CI 34–69 %), excluding neither chance nor moderate reproducibility.

In FunAqua, animals and plants are most compartment-bound (0.294, 0.227) while fungi (0.150) and
Alveolata (0.157) are indistinguishable; this survives five pre-registered checks but could not be tested
externally.

Per-taxon specialisation estimates should not be assumed transferable; matched read depth and taxonomic
rank are prerequisites for cross-dataset comparison.

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

One consequence of that last point turns out to govern everything else, and it is the main result of this
paper. Every such analysis needs a minimum-evidence rule — a taxon must appear in enough samples, or
enough discordant pairs, before it can be tested at all. That rule is not a neutral filter. A taxon
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
sample depth is 125 reads.

The fungal and non-fungal tables derive from **one amplicon, one reaction and one DNA extract**, split only
by taxonomy: the project's own pipeline names them `fungiotu` and `eukotu` ("non-fungal eukaryotes OTU
table"). This is the basis of the within-assay comparison in § 3.5 and also its principal limitation,
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

### 3.3 Directional agreement with an independent dataset is not demonstrated

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

### 3.4 The testability threshold decides which taxa are measured

Section 3.3 is an empirical failure. Simulation under known truth turns it into a
mechanism (`analysis_simulation.py`; Figure 4).

A taxon at relative abundance *a*, in a sample rarefied to depth *D*, is detected with probability
1 − (1 − *a*)^*D* ≈ 1 − e^(−λ), λ = *aD*. Writing *P*_w and *P*_s for the two compartments, the discordant
pairs — the only ones the test uses — split as *P*_w(1 − *P*_s) and *P*_s(1 − *P*_w), so

> E[index] = *P*_w(1 − *P*_s) / [ *P*_w(1 − *P*_s) + *P*_s(1 − *P*_w) ]

**which depends on *D* even when the true abundances are fixed** (Figure 4a). A taxon whose true
abundance ratio is a constant 4× has an expected index of 0.80 at FunAqua's depth and 1.00 at NEON's.
Presence/absence at 211 reads and at 31 131 reads are not noisy measurements of one quantity.

That algebra predicts that deeper data should show more extreme indices. **It does not.** NEON's median
|index − 0.5| is 0.151 against FunAqua's 0.149 — a ratio of 1.02, Mann–Whitney *p* = 0.69 — despite the
~150-fold difference in depth. We report the failed prediction because the reason it fails is the result.

**The ≥ 20-discordant-pairs requirement is a filter.** A taxon contributes discordant pairs only when its
detection probability lies away from 0 and 1, so the rule favours taxa whose λ is near 1 *at whatever
depth is used*. In a simulated community spanning six orders of magnitude in abundance, log-uniformly
(Figure 4b): across a 1000-fold change in depth the median λ of the testable set holds at **0.22–0.49**,
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
(Figure 4c). The observed spread is somewhat wider, as expected given that the simulation assumes
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

Figure 5a and 6b show the ladder: the abundance of the measurable set falling with depth, the set turning
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
one threshold value. Full verdicts, including the two narrow margins, are in
`DEPTH_LADDER_RESULTS.md`.

### 3.4b How general is this, and does the model explain the exponent?

A sweep over 240 parameter combinations — four abundance distributions, effect sizes 1× to 16×, thresholds
of 10, 20 and 40 discordant pairs, 102 and 239 pairs, with and without site heterogeneity — was
pre-registered in `SWEEP_PLAN.md` with five predictions. It answers two questions and the answers differ
in kind.

**The depth-invariance of apparent effect size is universal.** In **239 of 239** interpretable cells the
median |index − 0.5| moves by ≤ 0.10 across the 100-fold depth range; the largest range anywhere is
**0.058** and the median is **0.0109**, with the ladder's observed 0.024 inside that distribution
(S5; Figure 5c).
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
rarer taxa qualify at low depth. Two further predictions came out inconclusive and are in
`SWEEP_RESULTS.md`.

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

Measured on real data: *b* = 0.485, Jaccard = 0.281, 36 → 77 orders. Figure 5a plots the four against the
measurement and Figure 5d their exponents. All five pre-registered predictions
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
occupancy measured on deep samples is itself depth-dependent, are in `PREVALENCE_RESULTS.md`.

### 3.4d Which estimator to use, and why the question has no single answer

A warning is useful only if some estimator avoids the problem, so five were compared on the calibrated
generator: the single-draw presence/absence test used here, the draw-stable rule of § 3.2, a paired
Wilcoxon on proportions, a paired *t*-test on centred log-ratios, and an intercept-only overdispersed
binomial on read counts. Each saw identical simulated counts (`BENCHMARK2_PLAN.md`, frozen beforehand).

**The ranking depends on how the preference is encoded, which settles the question in the negative.**
When a preferring taxon is only more *abundant* in water, the read-count estimator wins by 0.852 in power
and presence/absence detects nothing at all. When it is instead more often *present* in water, that
advantage is **−0.002** — it vanishes. Any benchmark built on one generator recommends whichever
estimator shares its parameterisation; our own first attempt did exactly that and was discarded for it
(`BENCHMARK_RESULTS.md`).

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

### 3.6 The result is not an artefact of our choice of metric

Published niche-breadth work reports the opposite ordering — fungi *narrower* than protists in pond
sediment (Hou et al., 2024) and smallest of four microbial groups in alpine soils (Malard et al., 2022).
Neither study measures what we measure. Malard et al. compute niche breadth along topoclimatic and edaphic
gradients, not between paired compartments; Hou et al. report habitat niche breadth following Xiong et al.
and use Levins' index only to split generalists from specialists, and their system is an indoor glass-tank
experiment with pond-sourced sediment rather than a field survey. The first question is therefore whether
a breadth measure and a paired-preference measure order taxa the same way on our own data. **They do**: Spearman between Levins' *B* and the preference index is **−0.590** (p = 1.7 × 10⁻¹⁰)
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
(51.9 %, 95 % CI 34.0–69.3 %, p = 1.0 — an interval wide enough to contain both chance and moderate
agreement, so the finding is a failure to demonstrate replication rather than a demonstration of its
absence) even though both datasets independently detect habitat partitioning at similar rates (49/100 and
37/103). Those two facts are not in tension once the threshold is taken into account — the first is a
statement about named taxa, the second about the marginal distribution, and the mechanism says the second
can hold while the first does not. A per-order preference table from a single survey describes that survey.

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
coarser aggregation weakens apparent preference (§ 3.5), such a test would have confirmed our hypothesis
for the wrong reason. We report the rank-parity audit (fungi 0.96× non-fungal in FunAqua at order level)
as a check that should accompany any between-group comparison.

### 4.3 Single-draw rarefaction is not a reporting unit

The count of significant orders moved between 43 and 57 on the same data under nothing but a different
random draw, and our own two implementations landed at 49 and 51. Where rarefaction is used — and for
paired designs with unequal depth it is the cleanest available fix — the result should be reported over
many draws, and the reporting unit should be the set of taxa that are stable across them. We found 35 of
100 orders stable at the 95 % level.

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

**Non-independence.** Orders are not independent observations. The checks address power and aggregation,
not phylogenetic non-independence.

**Thin groups.** Rhizaria (7 orders) and Cryptista (24 OTUs) positions should not be quoted.

### 4.5 Why the conservatism result is taxonomic and not phylogenetic

Three routes to a resolved fungal phylogeny were attempted and all are closed. ITSx (Bengtsson-Palme
et al., 2013) shows the amplicon carries a median of only ~61 bp of LSU, so neither the authors'
representative sequences nor the raw reads (downloaded, 5.7 GB) can resolve deep fungal relationships;
deep fungal phylogeny needs D1/D2, ~600 bp. Mapping orders onto the Open Tree of Life synthetic tree
(synthesis `opentree16.1`, taxonomy 3.7draft3; Hinchliff et al., 2015; Redelings and Holder, 2017)
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
in each quadrant. The dashed line is 1:1. Agreement is 14 of 27 (51.9 %, 95 % CI 34.0–69.3 %, two-sided
exact binomial *p* = 1.0); Spearman +0.28 (*p* = 0.16). FunAqua 239 pairs, NEON 149 pairs across 23 US
streams.

**Figure 3. Animals and plants are the most compartment-bound; fungi are not distinguishable from
Alveolata.** (**a**) Median effect magnitude |preference − 0.5| and (**b**) median Levins' niche breadth
*B* across the two compartments, by kingdom, on the same 204 site pairs at the same per-pair rarefaction
depth. Lower preference and higher *B* both indicate a more generalist distribution. Fungi (blue) are
lowest on preference and third on *B*; the gap against non-fungal eukaryotes as a pooled group is driven
by Metazoa and Straminipila. Rhizaria rests on 7 orders and should not be quoted.

**Figure 4. The testability threshold, not the depth itself, decides what is measured.**
(**a**) Expected preference index against rarefaction depth, from
E[index] = *P*_w(1 − *P*_s)/[*P*_w(1 − *P*_s) + *P*_s(1 − *P*_w)] with *P* = 1 − e^(−*aD*), for three fixed
true abundance ratios; the sediment abundance is 10⁻⁴ throughout. Dotted lines mark the two datasets'
depths. The index a study reports depends on how deeply it sequenced, at fixed biological truth.
(**b**) Median relative abundance of the orders that clear the ≥ 20-discordant-pairs requirement, against
depth, **in a simulated community** spanning 10⁻⁶ to 10⁻¹ drawn log-uniformly. In *this* community the
threshold holds the testable set at a fixed position on the detection curve (median λ 0.22–0.49) while the
abundance of the taxa occupying it falls 454-fold. Real data behave differently, for a reason that turns
out to be instructive: see § 3.4a and Figure 5. (**c**) Significant orders per draw in 200 simulated draws
of a FunAqua-like community (239 pairs, depth 211, observed abundance distribution), against the range
observed in the real data (shaded, median marked).

**Figure 5. The mechanism in real data, and what a model has to contain to reproduce it.**
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
listed in `PLAN_HASHES.md`, and all result tables: `github.com/kgurgen97/aquatic-fungi-habitat`.

---

## 6 References

Every entry below was checked against Crossref or the publisher's own record on 22 October 2026; author
lists, volumes, article numbers and DOIs are taken from those records rather than from memory. Where a
paper's finding is stated in the text, that finding was read in the paper's own abstract or text, not
inferred from its title. Database entries carry the version actually used.

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

Vaulot, D., Sim, C.W.H., Ong, D., Teo, B., Biwer, C., Jamy, M. and Lopes dos Santos, A. (2022) metaPR²: a
database of eukaryotic 18S rRNA metabarcodes with an emphasis on protists. *Molecular Ecology Resources*
**22**, 3188–3201. doi:10.1111/1755-0998.13674  *(version 4.0 used here,
doi:10.5281/zenodo.21106507)*

**Removed from the earlier draft:** Madin et al. (2020), *Scientific Data* **7**, 170 — it belongs to the
companion pseudogene study and is cited nowhere in this manuscript.

---

## Author's notes on what still has to happen before submission

1. **References — done, 2026-10-22.** All 14 entries verified against Crossref or the publisher's record;
   author lists, volumes, article numbers and DOIs come from those records. Four corrections came out of
   the check and are now in the text: (a) the Malfertheiner paper is *Community conservatism is widespread
   across microbial phyla and environments*, not an aquatic-specific paper, and the draft's gloss
   "leaving the fungal question open" was replaced by their own sentence, since they never mention fungi;
   (b) § 3.6 had attributed Levins' *B* to both niche-breadth studies — Malard et al. use
   environmental-gradient breadth, and Hou et al. use Levins only to split generalists from specialists;
   (c) Hou et al. is an indoor glass-tank experiment, not a field survey, and is now described as such;
   (d) Reñé et al. was in the list but cited nowhere, and reports the *opposite* vertical pattern, so it
   is now cited in § 3.1 as the disagreement it is. Madin et al. (2020) was removed as belonging to the
   companion study. The Open Tree synthesis version was never recorded, so the query was re-run: synthesis
   `opentree16.1`, all 77 OTT ids identical to the stored file, 70 placed in the tree — the number the
   text reports. NEON was pulled without pinning a release; the data availability statement now says so
   and shows that no provisional month entered the analysis.
2. **Framing — decided, 2026-10-22: mechanism-led.** The biology-led alternative was rejected on the
   paper's own evidence. Its headline would have been that fungi are generalists among eukaryotes, but
   § 3.5 shows fungi (0.150) are indistinguishable from Alveolata (0.157) and that the pooled gap against
   non-fungal eukaryotes is driven by Metazoa (0.294), Straminipila (0.271) and Viridiplantae (0.227); on
   Levins' *B* pooled against protists, p = 0.163. The title would have been contradicted by the paper's
   own table. § 3.5 also has no external validation available — NEON carries only 16S and ITS, and
   metaPR2 failed rank-comparability — whereas the mechanism is derived under known truth and confirmed in
   two independent datasets. The depth objection to § 3.5 does *not* hold, and that is worth keeping in
   mind if a reviewer raises it: the comparison uses the same 204 pairs at the same per-pair depth, the
   ≥ 20 discordant-pair rule applies identically to both groups, and robustness check (ii) leaves the
   difference intact in 4 of 4 discordant-count bins. The mechanism in § 3.4 reinforces this rather than
   undermining it, since both groups are filtered to the same position on the detection curve.
   Consequently: the title, abstract, introduction and § 4.1 now lead with the mechanism; the mechanism
   section was promoted from the interim label 3.3b to § 3.4, shifting the two sections after it; and the
   phrase "fungi sit at the generalist end of the eukaryotic range" was replaced in the abstract, the
   § 3.5 heading and the Figure 3 caption by what the table actually shows.
2a. **External review, 2026-10-23 — what it changed and what it leaves open.** A reviewer read the draft,
   the results files and the figures. Four criticisms were correct and are now fixed in the text:
   (a) § 3.5's heading read "What does transfer", which § 4.4 contradicts — no external validation of that
   section was possible, so nothing showed it transfers; it is now "A group-level contrast within
   FunAqua". (b) § 3.3's heading asserted non-transfer as established, and § 4.1 called the agreement
   "at chance"; with 14 of 27 and a 95 % CI of 34.0–69.3 %, the interval contains both chance and moderate
   agreement, so both now say replication was *not demonstrated* rather than *absent*. (c) The Jaccard
   0.51 turnover between the two datasets is **confounded** — they differ in geography, water-body type
   and in what "sediment" means, so orders testable in one and not the other may be ecologically absent
   rather than filtered out by depth. § 3.4 now says the simulation shows the mechanism is *sufficient* to
   produce such turnover, not that it produced this one. (d) The λ ≈ 1 selection was derived at a single
   effect size (3×) and one abundance distribution; a taxon with a large effect accumulates discordant
   pairs at λ ≫ 1, so the statement holds in the small-to-moderate effect regime, which § 3.4 now says.
   The abstract was rewritten to 241 words, inside Molecular Ecology Resources' 250-word cap.

   Still open, in the order I would do them: (i) the controlled within-NEON depth experiment —
   `rarefied_depth` runs 232 to 425 195, so a ladder of 100/300/1 000/3 000/10 000 keeps 102–149 pairs and
   holds ecosystem, samples and reference database fixed while varying only depth. This is the test (c)
   requires and the one that would turn the mechanism from sufficient to demonstrated. (ii) the parameter
   sweep for (d). (iii) a benchmark of alternative estimators — single and repeated rarefaction,
   abundance-based paired models, hierarchical occupancy, beta-binomial — on simulated truth, comparing
   bias, RMSE, FDR, power and coverage. (iv) sensitivity to site clustering (149 NEON pairs come from
   23 sites) and to taxonomic non-independence among orders, neither of which the agreement test accounts
   for. Also: the reviewer's literature check and ours agree that depth and rarefaction effects on
   diversity metrics are known (Shirazi et al., 2021, *Ecol Evol* 11(22)), and that the closest result to
   ours is García-Machado et al. (2023, *Environmental DNA*), who found depth strongly affected occupancy
   *before* thresholding and not after. They read that as the threshold curing depth sensitivity; we read
   the same phenomenon as the threshold concealing it by changing which taxa are measured. That contrast
   is the paper's clearest novelty claim and should be stated explicitly in the introduction.
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
   | 5 | the algebra and the idealised simulation, three panels | `fig4_simulation` |
   | 6 | the mechanism in real data: ladder, sweep, four abundance variants | `fig5_mechanism_in_data` |

3a. **Why § 3.4a–3.4c keep letters when § 3.3b did not.** These look like the interim label that was
   removed earlier, and they are not the same case. "3.3b" was an independent section wedged after § 3.3
   while it waited for a number. § 3.4a–3.4c are three parts of the single argument § 3.4 opens — the
   algebra, then the same question in real data, then its generality and the model that reproduces it —
   which is an ordinary convention. They are also referenced by `SWEEP_PLAN.md`, `PREVALENCE_PLAN.md` and
   `BENCHMARK_PLAN.md`, which are frozen pre-registrations and must not be edited; renumbering the
   manuscript would leave those references pointing at sections that no longer exist and weaken the
   provenance trail they exist to provide. A journal copy-editor may renumber at proof stage, which is
   the right moment for it.

   Still to decide: whether Fig 4 belongs in the main text or the supplement, given that it is post hoc
   and reaches significance at no cap.
4. **Abstract word count** is 248; check the target journal's limit.
5. **Statement on the corrected numbers.** § 3.2 reports that our own two implementations disagreed. Decide
   whether that belongs in the Results as written, or in a methods note. It is honest either way and a
   reviewer will respect it, but it does draw attention.
