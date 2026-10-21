# NEON external validation — result

2026-10-21. Analysis as fixed in `NEON_VALIDATION_PLAN.md`
(sha256 `85fd8aa2c9f8a14dc30c9bded5858d4f39dd0d6801f01753531fda0ef75fb939`, committed `bbc9892`
before the data finished downloading; implementation committed `69ce33d`). Run once, no tuning.

## Primary endpoint: NOT REPLICATED

Of the 49 orders significant in FunAqua's pre-registered primary, 27 were testable in NEON.

> **Sign agreement: 14 / 27 = 51.9 %**, 95 % CI [34.0 %, 69.3 %], two-sided exact binomial **p = 1.0**.

The plan declared replication only if agreement > 50 % **and** p < 0.05. Agreement is at chance. By the
pre-committed reading this is **not replicated**, and that is the headline.

Per the plan, FunAqua is not revised on the strength of this. The water-vs-sediment partitioning stands as
a finding in FunAqua; what this test says is that it did **not** carry over to NEON's streams.

## What NEON does show

NEON is not a null dataset. On its own it reproduces the *phenomenon* strongly:

| | FunAqua | NEON |
|---|---|---|
| qualifying pairs | 239 | 149 |
| orders tested (≥ 20 discordant pairs) | 100 | 103 |
| orders with a significant preference (BH q < 0.05) | 49 | **37** |
| split | 30 water / 19 sediment | 18 water / 19 benthic |

So habitat partitioning of fungal orders is real in both datasets. What fails to transfer is **which order
prefers which habitat**.

Secondary 1 — the preference index correlates weakly and non-significantly:
**Spearman +0.281 (p = 0.156)** over the 27 compared orders, **+0.186 (p = 0.126)** over all 69 orders
testable in both. Positive in sign, but indistinguishable from no relationship.

## Secondary 3 (pre-specified): the zoosporic signal *does* replicate

Section 7.3 named these lineages in advance. This is not a post-hoc rescue; it was written down before the
result was known.

| lineage | FunAqua water fraction | NEON | NEON q | verdict |
|---|---|---|---|---|
| Zygophlyctidales | 0.870 | **0.875** | 1.8e−04 | agrees, significant |
| Chytridiales | 0.752 | **0.933** | 2.3e−11 | agrees, significant |
| Zygorhizidiales | 0.984 | — | — | not testable in NEON |
| Rozellales | 0.725 | — | — | not testable in NEON |
| Microascales | 0.150 | 0.323 | 0.022 | agrees, significant |
| Sordariales | 0.168 | 0.491 | 1.0 | no preference in NEON |
| Zoopagales | 0.198 | — | — | not testable in NEON |

Both testable zoosporic orders replicate, with the two strongest q-values in the whole NEON table. With
n = 2 there is no formal test, so this is a consistent observation, not a demonstrated effect — but the
mechanistically specific part of FunAqua's claim (flagellate dispersal stage → water column) is the part
that survived, while the broad order-level pattern did not.

Two orders flip significantly between datasets and should be reported as such, not buried:
**Filobasidiales** 0.646 → 0.167 and **Amphisphaeriales** 0.676 → 0.292 (both water in FunAqua, benthic in
NEON); **GS11** 0.277 → 0.721 (sediment in FunAqua, water in NEON).

## Secondary 4 (pre-specified): sensitivity

Agreement is flat at chance across every threshold — this is not a knife-edge result.

| min reads | pairs | tested | significant | agreement |
|---|---|---|---|---|
| 50 | 150 | 103 | 35 | 14/27 = 51.9 % |
| 100 (primary) | 149 | 103 | 37 | 14/27 = 51.9 % |
| 200 | 149 | 103 | 37 | 14/27 = 51.9 % |
| 400 | 146 | 103 | 39 | 12/27 = 44.4 % |

| orderConfidence | pairs | tested | significant | agreement |
|---|---|---|---|---|
| 0.50 | 150 | 127 | 46 | 14/30 = 46.7 % |
| 0.80 (primary) | 149 | 103 | 37 | 14/27 = 51.9 % |
| 0.90 | 149 | 89 | 27 | 11/25 = 44.0 % |

## Which pairing rule fired, and why it matters

Same-site **same-date** pairing gave only **20** qualifying pairs — below the trigger of 60 fixed in plan
section 4 — so the primary became the declared **±7-day** fallback: **149 pairs**, 23 sites, 128 dates.
Offsets: 20 pairs at 0 days, 48 at 1–2, 81 at 3–7. NEON samples water and benthos in the same field bout
but not always on the same day.

This is the rule the plan said to use in exactly this situation. It was not chosen after seeing an outcome.

## Why the replication may have failed — read as hypotheses, not findings

**1. Sequencing depth differs by ~200×, and that changes what "present" means.**
FunAqua's pairs carry a median of 131–179 reads per sample; NEON's within-pair rarefaction depth is a
median of **31,131** reads (IQR 2,672–93,876). Within-pair rarefaction equalises detection power *inside*
each pair — which is what the design needs — but it cannot make presence/absence mean the same thing
across the two datasets. At ~130 reads, "present" means "among the most abundant sequences"; at 31,000 it
means "detectable at all". The discordant pairs that drive the test are therefore a different slice of the
community in each dataset.

`exploratory_neon_depth.py` (**post hoc, not in the plan**) probes this by capping NEON's pairs at
FunAqua-like depths:

| NEON depth cap | orders tested | agreement, pre-registered set | p | agreement, draw-stable set |
|---|---|---|---|---|
| 130 | 50 | 14/23 = 60.9 % | 0.405 | 11/16 = 68.8 % |
| 300 | 65 | 18/26 = 69.2 % | 0.0755 | 13/19 = 68.4 % |
| 1 000 | 79 | 16/26 = 61.5 % | 0.327 | 12/21 = 57.1 % |
| 3 000 | 85 | 15/26 = 57.7 % | 0.557 | 12/21 = 57.1 % |
| 10 000 | 93 | 15/27 = 55.6 % | 0.701 | 13/22 = 59.1 % |
| none (31 131) | 103 | 14/27 = 51.9 % | 1.0 | 13/22 = 59.1 % |

Agreement rises monotonically as NEON is brought toward FunAqua's depth — 51.9 % → 69 % — which is the
direction the hypothesis predicts. **But it never reaches significance** (best p = 0.0755), and capping
shrinks the comparison to 16–26 orders, so the probe is underpowered to confirm its own suggestion. It
makes depth mismatch a credible explanation worth designing for; it does not establish it, and it does not
change the pre-registered answer.

**2. The second habitat is not the same habitat.** Stated in plan section 8 before the result: FunAqua's
contrast is water vs **sediment**; NEON's is water vs **benthic biofilm** (epipelon, epilithon,
epipsammon, epixylon, epiphyton). A sediment-dwelling order and a biofilm-dwelling order need not be the
same order.

**3. Geography and water body.** NEON here is 23 US wadeable streams. No lakes, one continent. FunAqua's
239 pairs are dominated by Europe (152) with North America second (37). The comparison is not
like-for-like in either biome or region.

**4. FunAqua's own comparison set is soft.** See below — the count of significant orders in FunAqua moves
between 43 and 57 depending only on the rarefaction draw. Repeating the endpoint against FunAqua's
draw-stable set instead of its single-draw set raises agreement from 51.9 % to 59.1 % (13/22, p = 0.52) —
still not significant.

## A problem found in FunAqua's own results while running this

`analysis_primary.py` and `analysis_secondary.py` implement the same declared primary endpoint at the same
read threshold with the same seed, and disagree: **100 vs 99 orders tested, 49 vs 51 significant**. Neither
is miscoded. `analysis_secondary.py` consumes one module-level RNG across its four thresholds, so by the
≥ 100 row its rarefaction draw differs from the primary's fresh RNG. `RESULTS.md` had quoted the
*secondary* script's row as the primary result.

`analysis_draw_stability.py` (**post hoc**) repeats the primary over 200 draws:

| quantity | median | range |
|---|---|---|
| orders tested | 100 | 96 – 103 |
| orders significant (BH q < 0.05) | 51 | **43 – 57** (mean 50.6, sd 2.5) |

Both stored runs sit inside that range, so neither was wrong — **the declared single-draw design is
draw-dependent at roughly ±5 orders.** Per-order:

- **35 orders significant in ≥ 95 % of draws** (21 water, 14 sediment); 28 in every single draw.
- 37 orders significant in some draws but under 95 %.
- Direction is never ambiguous within the stable set.

**The defensible claim for the paper is 35 orders with a stable habitat preference, not 51.**
`RESULTS.md` has been corrected accordingly.

One consequence worth noting: the three orders `RESULTS.md` reported as disagreeing between Europe and
North America — **Auriculariales, Pezizales, Trechisporales** — are all in the unstable set
(27/200, 136/200, 185/200 draws significant). That reported disagreement was in substantial part draw
noise, not biogeography.

## What this means for the paper

The honest framing, and the one the evidence supports:

1. Fungal orders partition between the water column and the benthos. This is now shown in **two
   independent datasets** — FunAqua (49 of 100 orders, single draw; 35 stable) and NEON (37 of 103) —
   with different operators, primers, pipelines and substrate definitions. That is the robust claim.
2. The **identity** of which order prefers which habitat does not transfer between the two. Reporting a
   per-order preference table as a general result is not supportable on this evidence.
3. The **zoosporic lineages do transfer**, pre-specified and with the strongest statistics in NEON. That
   is the mechanistic hypothesis worth taking forward: a motile flagellate stage predicts water-column
   preference.
4. Depth mismatch is a credible but unproven explanation for (2), and it is a design requirement for any
   future comparison: datasets must be compared at matched read depth, not merely rarefied within pairs.

## Files

`NEON_VALIDATION_PLAN.md` (frozen) · `analysis_neon_validation.py` · `results/neon_order.tsv` ·
`results/neon_pairs.tsv` · `results/neon_vs_funaqua.tsv` · `results/neon_validation.json` ·
post hoc: `analysis_draw_stability.py` → `results/draw_stability.{tsv,json}`,
`exploratory_neon_depth.py` → `results/neon_depth_probe.{tsv,json}`.
Data: `~/research/datasets/neon2/` (324 water + 2 131 benthic ITS tables, 23 sites, 298 frozen
site-months, 0 download failures).
