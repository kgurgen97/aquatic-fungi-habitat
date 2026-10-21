# The flagellum hypothesis — result

2026-10-21. Analysis as fixed in `ZOOSPORIC_PLAN.md`
(sha256 `e02627db39a9934cc2100f360ec534ba912c302c79864edf1b19b64fd1c8090f`, committed `a12b677`
before any group comparison was computed). Run once.

## Primary endpoint: NOT REPLICATED

The plan required the predicted direction with p < 0.05 **in both datasets independently**.

| dataset | flagellated | aflagellate | Mann–Whitney p | rank-biserial | verdict |
|---|---|---|---|---|---|
| FunAqua (water vs sediment) | n = 22, median **0.439** | n = 75, median **0.548** | 0.179 | −0.179 | NOT SUPPORTED, direction reversed |
| NEON (water vs benthic) | n = 18, median **0.858** | n = 83, median **0.379** | **5.1e−08** | **+0.823** | SUPPORTED |

NEON supports the hypothesis emphatically — 16 of 18 flagellated orders prefer the water column against
18 of 83 aflagellate ones (Fisher OR 28.9, p = 1.3e−07), permutation p = 1.0e−04. FunAqua does not
support it at all, and what little direction it has runs the other way.

**One dataset out of two is not a replication.** The pre-registered answer is no.

## Why they disagree — the two flagellated phyla behave oppositely

Post hoc (`exploratory_zoosporic_byphylum.py`); this explains the split but does not change the answer
above. `ZOOSPORIC_PLAN.md` had already flagged "phylum is a coarse proxy for the trait" as a limitation,
and this is that limitation biting.

Median water preference by phylum:

| phylum | trait as assigned | FunAqua | NEON |
|---|---|---|---|
| Chytridiomycota | flagellated | 0.590 (n = 11) | **0.893** (n = 10) |
| Rozellomycota | flagellated | **0.309** (n = 9) | 0.519 (n = 5) |
| Blastocladiomycota | flagellated | 0.455 (n = 1) | 0.833 (n = 1) |
| Monoblepharomycota | flagellated | 0.360 (n = 1) | 0.755 (n = 2) |
| Ascomycota | aflagellate | 0.474 (n = 40) | 0.384 (n = 46) |
| Basidiomycota | aflagellate | 0.648 (n = 30) | 0.378 (n = 33) |

The one thing that **does** replicate is the contrast *within* the flagellated group:
Chytridiomycota prefer water more than Rozellomycota in **both** datasets —
FunAqua p = 0.012, NEON p = 0.0013. Pooling them as one trait group cancels the signal in FunAqua, where
Rozellomycota are strongly sediment-associated (0.309 vs the aflagellate 0.548, p = 0.0035).

### Chytridiomycota alone still does not replicate

The obvious move is to narrow the trait to chytrids. It does not rescue the result:

| dataset | chytrids | aflagellate | p | rank-biserial |
|---|---|---|---|---|
| FunAqua | n = 11, median 0.590 | n = 75, median 0.548 | 0.269 | +0.208 |
| NEON | n = 10, median 0.893 | n = 83, median 0.379 | **6.0e−07** | **+0.971** |

In NEON the separation is almost perfect. In FunAqua chytrids sit at 0.590 against an aflagellate 0.548 —
a difference of 0.04. That is not a power problem hiding a real effect; the effect itself is small in
FunAqua and enormous in NEON. The two datasets genuinely disagree about chytrids.

Narrowing the trait after seeing the result would be fitting the hypothesis to the data, and it does not
work anyway. It is recorded here as a post-hoc observation only.

### A problem with the trait assignment itself

| dataset | phylum | orders | of which are undescribed environmental clade labels |
|---|---|---|---|
| FunAqua | Chytridiomycota | 11 | 0 (0 %) |
| FunAqua | Rozellomycota | 9 | 8 (89 %) |
| NEON | Chytridiomycota | 10 | 3 (30 %) |
| NEON | Rozellomycota | 5 | 5 (100 %) |

Almost every Rozellomycota "order" in these datasets is an environmental clade label — GS05, GS07, GS10,
GS11, Branch02, Branch03 — not a described taxon. The flagellum was assigned to Rozellomycota on the
strength of *Rozella* itself, but nobody has observed the dispersal stage of these clades. For most of
the phylum the trait assignment is an extrapolation from one described genus, not a fact, and it should
not have been treated as equivalent to the Chytridiomycota assignment.

## What this means

Three things are now established across two independent datasets, and one is retracted.

1. **Habitat partitioning of fungal orders is real in both datasets.** Unchanged.
2. **Chytridiomycota prefer the water column more than Rozellomycota does.** Replicated, p = 0.012 and
   p = 0.0013. This is a contrast between two groups, not a trait claim, and it is the strongest
   cross-dataset signal in the project.
3. **In NEON, flagellated fungi separate from aflagellate fungi almost perfectly** (rank-biserial 0.823
   overall, 0.971 for chytrids). A striking result about NEON's streams, not yet a general one.
4. **Retracted: the flagellum hypothesis is not supported as a general claim.** It was promoted earlier
   on the strength of two orders agreeing between datasets, which was too small a basis. Tested properly
   on 22 and 18 orders it fails in one of the two datasets. The earlier framing was wrong and is
   withdrawn.

The honest next step is a third dataset with paired water/benthic sampling and comparable depth, with the
trait restricted to phyla whose zoospores have actually been described, pre-registered before it is run.
Testing it again on FunAqua or NEON would not be a test.

Also worth confronting rather than explaining away: Reñé et al. 2023 (Gulf of Naples) report *higher*
chytrid diversity and relative abundance in sediment than in the water column. That is a different
measurement at one site, but it points the same way as FunAqua and against NEON.

## Files

`ZOOSPORIC_PLAN.md` (frozen) · `analysis_zoosporic.py` · `results/zoosporic.tsv` ·
`results/zoosporic.json` · post hoc: `exploratory_zoosporic_byphylum.py` →
`results/zoosporic_byphylum.json` · `results/order_phylum_map.json` (order → phylum in both datasets;
216 orders shared, zero phylum disagreements).
