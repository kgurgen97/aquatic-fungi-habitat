# Supplementary information

Companion to `MANUSCRIPT.md`. Every table here is a copy of the authoritative file in `results/` under
its S-number, with provenance recorded in `supplementary/table_provenance.tsv`; nothing is recomputed for
a table, so a supplementary table cannot drift from the result it came from. The three supplementary
figures are produced by `make_supplementary.py` in the same style as the main figures.

---

## Supplementary figures

**Figure S1. The preference index and Levins' niche breadth measure the same thing on these data.**
Median effect magnitude |preference − 0.5| against Levins' *B* across the two compartments, per order,
for (**a**) fungi and (**b**) non-fungal eukaryotes, on the 204 matched pairs at matched per-pair
rarefaction depth. Spearman ρ = −0.59 (*p* = 1.7 × 10⁻¹⁰, *n* = 98) and ρ = −0.78 (*p* = 2.2 × 10⁻²⁶,
*n* = 122). This matters because the published niche-breadth literature uses Levins' *B* and reports an
ordering opposite to ours; the two metrics are not measuring different quantities here, so that
explanation is not available.

**Figure S2. The flagellum trait test replicated in NEON and failed in FunAqua.** Water preference index
for orders in flagellated (zoosporic) versus aflagellate fungal phyla, assigned by phylum membership
before any group value was computed, in (**a**) FunAqua and (**b**) NEON. Horizontal bars are group
medians; the dashed line is no preference. NEON supports the predicted direction emphatically
(*p* = 5.1 × 10⁻⁸, rank-biserial +0.823); FunAqua does not, and what direction it has runs the other way
(*p* = 0.179). Under the pre-declared rule — predicted direction with *p* < 0.05 in both datasets
independently — the hypothesis is not replicated.

**Figure S3. Why metaPR2 could not serve as a second external dataset.** Sequences per taxonomic label,
fungi divided by non-fungal eukaryotes, at four ranks in FunAqua (ITS) and metaPR2 (18S, PR2 5.1.0). A
value of 1 means a fungal label aggregates as much diversity as a non-fungal one. FunAqua is at parity
where it matters (0.96 at order, 0.95 at class); metaPR2 assigns fungi 1.5–2.2× coarser labels at every
rank coarse enough for the 44 paired coastal samples available, because Ascomycota and Basidiomycota
collapse into six subphylum labels while chytrids keep real orders. Since coarser aggregation
demonstrably weakens apparent preference (main text § 3.5), a test on those data would have confirmed our
hypothesis for the wrong reason.

---

## Supplementary tables

| | contents | rows | from |
|---|---|---|---|
| **S1** | FunAqua primary result, all 100 tested orders: discordant counts, preference index, Wilson CI, *p*, *q* | 100 | `analysis_primary.py` |
| **S2** | Per-order significance across 200 rarefaction draws, with times tested and mean preference | 107 | `analysis_draw_stability.py` |
| **S3** | Fungal orders on the 204 matched pairs | 98 | `analysis_euk_matched.py` |
| **S4** | Non-fungal eukaryote orders on the same 204 pairs at the same depth | 122 | `analysis_euk_matched.py` |
| **S5** | Per-kingdom metrics on one consistent footing: effect, Levins *B* per OTU and per order | 6 | `analysis_kingdom_consistent.py` |
| **S6** | Levins' *B* and the preference index, per order, both groups | 220 | `analysis_levins.py` |
| **S7** | NEON per-order result, all 103 tested orders | 103 | `analysis_neon_validation.py` |
| **S8** | The 27 orders comparable between FunAqua and NEON, with the agreement flag | 27 | `analysis_neon_validation.py` |
| **S9** | The 149 NEON water/benthic pairs: site, date, day offset, rarefaction depth | 149 | `analysis_neon_validation.py` |
| **S10** | Flagellum trait assignment and preference, both datasets | 203 | `analysis_zoosporic.py` |
| **S11** | Agreement under capping of NEON depth, both comparison sets | 14 | `exploratory_neon_depth.py` |
| **S12** | Per-order values behind the five robustness checks | 221 | `analysis_euk_robustness.py` |
| **S13** | The five pre-registered robustness checks: statistic, declared threshold, outcome | 5 | `analysis_euk_robustness.py` |
| **S14** | Rank granularity audit behind Figure S3 | 8 | `make_supplementary.py` |

---

## Supplementary methods

### S.1 Why the `.001` NEON products were rejected

NEON's superseded surface water (`DP1.20141.001`) and benthic (`DP1.20086.001`) community composition
products assign a taxonomic order to only **5.2 % of reads** — 104 837 of 2 032 229 across 123 files we
measured before abandoning them. Order-level analysis is not possible on that basis. The current
products, `DP1.20141.002` and `DP1.20086.002`, emit a rank name for every sequence regardless of support
(one probed sample: median `orderConfidence` 0.69, minimum 0), so the 0.80 confidence cutoff is a
requirement and not a refinement.

### S.2 The pairing rule that fired, and why

Same-site same-date pairing of the NEON products gave only 20 qualifying pairs, below the trigger of 60
fixed in `NEON_VALIDATION_PLAN.md` before the data were examined, so the pre-declared ±7-day fallback
applied: 149 pairs across 23 sites and 128 dates, with 20 pairs at 0 days offset, 48 at 1–2 days and 81
at 3–7. NEON samples water and benthos in one field bout but not always on one day.

### S.3 Pre-registration

Analysis plans were committed to version control before the corresponding outcomes were computed, so the
order of events is checkable rather than asserted:

| plan | frozen at | result at |
|---|---|---|
| NEON validation | `bbc9892` | `3d130c3` |
| NEON implementation | `69ce33d` (before the download finished) | — |
| Flagellum trait test | `a12b677` | `9dec095` |
| Eukaryote robustness checks | `a5c9fe6` | `ab5a24e` |
| Levins' *B* | `9e01ec4` (added before those outcomes were seen) | — |
| metaPR2 replication | `eba1edf` | `b298db6` (not executed; see § S.4) |

Each plan's SHA-256 is recorded in `PLAN_HASHES.md`, alongside the commit that introduced it, and can be
checked with `sha256sum *PLAN*.md`. The hashes are deliberately kept outside the plan files: writing a
file's own hash into the file changes the file, so the recorded value could never match, and a frozen plan
must not be edited for any reason. `verify_plan_hashes.py` asserts that every hash quoted anywhere in the
repository matches the file on disk; all five plans pass.

### S.4 The metaPR2 assessment

metaPR2 v4.0 (Vaulot et al., 2022; Zenodo 10.5281/zenodo.21106507, 187 datasets, 18S rRNA, annotated
against PR² 5.1.0, Guillou et al., 2013) was obtained and assessed as a second external dataset, attractive because 18S carries the opposite primer bias to ITS and
would therefore have tested whether the result is an artefact of the assay. It was not used, for reasons
recorded in `METAPR2_RESULTS.md` and summarised here:

- The pre-registered primary scope, paired freshwater water and sediment, has **zero** pairs. metaPR2
  holds 870 freshwater water-or-sediment samples of which only 118 are sediment, and none of its 9
  datasets carrying both substrates is freshwater.
- The declared secondary scope, all aquatic ecosystems, gives 54 paired coastal sites, 44 after
  restricting to one amplicon region.
- Decisively, PR2 assigns fungi systematically coarser ranks than protists (Figure S3). Running the test
  would have produced a spurious confirmation.

### S.5 A correction made during this work

`analysis_primary.py` and `analysis_secondary.py` implement the same declared primary endpoint at the
same read threshold with the same seed and returned 100/49 and 99/51 respectively. Neither is miscoded:
`analysis_secondary.py` consumes one module-level random generator across its four thresholds, so by the
≥ 100 row its rarefaction draw differs from the primary's fresh generator. An earlier version of
`RESULTS.md` quoted the secondary script's row as the primary result. The pre-registered primary is
`analysis_primary.py`, which reproduces `results/primary_order.tsv` bit-identically on rerun. Both values
fall inside the 43–57 range measured over 200 draws (Figure 1, Table S2), which is the substantive point:
the single-draw design is unstable at roughly ±5 orders, and the draw-stable set of 35 is the defensible
unit of reporting.

A consequence worth recording: three orders previously reported as disagreeing in direction between
Europe and North America — Auriculariales, Pezizales, Trechisporales — are significant in only 27, 136
and 185 of 200 draws, so that reported disagreement was substantially draw noise.

---

## Data and code

All analysis code, every frozen plan with its SHA-256, and all result tables:
`github.com/kgurgen97/aquatic-fungi-habitat`. Supplementary tables and figures are regenerated by
`make_supplementary.py`; main figures by `make_figures.py`.


**Figure S4. Agreement with FunAqua under capping of NEON's sequencing depth.** Direction agreement for
the 49 FunAqua-significant orders (blue) and for the 35 draw-stable orders (orange), as NEON's within-pair
rarefaction depth is capped at successively lower values; "none" is the uncapped median of 31 131 reads.
FunAqua's own pairs carry a median of 131–179 reads per sample, the regime at the left-hand end. Agreement
rises toward matched depth but reaches significance at no cap (best *p* = 0.076). Post hoc; not part of
the frozen validation plan.

Moved here from the main text: this probe is post hoc, outside the frozen validation plan, reaches significance at no cap (best *p* = 0.076), and was never cited in the body. It is kept because it is the only direct look at depth matching between the two datasets, and because § 3.4a now asks the same question inside NEON, where the comparison is clean.