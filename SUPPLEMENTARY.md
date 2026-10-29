# Supplementary information

Companion to `MANUSCRIPT.md`. Every table here is a copy of the authoritative file in `results/` under
its S-number, with provenance recorded in `supplementary/table_provenance.tsv`; nothing is recomputed for
a table, so a supplementary table cannot drift from the result it came from. Figures S1–S3 are produced by
`make_supplementary.py` and Figure S4 by `analysis_agreement_power.py`, in the same style as the main
figures. Tables S15–S17 are the only tables here not copied from `results/`: they are written directly by
`analysis_agreement_power.py` from the other supplementary tables, which is recorded as their
provenance.

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

**Figure S4. What the agreement endpoint could detect, and a better-powered version of it.**
(**a**) Direction agreement with 95 % intervals for four nested selection sets, from all 65 orders
testable in both datasets down to the nine significant in both. Blue, Wilson intervals; orange, intervals
from bootstrapping phyla rather than orders, so that related orders are resampled together. The vertical
line is chance. Agreement rises as the set is narrowed, which is what conditioning on significance in both
datasets produces mechanically; no selection makes it significant.
(**b**) Power of a two-sided exact binomial test against true direction agreement, for the pre-registered
26-order endpoint and for the 65-order version. Dotted vertical lines mark the agreement each endpoint
would need to reach *p* < 0.05 — 73.1 % and 63.1 %. The pre-registered endpoint has 46 % power against a
true agreement of 70 %; the 65-order version has 91 %. Post hoc, `analysis_agreement_power.py`.

**Figure S3. Why metaPR2 could not serve as a second external dataset.** Sequences per taxonomic label,
fungi divided by non-fungal eukaryotes, at four ranks in FunAqua (ITS) and metaPR2 (18S, PR2 5.1.0). A
value of 1 means a fungal label aggregates as much diversity as a non-fungal one. FunAqua is at parity
where it matters (0.96 at order, 0.95 at class); metaPR2 assigns fungi 1.5–2.2× coarser labels at every
rank coarse enough for the 44 paired coastal samples available, because Ascomycota and Basidiomycota
collapse into six subphylum labels while chytrids keep real orders. The axis carries FunAqua's rank names;
PR2's equivalents are *division*, *subdivision*, *class* and *order*, so the coarsest pair compares
FunAqua's Phylum with PR2's division. Since coarser aggregation
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
| **S15** | Direction agreement across four nested selection sets, with Wilson and phylum-clustered intervals | 4 | `analysis_agreement_power.py` |
| **S16** | Power of the agreement endpoint: minimum detectable agreement, power curve, *n* for 80 % power | 10 | `analysis_agreement_power.py` |
| **S17** | The 65 orders testable in both datasets, both indices, both *q*-values and phylum | 65 | `analysis_agreement_power.py` |
| **S18** | Per-order between-dataset heterogeneity: both count pairs, both indices, Fisher *p* and BH *q* | 69 | `analysis_heterogeneity.py` |
| **S19** | Endpoint recomputed at discordant-pair thresholds 20–60, both datasets | 10 | `analysis_heterogeneity.py` |
| **S20** | Sensitivity of heterogeneity counts to assumed site-level intraclass correlation; exploratory scenarios, not measured ICC | 6 | Supplementary sensitivity table (unreproduced calculation) |

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
| Depth ladder | `464e481` (before the analysis script existed) | `5058008` |
| Parameter sweep | `7a0d398` (before `analysis_sweep.py` existed) | `fa995db` |
| Abundance representation | `4ac0440` (before `analysis_prevalence.py` existed) | `96be558` |
| Estimator benchmark | `828f8da` (before `analysis_benchmark.py` existed) | `6a94646` |
| Estimator benchmark, corrected generator | `6a94646` | `8c6df32` |
| Occupancy model | `be9490c` (before `analysis_occupancy.py` existed) | `6fd306a` |

Each plan's SHA-256 is recorded in `PLAN_HASHES.md`, alongside the commit that introduced it, and can be
checked with `sha256sum *PLAN*.md`. The hashes are deliberately kept outside the plan files: writing a
file's own hash into the file changes the file, so the recorded value could never match, and a frozen plan
must not be edited for any reason. `verify_plan_hashes.py` asserts that every hash quoted anywhere in the
repository matches the file on disk; all **eleven** plans pass. Two frozen plans are known to contain
errors that have deliberately not been corrected, since correcting them would defeat the purpose of
freezing: `BENCHMARK2_PLAN.md` § 3 states "81 cells" where its own axes define 27, and
`METAPR2_REPLICATION_PLAN.md` retains the superseded "generalist end" wording.

Analyses labelled **post hoc** are covered by no plan: `analysis_draw_stability.py`,
`exploratory_neon_depth.py`, `exploratory_zoosporic_byphylum.py`, the log–log fits of main-text § 3.4a,
and `analysis_agreement_power.py` (§ S.6 below).

### S.4 The metaPR2 assessment

metaPR2 v4.0 (Vaulot et al., 2022; Zenodo 10.5281/zenodo.21106507, 187 datasets, 18S rRNA, annotated
against PR² 5.1.0, Guillou et al., 2013) was obtained and assessed as a second external dataset, attractive because 18S carries the opposite primer bias to ITS and
would therefore have tested whether the result is an artefact of the assay. It was not used, for reasons
recorded in the repository (§ 5 of the main text) and summarised here:

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

### S.6 The post hoc agreement analyses

`analysis_agreement_power.py` reads only the supplementary tables and does three things not covered by
`NEON_VALIDATION_PLAN.md`, which is unchanged and remains the pre-registered record.

- **A tie correction.** `results/neon_vs_funaqua.tsv` flagged Coniochaetales as agreeing in direction
  while recording a NEON index of exactly 0.500 and no NEON preference. A tie has no direction, so it is
  excluded: the endpoint becomes 13 of 26 = 50.0 % rather than 14 of 27 = 51.9 %. Four orders are ties
  across the full set testable in both datasets, so 69 shared testable orders give 65 directional
  comparisons.
- **A better-powered version.** The pre-registered endpoint compares only orders significant in FunAqua,
  discarding 39 of the 65 comparable orders and conditioning on significance in a draw-dependent design.
  Agreement across all 65 is 35 = 53.8 %.
- **Power and clustering.** Exact binomial power, the minimum detectable agreement, the *n* required for
  80 % power, and intervals bootstrapped over phyla (20 000 resamples, seed 11) rather than over orders.

Full account, including the limits of using phylum as the clustering unit, in the agreement-power
results report in the repository.

### S.7 Between-dataset heterogeneity and threshold sensitivity

`analysis_heterogeneity.py`, post hoc, reads only the supplementary tables. Full account, including the
limits of the test and what it does and does not attribute, in the heterogeneity results report in the
repository.

For each of the 69 orders testable in both datasets a 2 × 2 Fisher exact test compares (water-only,
sediment-only) in FunAqua with (water-only, benthic-only) in NEON, BH-corrected across the 69. This asks
whether the two datasets give the same index *value*, where the pre-registered endpoint asked only whether
they fall on the same side of 0.5. Twenty-six differ significantly; 16 of the 30 direction-disagreements
are significant, as are 11 of the 13 inside the pre-registered 26-order set, and a further 7 orders agree
in direction while differing significantly in magnitude.

Because the test has more power where orders carry more discordant pairs, and heterogeneous orders carry
marginally more (median minimum discordant count 38 against 33, Mann–Whitney *p* = 0.081), 37.7 % is a
lower bound on disagreement among well-measured orders rather than an unbiased estimate over all of them.

Threshold sensitivity recomputes each dataset's exact binomial test and BH correction at discordant-pair
thresholds of 20, 30, 40, 50 and 60. Thresholds below 20 are unreachable: orders failing the original rule
were never tabulated. At 50 and 60 the NEON comparison falls to 3 and 2 orders and is reported but not
interpreted.

---

### S.8 Assumed site-level correlation sensitivity

Table S20 uses six assumed intraclass correlations (ICC = 0–0.50) and their corresponding design effects. The ICC values are scenarios rather than estimates from the sampled sites. The calculation has not been independently reproduced, so this table is not a correction for site-level dependence.

## Data and code

All analysis code, every frozen plan with its SHA-256, and all result tables:
`github.com/kgurgen97/aquatic-fungi-habitat`. Supplementary tables and figures are regenerated by
`make_supplementary.py`; main figures by `make_figures.py`.