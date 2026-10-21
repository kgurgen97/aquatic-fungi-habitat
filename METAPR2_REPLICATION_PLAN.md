# External replication of the generalist result in metaPR2 (18S)

**Written 2026-10-21, before any metaPR2 water-vs-sediment outcome was computed**, and before the
download finished. Sample counts in section 2 are filled in from metadata only — how many paired sites
exist is a feasibility fact, not an outcome. Nothing else here may be changed because of a result.

## 1. Why this dataset, and why it is the right test

`EUK_RESULTS.md` reports that fungi sit at the generalist end of the eukaryotic range across the
water–sediment boundary, and names the limitation that will decide the paper's fate: **FunAqua's amplicon
is ITS, a fungal barcode.** Non-fungal eukaryotes there are off-target, with thinner reference coverage
and taxon-dependent recovery. Every check in `EUK_ROBUSTNESS_PLAN.md` passed, but none of them can fix
that, because the bias is in the assay rather than in the analysis.

metaPR2 fixes exactly that, because **its bias runs the other way**. It is built from 18S rRNA studies —
a eukaryote-universal marker curated with an emphasis on protists, taxonomy from PR2. Under 18S the
protists are the on-target group and **fungi are the comparatively under-served ones**: 18S resolves fungi
poorly relative to ITS, and the reference situation is reversed.

So the two datasets have opposite primer bias with respect to the very contrast being tested. If fungi
come out at the generalist end under a marker that disfavours them, the finding cannot be an artefact of
primer choice. If they do not, the FunAqua result is most likely an ITS artefact and must be reported as
such. Either way this is decisive in a way no further internal check could be.

Second independent difference: metaPR2 is a compilation of many separate studies with different operators,
depths and sampling protocols, whereas FunAqua is one coordinated survey. Agreement across that as well
would be worth more than agreement within one survey.

| | FunAqua | metaPR2 |
|---|---|---|
| marker | ITS (fungal) | 18S rRNA (eukaryote-universal) |
| favours | fungi | protists |
| taxonomy | UNITE-style, fungal | PR2 5.1.0 |
| structure | one coordinated survey | 187 datasets, many studies |
| unit here | Order | Order (PR2 rank) |

## 2. Data (frozen)

metaPR2 version 4.0, downloaded from Zenodo record **21106507**
(`metapr2_4.0_datasets_samples_asv.zip`, `metapr2_4.0_asv_abundance.tsv.gz`).
Filter at source: `gene = 18S rRNA`, `substrate ∈ {water, sediment}`, all gene regions, DNA and RNA,
all size fractions, all depth levels. The app reports 19 306 samples under that filter before pairing.

**Primary ecosystem scope: freshwater (lake + river)**, to match where FunAqua's 204 matched pairs
actually come from (95 river, 116 lake). All aquatic ecosystems together is a declared secondary.

Clustered ASVs are used, as the app serves them by default.

## 3. Pairing

A pair is one **water** sample and one **sediment** sample from the same site. metaPR2 carries
`dataset_id`, latitude, longitude and, where reported, date. Site identity is therefore defined as:

- same `dataset_id`, and
- coordinates agreeing to **3 decimal places** (~100 m), and
- where both samples carry a date, the same date; where a date is missing for either, the pair is still
  allowed and the count of such pairs is reported separately.

Where a site offers more than one candidate on a side, take the **single deepest** sample on that side,
as in the NEON work. A pair qualifies when both samples have ≥ 100 reads after the taxonomic filter.

**Declared fallback, fixed now:** if coordinate-based pairing yields fewer than 40 qualifying pairs,
the primary becomes pairing on `dataset_id` + nearest-coordinate match within **1 km**, and the
100 m analysis is reported alongside it. Whichever rule fires is stated with its pair count.

## 4. Measurement — identical to FunAqua

- Unit: taxonomic **Order** (PR2 `order`). Rows with no named order are dropped.
- Groups: **Fungi** (PR2 division/subdivision Fungi) versus **all other eukaryotes**, excluding
  Metazoa-only and Streptophyta-only considerations — i.e. other eukaryotes means every non-fungal
  eukaryote, exactly as in FunAqua, with the per-supergroup breakdown reported separately.
- Within each pair, rarefy both samples to that pair's minimum depth
  (`multivariate_hypergeometric`, seed **20261021**), then reduce to presence/absence per Order.
- Per order, count water-only and sediment-only discordant pairs; test orders with ≥ **20** discordant
  pairs; `water_fraction` = water-only / discordant; effect = |water_fraction − 0.5|.

## 5. Primary endpoint — one number, declared now

> **Endpoint: the median effect magnitude of fungal orders versus that of non-fungal eukaryote orders.**
> **Test: two-sided Mann–Whitney U.**
> **Replication is declared if fungi are LOWER (more generalist) with p < 0.05.**

Pre-committed readings:

| outcome | what it means |
|---|---|
| fungi lower, p < 0.05 | replicated under opposite primer bias; the claim is about organisms, not the assay |
| not significant | not replicated; report with achieved N and power, and say the FunAqua result stands on one dataset |
| fungi higher, p < 0.05 | contradicted; the FunAqua result is most likely an ITS artefact and must be reported that way |

The FunAqua result is not revised on the strength of this, and this test is run **once**.

## 6. Secondary endpoints (declared)

1. **Draw stability**: 200 rarefaction draws; report how often the ordering reverses. The gate used in
   `EUK_ROBUSTNESS_PLAN.md` was ≤ 5 %.
2. **Per supergroup/division**, as in `results/kingdom_consistent.tsv`: is Metazoa again the most
   compartment-bound and are the protist groups again at the generalist end? A reproduced *ordering*
   would be stronger evidence than the fungi-vs-rest contrast alone.
3. **Matched on power**: the same discordant-count stratification that was decisive in FunAqua.
4. **All aquatic ecosystems**, not just freshwater.
5. **Levins' B** across the two compartments and across all samples, for comparability with
   `results/levins.json` and with the published niche-breadth literature.

## 7. Known limitations, stated before the result

- **18S resolves fungi poorly.** That is the point of using it, but it means fungal orders here will be
  fewer and coarser than in FunAqua. If too few fungal orders clear the 20-discordant bar the test may
  simply be underpowered, and that must be reported as underpowered rather than as a null.
- **metaPR2 is a compilation.** Sequencing depth, region (V4 / V6-V8 / V9), platform and protocol vary
  between its datasets. Within-pair rarefaction handles depth inside a pair; it does not make datasets
  comparable to each other, and no cross-dataset comparison will be made.
- **Sediment coverage may be thin.** metaPR2 emphasises protist plankton, so sediment samples are likely
  a minority. The achieved pair count decides whether this test is informative at all.
- **Different sites from FunAqua.** This is replication of a pattern, not of particular orders.
- Orders are not independent observations.

## 8. Deliverable

`analysis_metapr2.py`, `results/metapr2_*.{tsv,json}`, and a result section reporting the primary first,
whatever it says. If the pair count is too low for the test to be informative, that is the reported
outcome and no test is dressed up as one.
