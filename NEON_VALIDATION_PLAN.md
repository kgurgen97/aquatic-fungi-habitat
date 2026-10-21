# External validation of the FunAqua habitat-partitioning result in NEON

**Written 2026-10-21, before any NEON water-vs-benthic outcome was computed.** Nothing below may be
changed because of a result. If a decision here turns out to be wrong, the fix is a new, separately
labelled analysis — not an edit to this file.

## 1. Why

`RESULTS.md` reports that fungal orders partition between the water column and the sediment in FunAqua:
239 paired sites, 99 orders testable, 51 with a significant preference, and a biologically coherent split
(zoosporic lineages in the water, terrestrial saprotrophs in the sediment). Everything supporting that
claim so far comes from **one dataset processed by one pipeline**. Rivers-vs-lakes and
Europe-vs-North-America agreement are internal splits of that same dataset; they cannot rule out a
pipeline artefact common to all of it.

This plan tests the claim in data that shares nothing with FunAqua: different sampling programme,
different operators, different primers, different bioinformatics, different substrate definition.

## 2. Data (frozen)

| | |
|---|---|
| water | NEON **DP1.20141.002** "Surface water microbe community taxonomy" |
| benthic | NEON **DP1.20086.002** "Benthic microbe community taxonomy" |
| files | the per-run ITS tables, `*_ITS*_map.csv` |
| universe | the **298 site-months where both products have data**: 23 wadeable-stream sites, 2018-03 … 2024-10 |
| sites | ARIK BIGC BLDE BLUE CARI COMO CUPE GUIL HOPB KING LECO LEWI MART MCDI MCRA OKSR POSE PRIN REDB SYCA TECR WALK WLOU |
| frozen list | `~/research/datasets/neon2/universe.json` |

The superseded `.001` products are **not** used: they assign only ~5 % of reads to an order
(measured: 104 837 / 2 032 229 reads over 123 files), which cannot support an order-level test.

## 3. Row inclusion

1. `kingdom == "Fungi"`.
2. `orderConfidence >= 0.80`.

Rank 2 is not optional. The `.002` products emit a name at every rank for every sequence, but the
confidence attached to it ranges down to 0 (one probed sample: median `orderConfidence` 0.69, 46.4 % of
rows ≥ 0.80). Without a cutoff the order column is partly noise. 0.80 is the long-standing
RDP/DADA2 convention and is chosen here for that reason, not by inspection of any outcome.

## 4. Sample and pair construction

- A **sample** is one `dnaSampleID`. Its depth is the sum of `individualCount` over included rows.
- Water IDs look like `LEWI.SS.20190508.DNA-DNA1`; benthic IDs like `LEWI.20241022.EPIPELON.4.DNA-DNA1`.
  Site and collection date are taken from the ID.
- Benthic sampling has replicates and several substrate types per date. Where a site-date offers more than
  one candidate on a side, take the **single deepest** sample on that side. Deterministic, and it
  maximises power without reference to any result.
- **Primary pairing: same site and the same calendar date.** Declared fallback, fixed now: if exact-date
  pairing yields fewer than 60 qualifying pairs, the primary becomes same site within **±7 days**, and the
  exact-date analysis is reported alongside it as the secondary. Whichever is primary is stated with the
  pair count, so the reader can see which rule fired.
- A pair qualifies when **both** samples have ≥ 100 included reads.

## 5. Measurement — identical to the FunAqua primary

Copied from `ANALYSIS_PLAN.md` / `analysis_primary.py`, unchanged:

- Within each pair, rarefy both samples to that pair's minimum depth
  (`numpy.random.Generator.multivariate_hypergeometric`, seed **20261021**). Detection power is then
  equal within every pair by construction, so a difference cannot be sequencing effort.
- Reduce to **presence/absence per Order** after rarefaction.
- Per order, count water-only and benthic-only discordant pairs. Test orders with ≥ **20** discordant
  pairs by two-sided exact binomial against 0.5. Benjamini–Hochberg across tested orders, significant at
  **q < 0.05**.
- `water_fraction` = water-only / discordant — the same preference index as in FunAqua.

## 6. Primary validation endpoint — one number, declared now

Take the orders that are **significant in FunAqua** (`results/primary_order.tsv`, q < 0.05) **and**
testable in NEON (≥ 20 discordant pairs). For each, ask whether NEON's `water_fraction` falls on the same
side of 0.5 as FunAqua's.

> **Endpoint: the fraction of those orders that agree in sign.**
> **Test: two-sided exact binomial against 0.5.**
> **Replication is declared if agreement > 50 % and p < 0.05.**

Pre-committed readings:

| outcome | what it means |
|---|---|
| agreement > 50 %, p < 0.05 | the FunAqua result replicates in independent data |
| agreement not significant | not replicated; reported as such, with the achieved N and power |
| agreement < 50 % | contradicted; reported as such |

Any of the three gets written up. The FunAqua result is not revised on the strength of this test — it
either carries external support or it does not.

## 7. Secondary endpoints (declared, not primary)

1. Spearman correlation of `water_fraction` between the two datasets, over orders testable in both.
2. NEON's own independent count: orders tested, orders significant, and the water/benthic split.
3. Direction check on the zoosporic lineages named in `RESULTS.md` (Zygophlyctidales, Chytridiales,
   Zygorhizidiales, Rozellales) and on the terrestrial saprotrophs (Sordariales, Microascales,
   Zoopagales), where testable.
4. Sensitivity to the read threshold at 50 / 100 / 200 / 400, and to `orderConfidence` at 0.5 / 0.8 / 0.9.
   Sensitivity analyses describe robustness; they never replace the primary.

## 8. Known differences from FunAqua — stated before the result, not after

These are real and are not to be presented later as if they had been anticipated only when convenient:

- **Substrate.** FunAqua's second habitat is sediment. NEON's is benthic biofilm — epipelon, epilithon,
  epipsammon, epixylon, epiphyton. Related, not identical. A replication across this difference is
  stronger evidence than one within a single substrate definition; a failure is correspondingly weaker
  evidence against, and will be reported that way.
- **Geography.** 23 US wadeable streams. No lakes, one continent. FunAqua's lake and river agreement
  cannot be re-tested here.
- **Pipeline.** NEON's classifier and FunAqua's are different, as are the primers and the amplicon.
  This is the point of the exercise; it also means order names must be matched as plain strings and
  non-overlapping orders simply drop out.
- **Depth.** NEON ITS libraries are deeper than FunAqua's. Within-pair rarefaction handles this inside
  each pair; it does not make the two datasets' absolute richness comparable, and no such comparison
  will be made.

## 9. Deliverables

`analysis_neon_validation.py` (one run, no tuning), `results/neon_order.tsv`,
`results/neon_validation.json`, and a `NEON_RESULTS.md` reporting the primary endpoint first,
whatever it says.
