# metaPR2 replication: not executable, and why that is the finding

2026-10-21. Plan frozen in `METAPR2_REPLICATION_PLAN.md`
(sha256 `52d75c02f3091fdbeab40148e13559390f5f742983d35192e84474bca86d77b2`, commit `eba1edf`)
before the data finished downloading. **No test was run.** Section 8 of that plan said: "If the pair
count is too low for the test to be informative, that is the reported outcome and no test is dressed up
as one." This is that outcome, for a reason worse than pair count.

## Verdict

> The pre-registered primary, paired freshwater water/sediment, has **zero** pairs in metaPR2.
> The declared secondary, all aquatic ecosystems, has 54 paired coastal sites, 44 of them on a single
> amplicon region. But PR2's 18S taxonomy applies **systematically coarser ranks to fungi than to
> protists**, in exactly the direction that would manufacture the result being tested. Running it would
> have produced a spurious confirmation.

So metaPR2 is not a valid replication set for this contrast. Recorded here so the next person does not
spend the same day on it.

## What is in metaPR2 (v4.0, Zenodo 21106507, MIT)

187 datasets, 26 315 samples, 229 313 ASVs, taxonomy PR2 5.1.0, gene 18S rRNA throughout
(V4 24 839, V9 1 402, V6-V8 74).

| substrate | samples |
|---|---|
| water | 17 476 |
| soil | 4 598 |
| **sediment** | **1 830** |
| biofilm | 646 |
| everything else | < 600 each |

| ecosystem | samples |
|---|---|
| coastal | 10 889 |
| oceanic | 9 646 |
| terrestrial | 4 794 |
| **freshwater lake** | **760** |
| **freshwater river** | **198** |

Freshwater totals 870 water-or-sediment samples, of which only **118 are sediment**.

## Pairing: 54 sites, none of them freshwater

Only **9 of 138** 18S datasets contain both water and sediment, and not one is freshwater:

| dataset | water | sediment | ecosystem |
|---|---|---|---|
| Australian marine survey | 2 946 | 326 | coastal, oceanic |
| Biomarks project 2009 | 100 | 20 | coastal |
| Singapore plastic study 2021 | 61 | 54 | coastal |
| Beaches British Columbia | 15 | 75 | coastal |
| Indonesia coral reef 2018 | 17 | 24 | coastal |
| Intertidal Xiamen | 21 | 15 | coastal |
| Seaweed communities | 26 | 28 | coastal, oceanic |
| Tsunami deposits SE Asia | 16 | 64 | terrestrial, coastal |
| Arctic ARK-XXVII/3 2012 | 8 | 3 | oceanic |

Same dataset plus coordinates to 3 decimal places (~100 m) gives **54 paired sites, all coastal**;
at ~1 km it gives 52, so the pairing is not sensitive to that threshold and the declared coordinate
fallback was not needed. 53 of 54 sites have a water and a sediment sample on the same date, and the
largest single dataset holds 33 % of the sites, so study independence would have been acceptable.
Gene region splits V4 917 / V9 35 among those samples, and mixing V4 with V9 would compare different
amplicons, so the usable set is ~44 V4 coastal pairs.

That alone would have made the test thin, an order must be discordant in ≥ 20 of 44 pairs, i.e. 45 % of
them, to be testable at all. But it is not the reason the test was abandoned.

## The reason: PR2 does not give fungi and protists the same rank

PR2's `order` column holds 43 distinct labels for all fungi. Of those, only **19 are genuinely
order-rank** names. Twelve are **subphyla**, Agaricomycotina, Pezizomycotina, Saccharomycotina,
Pucciniomycotina, Ustilaginomycotina, Taphrinomycotina, Glomeromycotina, Mucoromycotina,
Blastocladiomycotina, Entomophthoromycotina, Kickxellomycotina, Zoopagomycotina, eight are placeholders
(`Rozellomycota_X`, `Ascomycota_X`, `Fungi_XX`, …), one is a family (Mortierellaceae) and two are classes.

The effect is concentrated where it does most damage. **Ascomycota and Basidiomycota, 263 of FunAqua's
346 named fungal orders, collapse into six subphylum labels plus two placeholders.** Meanwhile the
chytrids keep real orders: Chytridiales, Rhizophydiales, Spizellomycetales, Lobulomycetales,
Zygophlyctidales, Zygorhizidiales.

Granularity, measured as ASVs per label:

| rank | fungal labels | ASVs/label | other labels | ASVs/label | ratio |
|---|---|---|---|---|---|
| division | 1 | 23 987 | 33 | 6 222 | 3.86 |
| subdivision | 1 | 23 987 | 60 | 3 422 | 7.01 |
| class | 13 | 1 845 | 246 | 835 | **2.21** |
| order | 43 | 558 | 564 | 364 | **1.53** |
| family | 96 | 250 | 1 276 | 161 | **1.55** |
| genus | 1 197 | 20 | 6 760 | 30 | 0.66 |

At every rank coarse enough for 44 pairs, a fungal label aggregates 1.5–2.2× more diversity than a
non-fungal one. At the one rank with rough parity (genus, 0.66) the units are far too fine for 44 pairs.
There is no usable rank-matched option. (An automatic search for "closest to parity" returns
*supergroup* at 1.17, but there is exactly **one** fungal supergroup, so that is an artefact of ranking
by ratio alone, not an answer.)

**Why this is disqualifying rather than merely awkward.** `EUK_RESULTS.md` establishes, by measurement,
that coarser aggregation weakens apparent habitat preference: within FunAqua's fungi, orders built from
more OTUs show weaker preference (Spearman −0.246, p = 0.0144). PR2 aggregates fungi more coarsely than
protists. A test on this data would therefore push fungi toward looking generalist **by construction**,
reproducing the claim for the wrong reason. That is worse than no test, because it would read as
independent confirmation.

## A side finding that strengthens the original analysis

The same audit run on FunAqua shows **near-perfect rank parity** between the two groups:

| rank | fungi labels | OTUs/label | other euk labels | OTUs/label | ratio |
|---|---|---|---|---|---|
| Phylum | 37 | 1 355 | 135 | 649 | 2.09 |
| **Class** | 142 | 353 | 236 | 371 | **0.95** |
| **Order** | 346 | 145 | 580 | 151 | **0.96** |
| Family | 814 | 62 | 1 083 | 81 | 0.76 |

At the rank the analysis actually uses, a fungal order aggregates 0.96× as much diversity as a
non-fungal one. This answers directly the reviewer question "were your two groups treated at comparable
taxonomic resolution?"; they were, to within 4 %. It is new evidence in favour of
`EUK_RESULTS.md` and was not available before this check.

## Where external validation could still come from

Honestly, not from anything already assembled:

1. **NEON: ruled out.** No 18S anywhere, not in the aquatic marker-gene products, not in soil, not in
 any of the 203 catalogue entries. Only 16S and ITS. Its metagenomes (302 benthic + 2 568 water
 site-months) do contain SSU reads, so extracting eukaryotic SSU from shotgun data is a genuine
 possibility, but that is a project on terabytes, not a validation run.
2. **metaPR2: ruled out above.**
3. **Reprocessing raw 18S reads** from the paired coastal studies with a rank-consistent reference, so
 fungi and protists are assigned at the same level. This is the only route that would give a clean
 test on existing public data, and it is a substantial piece of work.
4. **A new paired survey with two markers** (ITS and 18S on the same DNA) would settle it properly and
 is the right thing to propose rather than to improvise.

Until one of those exists, `EUK_RESULTS.md` rests on one dataset, and the paper must say so.

## Files

`METAPR2_REPLICATION_PLAN.md` (frozen, not executed) · data at `~/research/datasets/metapr2/`
(Zenodo 21106507: `datasets.xlsx`, `samples.xlsx`, `asv.xlsx`, `asv.fasta`,
`metapr2_4.0_asv_abundance.tsv.gz`).
