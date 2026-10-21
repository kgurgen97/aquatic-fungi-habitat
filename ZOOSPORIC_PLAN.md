# The flagellum hypothesis: does a motile zoospore predict water-column preference?

**Written 2026-10-21, before computing any group-level comparison.** Nothing here may be changed because
of a result. If a decision turns out to be wrong, the fix is a new, separately labelled analysis.

## Why this test and not the last one

`NEON_RESULTS.md` reports that the per-order preference table does not transfer between FunAqua and NEON
(sign agreement 14/27, p = 1.0), but that the two zoosporic orders testable in both did transfer, with
NEON's strongest two q-values. With n = 2 there was no test to run.

That was the wrong unit. The hypothesis was never about Zygophlyctidales in particular; it is about a
**trait**. A fungus with a flagellated, swimming dispersal stage should be found in the water column more
often than one that disperses as a non-motile spore. Tested as a trait, the comparison uses every testable
order instead of four named ones, and it can be run independently in each dataset.

## Hypothesis

> **H:** among orders testable in a paired water/benthic design, those in flagellated (zoosporic) fungal
> phyla have a higher water preference index than those in aflagellate phyla.
>
> Direction is predicted: **flagellated > aflagellate**. A result in the opposite direction is a
> contradiction, not a success.

## The trait assignment — fixed now, from phylum membership, not from the data

Assignment is by phylum, from fungal biology, and was written down before any group means were computed.
Phylum labels come from each dataset's own taxonomy. Both datasets were checked first for label
compatibility: **216 orders are named in both and the two disagree on phylum for none of them**, so the
assignment transfers without a crosswalk.

**Flagellated (zoosporic) — motile zoospore with a functional flagellum:**
Chytridiomycota, Monoblepharomycota, Neocallimastigomycota, Blastocladiomycota, Olpidiomycota,
Rozellomycota.

**Aflagellate — no flagellated stage:**
Ascomycota, Basidiomycota, Mucoromycota, Glomeromycota, Kickxellomycota, Zoopagomycota,
Entomophthoromycota, Basidiobolomycota, Mortierellomycota, Calcarisporiellomycota, Entorrhizomycota.

**Excluded as biologically ambiguous, decided before seeing their values:**

| phylum | why excluded |
|---|---|
| Sanchytriomycota | sanchytrids lost the flagellum; their zoospores are amoeboid, so they are neither cleanly flagellated nor a typical aflagellate fungus |
| Aphelidiomycota | aphelids vary by genus between flagellated and amoeboid zoospores, so the phylum cannot be assigned as a whole |
| `Fungi_phy_Incertae_sedis` and any unplaced label | no phylum to assign |

These three are reported separately as a declared aside. They are **not** silently folded into whichever
group would help.

## Measurement

The preference index is the one already defined and used in both datasets: for each order,
`water_fraction` = water-only discordant pairs / all discordant pairs, computed after within-pair
rarefaction to the pair minimum, over orders with at least 20 discordant pairs. No new measurement is
introduced. FunAqua values come from `results/primary_order.tsv`; NEON values from
`results/neon_order.tsv`.

## Tests

Run **separately and independently in each dataset**, so that NEON is a genuine replication and not a
pooled re-description:

1. **Primary:** Mann–Whitney U on `water_fraction`, flagellated vs aflagellate, two-sided.
   Reported with the median of each group and the rank-biserial correlation as effect size.
2. **Secondary:** proportion of each group whose `water_fraction` exceeds 0.5, compared by Fisher exact.
3. **Secondary:** the same Mann–Whitney restricted to FunAqua's draw-stable orders
   (significant in ≥ 95 % of 200 draws), since the full set carries draw noise.
4. **Secondary:** a permutation test shuffling the flagellated/aflagellate label across orders
   (10 000 shuffles), as a distribution-free check on the primary.

**Replication is declared only if the effect is in the predicted direction with p < 0.05 in *both*
datasets independently.** One dataset alone is a finding about that dataset, not a replicated result.

## What would falsify this

A null or an opposite-direction result in either dataset. In particular, Reñé et al. 2023 (Gulf of Naples)
report *higher* chytrid diversity and relative abundance in sediment than in the water column. That is a
different measurement — diversity at one site over time, not paired presence/absence preference — but it
is a published result pointing the other way, and if this test comes out null it should be read alongside
that paper rather than explained away.

## Known limitations, stated before the result

- Phylum is a coarse proxy for the trait. Within Chytridiomycota some lineages are sediment-associated;
  the test asks whether the trait shifts the distribution, not whether it determines each order.
- The flagellated phyla are also, on the whole, the early-diverging ones, so phylum-level phylogeny and
  the trait are confounded. This test cannot separate "has a flagellum" from "is an early-diverging
  lineage". That limitation is inherent and must be stated in any write-up, not buried.
- NEON's second compartment is benthic biofilm, FunAqua's is sediment (see `NEON_VALIDATION_PLAN.md`
  section 8). The contrast is related, not identical.
- Orders are not independent observations; the permutation test in secondary 4 addresses the significance
  question but not the underlying non-independence.

## Deliverables

`analysis_zoosporic.py`, `results/zoosporic.tsv`, `results/zoosporic.json`, and a result section appended
to `NEON_RESULTS.md` reporting the primary first, whatever it says.
