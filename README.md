# Do aquatic fungi specialise between the water column and the sediment?

FunAqua is a global paired survey — the same sites sampled in both water and sediment — and nobody had
asked that question of the fungi in it. The closest work I could find, Malfertheiner et al. 2026 on
community conservatism, explicitly threw out all the eukaryotic reads, so the fungi were genuinely open.

I wrote the analysis plan before computing any water-vs-sediment outcome (`ANALYSIS_PLAN.md`), got a
positive result, and then tried to replicate it in NEON. It did not replicate. Both halves are in here,
because the second half is the more useful one.

## FunAqua

239 sites have both a water and a sediment sample deep enough to use. Sediment samples are systematically
deeper than water samples (median 179 reads against 131, p 1.5e-13) and richness tracks depth at Spearman
0.84, so any naive comparison would mostly measure sequencing effort. The fix is in the design: both
samples of a pair are rarefied down to that pair's own minimum before anything is counted, which makes
detection power equal inside every pair by construction.

At order level, 100 orders have enough discordant pairs to test and 49 show a significant habitat
preference, 30 toward water and 19 toward sediment. The split is biologically coherent, which matters
because nothing in the method knows any fungal biology: the zoosporic lineages with a motile flagellate
stage — Zygophlyctidales, Chytridiales, Zygorhizidiales, Rozellales — go to the water column, while
Sordariales, Microascales and Zoopagales, all soil and dung organisms, go to the sediment.

Preference is also conserved within classes more than chance allows (permutation p 0.0041). I am careful
to call this *taxonomic* and not phylogenetic conservatism. I tried three routes to a real phylogeny and
all three are closed: the amplicon carries only about 61 bp of LSU by ITSx, so neither the authors'
representative sequences nor the raw reads can resolve deep fungal relationships, and the Open Tree
synthetic fungal tree is almost entirely polytomous — 4 resolved sister pairs out of 70 orders.

One result I did not expect: fungi partition *less* than other eukaryotes. On the same 204 pairs at the
same depth, 43 % of fungal orders show a preference against 63 % of non-fungal eukaryote orders
(Fisher p 0.0029), and the median effect is smaller too. My first attempt at this compared 239 fungal
pairs against 580 eukaryote pairs and I threw it away — eukaryotes carry more reads, so more sites passed
the filter and the "difference" was power. The matched version survives.

## The count depends on the random draw, and I nearly published the wrong number

While setting up the NEON comparison I noticed that `analysis_primary.py` and `analysis_secondary.py`
implement the same declared endpoint, at the same threshold, with the same seed, and disagree: 100 orders
tested and 49 significant against 99 and 51. Neither is miscoded. `analysis_secondary.py` creates one RNG
at module level and consumes it across its four read thresholds, so by the time it reaches the ≥100 row
its rarefaction draw differs from the primary's fresh RNG. `RESULTS.md` had been quoting the secondary
script's row as the primary result.

So I ran the primary over 200 draws (`analysis_draw_stability.py`). The number of significant orders
ranges from 43 to 57, median 51, sd 2.5. Both stored runs sit inside that range, which means neither was
wrong — the single-draw design is just unstable at about ±5 orders. What is stable is the set:
**35 orders are significant in at least 95 % of draws** (21 water, 14 sediment) and 28 in every single
draw. That set is what belongs in a paper, not any one draw's count.

This also cost me a result I had already written down. The three orders I reported as disagreeing between
Europe and North America — Auriculariales, Pezizales, Trechisporales — turn out to be significant in only
27, 136 and 185 of 200 draws. That disagreement was substantially draw noise, not biogeography.

`results/quarantine/` holds a third results file whose provenance I can't establish; it contradicts the
deterministic primary and no script in here regenerates it. Kept so the discrepancy stays on record.

## NEON: not replicated

Plan frozen in `NEON_VALIDATION_PLAN.md` and committed before the data finished downloading, with one
pre-committed number: of the FunAqua-significant orders also testable in NEON, the fraction agreeing in
sign, tested by exact binomial against 0.5. Replication only if agreement beat 50 % with p < 0.05.

Agreement came out at **14 of 27, 51.9 %, p = 1.0**. That is a coin flip. By the rule I wrote down in
advance, this is not replicated, and that is the headline.

Worth being precise about what failed, though, because NEON is not a null dataset. On its own it finds 37
of 103 orders with a significant preference, 18 water and 19 benthic, against FunAqua's 49 of 100. So
habitat partitioning of fungal orders is real in both datasets, measured by different people with
different primers and a different pipeline. What does not transfer is *which* order prefers *which*
habitat. A per-order preference table cannot be presented as a general result on this evidence.

The part that did replicate was pre-specified, in section 7.3 of the plan, before I saw anything:
Zygophlyctidales goes 0.870 → 0.875 (q 1.8e-4) and Chytridiales 0.752 → 0.933 (q 2.3e-11) — the two
strongest q-values in the whole NEON table. Microascales agrees too (q 0.022). Sordariales doesn't
(q 1.0), and the other named lineages weren't testable. With two zoosporic orders there's no formal test
to run, so this is a consistent observation rather than a demonstrated effect — but the mechanistically
specific part of the claim is the part that survived, while the broad pattern didn't.

Why it might have failed, as hypotheses rather than findings. The sequencing depth differs by about 200×
(FunAqua 131–179 reads per sample, NEON a median of 31,131), and within-pair rarefaction equalises power
inside a pair without making "present" mean the same thing in both datasets — at 130 reads it means
"abundant", at 31,000 it means "detectable". Capping NEON toward FunAqua's depth does raise agreement from
51.9 % to around 69 % (`exploratory_neon_depth.py`), which is the direction the hypothesis predicts, but
it never reaches significance and the overlap shrinks to 16–26 orders, so the probe is too underpowered to
confirm its own suggestion. Separately, NEON's second habitat is benthic biofilm, not sediment, and it is
23 US wadeable streams with no lakes. I wrote both of those down in the plan before the result, so they
are limitations rather than excuses.

## Where I think this goes

The defensible claim is that fungal orders partition between water column and benthos, now shown in two
independent datasets. The per-order table goes in a supplement with an explicit note that it doesn't
transfer. And the thing actually worth building a paper around is the zoosporic hypothesis — a motile
flagellate stage predicts water-column preference — because it was predicted in advance, it replicated,
and it has a mechanism.

Any future cross-dataset comparison of this kind needs matched read depth, not just within-pair
rarefaction. That is the methodological lesson and it was not obvious to me beforehand.

## Running it

Analyses are plain Python, needing numpy, pandas, scipy and statsmodels. Each script is standalone and
writes into `results/`.

    analysis_primary.py           the pre-registered FunAqua primary (seed 20261021, deterministic)
    analysis_secondary.py         thresholds, conservatism, regions, FungalTraits, eukaryote benchmark
    analysis_euk_matched.py       the depth- and pair-matched eukaryote comparison
    analysis_biome_kingdom.py     rivers vs lakes, and the per-kingdom breakdown
    analysis_draw_stability.py    200 draws; post hoc, added after finding the RNG problem
    analysis_neon_validation.py   the pre-registered NEON test
    exploratory_neon_depth.py     the depth probe; post hoc

FunAqua's OTU tables go in `~/research/datasets/funaqua/`, FungalTraits in
`~/research/datasets/fungaltraits/`. The NEON tables come from the API and need a token in
`~/.neon_token` — DP1.20141.002 for surface water and DP1.20086.002 for benthic. Do not use the `.001`
versions of those products; they are superseded and assign only about 5 % of reads to an order, which I
found out the slow way.

Nothing larger than the metadata is committed here.
