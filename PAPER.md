# Aquatic fungi: what a minimum-evidence threshold decides a survey can measure

**Status, 2026-10-22: draft complete, figures done, references verified, framing decided.** No blocker
remains that is mine to clear. What is left is editorial — see the last section.

Title: *Minimum-evidence thresholds decide which taxa a metabarcoding survey can measure, not how
specialised they look: why per-taxon habitat preferences do not transfer between datasets.*
Target journal: Molecular Ecology Resources (Q1, methods-leaning); alternatives Environmental
Microbiology, ISME Communications.

---

## What the paper says

Four findings, in the order the manuscript makes them, which is the order in which they constrain each
other.

1. **The mechanism, and the paper's main result (§ 3.4).** Any analysis that tests taxa for habitat
   preference must first decide which taxa carry enough evidence to test. That threshold is not a neutral
   filter: a taxon contributes evidence only while its detection probability sits away from both 0 and 1,
   so the rule admits taxa whose abundance × depth is near one, at whatever depth is used. Simulated
   across a 1000-fold depth range, the testable set's median detection parameter holds at **0.22–0.49**
   while its median abundance falls **454-fold**. Depth changes *which* taxa are measurable, not how
   specialised they appear.
2. **The headline count is unstable to the rarefaction draw (§ 3.2).** Over 200 draws the number of
   fungal orders with a significant water/sediment preference ranges **43–57**. Only **35 orders** are
   significant in at least 95 % of draws, and that set — not any single draw's count — is the defensible
   unit.
3. **Per-order preferences do not transfer (§ 3.3).** In NEON, an independent paired dataset ~150-fold
   deeper, habitat partitioning reproduces (37 of 103 orders) but the *direction* agrees with FunAqua for
   only **14 of 27** comparable orders (51.9 %, *p* = 1.0). Chance. The mechanism accounts for this: the
   two datasets share only **69 of 134** testable orders, so they are largely measuring different taxa.
4. **What does transfer is a contrast between large groups (§ 3.5).** Animals and plants are the most
   compartment-bound (median effect **0.294** and **0.227**); fungi and the protist groups are not.
   Non-overlapping ranges over 200 draws, zero reversals, all five pre-registered robustness checks
   passed.

**Stated narrowly on purpose.** Fungi (0.150) are **indistinguishable from Alveolata** (0.157), and the
gap against "non-fungal eukaryotes" as a pooled group is driven by Metazoa, Straminipila and
Viridiplantae. On Levins' *B* pooled against protists, *p* = 0.163. The paper does **not** claim fungi are
generalists relative to eukaryotes generally — an earlier draft did, and § 3.5 contradicts it.

## The mechanism, in four pre-registered steps

Added after the first draft, each with its plan frozen and hashed before its analysis existed. Together
they are now the paper's core, and two of them failed in ways that changed what it claims.

1. **Depth ladder (§ 3.4a).** NEON's own libraries subsampled to 100/300/1 000/3 000/10 000 reads on a
   fixed set of 102 pairs, so ecology is held constant and only depth moves. The testable set turns over
   hard (Jaccard 0.281) while apparent effect size barely shifts (0.122–0.145). Four of five predictions
   confirmed; the fifth showed the set does *not* hold a fixed position on the detection curve.
2. **Parameter sweep (§ 3.4b).** 240 cells. Depth-invariance of apparent effect holds in **239 of 239**
   interpretable cells — the paper's firmest statement. But two predictions failed, and they carried my
   explanation of the ladder's exponent, which was withdrawn.
3. **Per-site abundance (§ 3.4c).** The failure was pooling. Representing abundance per sample moves the
   exponent from 0.958 to 0.484 against the measured 0.485; occupancy plus mean abundance — two numbers
   per taxon — recovers 0.432. The median order is present in only 9.8 % of samples.
4. **Estimator benchmark (§ 3.4d).** The first attempt was discarded: its generator encoded preference as
   a read ratio and so handed one estimator its own estimand. The replacement shows the ranking **flips**
   with the encoding, so there is no best estimator. Under realistic encoding a single draw runs at
   FDR 0.125, the draw-stable rule at 0.037.

## What is in this folder

| item | contents |
|---|---|
| `MANUSCRIPT.md` | the full draft, every number quoted from a file in `results/` |
| `SUPPLEMENTARY.md` | supplementary figures S1–S3, tables S1–S14, supplementary methods |
| `figures/` | the **six** main figures, vector PDF and 600 dpi PNG |
| `supplementary/` | Table_S1…S14, Figure_S1…S3, and `table_provenance.tsv` |
| `PLAN_HASHES.md` | every frozen plan's SHA-256 and the commit that introduced it |
| `*_RESULTS.md` | the long-form write-up behind each section. The four above are current; `RESULTS.md`
  and `EUK_RESULTS.md` predate the simulation and carry a banner saying so. Where any of them differs
  from `MANUSCRIPT.md`, the manuscript is authoritative. |

Not copied here: the frozen pre-registration plans, the analysis code and the 33 raw result tables. All
are in the repository, which is the reproducible record.

## Links

- **Code, frozen plans and full history:** https://github.com/kgurgen97/aquatic-fungi-habitat (private)
- **FunAqua:** Prins, V., Tedersoo, L., Mikryukov, V. *et al.* (2026) The FunAqua dataset of global fungal
  biodiversity in aquatic ecosystems. *Scientific Data* **13**, 1346. doi:10.1038/s41597-026-07729-8
  (240 authors). Data at doi:10.5281/zenodo.17855266, CC-BY. ITS long-read metabarcoding, 50 152 fungal
  OTUs, 2 375 samples, 1 452 sites, 87 countries.
- **NEON:** *Surface water microbe community taxonomy* `DP1.20141.002`, doi:10.48443/4xwg-p967, and
  *Benthic microbe community taxonomy* `DP1.20086.002`, doi:10.48443/fe8q-qh09; RELEASE-2026, retrieved
  through the NEON Data API on 21 October 2026. The superseded `.001` products assign only 5.2 % of reads
  to an order and must not be used.
- **metaPR2 v4.0:** doi:10.5281/zenodo.21106507 (MIT). Assessed and deliberately **not** used — see
  `METAPR2_RESULTS.md`. Note the frozen plan for it still carries the old "generalist end" wording; that
  is pre-registration and is left untouched on purpose.
- **Open Tree of Life:** synthesis `opentree16.1`, taxonomy 3.7draft3, queried 22 October 2026.

## Where the data actually lives

Not in this folder; it is gigabytes. All of it is under WSL at `~/research/datasets/`:
`funaqua/` (656 MB), `neon2/` (624 MB, 324 water + 2 131 benthic ITS tables),
`metapr2/` (551 MB), `fungaltraits/`.

## What still has to happen before submission

Both are editorial decisions, not analysis.

1. **Convert to the journal's format.** The draft is Markdown. Molecular Ecology Resources wants Word or
   LaTeX, with figures supplied separately — the PDFs are already at 89 mm and 183 mm.
2. **Decide whether Figure 4 belongs in the main text or the supplement.** It is post hoc, outside the
   frozen validation plan, and reaches significance at no depth cap (best *p* = 0.076). It is honest
   either way; it is weaker evidence than the rest of the main text.

Also outstanding and not mine to supply: the author affiliation, currently "to be completed".

## The honest limitation, which cannot be fixed with existing data

ITS is a fungal barcode, so the non-fungal eukaryotes it captures are off-target amplification. The design
controls sample, day, assay and depth — the same 204 pairs, the same per-pair depth, the same ≥ 20
discordant-pair rule for both groups, and the difference survives stratification on discordant count in
4 of 4 bins — but it cannot make the two groups equally well served by the primers or the reference
database, and which non-fungal clades ITS amplifies is not a random sample of eukaryotes.

**No external validation of § 3.5 is possible** with existing data: NEON carries no eukaryote marker in
its paired products, and metaPR2 assigns fungi systematically coarser ranks than protists, so a test there
would have confirmed the hypothesis for the wrong reason. Closing this needs a new paired survey carrying
ITS and 18S on the same DNA — a grant, not a revision.

Worth knowing if a reviewer raises depth as an objection to § 3.5: it does not hold, and § 3.4 is the
reason. Both groups pass through the same threshold and are therefore filtered to the same position on the
detection curve, so the remaining difference between them cannot be a depth artefact.
