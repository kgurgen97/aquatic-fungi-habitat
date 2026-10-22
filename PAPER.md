# Aquatic fungi: habitat specialisation and what transfers between datasets

**Status: draft written, figures done, not submittable yet.** The blocker is the reference list.

Working title: *What transfers between amplicon datasets and what does not: habitat specialisation of
aquatic fungi measured twice.* Target journal: Molecular Ecology Resources (Q1, methods-leaning);
alternatives Environmental Microbiology, ISME Communications.

---

## What the paper says

Three findings, in the order the manuscript makes them.

1. **The headline count is unstable to the rarefaction draw.** Over 200 draws the number of fungal orders
   with a significant water/sediment preference ranges 43–57. Only **35 orders** are significant in at
   least 95 % of draws, and that set — not any single draw's count — is the defensible unit.
2. **Per-order preferences do not transfer.** In NEON, an independent paired dataset, habitat
   partitioning reproduces (37 of 103 orders) but the *direction* agrees with FunAqua for only
   **14 of 27** comparable orders (51.9 %, *p* = 1.0). Chance.
3. **What does transfer is a group-level contrast.** Fungi sit at the generalist end of the eukaryotic
   range: median effect **0.156** against **0.228** for non-fungal eukaryotes in the same assay and
   samples, with non-overlapping ranges over 200 draws, zero reversals, and all five pre-registered
   robustness checks passed.

The claim is stated narrowly on purpose: fungi are **indistinguishable from Alveolata**, and the gap
against "non-fungal eukaryotes" as a pooled group is driven by Metazoa and Straminipila.

## What is in this folder

| folder | contents |
|---|---|
| `MANUSCRIPT.md` | the full draft, every number quoted from a file in `results/` |
| `SUPPLEMENTARY.md` | supplementary figures S1–S3, tables S1–S14, supplementary methods |
| `figures/` | the four main figures, vector PDF and 600 dpi PNG |
| `supplementary/` | Table_S1…S14 and Figure_S1…S3 |
| `plans/` | the frozen pre-registration plans, each with its SHA-256 in its header |
| `code/` | every analysis script that produced a number in the paper |
| `results/` | the result tables the manuscript cites |
| `*_RESULTS.md` | the long-form write-up behind each section |
| `SNAPSHOT.txt` | which commit this copy was taken from, and the full history |

## Links

- **Code and full history:** https://github.com/kgurgen97/aquatic-fungi-habitat (private)
- **FunAqua dataset:** CC-BY; ITS long-read metabarcoding, 50 152 fungal OTUs, 2 375 samples,
  1 452 sites, 87 countries. *Exact citation and Zenodo DOI still to be filled in.*
- **NEON:** products `DP1.20141.002` (surface water) and `DP1.20086.002` (benthic), obtained through the
  NEON API. The superseded `.001` products assign only 5.2 % of reads to an order and must not be used.
- **metaPR2 v4.0:** https://doi.org/10.5281/zenodo.21106507 (MIT). Assessed and deliberately **not**
  used — see `METAPR2_RESULTS.md`.

## Where the data actually lives

Not in this folder; it is gigabytes. All of it is under WSL at `~/research/datasets/`:
`funaqua/` (656 MB), `neon2/` (624 MB, 324 water + 2 131 benthic ITS tables),
`metapr2/` (551 MB), `fungaltraits/`.

## What still has to happen before submission

1. **Complete and verify every reference.** The list in `MANUSCRIPT.md` § 6 is incomplete and several
   entries are recorded by topic and first author only. Nothing there is safe to cite as written. This is
   the only thing standing between the draft and a submission.
2. **Decide the framing.** The draft leads with transferability, which makes the two failed external
   validations contributions rather than weaknesses. The alternative leads with the fungal generalist
   result; it reads better but is more exposed, since that headline then rests on one dataset with an
   off-target amplicon.
3. **Decide whether Figure 4 belongs in the main text or the supplement** — it is post hoc and reaches
   significance at no depth cap.

## The honest limitation, which cannot be fixed with existing data

ITS is a fungal barcode, so the non-fungal eukaryotes it captures are off-target amplification. The
design controls sample, day, assay and depth; it cannot make the two groups equally well served by the
primers or the reference database. **No external validation of the headline is possible**: NEON carries
no eukaryote marker in its paired products, and metaPR2 assigns fungi systematically coarser ranks than
protists, so a test there would have confirmed the hypothesis for the wrong reason. Closing this needs a
new paired survey carrying ITS and 18S on the same DNA — a grant, not a revision.
