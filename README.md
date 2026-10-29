# Sequencing depth and fungal habitat preferences

This repository contains the data products, code, analysis plans and manuscript for a methodological study of habitat-preference inference in aquatic metabarcoding surveys. The manuscript is being prepared for *Journal of Microbiological Methods*.

**Authors:** Gurgen Karapetyan and Astghik Z. Pepoyan.

## Study summary

A minimum discordant-pair threshold determines which taxonomic orders can be tested for water or sediment preference. In the FunAqua analysis, 49 of 100 testable orders were significant in the primary draw; the count varied from 43 to 57 across 200 rarefaction draws, and 35 orders were significant in at least 95% of draws.

In a controlled NEON sequencing-depth experiment using the same 102 pairs, the sets of testable orders at 100 and 10,000 reads overlapped by only 28.1% (Jaccard), while median absolute preference changed by 0.024. A post hoc count-based comparison found heterogeneity for 26 of 69 shared testable orders after BH correction. These datasets also differ in habitat definition and taxonomy, so depth alone cannot explain their between-study disagreement.

The study reports unsuccessful preregistered predictions and sensitivity analyses as well as positive findings.

## Files

| Location | Contents |
|---|---|
| [MANUSCRIPT.md](MANUSCRIPT.md) | Current scientific manuscript, including declarations |
| [SUPPLEMENTARY.md](SUPPLEMENTARY.md) | Supplementary methods, figure captions and table index |
| [figures/](figures/) | Six main figures, PDF and PNG |
| [supplementary/](supplementary/) | Tables S1–S20, Figures S1–S3, generated Figure S4 |
| [results/](results/) | Derived tabular and JSON results |
| [docs/analysis/](docs/analysis/) | Individual analysis reports and negative findings |
| [docs/history/](docs/history/) | Older drafts and editorial decisions, kept for provenance |
| `*PLAN*.md`, [PLAN_HASHES.md](PLAN_HASHES.md) | Frozen plans and SHA-256 checks |
| [submission/JMM/](submission/JMM/) | Cover letter, highlights and journal-specific notes |

## Reproduce the analyses

Python 3.11 or later is recommended.

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python verify_plan_hashes.py
python analysis_heterogeneity.py
python analysis_agreement_power.py
python tests/test_publication.py
```

The last two analyses use only the committed supplementary TSV files. Complete reprocessing of FunAqua and NEON requires the source datasets and the environment variables described in [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md). This repository does **not** include raw sequencing reads or every intermediate matrix.

Regenerate the main figures from the committed derived results with `python make_figures.py`. Regenerate supplementary Figure S4 from the tables with `python analysis_agreement_power.py`.

## Data sources

- FunAqua: Prins et al., *Scientific Data* (2026), DOI [10.1038/s41597-026-07729-8](https://doi.org/10.1038/s41597-026-07729-8); dataset [10.5281/zenodo.17855266](https://doi.org/10.5281/zenodo.17855266).
- NEON surface water microbial taxonomy: DP1.20141.002, [10.48443/4xwg-p967](https://doi.org/10.48443/4xwg-p967).
- NEON benthic microbial taxonomy: DP1.20086.002, [10.48443/fe8q-qh09](https://doi.org/10.48443/fe8q-qh09).

See the manuscript for the data-processing and eligibility rules. All eleven preregistration documents remain unchanged.

## Reuse and publication status

The repository is public, but the manuscript has not been accepted or published in a journal. The final version and submission declarations require approval from both authors. Details of the remaining verification scope appear in [docs/RELEASE_CHECKS.md](docs/RELEASE_CHECKS.md).

Code is available under the MIT license. Original manuscript text, figures and derived tables are shared under CC BY 4.0; third-party input datasets retain their own terms.
