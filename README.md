# Aquatic metabarcoding: depth-dependent testability and habitat preference

**Public research code and data repository** for analyses of habitat-preference inference in aquatic fungal metabarcoding data.

**Authors:** Gurgen Karapetyan and Astghik Z. Pepoyan.

The associated research article is being prepared for submission to *Journal of Microbiological Methods*. **The manuscript and editorial files are not distributed in the current repository tree.** The repository provides analysis code, frozen study plans, derived results, figures, and reproducibility information.

## Scientific scope

In a FunAqua water–sediment analysis, 49 of 100 testable fungal orders were significant in the primary subsample, while repeated rarefaction yielded 43–57 significant orders across draws. In an independent NEON controlled depth ladder (102 paired observations), the Jaccard overlap between testable sets at 100 and 10,000 reads was 0.281. The median absolute preference effect changed by 0.024 over the same depth range. Count-based post hoc heterogeneity testing found differences for 26 of 69 comparable orders after BH correction.

FunAqua and NEON differ in habitat definitions and taxonomic resolution; cross-dataset differences must not be attributed exclusively to sequencing depth. These are results of the accompanying analysis, not claims of biological causation.

## Repository map

| Path | Description |
|---|---|
| [scripts/](scripts/) | Analysis, plotting and integrity-check code |
| [plans/](plans/) | Original, unmodified preregistration plans and SHA-256 manifest |
| [results/](results/) | Derived CSV/TSV/JSON and model outputs |
| [figures/](figures/) | Six numbered main figures, PNG and vector PDF |
| [supplementary/](supplementary/) | Tables S1–S19 and supplementary figures |
| [docs/analysis/](docs/analysis/) | Detailed findings, unsuccessful tests and methodological decisions |
| [docs/REPRODUCIBILITY.md](docs/REPRODUCIBILITY.md) | Required data, dependencies, commands and limitations |
| [docs/RELEASE_CHECKS.md](docs/RELEASE_CHECKS.md) | What has and has not been independently verified |
| [requirements.txt](requirements.txt) | Python dependencies |

## Quick verification

Use Python 3.11 or newer.

```bash
git clone https://github.com/kgurgen97/aquatic-fungi-habitat.git
cd aquatic-fungi-habitat
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/verify_plan_hashes.py
python scripts/analysis_heterogeneity.py
python scripts/analysis_agreement_power.py
python tests/test_publication.py
```

On Windows PowerShell activate with `.venv\Scripts\Activate.ps1`. The last three commands rely on committed supplementary input tables and do not need raw reads. Run `python scripts/make_figures.py` to regenerate the six primary figures from derived result files. Full source-data reprocessing requires the original FunAqua, NEON and related reference datasets, as detailed in the reproducibility guide.

## Data provenance

- FunAqua: Prins et al. (2026), *Scientific Data*, [10.1038/s41597-026-07729-8](https://doi.org/10.1038/s41597-026-07729-8); archived dataset [10.5281/zenodo.17855266](https://doi.org/10.5281/zenodo.17855266).
- NEON surface-water taxonomy: DP1.20141.002, [10.48443/4xwg-p967](https://doi.org/10.48443/4xwg-p967).
- NEON benthic taxonomy: DP1.20086.002, [10.48443/fe8q-qh09](https://doi.org/10.48443/fe8q-qh09).

No raw sequence datasets are redistributed. Frozen analysis plans are retained byte-for-byte and can be checked using `scripts/verify_plan_hashes.py`.

## Citation and licence

Use [CITATION.cff](CITATION.cff) to cite the code. Cite the underlying datasets separately when using their data. The code is MIT licensed. Authored text, figures and derived tables are available under CC BY 4.0; the upstream datasets remain under their own terms.

The accompanying journal article has not yet been peer-reviewed or published. The results here are provided to support inspection and reproducibility.
