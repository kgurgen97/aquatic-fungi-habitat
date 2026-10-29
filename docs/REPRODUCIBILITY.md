# Reproducibility

## Scope

This repository contains derived tables and code. The committed supplementary tables allow the agreement/power and heterogeneity analyses to be checked without a complete source-data download. It is **not** a self-contained raw-read archive or a verified clean-room run of the full pipeline.

## Installation

Python 3.11 or newer:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The analysis programs are in `scripts/`. Their shared `scripts/project_paths.py` resolves the repository root independently of the working directory.

Set local source data directories using these environment variables:

| Variable | Dataset directory |
|---|---|
| `FUNAQUA_DATA_DIR` | FunAqua OTU and sample-level tables |
| `NEON_DATA_DIR` | NEON processed taxonomy with `water/` and `benthic/` |
| `FUNGALTRAITS_DATA_DIR` | Fungal traits |
| `METAPR2_DATA_DIR` | metaPR2 extracts including `asv.xlsx` |

Defaults are the original author's `~/research/datasets/...` directories. Users should set their own paths and never commit credentials, reads or external downloads. NEON acquisition may require a personal token.

## Lightweight checks using committed results

```bash
python scripts/verify_plan_hashes.py
python scripts/analysis_heterogeneity.py
python scripts/analysis_agreement_power.py
python tests/test_publication.py
```

The two analyses update Tables S15–S19 and supplementary Figure S4; CI checks that the regenerated tables are identical to those committed. Those outputs depend on upstream filtered results S1/S7/S8 and therefore do **not** validate upstream processing from raw reads.

## Full analyses and graphics

The programs for primary and validation work include:

```bash
python scripts/analysis_primary.py
python scripts/analysis_draw_stability.py
python scripts/analysis_neon_validation.py
python scripts/analysis_depth_ladder.py
python scripts/analysis_sweep.py
python scripts/analysis_benchmark2.py
python scripts/analysis_occupancy.py
python scripts/make_figures.py
```

These are **entry points, not an independently tested end-to-end execution sequence**. Each script may require generated intermediates and raw downloads from preceding steps. In particular, `results/presence_matrices.npz` is not committed. Consult the original frozen protocols in `plans/` before running and do not change preregistered settings to produce different outcomes.

The main plotting script uses committed derived results. The optional `scripts/make_supplementary.py` regenerates Tables S1–S14 and Figures S1–S3 but requires source material not all included here.

## Integrity and limits

The 11 original preregistration files are preserved verbatim, with SHA-256 values and originating commits in `plans/PLAN_HASHES.md`. The integrity test only checks agreement with those recorded hashes. Historical negative and superseded analysis reports remain in `docs/analysis/` to preserve the scientific record.

Table S20 shows scenario sensitivity under *assumed* intraclass correlations, not estimated site-level ICC. A complete independent raw-data reproduction and a separate independent audit of S20 remain outstanding.
