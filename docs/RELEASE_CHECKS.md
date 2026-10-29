# Public code-and-data release checks

Updated: 29 October 2026

## Verified in the automated checks
- Original frozen plans have unchanged contents and match the SHA-256 manifest.
- Python sources compile under the configured test environment.
- Derived checks recover the recorded main counts: 49/100 FunAqua, 37/103 NEON and 26/69 heterogeneous shared orders.
- The S15–S19 source tables reproduce byte-for-byte using committed inputs.
- Main figure filenames and formats are present, along with all 20 supplementary tables.

## Publication boundaries
- The current branch does **not** contain the manuscript, author cover letter, editorial highlights or historical manuscript drafts.
- Public Git history **may still contain earlier versions**. Removing a tracked file in a new commit does not expunge its previous commits or cached copies.
- The full source-data pipeline has not been independently rerun from the original FunAqua and NEON files.
- Table S20 is a scenario sensitivity exercise and has not been independently regenerated.
- Figure rendering and visual comparison against a submission-ready manuscript still require the authors' final check.
- Journal declarations, funding wording and both authors' approval are separate from this code repository.

Preserve frozen protocols, sources and negative analyses. Do not rewrite the git history without a separately planned provenance-preserving migration.
