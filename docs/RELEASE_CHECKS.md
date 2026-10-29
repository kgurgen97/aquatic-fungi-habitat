# Public release and submission checks (29 October 2026)

## Checked

- The public-facing layout separates current manuscript files, frozen plans, analysis reports and historical drafts.
- The six principal figure files are named consistently with the figure captions.
- The new S15–S19 scripts were run against supplied supplementary inputs and their resulting TSV contents matched the supplied tables numerically and structurally.
- The analysis scripts now use repository-relative project paths and configurable external dataset directories.
- Frozen plan documents were not modified.

## Remaining author and data checks

- A full clean-clone execution from raw FunAqua and NEON sequencing products has **not** been verified. Some upstream intermediate matrices are not distributed.
- Table S20 evaluates assumed site-level ICC scenarios and is not a measured ICC analysis; there is no independently verified generator for S20 in this repository.
- Review all declarations (grant support and funder roles, author contributions, competing interests, use of generative AI) and obtain both authors' approval before journal submission.
- Verify the exact output of `make_figures.py` against the submitted images on a machine with complete derived results and supported fonts.
- The submission PDF/Word files live in a separate local submission package, not in this GitHub source repository.

Deleted or moved paths from a public branch remain in the Git history. Do not treat moving historical files as erasing their contents.
