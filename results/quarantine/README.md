# Quarantined

`primary_order_rerun.tsv` — 99 orders, 54 significant, different column names, and 94 of its 98 shared
q-values differ from `results/primary_order.tsv`. No script in this repository regenerates it, and
`analysis_primary.py` is deterministic (verified: rerunning reproduces `primary_order.tsv` bit-identically).
Its provenance is therefore unknown and it must not be cited. Kept only so the discrepancy stays on record.
Superseded by `results/draw_stability.tsv`, which measures draw-to-draw variation properly.
