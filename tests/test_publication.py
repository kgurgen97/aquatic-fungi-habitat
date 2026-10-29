"""Minimal derived-result checks and public-release file inventory."""
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
sup = root / "supplementary"
f = pd.read_csv(sup / "Table_S1.tsv", sep="\t")
n = pd.read_csv(sup / "Table_S7.tsv", sep="\t")
a = pd.read_csv(sup / "Table_S15.tsv", sep="\t")
h = pd.read_csv(sup / "Table_S18.tsv", sep="\t")
assert len(f) == 100 and int((f.q < 0.05).sum()) == 49
assert len(n) == 103 and int((n.q < 0.05).sum()) == 37
assert len(h) == 69 and int((h.q_het < 0.05).sum()) == 26
assert int(a.loc[a.selection_set == "all orders testable in both", "n_orders"].iloc[0]) == 65
for i in range(1, 20):
    assert (sup / f"Table_S{i}.tsv").is_file(), f"Missing Table S{i}"
assert not (sup / "Table_S20.tsv").exists(), "Excluded unverified supplementary Table S20 should not be released"
for i, stem in enumerate(("draw_stability", "transfer_failure", "kingdoms", "depth_probe", "simulation", "mechanism_in_data"), 1):
    for extension in ("pdf", "png"):
        assert (root / "figures" / f"fig{i}_{stem}.{extension}").is_file()
for forbidden in ("MANUSCRIPT.md", "SUPPLEMENTARY.md", "SUBMISSION.md", "COVER_LETTER.md"):
    assert not (root / forbidden).exists(), f"Editorial text in public root: {forbidden}"
assert (root / "plans" / "PLAN_HASHES.md").is_file()
assert (root / "scripts" / "project_paths.py").is_file()
print("Public code/data checks passed")
