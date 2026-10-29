"""Lightweight checks on committed derived data. Not a raw-data rerun."""
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
sup = root / "supplementary"
s1 = pd.read_csv(sup / "Table_S1.tsv", sep="\t")
s7 = pd.read_csv(sup / "Table_S7.tsv", sep="\t")
s15 = pd.read_csv(sup / "Table_S15.tsv", sep="\t")
s18 = pd.read_csv(sup / "Table_S18.tsv", sep="\t")
assert len(s1) == 100 and int((s1.q < 0.05).sum()) == 49
assert len(s7) == 103 and int((s7.q < 0.05).sum()) == 37
assert len(s18) == 69 and int((s18.q_het < 0.05).sum()) == 26
assert int(s15.loc[s15.selection_set == "all orders testable in both", "n_orders"].iloc[0]) == 65
for i in range(1, 21):
    assert (sup / f"Table_S{i}.tsv").is_file(), f"Missing Table S{i}"
for i, stem in enumerate(("draw_stability", "transfer_failure", "kingdoms", "depth_probe", "simulation", "mechanism_in_data"), 1):
    assert (root / "figures" / f"fig{i}_{stem}.pdf").is_file()
    assert (root / "figures" / f"fig{i}_{stem}.png").is_file()
print("Derived-results and file-manifest checks passed")
