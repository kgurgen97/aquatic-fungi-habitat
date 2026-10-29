#!/usr/bin/env python3
"""
analysis_heterogeneity.py  --  POST HOC, 2026-10-23.

Not covered by any frozen plan; labelled post hoc throughout. Reads only the supplementary tables.

Two questions the pre-registered sign-agreement endpoint could not answer.

  1. Is the per-order disagreement between FunAqua and NEON REAL, or just imprecision? The sign test
     throws away the discordant-pair counts and keeps only their sign, which is why it is underpowered
     (section 3.3a). Testing the counts directly -- a 2x2 Fisher exact test of (water-only, sediment-only)
     in FunAqua against (water-only, benthic-only) in NEON, BH-corrected across the 69 shared testable
     orders -- asks whether the two datasets give significantly DIFFERENT index values.

  2. Does the non-replication depend on the >= 20-discordant-pair threshold? DEPTH_LADDER_RESULTS.md
     lists "one threshold" as a limitation and the sweep varied it only in simulation. Both datasets'
     per-order tables carry discordant counts, so the endpoint can be recomputed at higher thresholds on
     real data (not lower: orders below 20 were never tabulated).

Outputs supplementary/Table_S18.tsv and Table_S19.tsv.
"""

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

import pathlib
SUP = str(pathlib.Path(__file__).resolve().parent.parent / "supplementary")


def heterogeneity(s1, s7):
    shared = sorted(set(s1["order"]) & set(s7["order"]))
    f = s1.set_index("order").loc[shared]
    n = s7.set_index("order").loc[shared]
    d = pd.DataFrame({
        "order": shared,
        "w_fun": f.water_only.values.astype(int), "s_fun": f.sediment_only.values.astype(int),
        "w_neon": n.water_only.values.astype(int), "s_neon": n.benthic_only.values.astype(int),
        "idx_fun": f.water_fraction.values, "idx_neon": n.water_fraction.values,
        "q_fun": f.q.values, "q_neon": n.q.values,
    })
    d["p_het"] = [stats.fisher_exact([[r.w_fun, r.s_fun], [r.w_neon, r.s_neon]])[1]
                  for r in d.itertuples()]
    d["q_het"] = multipletests(d.p_het, method="fdr_bh")[1]
    d["agree"] = np.where((d.idx_fun == 0.5) | (d.idx_neon == 0.5), np.nan,
                          (d.idx_fun > 0.5) == (d.idx_neon > 0.5))
    d["disc_min"] = np.minimum(d.w_fun + d.s_fun, d.w_neon + d.s_neon)
    return d


def threshold_sensitivity(s1, s7, thresholds=(20, 30, 40, 50, 60)):
    rows = []
    for thr in thresholds:
        per = {}
        for label, df, wcol, scol in [("FunAqua", s1, "water_only", "sediment_only"),
                                      ("NEON", s7, "water_only", "benthic_only")]:
            d = df[df.discordant >= thr].copy()
            p = [stats.binomtest(int(r[wcol]), int(r[wcol] + r[scol]), 0.5).pvalue
                 for _, r in d.iterrows()]
            d["q2"] = multipletests(p, method="fdr_bh")[1] if len(p) else np.array([])
            per[label] = d
            rows.append({"threshold": thr, "dataset": label, "tested": len(d),
                         "significant": int((d.q2 < 0.05).sum()),
                         "median_effect": float((d.water_fraction - 0.5).abs().median()),
                         "n_comparable": np.nan, "n_agree": np.nan, "agreement": np.nan})
        sig = per["FunAqua"][per["FunAqua"].q2 < 0.05]
        c = sig.merge(per["NEON"], on="order", suffixes=("_f", "_n"))
        c = c[(c.water_fraction_f != 0.5) & (c.water_fraction_n != 0.5)]
        k = int(((c.water_fraction_f > 0.5) == (c.water_fraction_n > 0.5)).sum())
        rows[-1].update(n_comparable=len(c), n_agree=k,
                        agreement=k / len(c) if len(c) else np.nan)
    return pd.DataFrame(rows)


def main():
    s1 = pd.read_csv(f"{SUP}/Table_S1.tsv", sep="\t")
    s7 = pd.read_csv(f"{SUP}/Table_S7.tsv", sep="\t")

    het = heterogeneity(s1, s7)
    het.to_csv(f"{SUP}/Table_S18.tsv", sep="\t", index=False, float_format="%.6g")
    sig = het.q_het < 0.05
    dis = het.agree == False                                          # noqa: E712
    pre = het[(het.q_fun < 0.05) & het.agree.notna()]
    print(f"{len(het)} shared testable orders; {int(sig.sum())} significantly heterogeneous "
          f"({sig.mean():.1%})")
    print(f"direction-disagreements: {int(dis.sum())}, of which {int((dis & sig).sum())} significant")
    print(f"agree in direction yet significantly different in magnitude: "
          f"{int(((het.agree == True) & sig).sum())}")                # noqa: E712
    print(f"pre-registered 26: {int((pre.q_het < 0.05).sum())} heterogeneous; of its "
          f"{int((pre.agree == False).sum())} disagreements "         # noqa: E712
          f"{int(((pre.agree == False) & (pre.q_het < 0.05)).sum())} significant")  # noqa: E712
    a, b = het.loc[sig, "disc_min"], het.loc[~sig, "disc_min"]
    print(f"power check -- median min-discordant {a.median():.0f} vs {b.median():.0f}, "
          f"Mann-Whitney p = {stats.mannwhitneyu(a, b)[1]:.3f}")

    ts = threshold_sensitivity(s1, s7)
    ts.to_csv(f"{SUP}/Table_S19.tsv", sep="\t", index=False, float_format="%.6g")
    print(ts.to_string(index=False, float_format=lambda x: f"{x:.4g}"))


if __name__ == "__main__":
    main()