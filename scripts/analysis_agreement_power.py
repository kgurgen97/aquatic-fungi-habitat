#!/usr/bin/env python3
"""
analysis_agreement_power.py  --  POST HOC, 2026-10-23.

Not covered by any frozen plan. Written after NEON_VALIDATION_PLAN.md's primary endpoint had
been run and reported, and labelled post hoc throughout, in the same way as
exploratory_neon_depth.py. It does three things the pre-registered endpoint did not:

  1. Corrects a tie in the published comparison set. Coniochaetales has a NEON water-preference
     index of exactly 0.500 (prefers_neon = "-"), i.e. no direction at all, yet Table_S8 flags it
     agree = True. A tie cannot agree in direction, so it is excluded: 14/27 becomes 13/26.
  2. Reports the agreement endpoint on ALL orders testable in both datasets, not only those
     significant in FunAqua. The pre-registered endpoint conditions on significance in one
     dataset, which both discards comparable orders and selects on noise.
  3. Quantifies what the endpoint could have detected (exact binomial power, minimum detectable
     agreement, n required) and whether the result is equivalent to chance, and checks both
     endpoints against taxonomic non-independence by bootstrapping phyla rather than orders.

Reads only the supplementary tables shipped with the manuscript; no raw data needed.
Run from the repository root or the paper folder.
"""

import math
import numpy as np
import pandas as pd
from scipy import stats

import pathlib
SUP = str(pathlib.Path(__file__).resolve().parents[1] / "supplementary")
Z95, Z90 = 1.959963985, 1.6448536269


def wilson(k, n, z=Z95):
    """Wilson score interval, the convention used for every CI in this project."""
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h) * 100, (c + h) * 100


def min_k(n, alpha=0.05):
    """Upper rejection threshold for the two-sided exact binomial test."""
    return int(stats.binom.isf(alpha / 2, n, 0.5)) + 1


def exact_power(n, pi, alpha=0.05):
    """Exact two-sided power, evaluated from the upper and lower binomial tails."""
    k = min_k(n, alpha)
    return float(stats.binom.sf(k - 1, n, pi) + stats.binom.cdf(n - k, n, pi))


def cluster_ci(df, pcts, B=20000, seed=11):
    """Bootstrap whole phyla rather than individual orders."""
    rng = np.random.default_rng(seed)
    g = df.groupby("phylum")["agree"].agg(["sum", "count"])
    successes = g["sum"].to_numpy(dtype=float)
    sizes = g["count"].to_numpy(dtype=float)
    draw = rng.integers(0, len(g), size=(B, len(g)))
    out = successes[draw].sum(axis=1) / sizes[draw].sum(axis=1)
    return tuple(np.percentile(out, pcts) * 100)


def main():
    s1 = pd.read_csv(f"{SUP}/Table_S1.tsv", sep="\t")    # FunAqua, 100 tested orders
    s2 = pd.read_csv(f"{SUP}/Table_S2.tsv", sep="\t")    # draw stability, 107 orders
    s7 = pd.read_csv(f"{SUP}/Table_S7.tsv", sep="\t")    # NEON, 103 tested orders
    s10 = pd.read_csv(f"{SUP}/Table_S10.tsv", sep="\t")  # order -> phylum

    shared = sorted(set(s1["order"]) & set(s7["order"]))
    f = s1.set_index("order").loc[shared]
    n = s7.set_index("order").loc[shared]
    both = pd.DataFrame({
        "order": shared,
        "frac_fun": f.water_fraction.values, "q_fun": f.q.values, "disc_fun": f.discordant.values,
        "frac_neon": n.water_fraction.values, "q_neon": n.q.values, "disc_neon": n.discordant.values,
    })
    n_shared = len(both)
    ties = both[(both.frac_fun == 0.5) | (both.frac_neon == 0.5)]
    both = both.drop(ties.index)                       # (1) a tie has no direction to agree about
    both["agree"] = (both.frac_fun > 0.5) == (both.frac_neon > 0.5)
    both["phylum"] = both["order"].map(
        s10.query("dataset == 'funaqua'").set_index("order")["phylum"].to_dict())
    print(f"{n_shared} orders testable in both; {len(ties)} dropped as ties "
          f"({', '.join(ties['order'])}); {len(both)} directional comparisons")

    stable = set(s2.loc[s2.times_significant / 200 >= 0.95, "order"])
    sets = {
        "all orders testable in both": both,
        "FunAqua-significant (published endpoint, tie-corrected)": both[both.q_fun < 0.05],
        "FunAqua draw-stable (>= 95 % of draws)": both[both["order"].isin(stable)],
        "significant in both datasets": both[(both.q_fun < 0.05) & (both.q_neon < 0.05)],
    }

    rows = []
    for label, sub in sets.items():
        k, m = int(sub.agree.sum()), len(sub)
        lo, hi = wilson(k, m)
        clo, chi = cluster_ci(sub, [2.5, 97.5])
        rows.append([label, m, k, k / m, stats.binomtest(k, m, 0.5).pvalue,
                     lo, hi, clo, chi])
    s15 = pd.DataFrame(rows, columns=[
        "selection_set", "n_orders", "n_agree", "agreement", "binom_p",
        "wilson_lo_pct", "wilson_hi_pct", "phylum_clustered_lo_pct", "phylum_clustered_hi_pct"])
    s15.to_csv(f"{SUP}/Table_S15.tsv", sep="\t", index=False, float_format="%.6g")

    rows = []
    for label, sub in [("published endpoint (tie-corrected)", sets[
                            "FunAqua-significant (published endpoint, tie-corrected)"]),
                       ("all orders testable in both", sets["all orders testable in both"])]:
        m = len(sub)
        kmin = min_k(m)
        for pi in (0.60, 0.65, 0.70, 0.75, 0.80):
            rows.append([label, m, kmin, kmin / m, pi, exact_power(m, pi),
                         next(x for x in range(10, 2000) if exact_power(x, pi) >= 0.80)])
    s16 = pd.DataFrame(rows, columns=[
        "endpoint", "n_orders", "min_k_for_p05", "min_agreement_for_p05",
        "true_agreement", "power", "n_for_80pct_power"])
    s16.to_csv(f"{SUP}/Table_S16.tsv", sep="\t", index=False, float_format="%.6g")

    both.to_csv(f"{SUP}/Table_S17.tsv", sep="\t", index=False, float_format="%.6g")

    print(s15.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    print(s16.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
    for label, sub in sets.items():
        k, m = int(sub.agree.sum()), len(sub)
        print(f"90 % Wilson CI, {label}: {wilson(k, m, Z90)[0]:.1f}-{wilson(k, m, Z90)[1]:.1f} %")



# ---------------------------------------------------------------- Figure S4
# Reads the tables written above, so the figure cannot drift from them.
# Palette and sizing match make_figures.py: blue #2a78d6, orange #eb6834,
# 183 mm double-column, vector PDF + 600 dpi PNG.

def make_figure_s4():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    BLUE, ORANGE = "#2a78d6", "#eb6834"
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.spines.top": False, "axes.spines.right": False})

    s15 = pd.read_csv(f"{SUP}/Table_S15.tsv", sep="\t")
    s16 = pd.read_csv(f"{SUP}/Table_S16.tsv", sep="\t")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(183 / 25.4, 75 / 25.4))

    order = ["all orders testable in both",
             "FunAqua-significant (published endpoint, tie-corrected)",
             "FunAqua draw-stable (>= 95 % of draws)",
             "significant in both datasets"]
    short = ["all testable\nin both (65)", "FunAqua-significant\n(published, 26)",
             "draw-stable\n(21)", "significant in\nboth (9)"]
    d = s15.set_index("selection_set").loc[order]
    y = np.arange(len(d))[::-1]
    ax1.axvline(50, color="0.35", lw=1, zorder=1)
    for yi, (_, r) in zip(y, d.iterrows()):
        ax1.plot([r.wilson_lo_pct, r.wilson_hi_pct], [yi, yi], color=BLUE, lw=2.4,
                 solid_capstyle="butt", zorder=2)
        ax1.plot([r.phylum_clustered_lo_pct, r.phylum_clustered_hi_pct], [yi - 0.22, yi - 0.22],
                 color=ORANGE, lw=1.6, solid_capstyle="butt", zorder=2)
        ax1.plot(r.agreement * 100, yi, "o", color=BLUE, ms=5.5, zorder=3)
    ax1.set_yticks(y); ax1.set_yticklabels(short)
    ax1.set_xlim(0, 100); ax1.set_ylim(-1.7, len(d) - 0.5)
    ax1.set_xlabel("direction agreement (%)")
    ax1.plot([], [], color=BLUE, lw=2.4, label="95 % Wilson")
    ax1.plot([], [], color=ORANGE, lw=1.6, label="95 % phylum-clustered")
    ax1.legend(frameon=False, loc="lower center", fontsize=7, ncol=2,
               bbox_to_anchor=(0.5, -0.02), handlelength=1.6, columnspacing=1.4)
    ax1.text(51.5, len(d) - 0.65, "chance", color="0.35", fontsize=7, va="top")
    ax1.set_title("a", loc="left", fontweight="bold", fontsize=10)

    for lab, col, mk in [("published endpoint (tie-corrected)", BLUE, "o"),
                         ("all orders testable in both", ORANGE, "s")]:
        g = s16[s16.endpoint == lab]
        ax2.plot(g.true_agreement * 100, g.power, marker=mk, color=col, lw=1.6, ms=5,
                 label=f"{lab.split(' (')[0]} (n = {int(g.n_orders.iloc[0])})")
        ax2.axvline(g.min_agreement_for_p05.iloc[0] * 100, color=col, ls=":", lw=1.1)
    ax2.axhline(0.80, color="0.35", ls="--", lw=1)
    ax2.text(60.4, 0.815, "80 % power", color="0.35", fontsize=7)
    ax2.set_xlabel("true direction agreement (%)")
    ax2.set_ylabel("power, two-sided exact binomial")
    ax2.set_ylim(0, 1.02); ax2.set_xlim(59, 81)
    ax2.legend(frameon=False, loc="upper left", fontsize=7)
    ax2.set_title("b", loc="left", fontweight="bold", fontsize=10)

    fig.tight_layout()
    fig.savefig(f"{SUP}/Figure_S4_agreement_power.pdf")
    fig.savefig(f"{SUP}/Figure_S4_agreement_power.png", dpi=600)
    plt.close(fig)
    print("wrote Figure_S4_agreement_power.{pdf,png}")


if __name__ == "__main__":
    main()
    make_figure_s4()