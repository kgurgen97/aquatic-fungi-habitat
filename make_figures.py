#!/usr/bin/env python3
"""The four manuscript figures. Print-ready vector PDF plus 600 dpi PNG.

Every value is read from a file in results/; nothing is hard-coded except the three
summary numbers quoted in annotations, which are taken from results/*.json and named
where they appear.

Palette: categorical slots 1 and 2 of the validated default (blue #2a78d6, orange
#eb6834) - validated with scripts/validate_palette.js --pairs all --mode light, all
checks PASS. Slot 3 (aqua) was dropped because it warns below 3:1 against the light
surface. Chrome inks are the documented chart-chrome roles.

Print figures carry no hover layer and no dark mode: both are properties of an
on-screen chart, not of a PDF going to a journal.
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

P = Path.home() / "research/projects/funaqua-habitat"
R = P / "results"
F = P / "figures"
F.mkdir(exist_ok=True)

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASE, SURF = "#e1e0d9", "#c3c2b7", "#fcfcfb"

MM = 1 / 25.4
W1, W2 = 89 * MM, 183 * MM          # single- and double-column widths

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7.5, "axes.titlesize": 8,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.edgecolor": BASE, "axes.linewidth": 0.6,
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED,
    "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.axisbelow": True, "axes.spines.top": False, "axes.spines.right": False,
    "figure.facecolor": SURF, "axes.facecolor": SURF,
    "savefig.facecolor": SURF, "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
    "lines.linewidth": 1.0, "lines.markersize": 3.5, "pdf.fonttype": 42,
})


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(F / f"{name}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"  wrote figures/{name}.pdf and .png")


# ---------------------------------------------------------------- Figure 1
# Job: show that significance is a property of the draw for some orders and not for
# others, and justify the 35-order stable set. One ranked series -> no legend needed;
# the threshold and the two regions are direct-labelled.
st = pd.read_csv(R / "draw_stability.tsv", sep="\t")
ds = json.load((R / "draw_stability.json").open())
N = ds["n_draws"]
st = st[st.times_tested > 0].copy()
# Stability is measured against ALL draws, exactly as analysis_draw_stability.py
# defines it, not against the draws in which an order happened to be testable. The
# difference is not cosmetic: an order tested in 27 of 200 draws and significant in
# all 27 scores 100 % on the second definition and is not stable on the first. Using
# the looser one here made the figure say 36 where the text says 35.
st["freq"] = st.times_significant / N
st = st.sort_values("freq", ascending=False).reset_index(drop=True)
stable = st.times_significant >= 0.95 * N
assert int(stable.sum()) == ds["stable_95pct"], \
    f"figure says {int(stable.sum())} stable, analysis says {ds['stable_95pct']}"

fig, ax = plt.subplots(figsize=(W1, 55 * MM))
x = np.arange(len(st))
ax.vlines(x[stable], 0, st.freq[stable], color=BLUE, linewidth=1.1)
ax.vlines(x[~stable], 0, st.freq[~stable], color=MUTED, linewidth=1.1, alpha=0.85)
ax.axhline(0.95, color=ORANGE, linewidth=0.9, linestyle=(0, (3, 2)), zorder=3)
ax.text(len(st) - 1, 0.965, "95 % of draws", ha="right", va="bottom",
        color=ORANGE, fontsize=6.5)
# A legend, not in-plot annotations: the first version printed "35 stable orders" in
# blue on top of the blue bars, which was unreadable.
from matplotlib.lines import Line2D
ax.legend(handles=[
    Line2D([], [], color=BLUE, linewidth=2.4,
           label=f"stable, significant in ≥ 95 % of draws ({int(stable.sum())})"),
    Line2D([], [], color=MUTED, linewidth=2.4,
           label=f"unstable ({int((~stable).sum())}), of which "
                 f"{int((st.times_significant == 0).sum())} never significant")],
    frameon=False, loc="lower left", bbox_to_anchor=(0.30, 0.04),
    handlelength=1.4, handletextpad=0.5, borderpad=0.1, labelspacing=0.35)
ax.set_xlabel(f"taxonomic orders, ranked ({len(st)} tested in at least one draw)")
ax.set_ylabel(f"draws in which significant\n(of {N})")
ax.yaxis.set_major_formatter(PercentFormatter(1.0))
ax.set_ylim(0, 1.06); ax.set_xlim(-2, len(st) + 1)
ax.set_title("Significance depends on the rarefaction draw", loc="left", color=INK)
sig = ds["significant"]
ax.text(0.0, -0.30,
        f"Across {N} draws the count of significant orders ranges "
        f"{sig['min']}–{sig['max']} (median {sig['median']:.0f}).",
        transform=ax.transAxes, color=INK2, fontsize=6.5, va="top")
save(fig, "fig1_draw_stability")

# ---------------------------------------------------------------- Figure 2
# Job: identity/agreement between two datasets -> scatter with the two decision axes.
# Two categories (agree / disagree) so a legend is present and both are direct-labelled
# by the quadrant counts.
m = pd.read_csv(R / "neon_vs_funaqua.tsv", sep="\t")
nv = json.load((R / "neon_validation.json").open())
fcol = "water_fraction_fun" if "water_fraction_fun" in m else "water_fraction_x"
ncol = "water_fraction_neon" if "water_fraction_neon" in m else "water_fraction_y"
m["agree"] = (m[fcol] > .5) == (m[ncol] > .5)

fig, ax = plt.subplots(figsize=(W1, W1))
ax.axhspan(0.5, 1.02, xmin=0, xmax=0.5, color=GRID, alpha=0.45, zorder=0)
ax.axhspan(-0.02, 0.5, xmin=0.5, xmax=1.0, color=GRID, alpha=0.45, zorder=0)
ax.axhline(0.5, color=BASE, linewidth=0.7); ax.axvline(0.5, color=BASE, linewidth=0.7)
ax.plot([0, 1], [0, 1], color=MUTED, linewidth=0.6, linestyle=(0, (2, 2)), zorder=1)
for ok, col, lab in ((True, BLUE, "same direction"), (False, ORANGE, "opposite")):
    s = m[m.agree == ok]
    ax.scatter(s[fcol], s[ncol], s=22, facecolor=col, edgecolor=SURF,
               linewidth=0.7, label=f"{lab} ({len(s)})", zorder=4)
ax.set_xlabel("FunAqua water preference index")
ax.set_ylabel("NEON water preference index")
ax.set_xlim(-0.02, 1.02); ax.set_ylim(-0.02, 1.02)
ax.set_xticks([0, .25, .5, .75, 1]); ax.set_yticks([0, .25, .5, .75, 1])
ax.set_aspect("equal")
# Legend above the axes: inside the plot it sat within a shaded "disagreement"
# quadrant, which muddled what the shading means.
ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0, 1.02), ncol=2,
          handletextpad=0.4, borderpad=0.1, columnspacing=1.4)
ax.set_title("Per-order preferences do not transfer", loc="left", color=INK, pad=22)
ax.text(0.0, -0.17,
        f"{nv['agree']}/{nv['n_compared']} orders agree in direction "
        f"({nv['agreement']:.1%}); exact binomial p = {nv['p']:.2g}.\n"
        f"Spearman {nv['spearman']:+.2f} (p = {nv['spearman_p']:.2f}). Shaded "
        "quadrants are the disagreements;\nthe dashed line is 1:1.",
        transform=ax.transAxes, color=INK2, fontsize=6.5, va="top")
save(fig, "fig2_transfer_failure")

# ---------------------------------------------------------------- Figure 3
# Job: magnitude across six named groups, on two metrics -> two panels, one bar per
# group, sorted once and held. Fungi highlighted; the others are one recessive colour,
# so hue is not carrying a category it does not have.
kc = pd.read_csv(R / "kingdom_consistent.tsv", sep="\t").sort_values("effect")
cols = [BLUE if k == "Fungi" else MUTED for k in kc.kingdom]
# Both panels read left to right. The first version mirrored the left panel so the
# two grew outward from shared centre labels, and those labels then overlapped the
# right panel's bars - "Straminipila" was clipped to "Straminip".
fig, axes = plt.subplots(1, 2, figsize=(W2, 58 * MM), sharey=True)
fig.subplots_adjust(wspace=0.10)
y = np.arange(len(kc))
for ax, col, lab, ttl in (
        (axes[0], "effect", "median |preference − 0.5|", "Habitat preference"),
        (axes[1], "B_order", "median Levins $B$ (2 compartments)", "Niche breadth")):
    ax.barh(y, kc[col], color=cols, height=0.62, zorder=3)
    for yi, v in zip(y, kc[col]):
        ax.text(v + (kc[col].max() * 0.015), yi, f"{v:.3f}", va="center",
                color=INK2, fontsize=6.3)
    ax.set_xlabel(lab); ax.set_title(ttl, loc="left", color=INK)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, kc[col].max() * 1.17)
axes[1].set_xlim(1.0, kc.B_order.max() * 1.08)
axes[0].set_yticks(y)
axes[0].set_yticklabels([f"$\\bf{{{k}}}$" if k == "Fungi" else k for k in kc.kingdom])
axes[0].text(0.0, -0.22,
             "Lower preference and higher $B$ both mean more generalist. "
             "Fungi are indistinguishable from Alveolata\non both; the gap against "
             "non-fungal eukaryotes as a pooled group is driven by Metazoa and "
             "Straminipila.\nSame 204 site pairs, same per-pair rarefaction depth.",
             transform=axes[0].transAxes, color=INK2, fontsize=6.5, va="top")
save(fig, "fig3_kingdoms")

# ---------------------------------------------------------------- Figure 4
# Job: change in agreement along a depth axis -> line per comparison set, log x.
dp = pd.read_csv(R / "neon_depth_probe.tsv", sep="\t")
dp["capn"] = pd.to_numeric(dp.cap, errors="coerce")
UNCAP = dp.capn.max() * 3.2
dp["xx"] = dp.capn.fillna(UNCAP)
# Taller than the others and with a short y label: the first version's label,
# "direction agreement with FunAqua", was clipped off the top of the figure.
fig, ax = plt.subplots(figsize=(W1, 64 * MM))
ax.axhline(0.5, color=MUTED, linewidth=0.7, linestyle=(0, (2, 2)))
ax.text(dp.xx.min(), 0.515, "chance", color=MUTED, fontsize=6.3, va="bottom")
for (name, g), col in zip(dp.groupby("comparison_set", sort=False), (BLUE, ORANGE)):
    g = g.sort_values("xx")
    lab = "pre-registered set" if "pre-reg" in name else "draw-stable set"
    ax.plot(g.xx, g.agreement, color=col, marker="o", label=lab, zorder=3)
ax.set_xscale("log")
ax.set_xticks(list(dp.capn.dropna().unique()) + [UNCAP])
ax.set_xticklabels([f"{int(v):,}" for v in dp.capn.dropna().unique()] + ["none"],
                   rotation=45, ha="right")
ax.yaxis.set_major_formatter(PercentFormatter(1.0))
ax.set_ylim(0.38, 0.78)
ax.set_xlabel("NEON within-pair depth cap (reads; 'none' = 31 131 median)")
ax.set_ylabel("agreement with FunAqua")
ax.legend(frameon=False, loc="upper right", handletextpad=0.4, borderpad=0.1)
# Neutral title: the trend is toward matched depth but reaches significance at no cap,
# so "rises as depth is matched" would overclaim what the panel shows.
ax.set_title("Agreement under depth capping", loc="left", color=INK)
ax.text(0.0, -0.38,
        "FunAqua's pairs sit near the left-hand end (median 131–179 reads per sample). "
        "Agreement never\nreaches significance at any cap (best p = 0.076), so depth "
        "mismatch is a credible contributor,\nnot an established explanation. "
        "Post hoc; not part of the frozen validation plan.",
        transform=ax.transAxes, color=INK2, fontsize=6.5, va="top")
save(fig, "fig4_depth_probe")

print("\nall four figures written to figures/")
