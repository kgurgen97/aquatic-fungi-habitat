#!/usr/bin/env python3
"""The manuscript figures, in journal style. Vector PDF plus 600 dpi PNG.

Every value is read from a file in results/; nothing is hard-coded except summary
numbers quoted in a label, which are read from results/*.json and named where used.

Style decisions, all deliberate:
  * Liberation Sans, metric-compatible with Arial, at 6-7 pt - the size these are
    printed at. Not DejaVu, which reads as "matplotlib default" on sight.
  * No titles and no explanatory text inside the panels. A journal figure carries its
    caption in the text; captions live in MANUSCRIPT.md under "Figure captions".
  * Panel letters a, b in bold outside the axes, which is where a copy-editor expects
    to find them.
  * Left and bottom spines only, outward ticks, horizontal hairline grid or none.
  * Thin marks: 0.7 pt lines, 2.6 pt markers.

Palette: categorical slots 1 and 2 of the validated default (blue #2a78d6, orange
#eb6834), checked with validate_palette.js --pairs all --mode light, all checks PASS.
The aqua slot was dropped because it warns below 3:1 against the light surface.

Print output carries no hover layer and no dark mode: both belong to an on-screen
chart, not to a PDF going to a journal.
"""
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import PercentFormatter

P = Path.home() / "research/projects/funaqua-habitat"
R, F = P / "results", P / "figures"
F.mkdir(exist_ok=True)

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, MUTED = "#111111", "#4a4a4a", "#8a8a8a"
GRID, SPINE, SURF = "#e8e8e4", "#4a4a4a", "#ffffff"

MM = 1 / 25.4
W1, W2 = 89 * MM, 183 * MM

FAM = "Liberation Sans" if any(f.name == "Liberation Sans" for f in fm.fontManager.ttflist) \
      else "DejaVu Sans"
if FAM != "Liberation Sans":
    print("  WARNING: Liberation Sans not found, falling back to DejaVu Sans")

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": [FAM],
    "font.size": 6.5, "axes.labelsize": 7, "axes.titlesize": 7,
    "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
    "axes.edgecolor": SPINE, "axes.linewidth": 0.6,
    "axes.labelcolor": INK, "text.color": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "xtick.color": SPINE, "ytick.color": SPINE,
    "xtick.labelcolor": INK, "ytick.labelcolor": INK,
    "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.pad": 2.0, "ytick.major.pad": 2.0,
    "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.5,
    "axes.axisbelow": True,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "savefig.bbox": "tight", "savefig.pad_inches": 0.015,
    "lines.linewidth": 0.9, "lines.markersize": 2.8,
    "legend.frameon": False, "legend.handlelength": 1.2,
    "legend.handletextpad": 0.45, "legend.borderpad": 0.0,
    "legend.labelspacing": 0.3, "legend.columnspacing": 1.1,
    "pdf.fonttype": 42, "ps.fonttype": 42,
})


def panel(ax, letter, dx=-0.16, dy=1.04):
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=8,
            fontweight="bold", va="bottom", ha="left", color=INK)


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(F / f"{name}.{ext}", dpi=600 if ext == "png" else None)
    plt.close(fig)
    print(f"  {name}.pdf / .png")


# ---------------------------------------------------------------- Figure 1
st = pd.read_csv(R / "draw_stability.tsv", sep="\t")
ds = json.load((R / "draw_stability.json").open())
N = ds["n_draws"]
st = st[st.times_tested > 0].copy()
# Stability is measured against ALL draws, as analysis_draw_stability.py defines it,
# not against the draws in which an order happened to be testable. Using the looser
# definition made an earlier version say 36 where the manuscript says 35.
st["freq"] = st.times_significant / N
st = st.sort_values("freq", ascending=False).reset_index(drop=True)
n_stable = int((st.times_significant >= 0.95 * N).sum())
assert n_stable == ds["stable_95pct"], \
    f"figure says {n_stable} stable, analysis says {ds['stable_95pct']}"

fig, ax = plt.subplots(figsize=(W1, 46 * MM))
x = np.arange(len(st))
# A light wash under a crisp step, rather than two solid blocks: the solid version
# read as a slab of colour and swamped the curve it is meant to describe.
ax.fill_between(x[:n_stable], st.freq[:n_stable], step="post",
                color=BLUE, alpha=0.20, linewidth=0)
ax.fill_between(x[n_stable - 1:], st.freq[n_stable - 1:], step="post",
                color=MUTED, alpha=0.18, linewidth=0)
ax.step(x[:n_stable], st.freq[:n_stable], where="post", color=BLUE, linewidth=0.9)
ax.step(x[n_stable - 1:], st.freq[n_stable - 1:], where="post",
        color=MUTED, linewidth=0.9)
ax.axhline(0.95, color=ORANGE, linewidth=0.7, dashes=(3.5, 2.2), zorder=4)
ax.axvline(n_stable, color=INK, linewidth=0.5, dashes=(2, 2), zorder=4)
ax.text(len(st) - 2, 0.935, "95 %", ha="right", va="top", color=ORANGE, fontsize=6)
ax.annotate(f"{n_stable} stable", xy=(n_stable, 0.42), xytext=(n_stable + 7, 0.42),
            fontsize=6.5, color=INK, va="center",
            arrowprops=dict(arrowstyle="-", color=INK, linewidth=0.5,
                            shrinkA=0, shrinkB=1))
ax.set_xlabel("taxonomic order, ranked")
ax.set_ylabel(f"draws significant (of {N})")
ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
ax.set_ylim(0, 1.0); ax.set_xlim(0, len(st))
ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax.set_xticks([0, 25, 50, 75, 100])
save(fig, "fig1_draw_stability")

# ---------------------------------------------------------------- Figure 2
m = pd.read_csv(R / "neon_vs_funaqua.tsv", sep="\t")
nv = json.load((R / "neon_validation.json").open())
fcol = "water_fraction_fun" if "water_fraction_fun" in m else "water_fraction_x"
ncol = "water_fraction_neon" if "water_fraction_neon" in m else "water_fraction_y"
m["agree"] = (m[fcol] > .5) == (m[ncol] > .5)

fig, ax = plt.subplots(figsize=(62 * MM, 62 * MM))
ax.plot([0, 1], [0, 1], color=MUTED, linewidth=0.5, dashes=(2.5, 2), zorder=1)
ax.axhline(0.5, color=INK, linewidth=0.5, zorder=2)
ax.axvline(0.5, color=INK, linewidth=0.5, zorder=2)
ax.scatter(m.loc[m.agree, fcol], m.loc[m.agree, ncol], s=13, facecolor=BLUE,
           edgecolor="white", linewidth=0.45, zorder=5)
ax.scatter(m.loc[~m.agree, fcol], m.loc[~m.agree, ncol], s=13, facecolor=ORANGE,
           edgecolor="white", linewidth=0.45, zorder=5)
# Quadrant counts, direct-labelled - no legend box needed.
q = dict(ha="center", va="center", fontsize=6.5)
ax.text(0.25, 0.97, f"{int(((m[fcol] <= .5) & (m[ncol] > .5)).sum())}", color=ORANGE, **q)
ax.text(0.75, 0.97, f"{int(((m[fcol] > .5) & (m[ncol] > .5)).sum())}", color=BLUE, **q)
ax.text(0.25, 0.03, f"{int(((m[fcol] <= .5) & (m[ncol] <= .5)).sum())}", color=BLUE, **q)
ax.text(0.75, 0.03, f"{int(((m[fcol] > .5) & (m[ncol] <= .5)).sum())}", color=ORANGE, **q)
ax.set_xlabel("water preference, FunAqua")
ax.set_ylabel("water preference, NEON")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xticks([0, .25, .5, .75, 1]); ax.set_yticks([0, .25, .5, .75, 1])
ax.set_aspect("equal")
save(fig, "fig2_transfer_failure")

# ---------------------------------------------------------------- Figure 3
kc = pd.read_csv(R / "kingdom_consistent.tsv", sep="\t").sort_values("effect")
cols = [BLUE if k == "Fungi" else MUTED for k in kc.kingdom]
fig, axes = plt.subplots(1, 2, figsize=(W2 * 0.82, 44 * MM), sharey=True)
fig.subplots_adjust(wspace=0.08)
y = np.arange(len(kc))
for ax, col, lab, lo in ((axes[0], "effect", "median |preference − 0.5|", 0),
                         (axes[1], "B_order", "median Levins $B$", 1.0)):
    ax.barh(y, kc[col] - lo, left=lo, color=cols, height=0.58, linewidth=0, zorder=3)
    span = kc[col].max() - lo
    for yi, v in zip(y, kc[col]):
        ax.text(v + span * 0.025, yi, f"{v:.3f}", va="center", color=INK2, fontsize=5.8)
    ax.set_xlabel(lab)
    ax.set_xlim(lo, lo + span * 1.22)
    ax.xaxis.grid(True, linewidth=0.45, color=GRID)
    ax.set_ylim(-0.6, len(kc) - 0.4)
axes[0].set_yticks(y)
axes[0].set_yticklabels([f"$\\bf{{{k}}}$" if k == "Fungi" else k for k in kc.kingdom])
for a in axes:
    a.tick_params(axis="y", length=0)   # the shared y axis left tick stubs on panel b
panel(axes[0], "a", dx=-0.42); panel(axes[1], "b", dx=-0.06)
save(fig, "fig3_kingdoms")

# ---------------------------------------------------------------- Figure 4
dp = pd.read_csv(R / "neon_depth_probe.tsv", sep="\t")
dp["capn"] = pd.to_numeric(dp.cap, errors="coerce")
UNCAP = dp.capn.max() * 3.0
dp["xx"] = dp.capn.fillna(UNCAP)
fig, ax = plt.subplots(figsize=(W1, 46 * MM))
ax.axhline(0.5, color=MUTED, linewidth=0.5, dashes=(2.5, 2))
ax.text(dp.xx.min() * 0.92, 0.504, "chance", color=MUTED, fontsize=5.8, va="bottom")
for (name, g), col in zip(dp.groupby("comparison_set", sort=False), (BLUE, ORANGE)):
    g = g.sort_values("xx")
    lab = "pre-registered" if "pre-reg" in name else "draw-stable"
    ax.plot(g.xx, g.agreement, color=col, marker="o", markeredgecolor="white",
            markeredgewidth=0.4, label=lab, zorder=3)
ax.set_xscale("log")
ticks = list(dp.capn.dropna().unique()) + [UNCAP]
ax.set_xticks(ticks)
ax.set_xticklabels([f"{int(v):,}" for v in dp.capn.dropna().unique()] + ["none"])
ax.tick_params(axis="x", which="minor", length=0)
ax.yaxis.set_major_formatter(PercentFormatter(1.0, decimals=0))
ax.set_ylim(0.47, 0.72); ax.set_yticks([0.5, 0.6, 0.7])   # was 0.42-0.76, half empty
ax.set_xlabel("NEON within-pair depth cap (reads)")
ax.set_ylabel("direction agreement")
ax.legend(loc="upper right", bbox_to_anchor=(1.0, 1.06))
save(fig, "fig4_depth_probe")

# ---------------------------------------------------------------- Figure 5
sim = json.load((R / "simulation.json").open())
ss = pd.read_csv(R / "sim_self_selection.tsv", sep="\t")
sc = pd.read_csv(R / "sim_community.tsv", sep="\t")
ds2 = json.load((R / "draw_stability.json").open())["significant"]

fig, axes = plt.subplots(1, 3, figsize=(W2, 46 * MM))
fig.subplots_adjust(wspace=0.42)

# (a) what the index means depends on depth - the algebra
D = np.geomspace(20, 2e5, 300)
A_S = 1e-4
for ratio, col, dash in ((2.0, MUTED, (2, 2)), (4.0, BLUE, None), (8.0, ORANGE, None)):
    pw, ps = 1 - np.exp(-A_S * ratio * D), 1 - np.exp(-A_S * D)
    num = pw * (1 - ps)
    axes[0].plot(D, num / (num + ps * (1 - pw)), color=col, linewidth=1.0,
                 dashes=dash or (None, None), label=f"{ratio:.0f}×")
for d, lab in ((211, "FunAqua"), (31131, "NEON")):
    axes[0].axvline(d, color=INK, linewidth=0.5, dashes=(2, 2))
    axes[0].text(d, 0.52, lab, rotation=90, fontsize=5.6, color=INK2,
                 ha="right", va="bottom")
axes[0].set_xscale("log")
axes[0].set_xlabel("rarefaction depth (reads)")
axes[0].set_ylabel("expected preference index")
axes[0].set_ylim(0.45, 1.02)
# lower left, not lower right: the rotated depth labels sit on the right and the
# legend title collided with "NEON"
axes[0].legend(loc="lower left", title="true ratio", title_fontsize=5.8)

# (b) the testability threshold pins lambda and moves the abundance
axb = axes[1]
axb.plot(ss.depth, ss.med_abundance, color=BLUE, marker="o",
         markeredgecolor="white", markeredgewidth=0.4)
axb.set_xscale("log"); axb.set_yscale("log")
axb.set_xlabel("rarefaction depth (reads)")
axb.set_ylabel("median abundance of\nthe testable set")
axt = axb.twiny()      # a second x axis is fine; this is not a second y scale
axt.remove()
axb.text(0.04, 0.10,
         f"median λ holds at\n{ss.med_lambda.min():.2f}–{ss.med_lambda.max():.2f}",
         transform=axb.transAxes, fontsize=6, color=INK2, va="bottom")

# (c) simulated against observed spread
axc = axes[2]
bins = np.arange(min(sc.n_significant.min(), ds2["min"]) - 0.5,
                 max(sc.n_significant.max(), ds2["max"]) + 1.5, 1)
axc.hist(sc.n_significant, bins=bins, color=BLUE, alpha=0.75, linewidth=0,
         label="simulated")
axc.axvspan(ds2["min"], ds2["max"], color=ORANGE, alpha=0.22, linewidth=0,
            label="observed range")
axc.axvline(ds2["median"], color=ORANGE, linewidth=1.0)
axc.set_xlabel("significant orders per draw")
axc.set_ylabel("draws")
axc.legend(loc="upper left")

for a, L in zip(axes, "abc"):
    panel(a, L, dx=-0.30)
save(fig, "fig5_simulation")

# Counted, not typed: this message said "four" for a while after Figure 5 landed.
print(f"\n{len(list(F.glob('fig*.png')))} figures in figures/, set in {FAM}")
