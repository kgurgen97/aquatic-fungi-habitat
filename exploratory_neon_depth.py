#!/usr/bin/env python3
"""POST HOC / EXPLORATORY, added 2026-10-21 after the pre-registered primary returned
NOT REPLICATED. Not part of NEON_VALIDATION_PLAN.md and not a substitute for it.

The question it asks: FunAqua's pairs carry a median of ~131-179 reads per sample, NEON's
a median of 31,131 - about 200x deeper. Within-pair rarefaction equalises power *inside*
each pair, but it cannot make "present" mean the same thing in both datasets. At ~130
reads, presence means "among the most abundant sequences in the sample"; at 31,000 it
means "detectable at all". If that mismatch drives the failed replication, then capping
NEON's pairs at FunAqua-like depths should raise the sign agreement.

This is a mechanism probe. Whatever it shows, the pre-registered answer stays NOT REPLICATED.
"""
import csv, json, re, collections
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

D = Path.home() / "research/datasets/neon2"
P = Path.home() / "research/projects/funaqua-habitat"
MIN_CONF, MIN_READS, MIN_DISCORDANT, SEED = 0.80, 100, 20, 20261021
WINDOW = 7
CAPS = [130, 300, 1000, 3000, 10000, 30000, None]

DATE = re.compile(r"(?<!\d)(20\d{6})(?!\d)")
SITE = re.compile(r"^([A-Z]{4})\b")


def parse_id(sid):
    m, s = DATE.search(sid), SITE.match(sid)
    if not m or not s:
        return None
    d = m.group(1)
    try:
        return s.group(1), pd.Timestamp(f"{d[:4]}-{d[4:6]}-{d[6:]}")
    except ValueError:
        return None


def load(kind):
    per_run = collections.defaultdict(collections.Counter)
    meta = {}
    for fp in sorted((D / kind).glob("*.csv")):
        with fp.open(newline="", encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh):
                sid = (r.get("dnaSampleID") or "").strip()
                if not sid or (r.get("kingdom") or "").strip() != "Fungi":
                    continue
                try:
                    if float(r.get("orderConfidence") or 0) < MIN_CONF:
                        continue
                    c = int(float(r.get("individualCount") or 0))
                except ValueError:
                    continue
                order = (r.get("order") or "").strip()
                pid = parse_id(sid)
                if not order or c <= 0 or pid is None:
                    continue
                run = (r.get("sequencerRunID") or "").strip()
                per_run[(sid, run)][((r.get("sequenceName") or "").strip() or f"{order}|anon", order)] += c
                meta[(sid, run)] = pid
    by = collections.defaultdict(list)
    for k in per_run:
        by[k[0]].append(k)
    best = {}
    for sid, keys in by.items():
        k = max(keys, key=lambda z: sum(per_run[z].values()))
        site, date = meta[k]
        best[sid] = dict(site=site, date=date, otu=per_run[k])
    return best


print("loading NEON tables ...", flush=True)
water, benthic = load("water"), load("benthic")
depth = lambda s: sum(s["otu"].values())

bidx = collections.defaultdict(list)
for sid, s in benthic.items():
    bidx[s["site"]].append(sid)
wbest = {}
for sid, s in water.items():
    k = (s["site"], s["date"])
    if k not in wbest or depth(s) > depth(water[wbest[k]]):
        wbest[k] = sid
pairs = []
for (site, wd), wsid in sorted(wbest.items()):
    cand = [b for b in bidx.get(site, []) if abs((benthic[b]["date"] - wd).days) <= WINDOW]
    if cand:
        pairs.append((wsid, max(cand, key=lambda b: depth(benthic[b]))))
print(f"{len(pairs)} candidate pairs")

orders = sorted({o for s in list(water.values()) + list(benthic.values()) for _, o in s["otu"]})
oidx = {o: i for i, o in enumerate(orders)}


def vec(s):
    ks = list(s["otu"])
    return (np.array([oidx[o] for _, o in ks]),
            np.array([s["otu"][k] for k in ks], dtype=np.int64))


fun = pd.read_csv(P / "results/primary_order.tsv", sep="\t")
fsig = fun[fun.q < 0.05][["order", "water_fraction"]]
print(f"FunAqua-significant orders (pre-registered single draw): {len(fsig)}")

# the draw-stable set from analysis_draw_stability.py: significant in >=95 % of 200 draws
st = pd.read_csv(P / "results/draw_stability.tsv", sep="\t")
NS = json.load((P / "results/draw_stability.json").open())["n_draws"]
stable = st[st.times_significant >= 0.95 * NS][["order", "mean_water_fraction"]].rename(
    columns={"mean_water_fraction": "water_fraction"})
print(f"FunAqua draw-stable orders (>=95 % of {NS} draws)      : {len(stable)}\n")

SETS = {"pre-registered single draw": fsig, "draw-stable set": stable}

print(f"{'NEON cap':>10s} {'pairs':>6s} {'med depth':>10s} {'tested':>7s} {'sig':>5s} "
      f"| {'comparison set':<27s} {'agreement':>16s} {'binom p':>9s} {'Spearman':>9s}")
rows = []
for cap in CAPS:
    rng = np.random.default_rng(SEED)
    W, S, mins = [], [], []
    for wsid, bsid in pairs:
        kw, vw = vec(water[wsid]); kb, vb = vec(benthic[bsid])
        if vw.sum() < MIN_READS or vb.sum() < MIN_READS:
            continue
        n = int(min(vw.sum(), vb.sum()))
        if cap is not None:
            n = min(n, cap)
        if n < MIN_READS:
            continue
        mins.append(n)
        ra = vw if vw.sum() <= n else rng.multivariate_hypergeometric(vw, n)
        rb = vb if vb.sum() <= n else rng.multivariate_hypergeometric(vb, n)
        W.append(np.bincount(kw, weights=ra, minlength=len(orders)) > 0)
        S.append(np.bincount(kb, weights=rb, minlength=len(orders)) > 0)
    W, S = np.array(W), np.array(S)
    wo = (W & ~S).sum(0); bo = (~W & S).sum(0); disc = wo + bo
    sel = np.where(disc >= MIN_DISCORDANT)[0]
    p = np.array([stats.binomtest(int(wo[j]), int(disc[j]), 0.5).pvalue for j in sel])
    q = multipletests(p, method="fdr_bh")[1]
    nf = pd.DataFrame({"order": [orders[j] for j in sel],
                       "neon_frac": wo[sel] / disc[sel], "q": q})
    first = True
    for sname, sdf in SETS.items():
        m = sdf.merge(nf, on="order")
        k = int(((m.water_fraction > .5) == (m.neon_frac > .5)).sum()) if len(m) else 0
        bp = stats.binomtest(k, len(m), 0.5).pvalue if len(m) else float("nan")
        rho = stats.spearmanr(m.water_fraction, m.neon_frac)[0] if len(m) > 2 else float("nan")
        lead = (f"{str(cap) if cap else 'none':>10s} {len(mins):>6d} {int(np.median(mins)):>10,d} "
                f"{len(sel):>7d} {int((q<0.05).sum()):>5d}") if first else " " * 42
        print(f"{lead} | {sname:<27s} {k:>5d}/{len(m):<3d} "
              f"{k/len(m) if len(m) else 0:>6.1%} {bp:>9.3g} {rho:>+9.3f}")
        rows.append(dict(cap=cap, comparison_set=sname, pairs=len(mins),
                         median_depth=int(np.median(mins)), tested=len(sel),
                         significant=int((q < 0.05).sum()), overlap=len(m), agree=k,
                         agreement=k / len(m) if len(m) else None, binom_p=bp, spearman=rho))
        first = False

out = pd.DataFrame(rows)
out.to_csv(P / "results/neon_depth_probe.tsv", sep="\t", index=False)
json.dump(rows, open(P / "results/neon_depth_probe.json", "w"), indent=1, default=str)
print("\nFunAqua's own pairs sit at a median of ~131 reads per sample, i.e. the top row's regime.")
print("wrote results/neon_depth_probe.tsv")
