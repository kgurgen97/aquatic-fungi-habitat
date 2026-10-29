#!/usr/bin/env python3
"""NEON external validation, exactly as NEON_VALIDATION_PLAN.md specifies.

sha256 of the plan this implements:
  85fd8aa2c9f8a14dc30c9bded5858d4f39dd0d6801f01753531fda0ef75fb939
"""
import collections, csv, json, re, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
D = NEON_DIR
P = REPO_ROOT
(P / "results").mkdir(exist_ok=True)

MIN_CONF = 0.80          # plan 3.2
MIN_READS = 100          # plan 4
MIN_DISCORDANT = 20      # plan 5
SEED = 20261021          # plan 5
EXACT_DATE_MIN_PAIRS = 60   # plan 4, declared fallback trigger
WINDOW_DAYS = 7
CONF_SENS = (0.5, 0.8, 0.9)
READ_SENS = (50, 100, 200, 400)

DATE = re.compile(r"(?<!\d)(20\d{6})(?!\d)")
SITE = re.compile(r"^([A-Z]{4})\b")


def parse_id(sid):
    m, s = DATE.search(sid), SITE.match(sid)
    if not m or not s:
        return None
    d = m.group(1)
    try:
        ts = pd.Timestamp(f"{d[:4]}-{d[4:6]}-{d[6:]}")
    except ValueError:
        return None
    return s.group(1), ts


def load(kind, min_conf):
    """-> {sample: {'site','date','run','otu':Counter}}  keeping, per sample, the deepest run."""
    per_run = collections.defaultdict(lambda: collections.Counter())
    meta, dropped = {}, collections.Counter()
    files = sorted((D / kind).glob("*.csv"))
    for fp in files:
        with fp.open(newline="", encoding="utf-8", errors="replace") as fh:
            for r in csv.DictReader(fh):
                sid = (r.get("dnaSampleID") or "").strip()
                if not sid:
                    dropped["no sample id"] += 1; continue
                if (r.get("kingdom") or "").strip() != "Fungi":
                    dropped["not Fungi"] += 1; continue
                try:
                    if float(r.get("orderConfidence") or 0) < min_conf:
                        dropped["order confidence below cutoff"] += 1; continue
                except ValueError:
                    dropped["unparseable confidence"] += 1; continue
                order = (r.get("order") or "").strip()
                if not order:
                    dropped["no order name"] += 1; continue
                try:
                    c = int(float(r.get("individualCount") or 0))
                except ValueError:
                    dropped["unparseable count"] += 1; continue
                if c <= 0:
                    dropped["zero count"] += 1; continue
                pid = parse_id(sid)
                if pid is None:
                    dropped["unparseable sample id"] += 1; continue
                run = (r.get("sequencerRunID") or "").strip()
                otu = (r.get("sequenceName") or "").strip() or f"{order}|anon"
                per_run[(sid, run)][(otu, order)] += c
                meta[(sid, run)] = pid

    # one run per sample: the deepest (plan 4, same rule as for replicates)
    best, multi = {}, 0
    by_sample = collections.defaultdict(list)
    for key in per_run:
        by_sample[key[0]].append(key)
    for sid, keys in by_sample.items():
        if len(keys) > 1:
            multi += 1
        k = max(keys, key=lambda z: sum(per_run[z].values()))
        site, date = meta[k]
        best[sid] = dict(site=site, date=date, run=k[1], otu=per_run[k])
    print(f"  {kind:8s} files={len(files):5d}  samples={len(best):5d}  "
          f"(samples spanning >1 run, deepest kept: {multi})")
    for k, v in dropped.most_common():
        print(f"      dropped rows, {k}: {v:,}")
    return best


def build_pairs(water, benthic, window):
    """Same site; same date, or |Δ| <= window days. Deepest sample on each side."""
    bidx = collections.defaultdict(list)
    for sid, s in benthic.items():
        bidx[s["site"]].append(sid)
    depth = lambda s: sum(s["otu"].values())
    # candidate water sample per (site, date): the deepest
    wbest = {}
    for sid, s in water.items():
        k = (s["site"], s["date"])
        if k not in wbest or depth(s) > depth(water[wbest[k]]):
            wbest[k] = sid
    pairs = []
    for (site, wdate), wsid in sorted(wbest.items()):
        cand = [b for b in bidx.get(site, [])
                if abs((benthic[b]["date"] - wdate).days) <= window]
        if not cand:
            continue
        bsid = max(cand, key=lambda b: depth(benthic[b]))
        pairs.append((site, wdate, wsid, bsid,
                      abs((benthic[bsid]["date"] - wdate).days)))
    return pairs


def test_pairs(pairs, water, benthic, min_reads, seed=SEED):
    """Within-pair rarefaction at OTU level, presence per order, exact paired test."""
    orders = sorted({o for s in list(water.values()) + list(benthic.values())
                     for _, o in s["otu"]})
    oidx = {o: i for i, o in enumerate(orders)}
    rng = np.random.default_rng(seed)

    def vec(s):
        ks = list(s["otu"])
        return (np.array([oidx[o] for _, o in ks]),
                np.array([s["otu"][k] for k in ks], dtype=np.int64))

    def rarefy(c, n):
        return c if c.sum() <= n else rng.multivariate_hypergeometric(c, int(n))

    kept, W, S, mins = [], [], [], []
    for site, wdate, wsid, bsid, dd in pairs:
        cw, cb = water[wsid], benthic[bsid]
        kw, vw = vec(cw); kb, vb = vec(cb)
        if vw.sum() < min_reads or vb.sum() < min_reads:
            continue
        n = int(min(vw.sum(), vb.sum())); mins.append(n)
        rw, rb = rarefy(vw, n), rarefy(vb, n)
        W.append(np.bincount(kw, weights=rw, minlength=len(orders)) > 0)
        S.append(np.bincount(kb, weights=rb, minlength=len(orders)) > 0)
        kept.append((site, wdate, wsid, bsid, dd, n))
    if not kept:
        return pd.DataFrame(), kept, orders, 0
    W, S = np.array(W), np.array(S)

    rows = []
    for j, o in enumerate(orders):
        wo = int((W[:, j] & ~S[:, j]).sum())
        bo = int((~W[:, j] & S[:, j]).sum())
        disc = wo + bo
        if disc < MIN_DISCORDANT:
            continue
        bt = stats.binomtest(wo, disc, 0.5)
        lo, hi = bt.proportion_ci(0.95, method="wilson")
        rows.append(dict(order=o, water_only=wo, benthic_only=bo,
                         both=int((W[:, j] & S[:, j]).sum()), discordant=disc,
                         water_fraction=wo / disc, ci_lo=lo, ci_hi=hi, p=bt.pvalue))
    r = pd.DataFrame(rows)
    if len(r):
        r["q"] = multipletests(r.p, method="fdr_bh")[1]
        r["prefers"] = np.where(r.q >= 0.05, "—",
                                np.where(r.water_fraction > .5, "water", "benthic"))
        r = r.sort_values("q").reset_index(drop=True)
    return r, kept, orders, int(np.median(mins))


# ----------------------------------------------------------------- load
print("=== loading NEON ITS tables (kingdom=Fungi, orderConfidence >= %.2f)" % MIN_CONF)
water = load("water", MIN_CONF)
benthic = load("benthic", MIN_CONF)
if not water or not benthic:
    sys.exit("no samples loaded — is the download complete?")

# ----------------------------------------------------------------- pairing (plan 4)
p_exact = build_pairs(water, benthic, 0)
q_exact = [p for p in p_exact
           if sum(water[p[2]]["otu"].values()) >= MIN_READS
           and sum(benthic[p[3]]["otu"].values()) >= MIN_READS]
print(f"\nsame-site same-date pairs: {len(p_exact)}  qualifying (both >= {MIN_READS} reads): {len(q_exact)}")

if len(q_exact) >= EXACT_DATE_MIN_PAIRS:
    window, rule = 0, "same site, same calendar date"
else:
    window, rule = WINDOW_DAYS, f"same site, within +/-{WINDOW_DAYS} days"
    print(f"  below the declared trigger of {EXACT_DATE_MIN_PAIRS} -> primary becomes the "
          f"+/-{WINDOW_DAYS}-day rule, as pre-specified in plan section 4")
pairs = build_pairs(water, benthic, window)
print(f"PRIMARY pairing rule: {rule}  -> {len(pairs)} candidate pairs")

# ----------------------------------------------------------------- primary (plan 5)
neon, kept, orders, med = test_pairs(pairs, water, benthic, MIN_READS)
print(f"\nqualifying pairs: {len(kept)}   median within-pair rarefaction depth: {med:,} reads")
print(f"orders present: {len(orders)}   orders tested (>= {MIN_DISCORDANT} discordant): {len(neon)}")
sig = neon[neon.q < 0.05] if len(neon) else neon
print(f"orders significant in NEON alone (BH q<0.05): {len(sig)}"
      + (f"   water {int((sig.prefers=='water').sum())} / benthic {int((sig.prefers=='benthic').sum())}"
         if len(sig) else ""))
neon.to_csv(P / "results/neon_order.tsv", sep="\t", index=False)
pd.DataFrame(kept, columns=["site", "water_date", "water_sample", "benthic_sample",
                            "day_offset", "rarefied_depth"]).to_csv(
    P / "results/neon_pairs.tsv", sep="\t", index=False)

# ----------------------------------------------------------------- primary endpoint (plan 6)
fun = pd.read_csv(P / "results/primary_order.tsv", sep="\t")
fsig = fun[fun.q < 0.05][["order", "water_fraction", "prefers", "discordant"]]
m = fsig.merge(neon[["order", "water_fraction", "discordant", "q", "prefers"]],
               on="order", suffixes=("_fun", "_neon"))
print(f"\n=== PRIMARY VALIDATION ENDPOINT (plan section 6)")
print(f"FunAqua-significant orders: {len(fsig)}   also testable in NEON: {len(m)}")

res = {"pairing_rule": rule, "pairs": len(kept), "median_depth": med,
       "orders_tested": int(len(neon)), "neon_significant": int(len(sig)),
       "funaqua_significant": int(len(fsig)), "overlap_testable": int(len(m))}
if len(m):
    m["agree"] = (m.water_fraction_fun > .5) == (m.water_fraction_neon > .5)
    k, n = int(m.agree.sum()), len(m)
    bt = stats.binomtest(k, n, 0.5)
    lo, hi = bt.proportion_ci(0.95, method="wilson")
    print(f"\nsign agreement: {k}/{n} = {k/n:.1%}   95% CI [{lo:.1%}, {hi:.1%}]"
          f"   exact binomial p = {bt.pvalue:.3g}")
    verdict = ("REPLICATED" if k / n > .5 and bt.pvalue < .05 else
               "CONTRADICTED" if k / n < .5 and bt.pvalue < .05 else
               "NOT REPLICATED (agreement not significant)")
    print(f"PRE-COMMITTED READING: {verdict}")
    res.update(agree=k, n_compared=n, agreement=k / n, ci=[lo, hi],
               p=bt.pvalue, verdict=verdict)

    rho, prho = stats.spearmanr(m.water_fraction_fun, m.water_fraction_neon)
    print(f"\nsecondary 1 — Spearman of the preference index, FunAqua vs NEON, "
          f"over these {len(m)} orders: {rho:+.3f} (p = {prho:.3g})")
    res.update(spearman=rho, spearman_p=prho)

    both = fun[["order", "water_fraction"]].merge(
        neon[["order", "water_fraction"]], on="order", suffixes=("_fun", "_neon"))
    rho2, prho2 = stats.spearmanr(both.water_fraction_fun, both.water_fraction_neon)
    print(f"              over all {len(both)} orders testable in both: {rho2:+.3f} (p = {prho2:.3g})")
    res.update(spearman_all=rho2, spearman_all_p=prho2, n_all=len(both))

    m = m.sort_values("water_fraction_fun", ascending=False)
    print("\nper-order comparison (water_fraction; >0.5 = prefers water):")
    print(m.rename(columns={"water_fraction_fun": "funaqua", "water_fraction_neon": "neon",
                            "discordant_fun": "disc_fun", "discordant_neon": "disc_neon",
                            "prefers_fun": "fun_pref", "prefers_neon": "neon_pref"})
          [["order", "funaqua", "neon", "agree", "disc_fun", "disc_neon", "fun_pref", "neon_pref"]]
          .to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    m.to_csv(P / "results/neon_vs_funaqua.tsv", sep="\t", index=False)

    # plan 7.3 — the named lineages
    named = {"zoosporic": ["Zygophlyctidales", "Chytridiales", "Zygorhizidiales", "Rozellales"],
             "terrestrial saprotroph": ["Sordariales", "Microascales", "Zoopagales"]}
    print("\nsecondary 3 — the lineages named in RESULTS.md:")
    for grp, names in named.items():
        for o in names:
            row = neon[neon.order == o]
            if len(row):
                r0 = row.iloc[0]
                print(f"  {grp:22s} {o:20s} NEON water_fraction {r0.water_fraction:.3f} "
                      f"(disc {int(r0.discordant)}, q {r0.q:.3g}) -> {r0.prefers}")
            else:
                print(f"  {grp:22s} {o:20s} not testable in NEON")

# ----------------------------------------------------------------- sensitivity (plan 7.4)
print("\nsecondary 4 — sensitivity")
print(f"{'min_reads':>10s} {'pairs':>6s} {'tested':>7s} {'sig':>5s} {'agree':>12s}")
for mr in READ_SENS:
    r2, k2, _, _ = test_pairs(pairs, water, benthic, mr)
    if not len(r2):
        print(f"{mr:>10d} {len(k2):>6d}       -     -            -"); continue
    m2 = fsig.merge(r2[["order", "water_fraction"]], on="order", suffixes=("_fun", "_neon"))
    a = ((m2.water_fraction_fun > .5) == (m2.water_fraction_neon > .5)).sum() if len(m2) else 0
    print(f"{mr:>10d} {len(k2):>6d} {len(r2):>7d} {int((r2.q<0.05).sum()):>5d} "
          f"{a:>5d}/{len(m2):<3d} {a/len(m2) if len(m2) else 0:>5.1%}")

print(f"\n{'conf':>10s} {'pairs':>6s} {'tested':>7s} {'sig':>5s} {'agree':>12s}")
for cf in CONF_SENS:
    if cf == MIN_CONF:
        w2, b2 = water, benthic
    else:
        w2, b2 = load("water", cf), load("benthic", cf)
    p2 = build_pairs(w2, b2, window)
    r2, k2, _, _ = test_pairs(p2, w2, b2, MIN_READS)
    if not len(r2):
        print(f"{cf:>10.2f} {len(k2):>6d}       -     -            -"); continue
    m2 = fsig.merge(r2[["order", "water_fraction"]], on="order", suffixes=("_fun", "_neon"))
    a = ((m2.water_fraction_fun > .5) == (m2.water_fraction_neon > .5)).sum() if len(m2) else 0
    print(f"{cf:>10.2f} {len(k2):>6d} {len(r2):>7d} {int((r2.q<0.05).sum()):>5d} "
          f"{a:>5d}/{len(m2):<3d} {a/len(m2) if len(m2) else 0:>5.1%}")

json.dump(res, open(P / "results/neon_validation.json", "w"), indent=1, default=str)
print(f"\nwrote results/neon_order.tsv, results/neon_pairs.tsv, "
      f"results/neon_vs_funaqua.tsv, results/neon_validation.json")
