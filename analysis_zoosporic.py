#!/usr/bin/env python3
"""The flagellum hypothesis, exactly as ZOOSPORIC_PLAN.md specifies.

Plan sha256: e02627db39a9934cc2100f360ec534ba912c302c79864edf1b19b64fd1c8090f

Run once. Each dataset is tested independently so NEON is a replication, not a pooled
re-description.
"""
import json
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

P = Path.home() / "research/projects/funaqua-habitat"
SEED = 20261021
N_PERM = 10_000

# ---- trait assignment, fixed in the plan before any group value was computed
FLAGELLATED = {"Chytridiomycota", "Monoblepharomycota", "Neocallimastigomycota",
               "Blastocladiomycota", "Olpidiomycota", "Rozellomycota"}
AFLAGELLATE = {"Ascomycota", "Basidiomycota", "Mucoromycota", "Glomeromycota",
               "Kickxellomycota", "Zoopagomycota", "Entomophthoromycota",
               "Basidiobolomycota", "Mortierellomycota", "Calcarisporiellomycota",
               "Entorrhizomycota"}
AMBIGUOUS = {"Sanchytriomycota", "Aphelidiomycota"}

pm = json.load((P / "results/order_phylum_map.json").open())
stab = pd.read_csv(P / "results/draw_stability.tsv", sep="\t")
NS = json.load((P / "results/draw_stability.json").open())["n_draws"]
STABLE = set(stab[stab.times_significant >= 0.95 * NS].order)


def trait(ph):
    if ph in FLAGELLATED:
        return "flagellated"
    if ph in AFLAGELLATE:
        return "aflagellate"
    if ph in AMBIGUOUS:
        return "ambiguous"
    return "unplaced"


def build(tsv, phylum_map, frac_col):
    d = pd.read_csv(P / tsv, sep="\t")[["order", frac_col, "discordant", "q"]]
    d = d.rename(columns={frac_col: "water_fraction"})
    d["phylum"] = [phylum_map.get(o, "?") for o in d.order]
    d["trait"] = [trait(p) for p in d.phylum]
    return d


def run(name, d):
    out = {"dataset": name, "orders_tested": int(len(d))}
    print("=" * 76)
    print(f"{name}: {len(d)} testable orders")
    print("  trait groups: " + ", ".join(
        f"{k} {v}" for k, v in d.trait.value_counts().items()))

    f = d[d.trait == "flagellated"].water_fraction.to_numpy()
    a = d[d.trait == "aflagellate"].water_fraction.to_numpy()
    out.update(n_flag=len(f), n_aflag=len(a))
    if len(f) < 3 or len(a) < 3:
        print("  too few in one group to test"); return out

    # ---- primary: Mann-Whitney
    u = stats.mannwhitneyu(f, a, alternative="two-sided")
    rbc = 2 * u.statistic / (len(f) * len(a)) - 1          # rank-biserial
    print(f"\n  PRIMARY  Mann-Whitney U")
    print(f"    flagellated  n={len(f):3d}  median {np.median(f):.3f}  mean {f.mean():.3f}")
    print(f"    aflagellate  n={len(a):3d}  median {np.median(a):.3f}  mean {a.mean():.3f}")
    print(f"    U={u.statistic:.1f}  p={u.pvalue:.4g}  rank-biserial={rbc:+.3f}")
    direction = "flagellated higher" if np.median(f) > np.median(a) else "aflagellate higher"
    verdict = ("SUPPORTED" if u.pvalue < 0.05 and np.median(f) > np.median(a)
               else "CONTRARY" if u.pvalue < 0.05 else "NOT SUPPORTED")
    print(f"    direction: {direction}   ->  {verdict}")
    out.update(median_flag=float(np.median(f)), median_aflag=float(np.median(a)),
               mw_u=float(u.statistic), mw_p=float(u.pvalue), rank_biserial=float(rbc),
               verdict=verdict)

    # ---- secondary 2: proportion above 0.5
    fw, aw = int((f > .5).sum()), int((a > .5).sum())
    ft = stats.fisher_exact([[fw, len(f) - fw], [aw, len(a) - aw]])
    print(f"\n  secondary 2  share with water_fraction > 0.5")
    print(f"    flagellated {fw}/{len(f)} = {fw/len(f):.1%}   "
          f"aflagellate {aw}/{len(a)} = {aw/len(a):.1%}")
    print(f"    Fisher OR={ft.statistic:.3f}  p={ft.pvalue:.4g}")
    out.update(share_flag=fw / len(f), share_aflag=aw / len(a),
               fisher_or=float(ft.statistic), fisher_p=float(ft.pvalue))

    # ---- secondary 3: draw-stable orders only (FunAqua's stability set)
    ds = d[d.order.isin(STABLE)]
    f2 = ds[ds.trait == "flagellated"].water_fraction.to_numpy()
    a2 = ds[ds.trait == "aflagellate"].water_fraction.to_numpy()
    print(f"\n  secondary 3  restricted to FunAqua's draw-stable orders")
    if len(f2) >= 3 and len(a2) >= 3:
        u2 = stats.mannwhitneyu(f2, a2, alternative="two-sided")
        print(f"    flagellated n={len(f2)} median {np.median(f2):.3f}   "
              f"aflagellate n={len(a2)} median {np.median(a2):.3f}   p={u2.pvalue:.4g}")
        out.update(stable_n_flag=len(f2), stable_n_aflag=len(a2),
                   stable_p=float(u2.pvalue),
                   stable_median_flag=float(np.median(f2)),
                   stable_median_aflag=float(np.median(a2)))
    else:
        print(f"    too few ({len(f2)} flagellated, {len(a2)} aflagellate)")

    # ---- secondary 4: label permutation
    rng = np.random.default_rng(SEED)
    lab = np.r_[np.ones(len(f), bool), np.zeros(len(a), bool)]
    vals = np.r_[f, a]
    obs = np.median(f) - np.median(a)
    null = np.empty(N_PERM)
    for i in range(N_PERM):
        p = rng.permutation(lab)
        null[i] = np.median(vals[p]) - np.median(vals[~p])
    pperm = (np.sum(np.abs(null) >= abs(obs)) + 1) / (N_PERM + 1)
    print(f"\n  secondary 4  label permutation ({N_PERM:,} shuffles)")
    print(f"    observed median difference {obs:+.3f}   permutation p={pperm:.4g}")
    out.update(perm_obs=float(obs), perm_p=float(pperm))

    # ---- declared aside: the ambiguous phyla
    amb = d[d.trait == "ambiguous"]
    print(f"\n  declared aside  phyla excluded as ambiguous: {len(amb)} orders")
    for _, r in amb.iterrows():
        print(f"    {r.order:34s} {r.phylum:20s} water_fraction {r.water_fraction:.3f}")
    out["ambiguous"] = [{"order": r.order, "phylum": r.phylum,
                         "water_fraction": float(r.water_fraction)}
                        for _, r in amb.iterrows()]

    print(f"\n  flagellated orders, highest water preference first:")
    ff = d[d.trait == "flagellated"].sort_values("water_fraction", ascending=False)
    for _, r in ff.iterrows():
        print(f"    {r.water_fraction:.3f}  {r.order:34s} {r.phylum:20s} "
              f"disc {int(r.discordant):3d}  q {r.q:.3g}")
    return out


fa = build("results/primary_order.tsv", pm["funaqua"], "water_fraction")
ne = build("results/neon_order.tsv", pm["neon"], "water_fraction")
r_fa, r_ne = run("FunAqua (water vs sediment)", fa), run("NEON (water vs benthic)", ne)

pd.concat([fa.assign(dataset="funaqua"), ne.assign(dataset="neon")]).to_csv(
    P / "results/zoosporic.tsv", sep="\t", index=False)

print("=" * 76)
print("REPLICATION RULE (plan): predicted direction with p < 0.05 in BOTH datasets independently.")
ok = all(r.get("verdict") == "SUPPORTED" for r in (r_fa, r_ne))
for r in (r_fa, r_ne):
    print(f"  {r['dataset']:30s} p={r.get('mw_p', float('nan')):.4g}  "
          f"medians {r.get('median_flag', float('nan')):.3f} vs "
          f"{r.get('median_aflag', float('nan')):.3f}  -> {r.get('verdict')}")
print(f"\n  REPLICATED: {'YES' if ok else 'NO'}")

json.dump({"plan_sha256": "e02627db39a9934cc2100f360ec534ba912c302c79864edf1b19b64fd1c8090f",
           "replicated": bool(ok), "funaqua": r_fa, "neon": r_ne},
          open(P / "results/zoosporic.json", "w"), indent=1)
print("\nwrote results/zoosporic.tsv and results/zoosporic.json")
