# The occupancy model does not justify itself on power — but it nearly doubles coverage

2026-10-23. Analysis as fixed in `OCCUPANCY_PLAN.md` (sha256 `74178c173ef708dd…`, full value in
`PLAN_HASHES.md`, commit `be9490c`, frozen before `analysis_occupancy.py` existed).

**The prediction the plan called "the claim that justifies the extra machinery" FAILED.** Three of four
others were confirmed. The plan committed in advance to what a failure here means, and that reading is
applied below rather than argued around.

## Six estimators, same generator, same seeds

| arm | estimator | evaluable | power | FDR | direction |
|---|---|---|---|---|---|
| **abundance** | `betabinom` | 117 | **0.852** | 0.056 | 0.973 |
| | `clr_ttest` | 254 | 0.180 | 0.880 | 0.812 |
| | `prop_wilcoxon` | 153 | 0.212 | 0.819 | 0.813 |
| | `occupancy` | 151 | **0.016** | 0.000 | 0.500 |
| | `pa_single` / `pa_repeated` | 79 / 76 | 0.000 | 0.000 | 0.000 |
| **occupancy** | `pa_single` | 86 | **0.385** | 0.000 | 1.000 |
| | `betabinom` | 114 | 0.383 | 0.052 | 0.981 |
| | `pa_repeated` | 83 | 0.317 | 0.000 | 1.000 |
| | **`occupancy`** | 152 | **0.286** | 0.000 | 1.000 |
| | `clr_ttest` | 254 | 0.129 | 0.853 | 0.836 |
| | `prop_wilcoxon` | 154 | 0.062 | 0.615 | 0.486 |
| **both** | `pa_single` | 87 | 0.920 | 0.125 | 1.000 |
| | `betabinom` | 126 | 0.908 | 0.050 | 0.980 |
| | `pa_repeated` | 82 | 0.875 | 0.037 | 1.000 |
| | **`occupancy`** | **155** | 0.833 | **0.030** | 1.000 |
| | `clr_ttest` | 254 | 0.441 | 0.813 | 0.884 |
| | `prop_wilcoxon` | 161 | 0.376 | 0.755 | 0.742 |

| # | prediction | observed | verdict |
|---|---|---|---|
| O1 | highest power of all six in the occupancy arm | ranks **4 of 6** — 0.286 against `pa_single`'s 0.385 | **FAILED** |
| O2 | near-blind in the abundance arm, by design | power 0.016, FDR 0.000 | CONFIRMED |
| O3 | controls FDR better than a single presence/absence draw | 0.030 against 0.125 in the `both` arm | CONFIRMED |
| O4 | evaluates at least as many orders as `pa_single` | 151/79, 152/86, 155/87 across the three arms | CONFIRMED |

## O1 failed, and the plan said what that means

> "O1 is the claim that justifies the extra machinery. If the model is not best in the arm built for it,
> there is no case for recommending it over the simpler draw-stable rule, and the manuscript will say so."

It is not best in that arm. **On power there is no case for it over the draw-stable rule**, and the
manuscript says so.

**Why it loses, which is the useful part.** Power in the occupancy arm is 0.031 at *D* = 100 with *r* = 2
and 0.900 at *D* = 10 000 with *r* = 8. The model must estimate the detection probability *p* from the
detection histories themselves, and where detections are sparse that estimate is noisy, so the correction
costs more precision than it buys. The extra machinery helps where data are already plentiful and hurts
where they are thin — the opposite of where a detection correction would be most wanted. A simple
discordant-pair test makes no such demand on the data and wins at low depth.

## What it does earn, stated separately so it does not disguise the failure

Two results stand on their own and do not depend on O1.

**It nearly doubles coverage at better error (O3, O4).** Under the realistic `both` encoding it evaluates
**155 orders against the draw-stable rule's 82** — 1.9× — at power 0.833 against 0.875 and FDR 0.030
against 0.037. That speaks directly to § 3.4a's problem, which is not that preferences are mis-estimated
but that too few taxa can be assessed at all. Trading 0.042 of power for 73 more evaluable orders at no
cost in FDR is a real option, and unlike `clr_ttest`'s nominal 254 it is coverage with an evidence floor
rather than coverage by having none.

**It is correctly silent when occupancy does not differ (O2).** In the `abundance` arm, where a
preferring taxon is more abundant in water but no more often present, it calls almost nothing: power
0.016, FDR 0.000. Every other estimator either sees nothing at all (`pa_single`, `pa_repeated`) or fires
hard (`betabinom`, 0.852). The occupancy model is the only one whose answer tracks its own estimand —
it reports that occupancy does not differ, because it does not. Under this benchmark's definition of
truth that scores as missed power; the plan pre-committed to reading it as correct behaviour, and it is.

## So what should be used

No single recommendation survives all three encodings, which was already
`BENCHMARK2_RESULTS.md`'s conclusion and is unchanged. Within the realistic arm:

- **For the fewest false discoveries with the least machinery:** the draw-stable rule. FDR 0.037,
  power 0.875, 82 orders.
- **To assess roughly twice as many taxa at comparable error:** the occupancy model (155 orders,
  FDR 0.030) or the overdispersed binomial (126, FDR 0.050).
- **Never without an evidence floor:** `clr_ttest` reaches all 254 orders at FDR 0.813.

The occupancy model's case rests on coverage, not on power, and it fails below about 1 000 reads per
sample where its detection estimate becomes unstable.

## Limits

- ψ is occupancy **in the sequenced library**, not at the site, as the plan stated before the run. The
  model corrects subsampling, not sequencing, and cannot recover a taxon the library missed.
- *p* is constant across samples within a compartment, which is wrong — detectability depends on a
  taxon's abundance in each particular sample. A per-sample detection covariate is the obvious next
  step and would likely fix the low-depth weakness that sank O1.
- Subsamples are drawn with replacement from one library: exchangeable, not independent visits.
- Single-season, two-group, no covariates, no random effects, *K* = 8 fixed.
- The generator's limits from `BENCHMARK2_PLAN.md` § 7 are inherited unchanged.
- Nothing was run on real FunAqua or NEON data.

## Files

`OCCUPANCY_PLAN.md` (frozen) · `analysis_occupancy.py` → `results/occupancy_model.tsv`,
`results/occupancy_model.json`. Read with `BENCHMARK2_RESULTS.md`, whose five estimators it joins.
