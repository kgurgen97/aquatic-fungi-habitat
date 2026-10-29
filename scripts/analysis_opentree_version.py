#!/usr/bin/env python3
"""Recover the Open Tree synthesis version and re-check the §4.5 mapping numbers.

results/order_ott_ids.json exists but no script in the repo produced it, so the
synthesis release it came from was never recorded. A database cited without a
version is a reviewer's first question, so re-query the API, record the version,
and check whether the counts the manuscript reports (70 orders mapped, 47 with a
preference value) still hold. If they moved, the manuscript gets the new numbers.
"""
import json
import pathlib
import urllib.request

import pandas as pd

API = "https://api.opentreeoflife.org/v3"
from project_paths import REPO_ROOT, FUNAQUA_DIR, NEON_DIR, TRAITS_DIR, METAPR2_DIR
ROOT = REPO_ROOT


def post(path, payload):
    """POST and return the parsed body. The induced_subtree endpoint answers 400
    when any ott_id is absent from the synthesis, and names them in the error
    body, so that body is the answer rather than a failure."""
    req = urllib.request.Request(
        f"{API}{path}",
        data=json.dumps(payload).encode(),
        headers={"content-type": "application/json"},
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.load(r)
    except urllib.error.HTTPError as e:
        raw = e.read().decode(errors="replace")
        try:
            return {"_http_error": e.code, **json.loads(raw)}
        except ValueError:
            return {"_http_error": e.code, "_raw": raw[:600]}


about = post("/tree_of_life/about", {})
print("=== Open Tree synthetic tree, as served now")
for k in ("synth_id", "date_created", "taxonomy_version",
          "num_source_trees", "num_source_studies"):
    print(f"  {k:20s} {about.get(k, '(absent)')}")
print(f"  {'root ntips':20s} {about.get('root', {}).get('num_tips', '(absent)')}")

# The orders the manuscript mapped: every order with a preference value, plus the
# full tested set, exactly as §4.5 describes it.
fungi = pd.read_csv(ROOT / "results/primary_order.tsv", sep="\t")
names = sorted(set(fungi["order"].dropna()))
print(f"\n=== mapping {len(names)} FunAqua fungal orders through TNRS")

matched, unmatched = {}, []
for i in range(0, len(names), 250):
    chunk = names[i:i + 250]
    res = post("/tnrs/match_names", {"names": chunk, "do_approximate_matching": False})
    for r in res.get("results", []):
        if r["matches"]:
            m = r["matches"][0]
            matched[r["name"]] = m["taxon"]["ott_id"]
    unmatched += res.get("unmatched_names", [])

print(f"  matched   {len(matched)}")
print(f"  unmatched {len(unmatched)}")

# How many of those are actually placed in the synthetic tree? The first call
# tells us which ids the synthesis does not carry; the second call, with those
# removed, is the one that returns a tree.
ott_ids = sorted(set(matched.values()))
sub = post("/tree_of_life/induced_subtree",
           {"ott_ids": ott_ids, "label_format": "name"})
absent = set()
for key in ("unknown", "broken"):
    absent |= {int(str(k).removeprefix("ott")) for k in sub.get(key, {})}
if absent:
    print(f"  not in synthesis: {len(absent)}")
    sub = post("/tree_of_life/induced_subtree",
               {"ott_ids": [o for o in ott_ids if o not in absent],
                "label_format": "name"})
in_tree = len(ott_ids) - len(absent)
assert "_http_error" not in sub, f"induced_subtree still failing: {sub}"
print(f"  in synthetic tree {in_tree} of {len(ott_ids)} matched")

# Compare with the stored file, which no script in the repo wrote.
stored = json.loads((ROOT / "results/order_ott_ids.json").read_text())
print(f"\n=== against results/order_ott_ids.json ({len(stored)} entries)")
same = sum(1 for n, o in matched.items() if stored.get(n) == o)
moved = {n: (stored[n], o) for n, o in matched.items()
         if n in stored and stored[n] != o}
print(f"  identical ott_id {same}")
print(f"  changed ott_id   {len(moved)}")
for n, (old, new) in list(moved.items())[:10]:
    print(f"    {n}: {old} -> {new}")
print(f"  in file but not matched now: {len(set(stored) - set(matched))}")
print(f"  matched now but not in file: {len(set(matched) - set(stored))}")

out = {"synth_id": about.get("synth_id"),
       "date_created": about.get("date_created"),
       "taxonomy_version": about.get("taxonomy_version"),
       "n_orders_submitted": len(names),
       "n_matched": len(matched),
       "n_in_synthetic_tree": in_tree,
       "n_identical_to_stored": same,
       "n_changed_vs_stored": len(moved)}
(ROOT / "results/opentree_version.json").write_text(json.dumps(out, indent=1))
print("\nwrote results/opentree_version.json")
