#!/usr/bin/env python3
"""Every plan hash quoted anywhere in the repo must match the file on disk.

The plans are the pre-registration. If a quoted hash and the file disagree, either
the file was edited after freezing or a hash was mistyped, and both destroy the
only evidence that the plan predates the result.
"""
import hashlib
import pathlib
import re
import sys

ROOT = pathlib.Path.home() / "research/projects/funaqua-habitat"
# Plan filenames contain digits (METAPR2_...), so the character class needs \w.
PAT = re.compile(r"`?(\w*PLAN\w*\.md)`?[^\n]*?([0-9a-f]{64})")

quoted = {}
for p in sorted(ROOT.glob("*.md")):
    for m in PAT.finditer(p.read_text()):
        quoted.setdefault(m.group(1), {}).setdefault(m.group(2), []).append(p.name)

plans = sorted(f.name for f in ROOT.glob("*PLAN*.md"))
print(f"plans on disk: {len(plans)}")

bad = []
for plan in plans:
    actual = hashlib.sha256((ROOT / plan).read_bytes()).hexdigest()
    hs = quoted.get(plan, {})
    if not hs:
        print(f"  {plan:30s} no hash quoted anywhere")
        continue
    for h, where in sorted(hs.items()):
        if h == actual:
            print(f"  {plan:30s} {h[:12]}… match   (in {', '.join(sorted(set(where)))})")
        else:
            bad.append((plan, h, actual, where))
            print(f"  {plan:30s} {h[:12]}… MISMATCH, file is {actual[:12]}…")

unknown = sorted(set(quoted) - set(plans))
if unknown:
    print(f"  hashes quoted for files that do not exist: {unknown}")
    bad.append(("missing files", unknown, "", []))

print("\nall quoted plan hashes match the frozen files" if not bad
      else "\nINTEGRITY PROBLEM — see above")
sys.exit(1 if bad else 0)
