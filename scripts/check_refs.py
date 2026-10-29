#!/usr/bin/env python3
"""Cross-check in-text citations against the reference list, both directions.

A reference list is wrong in two ways: a citation with no entry, and an entry
nothing cites. Both are caught here so neither reaches a reviewer.
"""
import argparse
from pathlib import Path
import re
import sys

parser = argparse.ArgumentParser(description="Check reference consistency in an external Markdown manuscript")
parser.add_argument("manuscript", type=Path, help="Path to the author's local Markdown manuscript")
args = parser.parse_args()
t = args.manuscript.read_text(encoding="utf-8")
m = re.search(r"^## (?:6 )?References\\s*$", t, re.MULTILINE)
if not m:
    parser.error("Reference heading not found (expected '## References' or '## 6 References')")
head, body = t[:m.start()], t[m.end():]

# Entries: a line starting with Surname, Initial.  or an institutional author.
entries = {}
for line in body.splitlines():
    m = re.match(r"^([A-Z][A-Za-zñéè'-]+),\s+[A-Z]\.", line)
    if m:
        entries[m.group(1)] = line
        continue
    m = re.match(r"^(National Ecological Observatory Network|OpenTree)\s*\(", line)
    if m:
        entries[m.group(1)] = line

# Citations in the text. A surname may sit anywhere inside a parenthetical, not
# just at its opening bracket, so match the surname-year pair itself wherever it
# occurs: "Surname et al., YEAR", "Surname and Other, YEAR", "Surname (YEAR)".
SUR = r"[A-Z][A-Za-zñéè'-]+"
cited = set()
for pat in (rf"\b({SUR})\s+et\s+al\.,?\s*\(?(\d{{4}})",
            rf"\b({SUR})\s+and\s+{SUR},?\s*\(?(\d{{4}})",
            rf"\b({SUR}),\s*(\d{{4}})\b"):
    for m in re.finditer(pat, head):
        cited.add(m.group(1))
cited.discard("Hochberg")
cited.discard("Holder")
# "Benjamini and Hochberg, 1995" style, second surname is not an entry key
print(f"reference entries : {len(entries)}")
for k in sorted(entries):
    print(f"    {k}")
print(f"\nsurnames cited in text : {len(cited)}")
for k in sorted(cited):
    print(f"    {k}")

missing = sorted(c for c in cited if c not in entries)
# NEON and OpenTree are cited by product/version in prose, not author-year
uncited = sorted(e for e in entries
                 if e not in cited
                 and e not in {"National Ecological Observatory Network", "OpenTree"})

print("\n=== citations with no reference entry")
print("    " + (", ".join(missing) if missing else "none"))
print("=== reference entries nothing cites")
print("    " + (", ".join(uncited) if uncited else "none"))

print("\n=== leftover placeholders anywhere")
bad = [f"line {i}: {l.strip()[:90]}"
       for i, l in enumerate(t.splitlines(), 1)
       if "[ref" in l or "needs full citation" in l.lower()
       or re.search(r"\*Verify\b", l) or "DO NOT SUBMIT" in l]
print("    " + ("\n    ".join(bad) if bad else "none"))

print("\n=== stale counts")
# The abstract once said "three findings" and the introduction "we report three
# things" when the simulation made it four. Catch those two claims specifically.
# An earlier, looser pattern also flagged "Three things follow" in 3.4d, which
# genuinely introduces three conclusions - a count of something else entirely.
stale = [f"line {i}: {l.strip()[:90]}"
         for i, l in enumerate(t.splitlines(), 1)
         if re.search(r"\b[Tt]hree findings\b", l)
         or re.search(r"\b(report|reports|reporting)\s+three\s+things\b", l)]
print("    " + ("\n    ".join(stale) if stale else "none"))

sys.exit(1 if (missing or uncited or bad or stale) else 0)
