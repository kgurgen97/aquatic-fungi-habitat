#!/usr/bin/env python3
"""Check the frozen preregistration files against their recorded SHA-256 hashes."""
import hashlib
import re
from project_paths import REPO_ROOT

PLAN_DIR = REPO_ROOT / "plans"
manifest = (PLAN_DIR / "PLAN_HASHES.md").read_text(encoding="utf-8")
recorded = dict(re.findall(r"\`([A-Z0-9_]+PLAN\.md)\`\s*\|\s*\`([a-f0-9]{64})\`", manifest))
present = {p.name for p in PLAN_DIR.glob("*PLAN.md")}
assert present, "No frozen preregistration plans found"
assert set(recorded) == present, (
    f"Manifest mismatch: missing={sorted(present-set(recorded))}; "
    f"extra={sorted(set(recorded)-present)}"
)
for filename in sorted(present):
    digest = hashlib.sha256((PLAN_DIR / filename).read_bytes()).hexdigest()
    assert recorded[filename] == digest, f"Frozen plan modified: {filename}"
print(f"Frozen plans verified: {len(present)}")
