# Frozen plan hashes

Each pre-registration plan's SHA-256 and the commit that introduced it. The hashes are kept
HERE and not inside the plan files: writing a file's own hash into the file changes the file,
so the recorded value could never match, and editing a frozen plan for any reason defeats the
point of freezing it. Verify with:

    sha256sum *PLAN*.md

| plan | sha256 | first commit | commit date |
|---|---|---|---|
| `ANALYSIS_PLAN.md` | `68ea46ccef0ef50b74373507c84e23dde59d78dac4551bed5599df6f8a731cc1` | `bbc9892` | 2026-10-21 |
| `NEON_VALIDATION_PLAN.md` | `85fd8aa2c9f8a14dc30c9bded5858d4f39dd0d6801f01753531fda0ef75fb939` | `bbc9892` | 2026-10-21 |
| `ZOOSPORIC_PLAN.md` | `e02627db39a9934cc2100f360ec534ba912c302c79864edf1b19b64fd1c8090f` | `a12b677` | 2026-10-21 |
| `EUK_ROBUSTNESS_PLAN.md` | `a12794f91b0df75b433df5fe37afb346ab3a603b5ed53152ca1a3dee85f879f8` | `a5c9fe6` | 2026-10-21 |
| `METAPR2_REPLICATION_PLAN.md` | `52d75c02f3091fdbeab40148e13559390f5f742983d35192e84474bca86d77b2` | `eba1edf` | 2026-10-21 |
| `DEPTH_LADDER_PLAN.md` | `b53dd38de0be59e304aa0844545d54ed79d05f1a601409574f234979b49165a5` | `464e481` | 2026-10-23 |
| `SWEEP_PLAN.md` | `2af9a5a9e4002fc02e9ee74e7cf5aa6b7e11a4e457ab51141ed105017ccf8a0a` | `7a0d398` | 2026-10-23 |
