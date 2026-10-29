"""Paths used by the public reproducibility scripts.

The analysis repository does not redistribute the raw FunAqua, NEON, traits or
metaPR2 downloads. Set environment variables to the local directories before
running scripts that require those datasets.
"""
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent
FUNAQUA_DIR = Path(os.environ.get("FUNAQUA_DATA_DIR", Path.home() / "research/datasets/funaqua")).expanduser()
NEON_DIR = Path(os.environ.get("NEON_DATA_DIR", Path.home() / "research/datasets/neon2")).expanduser()
TRAITS_DIR = Path(os.environ.get("FUNGALTRAITS_DATA_DIR", Path.home() / "research/datasets/fungaltraits")).expanduser()
METAPR2_DIR = Path(os.environ.get("METAPR2_DATA_DIR", Path.home() / "research/datasets/metapr2/extract")).expanduser()
