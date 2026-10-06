"""Reproducible numerical settings and source identity for exported runs."""
from datetime import datetime, timezone
import hashlib
import importlib.metadata
from pathlib import Path
import platform
import subprocess
import sys

VERSION = "3.2-rc1"


def run_manifest(solver):
    root = Path(__file__).resolve().parent
    packages = {}
    for name in ("numpy", "numba", "llvmlite", "matplotlib"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    # File hashes identify source independently of Git availability.
    hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(root.glob("*.py"))}
    git = {"commit": None, "working_tree_dirty": None}
    try:
        options = dict(cwd=root, capture_output=True, text=True, timeout=3, check=True)
        git["commit"] = subprocess.run(["git", "rev-parse", "HEAD"], **options).stdout.strip()
        git["working_tree_dirty"] = bool(subprocess.run(
            ["git", "status", "--porcelain", "--untracked-files=normal"], **options).stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return dict(hemoflow_version=VERSION, created_utc=datetime.now(timezone.utc).isoformat(),
                python=sys.version, platform=platform.platform(), processor=platform.processor(),
                argv=sys.argv, packages=packages, git=git, source_sha256=hashes,
                numerics=solver.numerics(),
                particle_advection="GUI only; headless runs compute fluid fields and particle reference parameters")
