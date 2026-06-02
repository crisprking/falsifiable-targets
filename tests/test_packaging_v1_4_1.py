"""Packaging regression tests.

These tests catch the v1.4.0 → v1.4.1 packaging bug where YAML data
files in sentinels/ and claims/ were not shipped with the wheel.

If these tests fail, do not release: a clean install (Kaggle, Colab,
fresh CI) will fail with FileNotFoundError at runtime even though
`pip install` succeeded.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = ROOT  # back-compat alias for existing test bodies

def _has_isolated_build() -> bool:
    """Whether the environment can `python -m build` cleanly."""
    try:
        import build  # noqa: F401
        return True
    except ImportError:
        return False


@pytest.mark.skipif(not _has_isolated_build(), reason="requires `build` package")
def test_wheel_includes_sentinels_yaml(tmp_path: Path) -> None:
    """Building a wheel must include sentinels/v1_sentinels.yaml.

    Regression test for v1.4.0 bug where setuptools silently dropped
    sentinels/ because it wasn't a package (no __init__.py).
    """
    dist_dir = tmp_path / "dist"
    result = subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(dist_dir), str(REPO_ROOT)],
        capture_output=True, text=True, timeout=180,
    )
    assert result.returncode == 0, f"wheel build failed: {result.stderr[-1000:]}"

    wheels = list(dist_dir.glob("*.whl"))
    assert len(wheels) == 1, f"expected exactly one wheel, got {wheels}"

    # Unpack and inspect
    unpack = tmp_path / "unpack"
    subprocess.run(
        [sys.executable, "-m", "zipfile", "-e", str(wheels[0]), str(unpack)],
        check=True, capture_output=True,
    )
    sentinel = unpack / "sentinels" / "v1_sentinels.yaml"
    assert sentinel.exists(), (
        f"sentinels/v1_sentinels.yaml missing from wheel. "
        f"Wheel contents: {sorted(p.relative_to(unpack) for p in unpack.rglob('*'))[:30]}"
    )


@pytest.mark.skipif(not _has_isolated_build(), reason="requires `build` package")
def test_wheel_includes_claim_yamls(tmp_path: Path) -> None:
    """Building a wheel must include claims/*.yaml example files."""
    dist_dir = tmp_path / "dist"
    result = subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(dist_dir), str(REPO_ROOT)],
        capture_output=True, text=True, timeout=180,
    )
    assert result.returncode == 0, f"wheel build failed: {result.stderr[-1000:]}"

    wheels = list(dist_dir.glob("*.whl"))
    unpack = tmp_path / "unpack"
    subprocess.run(
        [sys.executable, "-m", "zipfile", "-e", str(wheels[0]), str(unpack)],
        check=True, capture_output=True,
    )
    claim_yamls = list((unpack / "claims").glob("*.yaml"))
    assert len(claim_yamls) >= 3, (
        f"expected ≥3 example claim YAMLs in wheel, got {len(claim_yamls)}: {claim_yamls}"
    )


def test_install_resolves_sentinel_path(tmp_path: Path) -> None:
    """Install into an isolated prefix and verify the runtime sentinel path resolves.

    This is the integration-level guarantee: not just "is the file in the
    wheel", but "after pip install, does smoke_test.py find it where it
    expects to".
    """
    # Build a wheel first
    dist_dir = tmp_path / "dist"
    if not _has_isolated_build():
        pytest.skip("requires `build` package for hermetic wheel construction")
    result = subprocess.run(
        [sys.executable, "-m", "build", "--wheel", "--outdir", str(dist_dir), str(REPO_ROOT)],
        capture_output=True, text=True, timeout=180,
    )
    if result.returncode != 0:
        pytest.skip(f"build failed in this env: {result.stderr[-500:]}")

    wheel = next(dist_dir.glob("*.whl"))

    # Install into an isolated prefix
    prefix = tmp_path / "prefix"
    install = subprocess.run(
        [sys.executable, "-m", "pip", "install",
         "--target", str(prefix), "--no-deps", "--quiet", str(wheel)],
        capture_output=True, text=True, timeout=120,
    )
    assert install.returncode == 0, f"pip install --target failed: {install.stderr[-500:]}"

    sentinel = prefix / "sentinels" / "v1_sentinels.yaml"
    assert sentinel.exists(), (
        f"After `pip install --target`, sentinel file not at {sentinel}. "
        f"Prefix contents: {sorted(p.relative_to(prefix) for p in prefix.iterdir())}"
    )
