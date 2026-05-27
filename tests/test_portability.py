"""
Portability tests.

These tests exist because v1.4.0 fixed a real bug: three test files
hardcoded `/kaggle/working/` as the project root, breaking the
test suite for any researcher cloning the repo to a non-Kaggle
machine. The fix was to use `Path(__file__).resolve().parent.parent`
in each test file.

If anyone ever puts a hardcoded /kaggle/, /home/USER/, or similar
absolute-path back in the codebase, these tests will catch it
before it ships.

Why hand-rolled and not pylint/ruff? ruff doesn't know that
/kaggle/working/ is wrong - it's syntactically valid Python.
The domain-specific guard belongs in the test suite.
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Patterns that should NEVER appear in source files. If the project is
# meant to run anywhere, none of these are acceptable.
FORBIDDEN_PATH_PATTERNS = [
    r"/kaggle/",                # Kaggle-specific
    r"/content/",               # Colab-specific
    r"C:\\\\Users\\\\",         # Windows user-specific
    r'"/home/[a-zA-Z]+/"',      # Linux home-dir
    r'"/Users/[a-zA-Z]+/"',     # macOS home-dir
]

# Files that should be scanned. Skip docs/markdown (which may legitimately
# mention these paths as examples) and the v1.4.0 changelog itself.
SCAN_GLOBS = ["*.py", "tests/*.py", "adapters/*.py"]

# Lines that mention the patterns in comments/docstrings ABOUT how the
# bug was fixed are allowed. We allow lines that contain the literal
# phrase "fixed" or "ALLOWED" or "regression-guard".
ALLOWLIST_MARKERS = (
    "# allowed:",
    "# ALLOWED:",
    "# fixed",
    "# Fixed",
    "regression-guard",
    "FORBIDDEN_PATH_PATTERNS",  # this test file itself
)


def test_no_hardcoded_machine_specific_paths():
    """Fail loudly if any source file contains a non-portable absolute path."""
    offenders: list[tuple[Path, int, str]] = []

    for pattern_glob in SCAN_GLOBS:
        for path in ROOT.glob(pattern_glob):
            # Skip this very file - it intentionally contains the patterns
            if path == Path(__file__):
                continue
            text = path.read_text(encoding="utf-8")
            for line_no, line in enumerate(text.splitlines(), start=1):
                # Skip lines flagged as allowed
                if any(marker in line for marker in ALLOWLIST_MARKERS):
                    continue
                for pat in FORBIDDEN_PATH_PATTERNS:
                    if re.search(pat, line):
                        offenders.append((path.relative_to(ROOT), line_no, line.strip()))

    if offenders:
        msg = ["Found non-portable hardcoded paths:"]
        for path, line_no, line in offenders:
            msg.append(f"  {path}:{line_no}  {line}")
        msg.append(
            "\nFix: replace with `Path(__file__).resolve().parent.parent` "
            "or similar relative-to-source-file resolution."
        )
        raise AssertionError("\n".join(msg))


def test_tests_directory_uses_path_resolution():
    """Each test file must use Path-based ROOT computation, not literals."""
    expected_pattern = "Path(__file__).resolve()"
    test_files = list((ROOT / "tests").glob("test_*.py"))
    assert test_files, "no test files found?"

    missing: list[str] = []
    for tf in test_files:
        if tf.name == Path(__file__).name:
            continue  # don't recurse into self
        text = tf.read_text(encoding="utf-8")
        # A test file is OK if it either does Path(__file__) or doesn't
        # define ROOT at all (some tests may not need it)
        if "ROOT = " in text and expected_pattern not in text:
            missing.append(tf.name)

    assert not missing, (
        f"Test files have ROOT but don't use Path(__file__): {missing}. "
        f"Use ROOT = Path(__file__).resolve().parent.parent for portability."
    )


def test_pyproject_toml_present_and_parseable():
    """v1.4.0 introduced packaging; pyproject.toml is now a load-bearing file."""
    pp = ROOT / "pyproject.toml"
    assert pp.exists(), "pyproject.toml missing - packaging broken"

    try:
        import tomllib  # Python 3.11+
        data = tomllib.loads(pp.read_text())
    except ImportError:
        # 3.10 fallback - just verify it's there and has [project] header
        text = pp.read_text()
        assert "[project]" in text
        assert 'name = "falsifiable-targets"' in text
        return

    assert data["project"]["name"] == "falsifiable-targets"
    assert "ft-audit" in data["project"]["scripts"]
    assert "ft-smoke" in data["project"]["scripts"]
    assert "ft-validate" in data["project"]["scripts"]


def test_version_single_source_of_truth():
    """The _version.py module is the only place tool version is declared."""
    import sys
    sys.path.insert(0, str(ROOT))
    if "_version" in sys.modules:
        del sys.modules["_version"]
    import _version

    # Tool version is a string like "1.4.0"
    assert isinstance(_version.__version__, str)
    parts = _version.__version__.split(".")
    assert len(parts) == 3, f"version not semver: {_version.__version__}"
    assert all(p.isdigit() for p in parts), f"non-numeric semver: {parts}"

    # Ruleset SHA matches what tests/test_sentinels.py expects
    assert _version.RULESET_SHA256 == (
        "35ef2b2ab5363298097962a0b6ae52c70d551a1edddc341054f75cb6e4fb7221"
    ), "_version.py SHA disagrees with test_sentinels.py SHA - drift!"
