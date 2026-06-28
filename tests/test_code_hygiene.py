"""Static code-hygiene checks that catch common maintenance hazards.

These tests scan the source and script files for patterns that have caused
bugs or dead code in the past:

- Magic numbers that should be shared constants
- sys.path hacks that bypass the installed package
- Unused imports
- Invalid NDArray type parameters
- Scripts that are strict subsets of more general scripts
"""

from __future__ import annotations

import ast
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC_DIR = ROOT / "src" / "segbench"
SCRIPTS_DIR = ROOT / "scripts"

ALL_PY = sorted(SRC_DIR.rglob("*.py")) + sorted(SCRIPTS_DIR.glob("*.py"))


def _read(p: Path) -> str:
    return p.read_text(encoding="utf-8")


# ── Magic number duplication ──────────────────────────────────────────────

def test_total_transcripts_not_hardcoded_outside_constants():
    """The ROI transcript count must only appear in constants.py."""
    pattern = re.compile(r"\b3[_,]?392[_,]?051\b")
    violations = []
    for p in ALL_PY:
        if p.name == "constants.py":
            continue
        for i, line in enumerate(_read(p).splitlines(), 1):
            if pattern.search(line) and not line.lstrip().startswith("#"):
                violations.append(f"{p.relative_to(ROOT)}:{i}")
    assert violations == [], (
        f"Hardcoded transcript count found outside constants.py:\n"
        + "\n".join(f"  {v}" for v in violations)
    )


# ── sys.path hacks ───────────────────────────────────────────────────────

def test_no_sys_path_insert():
    """Scripts should import from the installed package, not hack sys.path."""
    pattern = re.compile(r"sys\.path\.(insert|append)\(")
    violations = []
    for p in ALL_PY:
        for i, line in enumerate(_read(p).splitlines(), 1):
            if pattern.search(line) and not line.lstrip().startswith("#"):
                violations.append(f"{p.relative_to(ROOT)}:{i}")
    assert violations == [], (
        f"sys.path manipulation found (use pip install -e . instead):\n"
        + "\n".join(f"  {v}" for v in violations)
    )


# ── Invalid NDArray type parameters ──────────────────────────────────────

def test_no_ndarray_np_number():
    """NDArray[np.number] is invalid — use np.floating, np.integer, etc."""
    pattern = re.compile(r"NDArray\[np\.number\]")
    violations = []
    for p in ALL_PY:
        for i, line in enumerate(_read(p).splitlines(), 1):
            if pattern.search(line):
                violations.append(f"{p.relative_to(ROOT)}:{i}")
    assert violations == [], (
        f"Invalid NDArray[np.number] annotation (use np.floating or np.integer):\n"
        + "\n".join(f"  {v}" for v in violations)
    )


# ── Unused imports ────────────────────────────────────────────────────────

def _find_unused_imports(filepath: Path) -> list[str]:
    """Return list of 'name (line N)' for imports whose names aren't used."""
    source = _read(filepath)
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return []

    imports: list[tuple[str, int]] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            for alias in node.names:
                name = alias.asname or alias.name.split(".")[0]
                imports.append((name, node.lineno))
        elif isinstance(node, ast.ImportFrom):
            if node.module == "__future__":
                continue
            for alias in node.names:
                name = alias.asname or alias.name
                imports.append((name, node.lineno))

    unused = []
    for name, lineno in imports:
        if name.startswith("_"):
            continue
        count = len(re.findall(rf"\b{re.escape(name)}\b", source))
        if count <= 1:
            unused.append(f"{name} (line {lineno})")
    return unused


def test_no_unused_imports():
    """Every imported name should be referenced at least once beyond the import."""
    violations = {}
    for p in ALL_PY:
        unused = _find_unused_imports(p)
        if unused:
            violations[str(p.relative_to(ROOT))] = unused
    assert violations == {}, (
        "Unused imports found:\n"
        + "\n".join(
            f"  {f}: {', '.join(names)}" for f, names in violations.items()
        )
    )


# ── Redundant / superseded scripts ────────────────────────────────────────

def test_no_superseded_scripts():
    """Scripts that have been replaced by a more general version should not exist."""
    superseded = {
        "add_cellpose_prior.py": "add_nuclear_prior.py",
    }
    for old, replacement in superseded.items():
        assert not (SCRIPTS_DIR / old).exists(), (
            f"scripts/{old} is superseded by scripts/{replacement} — delete it"
        )
