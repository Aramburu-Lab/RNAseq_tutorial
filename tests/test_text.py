"""Page text lives in src/ctx/report/content/*.toml: read strictly, all of it used, editable.

The riboseq status report's split, made Python: its test_build.js checked that the page
assembled from the content files; these check the same from both sides of the files.
"""

import re
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

from ctx.report import text

ROOT = Path(__file__).resolve().parents[1]
CONTENT = ROOT / "src" / "ctx" / "report" / "content"
CTX = Path(shutil.which("ctx") or "ctx")


def _paths(table: dict, prefix: str = "") -> list[str]:
    """Dotted path of every value; arrays and their tables are one path each."""
    paths = []
    for key, value in table.items():
        path = f"{prefix}{key}"
        paths += _paths(value, f"{path}.") if isinstance(value, dict) else [path]
    return paths


def _code() -> str:
    """Every file that can read a key: the Python sources, the page templates and the scripts.

    The scripts count because the course preview draws the pages of modules not built yet
    from their content files (scripts/python/build_preview.py).
    """
    files = [*(ROOT / "src" / "ctx").rglob("*"), *(ROOT / "scripts" / "python").glob("*.py")]
    return "\n".join(p.read_text() for p in files if p.suffix in {".py", ".j2"})


def test_every_key_is_read_somewhere():
    """A key no code names is text no page shows: edited, it would change nothing.

    A key counts as read when its path, or a table above it that the code reads whole, is
    named as a string ("machine.title", f"sections.{id}") or in a template (ui.calls.fc).
    """
    code = _code()
    unused = []
    for file in sorted(CONTENT.glob("*.toml")):
        for path in _paths(tomllib.loads(file.read_text())):
            parts = path.split(".")
            prefixes = [".".join(parts[:i]) for i in range(len(parts), 0, -1)]
            named = any(
                re.search(rf"""(["']|\bui\.){re.escape(p)}(["'.{{]|\b)""", code) for p in prefixes
            )
            if not named:
                unused.append(f"{file.name}: {path}")
    assert unused == []


def test_placeholders_are_strict(monkeypatch, tmp_path):
    (tmp_path / "t.toml").write_text('line = "{n} of {total} runs"\nplain = "none"\n')
    monkeypatch.setenv(text.CONTENT_ENV, str(tmp_path))
    text.load.cache_clear()
    try:
        assert text.fill("t", "line", n=2, total=5, spare=1) == "2 of 5 runs"
        with pytest.raises(text.TextError, match=r"t\.toml 'line' needs \['total'\]"):
            text.fill("t", "line", n=2)
        with pytest.raises(text.TextError, match=r"t\.toml has no 'missing'"):
            text.word("t", "missing")
        assert text.word("t", "plain") == "none"
    finally:
        text.load.cache_clear()


def test_an_edited_copy_changes_the_page(tmp_path):
    """CTX_CONTENT is read with no fallback: the packaged files are out of reach."""
    copy = tmp_path / "content"
    shutil.copytree(CONTENT, copy)
    edited = (
        (copy / "setup.toml")
        .read_text()
        .replace("Is this computer ready for the course", "Is this laptop ready for class")
    )
    (copy / "setup.toml").write_text(edited)
    env = {"CTX_CONTENT": str(copy), "PATH": str(CTX.parent)}
    root = tmp_path / "root"
    root.mkdir()
    done = subprocess.run([str(CTX), "--root", str(root), "setup"], cwd=root, env=env,
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert done.returncode == 0, done.stderr
    (page,) = root.rglob("*_setup_environment_v1.html")
    html = page.read_text()
    assert "Is this laptop ready for class" in html
    assert "Is this computer ready for the course" not in html
    (copy / "ui.toml").unlink()
    gone = subprocess.run([str(CTX), "--root", str(root), "setup"], cwd=root, env=env,
                          capture_output=True, text=True, check=False)  # fmt: skip
    assert gone.returncode != 0 and "ui.toml" in gone.stderr


def test_no_ui_key_shadows_a_dict_method():
    """Templates read ``ui`` keys as attributes, and Jinja finds a dict's own method first: a
    key named ``clear`` printed "<built-in method clear of dict object>" on the brush button.
    The other files are read by dotted path in Python, where a key named ``keys`` is safe.
    """
    methods = set(dir(dict))
    ui = tomllib.loads((CONTENT / "ui.toml").read_text())
    assert [path for path in _paths(ui) if methods & set(path.split("."))] == []
