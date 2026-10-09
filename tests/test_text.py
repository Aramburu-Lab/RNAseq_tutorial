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


def _templates() -> str:
    """The page templates, where a content file is read under its own name (course.intro)."""
    return "\n".join(p.read_text() for p in (ROOT / "src/ctx/report/templates").glob("*.j2"))


NAME = re.compile(r"[\w.\-]+")
WORD = re.compile(r"\w")
# A quote, the run of name characters it opens, and the character after the run.
QUOTED = re.compile(r"""["']([\w.\-]*)(?=(.?))""", re.DOTALL)


def _quoted_names(code: str) -> set[str]:
    """Every name ``["']name(["'.{]|\\b)`` matches in ``code``, found in one pass.

    One regex per key over 1.9 MB of code took 57 s for 2,526 keys (M3.14); this index answers
    the same question for any key of name characters, which test_quoted_index_agrees checks.
    """
    names = set()
    for found in QUOTED.finditer(code):
        run, after = found.groups()
        for cut in range(1, len(run) + 1):
            following = run[cut] if cut < len(run) else after
            boundary = bool(WORD.fullmatch(run[cut - 1])) != bool(WORD.fullmatch(following))
            if (following and following in "\"'.{") or boundary:
                names.add(run[:cut])
    return names


def _named_in_code(prefix: str, code: str, names: set[str]) -> bool:
    if NAME.fullmatch(prefix):
        return prefix in names
    return re.search(rf"""["']{re.escape(prefix)}(["'.{{]|\b)""", code) is not None


def test_quoted_index_agrees():
    """The index equals the per-key regex on code built to hit each of its cases.

    Mutations caught: the character after a run consumed (a run opened by it is lost); the end
    of the code taken as no boundary; a dot taken as a word character."""
    code = """x = 'a.b' + "c{d}" + 'e'f' + "g-h" + 'i.' + "j_k(" + 'l"""
    names = _quoted_names(code)
    for prefix in ("a", "a.b", "a.", "c", "e", "f", "g", "g-h", "i", "i.", "j", "j_k", "l", "h"):
        expected = re.search(rf"""["']{re.escape(prefix)}(["'.{{]|\b)""", code) is not None
        assert (prefix in names) is expected, prefix


def test_every_key_is_read_somewhere():
    """A key no code names is text no page shows: edited, it would change nothing.

    A key counts as read when its path, or a table above it that the code reads whole, is
    named as a string ("machine.title", f"sections.{id}") anywhere, or in a template under its
    file's name (ui.calls.fc, course.intro). The name counts only in a template and never
    after a dot, so process.argv or ctx.qc.multiqc_io in code cannot vouch for a dead key.
    """
    code, templates = _code(), _templates()
    names = _quoted_names(code)
    unused = []
    for file in sorted(CONTENT.glob("*.toml")):
        for path in _paths(tomllib.loads(file.read_text())):
            parts = path.split(".")
            prefixes = [".".join(parts[:i]) for i in range(len(parts), 0, -1)]
            named = any(
                _named_in_code(p, code, names)
                or re.search(rf"(?<![\w.]){file.stem}\.{re.escape(p)}\b", templates)
                for p in prefixes
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
        .replace("before any analysis runs on it", "before class starts")
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
    assert "before class starts" in html
    assert "before any analysis runs on it" not in html
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
