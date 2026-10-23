"""Page text: one TOML file per page in ``ctx/report/content``, or in ``$CTX_CONTENT``.

Every sentence a reader sees lives in those files, so text is edited without touching code
(the riboseq status report's split, Rytis 2026-09-29); code names a key and supplies the
numbers. A key the file lacks, or a placeholder the code does not fill, raises an error
naming both: a page never ships with a hole in it. ``content/README.md`` says what is where.
"""

import os
import string
import tomllib
from functools import cache
from importlib.resources import files
from pathlib import Path
from typing import Any

# A folder holding edited copies of the content files; unset reads the packaged ones.
CONTENT_ENV = "CTX_CONTENT"


class TextError(LookupError):
    """A content file lacks a key, or a template needs a value the code did not give."""


@cache
def load(page: str) -> dict[str, Any]:
    """The parsed content file of one page, from ``$CTX_CONTENT`` when that is set.

    Args:
        page: The file's stem, such as ``setup`` or ``ui``.

    Returns:
        The file's tables and strings.
    """
    folder = os.environ.get(CONTENT_ENV)
    if folder:
        return tomllib.loads((Path(folder) / f"{page}.toml").read_text(encoding="utf-8"))
    packaged = files("ctx.report").joinpath("content", f"{page}.toml")
    return tomllib.loads(packaged.read_text(encoding="utf-8"))


def raw(page: str, key: str) -> Any:
    """The value at a dotted key, unfilled: a list, a table, or a template the page fills.

    Args:
        page: The content file's stem.
        key: Dotted path, such as ``machine.title``.

    Returns:
        The value as the file holds it.

    Raises:
        TextError: If the file has no such key.
    """
    node: Any = load(page)
    for part in key.split("."):
        if not isinstance(node, dict) or part not in node:
            raise TextError(f"content/{page}.toml has no {key!r}")
        node = node[part]
    return node


def word(page: str, key: str) -> str:
    """The text at a key, which has no placeholders: a label, heading or fixed sentence.

    Args:
        page: The content file's stem.
        key: Dotted path to a string.

    Returns:
        The text.

    Raises:
        TextError: If the key is missing, not text, or has placeholders to fill.
    """
    return fill(page, key)


def words(page: str, key: str) -> tuple[str, ...]:
    """The list of texts at a key, such as a table's column headings."""
    value = raw(page, key)
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise TextError(f"content/{page}.toml {key!r} is not a list of text")
    return tuple(value)


def placeholders(template: str) -> set[str]:
    """The names a template's ``{placeholders}`` need, format specs stripped."""
    fields = (field for _, field, _, _ in string.Formatter().parse(template) if field)
    return {field.split(".")[0].split("[")[0] for field in fields}


def fill(page: str, key: str, **values: Any) -> str:
    """The text at a key, its placeholders filled from ``values``.

    Args:
        page: The content file's stem.
        key: Dotted path to a string.
        values: One value per placeholder; extra values are allowed, for a key whose
            variants use fewer of them.

    Returns:
        The filled text.

    Raises:
        TextError: If the key is missing or not text, or a placeholder has no value.
    """
    template = raw(page, key)
    if not isinstance(template, str):
        raise TextError(f"content/{page}.toml {key!r} is not text")
    if missing := placeholders(template) - values.keys():
        raise TextError(f"content/{page}.toml {key!r} needs {sorted(missing)}")
    return template.format_map(values)


def p_said(p: float) -> str:
    """A p-value with its relation, as a caption states it: ``= 0.012``, or ``< 5e-324``.

    A test whose p underflowed to 0.0 has a p below the smallest positive double, 4.9e-324;
    "p = 0" says something no test can find (Rytis 2026-10-23, M8's filter caption).
    """
    return "< 5e-324" if p == 0 else f"= {p:.2g}"


def lines(page: str, key: str, **values: Any) -> tuple[str, ...]:
    """The list of texts at a key, each filled, as ``fill`` fills one.

    Raises:
        TextError: If the key is missing or not a list of text, or a placeholder has no value.
    """
    found = words(page, key)
    if missing := set().union(*map(placeholders, found)) - values.keys():
        raise TextError(f"content/{page}.toml {key!r} needs {sorted(missing)}")
    return tuple(line.format_map(values) for line in found)


def _filled(page: str, key: str, template: str, values: dict[str, Any]) -> str:
    """One section's text filled from ``values``; a placeholder with no value is refused."""
    if missing := placeholders(template) - values.keys():
        raise TextError(f"content/{page}.toml {key!r} needs {sorted(missing)}")
    return template.format_map(values)


def skeleton(page: str, variant: str | None = None, /, **values: Any) -> dict[str, Any]:
    """The fixed parts of a module page, as ``ModulePage`` takes them, numbers filled.

    Args:
        page: The content file's stem.
        variant: The dotted key of a table whose sections replace the page's own, for a
            design the page's text was not written for (M7's ``variants.tumour_normal``),
            named whole so the key test sees it read; None, or a section the table lacks,
            keeps the page's own.
        values: One value per placeholder in the sections, terms, questions and pitfalls.

    Returns:
        The page's purpose, input, method, terms, interpretation, questions and pitfalls,
        under the ``ModulePage`` field names, which are the section keys (P2P.3). A page whose
        questions are a table of variants, asked by what its data hold, gets none here. With
        no ``values`` the texts are left as written, for a page that fills its own sections.
    """
    from ctx.report.builder import SelfCheck  # noqa: PLC0415 - builder imports this module

    def one(key: str, template: str) -> str:
        """One text of the section at ``key``, filled when the caller gave values."""
        return _filled(page, key, template, values) if values else template

    own = raw(page, variant) if variant else {}

    def section(key: str) -> Any:
        """The variant's text of a section where it has one, else the page's."""
        return own[key] if key in own else raw(page, key)

    asked = section("questions")
    return {
        **{k: one(k, section(k)) for k in ("purpose", "input", "method", "interpretation")},
        "terms": tuple((t["term"], one("terms", t["meaning"])) for t in section("terms")),
        "questions": tuple(
            SelfCheck(one("questions", t["prompt"]), one("questions", t["answer"])) for t in asked
        )
        if isinstance(asked, list)
        else (),  # fmt: skip
        "pitfalls": tuple(one("pitfalls", m) for m in section("pitfalls")),
    }
