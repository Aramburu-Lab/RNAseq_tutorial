"""The teacher guide: each module's compute time and memory, read from the pages it built.

Every number comes from a module page's Reproduce block, the record of the run that wrote the
page, so the guide never carries a typed time (TODO.md Phase 5, teacher guide). A module no root
holds a page of is listed as not measured. Run ``python -m ctx.report.teacher OUT ROOT...``;
the page text is in ``content/teacher.toml``.
"""

import argparse
import functools
import html
import re
import sys
from collections.abc import Sequence
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter

from ctx import editions, runtime
from ctx.io.naming import MODULES
from ctx.report import text
from ctx.report.builder import (
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    ResultTable,
    write,
)
from ctx.report.course import label, newest
from ctx.report.hub import NOT_MODULES, TITLES

PAGE = "teacher"
# The cells module.html.j2 writes in every Reproduce block; tests/test_budget.py reads them too.
_WALL = re.compile(r'<td class="wall">([\d.]+) s</td>')
_RSS = re.compile(r'<td class="rss-(self|children)">[\d.]+ MiB \((\d+) bytes\)</td>')


@dataclass(frozen=True)
class Measured:
    """One module page's Reproduce values, as the page recorded them."""

    module: str
    dataset: str
    wall_seconds: float
    peak_rss_self: int
    peak_rss_children: int
    profile: str
    page: str


def measured(module: str, page: Path) -> Measured | None:
    """The Reproduce values of one page, or None for a page with no Reproduce block."""
    found = page.read_text(encoding="utf-8", errors="replace")
    wall, rss = _WALL.search(found), dict(_RSS.findall(found))
    if not wall or set(rss) != {"self", "children"}:
        return None
    heading = re.escape(html.escape(text.word("ui", "reproduce.profile")))
    profile = re.search(rf"<th>{heading}</th><td>(.*?)</td>", found)
    return Measured(module, label(page), float(wall.group(1)), int(rss["self"]),
                    int(rss["children"]), html.unescape(profile.group(1)) if profile else "",
                    page.name)  # fmt: skip


def gathered(roots: Sequence[Path]) -> list[tuple[str, list[Measured]]]:
    """Each module of the selected edition in order, with each dataset's newest overview page."""
    found, edition = [], editions.active()
    for module, number in MODULES:
        if module in NOT_MODULES or not edition.has_module(module):
            continue
        suffix = "environment" if module == "setup" else "overview"
        pages = [p for root in roots for p in newest(root / "outputs" / f"{number:02d}_{module}",
                                                    suffix)]  # fmt: skip
        found.append((module, [m for p in pages if (m := measured(module, p))]))
    return found


def track(module: str) -> str:
    """A module's track and its limits, as its own page states them (CONVENTIONS.md, Budgets)."""
    limit = runtime.budget(module)
    name = "A" if module in runtime.TRACK_A_MODULES else "B"
    return text.fill(PAGE, "budget", track=name, minutes=limit.minutes, ram_gb=limit.ram_gb)


def table(found: Sequence[tuple[str, list[Measured]]]) -> ResultTable:
    """One row per measured page, and one per module no root measured."""
    word = functools.partial(text.word, PAGE)
    none = word("times.none")

    def size(n: int) -> str:
        """Bytes as the page's MiB."""
        return text.fill(PAGE, "times.mib", mib=n / 1048576)

    def row(m: Measured, cap: str) -> tuple[str, ...]:
        """One measured page as a table row."""
        wall = text.fill(PAGE, "times.seconds", seconds=m.wall_seconds)
        return (TITLES[m.module], m.dataset, wall, size(m.peak_rss_self),
                size(m.peak_rss_children), m.profile, cap, m.page)  # fmt: skip

    rows = []
    for module, pages in found:
        cap = track(module)
        rows += [row(m, cap) for m in pages] or [(TITLES[module], *(none,) * 5, cap, none)]
    return ResultTable(word("times.title"), text.words(PAGE, "times.columns"), tuple(rows),
                       numeric=(2, 3, 4))  # fmt: skip


def page(roots: Sequence[Path], command: str, started: float) -> ModulePage:
    """The teacher guide over the given roots, its own cost measured at the end.

    Args:
        roots: Folders holding ``outputs/``: ctx roots, or a kit's dataset folders.
        command: The command as a reader retypes it.
        started: ``perf_counter()`` at the command's start.

    Returns:
        The page.
    """
    word = functools.partial(text.word, PAGE)
    found = gathered(roots)
    # M0 names the computer, not a dataset.
    datasets = sorted({m.dataset for module, pages in found if module != "setup" for m in pages})
    count = functools.partial(text.fill, PAGE, "keys.value")
    keys = (
        KeyNumber(text.word(PAGE, "keys.pages"), f"{sum(len(p) for _, p in found):,}"),
        KeyNumber(text.word(PAGE, "keys.modules"),
                  count(n=sum(bool(p) for _, p in found), total=len(found))),
    )  # fmt: skip
    cards = tuple(ParameterCard(c["name"], c["meaning"], c["default"], c["effect"])
                  for c in text.raw(PAGE, "cards"))  # fmt: skip
    group = Group(text.word(PAGE, "times.group"), (table(found),),
                  lead=text.word(PAGE, "times.lead"))  # fmt: skip
    usage = runtime.measure(started)
    return ModulePage(
        title=text.fill(PAGE, "title", dataset=", ".join(datasets) or word("times.none")),
        **text.skeleton(PAGE),
        parameters=cards,
        results=(group,),
        key_numbers=keys,
        kicker=text.word(PAGE, "kicker"),
        reproduce=Reproduction(
            command=command, wall_seconds=usage.wall_seconds, peak_rss_self=usage.peak_rss_self,
            peak_rss_children=usage.peak_rss_children, profile=text.word(PAGE, "times.profile"),
            budget=track("report"),
            versions=runtime.versions(),
        ),
    )  # fmt: skip


def main(argv: Sequence[str] | None = None) -> Path:
    """Write the teacher guide; the command line is ``OUT ROOT...``.

    Returns:
        The written page.
    """
    started = perf_counter()
    parser = argparse.ArgumentParser(prog="python -m ctx.report.teacher")
    parser.add_argument("out", type=Path)
    parser.add_argument("roots", type=Path, nargs="+")
    args = parser.parse_args(argv)
    if missing := [str(r) for r in args.roots if not (r / "outputs").is_dir()]:
        parser.error(f"no outputs/ in {', '.join(missing)}")
    given = list(sys.argv[1:] if argv is None else argv)
    command = " ".join(["python -m ctx.report.teacher", runtime.retyped(given)])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    return write(page(args.roots, command, started), args.out)


if __name__ == "__main__":
    print(main())
