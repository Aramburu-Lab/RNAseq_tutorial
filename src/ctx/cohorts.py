"""The one reader of a cohort bundle (``conf/cohorts/<cancer>.toml``, TODO.md Cohort bundles).

A bundle lists a cancer's dataset files by Step and module, with each one's laptop cap and
teaching note. Contrasts, covariates, sanity rules and worked-example genes stay in the dataset
files, so each is defined once. ``ctx --cohort <cancer>`` selects a bundle, ``CTX_COHORT`` by
default; a figure's stable ID begins with the cohort it was drawn for (``breast_m5_volcano``).
Every schema error is reported at once, as ``ctx.dataset.load`` reports a dataset file's.
"""

import os
import re
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ctx.io.naming import MODULES

# The bundles' folder: $CTX_COHORTS, else the repository's own.
FOLDER_ENV = "CTX_COHORTS"
COHORT_ENV = "CTX_COHORT"
# The cohort of a dataset no bundle lists, such as a student's own: its figures stay reusable.
OWN = "own"
# A cohort name is one word, as a figure ID's first part must be (``<cohort>_<module>_...``).
NAME = re.compile(r"[a-z][a-z0-9]*")
# Step-1 immune cells, Step-2 another blood cancer, Step-3 tissue (TODO.md Goal).
STEPS = (1, 2, 3)
# A laptop cap is the samples kept per group on the laptop profile, or every sample.
ALL = "all"
DATASET_KEYS = {"slug", "step", "modules", "laptop_cap", "note"}
COHORT_KEYS = {"name", "label", "notes"}
# Modules a dataset is run through; setup and report take every dataset, so a bundle names none.
RUN_MODULES = tuple(m for m, _ in MODULES if m not in {"setup", "report", "lite_build"})


class CohortError(ValueError):
    """A cohort bundle that cannot be used; the message lists every problem."""


@dataclass(frozen=True)
class Entry:
    """One dataset of a cohort: its file, Step, the modules it is run through, cap and note."""

    slug: str
    path: Path
    step: int
    modules: tuple[str, ...]
    laptop_cap: int | str
    note: str


@dataclass(frozen=True)
class Cohort:
    """A cancer's bundle: its name, label, teaching notes and dataset entries."""

    name: str
    label: str
    notes: str
    datasets: tuple[Entry, ...]

    def entry(self, slug: str) -> Entry | None:
        """The entry of one dataset slug, or None when this cohort does not list it."""
        return next((e for e in self.datasets if e.slug == slug), None)


def folder() -> Path:
    """Where the bundles are: ``$CTX_COHORTS``, else ``conf/cohorts`` of the repository."""
    env = os.environ.get(FOLDER_ENV)
    return Path(env) if env else Path(__file__).resolve().parents[2] / "conf" / "cohorts"


def _text(table: Mapping[str, Any], key: str, where: str, problems: list[str]) -> str:
    """A required, non-blank string, or a problem recorded."""
    value = table.get(key)
    if not isinstance(value, str) or not value.strip():
        problems.append(f"{where}.{key}: expected non-blank text, found {value!r}")
        return ""
    return value


def _modules(value: Any, where: str, problems: list[str]) -> tuple[str, ...]:
    """A non-empty list of run modules, each once."""
    if not isinstance(value, list) or not value or not all(isinstance(m, str) for m in value):
        problems.append(f"{where}.modules: expected a list of modules, found {value!r}")
        return ()
    problems += [f"{where}.modules: {m!r} is not one of {', '.join(RUN_MODULES)}"
                 for m in value if m not in RUN_MODULES]  # fmt: skip
    if len(set(value)) != len(value):
        problems.append(f"{where}.modules: a module listed twice")
    return tuple(value)


def _cap(value: Any, where: str, problems: list[str]) -> int | str:
    """``"all"`` or a positive whole number of samples per group."""
    whole = isinstance(value, int) and not isinstance(value, bool) and value > 0
    if not (whole or value == ALL):
        problems.append(f"{where}.laptop_cap: expected {ALL!r} or a positive integer, "
                        f"found {value!r}")  # fmt: skip
    return value if whole or value == ALL else ALL


def _entry(raw: Any, i: int, datasets: Path, problems: list[str]) -> Entry | None:
    """One ``[[datasets]]`` entry, its dataset file checked to exist and to carry its slug."""
    where = f"datasets[{i}]"
    if not isinstance(raw, dict):
        problems.append(f"{where}: expected a table")
        return None
    problems += [f"{where}.{k}: unknown key" for k in sorted(set(raw) - DATASET_KEYS)]
    slug = _text(raw, "slug", where, problems)
    step = raw.get("step")
    if step not in STEPS or isinstance(step, bool):
        problems.append(f"{where}.step: expected one of {STEPS}, found {step!r}")
    path = datasets / f"{slug}.toml"
    if slug and not path.is_file():
        problems.append(f"{where}.slug: no dataset file {path}")
    elif slug:
        found = tomllib.loads(path.read_text()).get("dataset", {}).get("slug")
        if found != slug:
            problems.append(f"{where}.slug: {path.name} declares slug {found!r}")
    modules = _modules(raw.get("modules"), where, problems)
    cap, note = _cap(raw.get("laptop_cap"), where, problems), _text(raw, "note", where, problems)
    return Entry(slug, path, step if isinstance(step, int) else 0, modules, cap, note)


def load(path: Path, datasets: Path | None = None) -> Cohort:
    """Read and check one bundle.

    Args:
        path: The TOML file; its stem is the cohort's name.
        datasets: Where the dataset files are; ``conf/datasets`` beside the bundles' folder.

    Returns:
        The parsed cohort.

    Raises:
        CohortError: Listing every problem found, not only the first.
    """
    datasets = datasets or path.resolve().parent.parent / "datasets"
    try:
        raw = tomllib.loads(path.read_text())
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise CohortError(f"{path}: {error}") from error
    problems = [f"unknown table [{t}]" for t in sorted(set(raw) - {"cohort", "datasets"})]
    head = raw.get("cohort", {})
    problems += [f"cohort.{k}: unknown key" for k in sorted(set(head) - COHORT_KEYS)]
    name = _text(head, "name", "cohort", problems)
    if name and (name != path.stem or not NAME.fullmatch(name)):
        problems.append(f"cohort.name: {name!r} must be one lower-case word, the file's stem")
    label, notes = (
        _text(head, "label", "cohort", problems),
        _text(head, "notes", "cohort", problems),
    )
    listed = raw.get("datasets", [])
    if not isinstance(listed, list) or not listed:
        problems.append("[[datasets]]: a bundle lists one dataset file or more")
        listed = []
    entries = [e for i, d in enumerate(listed) if (e := _entry(d, i, datasets, problems))]
    slugs = [e.slug for e in entries if e.slug]
    twice = sorted({s for s in slugs if slugs.count(s) > 1})
    problems += [f"datasets: {s!r} listed twice" for s in twice]
    if problems:
        raise CohortError(f"{path}:\n  " + "\n  ".join(problems))
    return Cohort(name, label, notes, tuple(entries))


def bundles() -> dict[str, Cohort]:
    """Every bundle in ``folder()``, by name; a bad one raises."""
    return {p.stem: load(p) for p in sorted(folder().glob("*.toml"))}


def choose(name: str | None) -> None:
    """Make ``name`` this process's cohort (``ctx --cohort``), checked now; None leaves it.

    It is set as ``$CTX_COHORT``, so the modules ``ctx run`` starts draw for the same cohort.

    Raises:
        CohortError: If no bundle has that name, or the bundle has a problem.
    """
    if name is None:
        return
    path = folder() / f"{name}.toml"
    if not path.is_file():
        known = ", ".join(p.stem for p in sorted(folder().glob("*.toml")))
        raise CohortError(f"no cohort {name!r} in {folder()}; known: {known}")
    load(path)
    os.environ[COHORT_ENV] = name


def active() -> str | None:
    """The cohort ``ctx --cohort`` or ``$CTX_COHORT`` chose, or None."""
    return os.environ.get(COHORT_ENV) or None


def _slugs(path: Path) -> set[str]:
    """The slugs a bundle lists, read without checking it, as a figure's ID needs them."""
    listed = tomllib.loads(path.read_text()).get("datasets", [])
    if not isinstance(listed, list):
        return set()
    return {str(d.get("slug")) for d in listed if isinstance(d, dict)}


def cohort_of(slug: str) -> str:
    """The cohort a dataset's figures are named for.

    The chosen cohort when there is one; else the one bundle that lists the slug; else
    ``OWN``, so a dataset no bundle lists still has stable figure IDs.

    Raises:
        CohortError: If two bundles list the slug and none was chosen.
    """
    if chosen := active():
        return chosen
    if not folder().is_dir():
        return OWN
    listing = [p.stem for p in sorted(folder().glob("*.toml")) if slug in _slugs(p)]
    if len(listing) > 1:
        raise CohortError(f"{slug}: listed by {', '.join(listing)}; choose one with --cohort")
    return listing[0] if listing else OWN
