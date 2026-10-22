"""The one reader of an edition (``conf/editions/<name>.toml``, TODO.md ED.1, D-110).

An edition is the course as one audience gets it: the modules it teaches and the dataset files
it runs, from one codebase. ``ctx --edition <name>`` selects one, ``CTX_EDITION`` by default,
else ``full``. ``"all"`` in an edition's ``modules`` or ``datasets`` means every course module
(``cohorts.RUN_MODULES``) or every dataset file, as ``full`` lists both, so a module or dataset
added to the course joins it with no second list to edit. An installed package away from the
repository has no editions folder; it runs ``full`` by that same rule (``named``) and refuses
any other name. Every caller that lists modules or datasets asks the selected edition
(``active``) through ``has_module``, ``holds``, ``has_dataset`` and ``refusal``, so the rule is
decided here and nowhere else.

Not a cohort bundle: a bundle is one cancer's datasets and names the figure IDs they draw
(``<cohort>_m5_volcano``), each dataset in exactly one bundle. An edition crosses cancers, the
core edition Step-1, breast and lung, so as a bundle it would rename those figures or make a
slug's cohort ambiguous. It reuses the bundles' module names, folder layout and selection.
"""

import os
import re
import tomllib
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ctx.cohorts import NAME, RUN_MODULES
from ctx.report import text

FOLDER_ENV = "CTX_EDITIONS"
EDITION_ENV = "CTX_EDITION"
DEFAULT = "full"
# Every course module, or every dataset file of the datasets folder, as full lists both.
ALL = "all"
KEYS = {"name", "label", "notes", "modules", "datasets"}
# Built with every edition's pages, never taught on their own (cohorts.RUN_MODULES leaves them).
ALWAYS = ("setup", "report", "lite_build")
WORDS = "edition"
# An output folder, ``NN_<module>``.
FOLDER = re.compile(r"\d{2}_(?P<module>[a-z_]+)")


class EditionError(ValueError):
    """An edition that cannot be used; the message lists every problem."""


@dataclass(frozen=True)
class Edition:
    """One edition: its name, label, notes, modules in course order and dataset slugs.

    ``datasets`` is None for every dataset file of the folder ``datasets_folder``.
    """

    name: str
    label: str
    notes: str
    modules: tuple[str, ...]
    datasets: tuple[str, ...] | None
    path: Path
    datasets_folder: Path

    def has_module(self, module: str) -> bool:
        """Whether this edition runs and shows a module; setup and report are in every one."""
        return module in self.modules or module in ALWAYS

    def holds(self, folder: str) -> bool:
        """Whether an output folder (``07_proteo``) is this edition's; any other name is."""
        found = FOLDER.fullmatch(folder)
        return found is None or self.has_module(found["module"])

    def has_dataset(self, slug: str) -> bool:
        """Whether a dataset runs here: listed, or the reader's own, with no course file."""
        if self.datasets is None or slug in self.datasets:
            return True
        return not (self.datasets_folder / f"{slug}.toml").is_file()

    def course_slug(self, dataset_file: Path | None) -> str | None:
        """The slug of a course dataset file, its stem; None for a reader's own file."""
        if dataset_file is None or dataset_file.resolve().parent != self.datasets_folder.resolve():
            return None
        return dataset_file.stem

    def dataset_files(self) -> tuple[Path, ...]:
        """The course dataset files of this edition, in its own order or by name for all."""
        if self.datasets is None:
            return tuple(sorted(self.datasets_folder.glob("*.toml")))
        return tuple(self.datasets_folder / f"{slug}.toml" for slug in self.datasets)

    def refusal(self, module: str | None = None, dataset: str | None = None) -> str | None:
        """Why a command may not run here, naming the module or dataset, or None when it may."""
        if module in RUN_MODULES and not self.has_module(module):
            return text.fill(WORDS, "refused.module", module=module, edition=self.name,
                             modules=", ".join(self.modules), default=DEFAULT)  # fmt: skip
        if dataset is not None and not self.has_dataset(dataset):
            return text.fill(WORDS, "refused.dataset", dataset=dataset, edition=self.name,
                             path=self.path)  # fmt: skip
        return None


def folder() -> Path:
    """Where the editions are: ``$CTX_EDITIONS``, else ``conf/editions`` of the repository."""
    env = os.environ.get(FOLDER_ENV)
    return Path(env) if env else Path(__file__).resolve().parents[2] / "conf" / "editions"


def _text(table: Mapping[str, Any], key: str, problems: list[str]) -> str:
    """A required, non-blank string, or a problem recorded."""
    value = table.get(key)
    if not isinstance(value, str) or not value.strip():
        problems.append(f"edition.{key}: expected non-blank text, found {value!r}")
        return ""
    return value


def _names(value: Any, key: str, problems: list[str]) -> tuple[str, ...]:
    """A non-empty list of distinct strings, or a problem recorded."""
    if not isinstance(value, list) or not value or not all(isinstance(v, str) for v in value):
        problems.append(f"edition.{key}: expected a list of names, found {value!r}")
        return ()
    if len(set(value)) != len(value):
        problems.append(f"edition.{key}: a name listed twice")
    return tuple(value)


def _modules(value: Any, problems: list[str]) -> tuple[str, ...]:
    """``"all"``, or course modules, each known, gather among them, in course order."""
    if value == ALL:
        return RUN_MODULES
    modules = _names(value, "modules", problems)
    problems += [f"edition.modules: {m!r} is not one of {', '.join(RUN_MODULES)}"
                 for m in modules if m not in RUN_MODULES]  # fmt: skip
    if modules and "gather" not in modules:
        problems.append("edition.modules: gather is missing; every dataset enters through M1")
    known = [m for m in modules if m in RUN_MODULES]
    if known != sorted(known, key=RUN_MODULES.index):
        problems.append(f"edition.modules: not in course order, {', '.join(RUN_MODULES)}")
    return modules


def _datasets(value: Any, datasets: Path, problems: list[str]) -> tuple[str, ...] | None:
    """``"all"``, or slugs each with a dataset file that declares that slug."""
    if value == ALL:
        return None
    slugs = _names(value, "datasets", problems)
    for slug in slugs:
        path = datasets / f"{slug}.toml"
        if not path.is_file():
            problems.append(f"edition.datasets: no dataset file {path}")
        elif (found := tomllib.loads(path.read_text()).get("dataset", {}).get("slug")) != slug:
            problems.append(f"edition.datasets: {path.name} declares slug {found!r}")
    return slugs


def load(path: Path, datasets: Path | None = None) -> Edition:
    """Read and check one edition.

    Args:
        path: The TOML file; its stem is the edition's name.
        datasets: Where the dataset files are; ``conf/datasets`` beside the editions' folder.

    Returns:
        The parsed edition.

    Raises:
        EditionError: Listing every problem found, not only the first.
    """
    datasets = datasets or path.resolve().parent.parent / "datasets"
    try:
        raw = tomllib.loads(path.read_text())
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise EditionError(f"{path}: {error}") from error
    problems = [f"unknown table [{t}]" for t in sorted(set(raw) - {"edition"})]
    head = raw.get("edition", {})
    if not isinstance(head, dict):
        raise EditionError(f"{path}:\n  [edition]: expected a table")
    problems += [f"edition.{k}: unknown key" for k in sorted(set(head) - KEYS)]
    name = _text(head, "name", problems)
    if name and (name != path.stem or not NAME.fullmatch(name)):
        problems.append(f"edition.name: {name!r} must be one lower-case word, the file's stem")
    label, notes = _text(head, "label", problems), _text(head, "notes", problems)
    modules = _modules(head.get("modules"), problems)
    slugs = _datasets(head.get("datasets"), datasets, problems)
    if problems:
        raise EditionError(f"{path}:\n  " + "\n  ".join(problems))
    return Edition(name, label, notes, modules, slugs, path, datasets)


def choose(name: str | None) -> None:
    """Make ``name`` this process's edition (``ctx --edition``), checked now; None leaves it.

    It is set as ``$CTX_EDITION``, so the modules ``ctx run`` starts run under the same one.

    Raises:
        EditionError: If no edition has that name, or the edition has a problem.
    """
    if name is None:
        return
    named(name)
    os.environ[EDITION_ENV] = name


def named(name: str) -> Edition:
    """The edition of that name in ``folder()``.

    Raises:
        EditionError: If there is none, naming the ones there are, or it has a problem.
    """
    path = folder() / f"{name}.toml"
    if name == DEFAULT and not folder().is_dir():
        return Edition(DEFAULT, DEFAULT, "", RUN_MODULES, None, path, path.parent)
    if not path.is_file():
        known = ", ".join(p.stem for p in sorted(folder().glob("*.toml")))
        raise EditionError(f"no edition {name!r} in {folder()}; known: {known}")
    return load(path)


def active() -> Edition:
    """The edition ``ctx --edition`` or ``$CTX_EDITION`` chose, else ``full``.

    Raises:
        EditionError: If the chosen one is not there or has a problem.
    """
    return named(os.environ.get(EDITION_ENV) or DEFAULT)
