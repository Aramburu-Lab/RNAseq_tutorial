"""A figure's ``FigureSpec``, built when a module writes it into a ctx root.

The module's own call to ``ctx.viz.save.write_figure`` names the figure by its stem alone, so
the spec is found from that: the stem's module and description, the dataset file this process
read for the dataset slug it begins with, and the kind whose figures that file makes. Each
input is the newest table of its module for the dataset, as the module itself read or wrote
it, hashed now. A figure with no kind, or an input that cannot be found, has no spec.
"""

import hashlib
import sys
from pathlib import Path

import polars as pl

import ctx
from ctx import cohorts, runtime
from ctx.dataset import Dataset, DatasetError, load, loaded
from ctx.figures import kinds
from ctx.figures.kind import DATASET_FILE, QC_REGISTRY, REGISTRY_OF, RESOURCE, Kind
from ctx.figures.spec import FigureSpec, Input, figure_id
from ctx.io.naming import MODULES
from ctx.report import text

FOLDERS = {module: f"{number:02d}_{module}" for module, number in MODULES}


def sha256(path: Path) -> str:
    """The SHA-256 of a file's bytes, as ``sha256sum`` prints it."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def root_of(directory: Path) -> Path | None:
    """The ctx root a module folder sits in, or None for any other folder."""
    return directory.parent.parent if directory.name in FOLDERS.values() else None


def checkout() -> Path | None:
    """The repository the running ctx was imported from, or None for an installed one."""
    if not runtime.source().startswith("git "):
        return None
    return Path(ctx.__file__).resolve().parents[2]


def _module_and_description(figure_stem: str) -> tuple[str, str] | None:
    """The module token and the description of a stem ``<day>_<module>_<desc>_v<N>``."""
    body = figure_stem[11:].rsplit("_v", 1)[0]
    for module in sorted(FOLDERS, key=len, reverse=True):
        if body.startswith(f"{module}_"):
            return module, body[len(module) + 1 :]
    return None


def _located(path: Path, root: Path, repo: Path | None) -> tuple[str, str]:
    """(source, base) of an input: under the root, under the repository, or a shown path."""
    resolved = path.resolve()
    for base, top in (("root", root.resolve()), ("repo", repo)):
        if top is not None and resolved.is_relative_to(top):
            return resolved.relative_to(top).as_posix(), base
    return runtime.shown_path(resolved), "path"


def input_paths(kind: Kind, params: object, spec: Dataset, root: Path) -> dict[str, Path]:
    """Every input of one figure, found as its module reads it.

    Raises:
        DatasetError: If a table is missing, or the dataset file was not read here.
    """
    from ctx.de.run import newest  # noqa: PLC0415 - de imports ctx.viz.save, which imports this
    from ctx.qc import registry  # noqa: PLC0415 - same import order

    found = {}
    for name, (module, what) in kind.inputs(params, spec).items():
        if (module, what) == DATASET_FILE:
            path = loaded().get(spec.slug)
            if path is None:
                raise DatasetError(f"{spec.slug}: no dataset file read in this process")
            found[name] = path
        elif (module, what) == QC_REGISTRY:
            found[name] = registry.resolve(None)
        elif module == REGISTRY_OF:
            options = newest(root / "outputs" / FOLDERS[what], spec.slug, "figure_options")
            found[name] = Path(pl.read_parquet(options)["registry"][0])
        elif module == RESOURCE:
            found[name] = _fetched(root, spec, what)
        else:
            found[name] = newest(root / "outputs" / FOLDERS[module], spec.slug, what)
    return found


def _fetched(root: Path, spec: Dataset, name: str) -> Path:
    """A resource of the dataset file where ``ctx setup`` put it under the root.

    Raises:
        DatasetError: If the file declares no such resource, or it was not fetched.
    """
    found = [r.dest for r in spec.resources if r.name == name]
    path = root / found[0] if found else None
    if path is None or not path.is_file():
        raise DatasetError(f"{spec.slug}: resource {name} not fetched to {path}")
    return path


def written(tables: tuple[tuple[str, str], ...], spec: Dataset, root: Path) -> bool:
    """Whether the module wrote any of these tables for the dataset, on any day."""
    from ctx.de.run import newest  # noqa: PLC0415 - de imports ctx.viz.save, which imports this

    for module, what in tables:
        try:
            newest(root / "outputs" / FOLDERS[module], spec.slug, what)
        except DatasetError:
            continue
        return True
    return False


def describe(figure_stem: str, directory: Path, data_csv: Path) -> FigureSpec | None:
    """The spec of a figure a module just wrote into ``directory``, or None when it has none.

    Args:
        figure_stem: The figure's stem.
        directory: The module folder it was written to.
        data_csv: Its data table, already written.
    """
    root, parts = root_of(directory), _module_and_description(figure_stem)
    if root is None or parts is None:
        return None
    module, description = parts
    slugs = [s for s in sorted(loaded(), key=len, reverse=True) if description.startswith(f"{s}_")]
    for slug in slugs:
        spec = load(loaded()[slug])
        if (match := kinds.matched(module, description[len(slug) + 1 :], spec)) is None:
            continue
        kind, params = match
        if written(kind.unless, spec, root):
            return None
        try:
            paths = input_paths(kind, params, spec, root)
        except DatasetError:
            return None
        repo = checkout()
        inputs = tuple(Input(name, path.name, *_located(path, root, repo), sha256(path))
                       for name, path in paths.items())  # fmt: skip
        function = kind.plotting.split(":")[0]
        return FigureSpec(
            figure_id=figure_id(cohorts.cohort_of(slug), module, kind.figure), dataset=slug,
            module=module, params=params, inputs=inputs, plotting=kind.plotting,
            plotting_file=f"src/{function.replace('.', '/')}.py", stem=figure_stem,
            data_sha256=sha256(data_csv), root=runtime.shown_path(root.resolve()),
            commit=runtime.source(), version=ctx.__version__,
            upstream=runtime.retypable(sys.argv[1:]),
        )  # fmt: skip
    return None


def readme(spec: FigureSpec) -> str:
    """The figure rows appended to the figure's ``_README.md``."""
    drawn = _readme_word("drawn_value", plotting=spec.plotting, file=spec.plotting_file)
    redraw = _readme_word("redraw_value", figure_id=spec.figure_id)
    rows = [("id", f"`{spec.figure_id}`"), ("drawn", drawn), ("redraw", redraw)]
    lines = [f"## {_readme_word('title')}", "", "| Field | Value |", "| --- | --- |"]
    lines += [f"| {_readme_word(key)} | {value} |" for key, value in rows]
    return "\n".join([*lines, ""])


def _readme_word(key: str, **values: str) -> str:
    """One word of the README's figure rows, from the content file."""
    return text.fill("figure", f"readme.{key}", **values)
