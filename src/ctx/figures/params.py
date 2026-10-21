"""The one reader of a figure's params file, listing every problem at once (D-020).

A params file names the figure by its ID, fills its kind's parameter dataclass, names each
input file in the reader's folder and, where the panel wrote them, each input's SHA-256. The
dataset file is itself an input; its contrasts, lessons and gene sets bound the parameters.
"""

import dataclasses
import hashlib
import json
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from ctx.dataset import Dataset, DatasetError, load
from ctx.figures import kinds
from ctx.figures.kind import Kind, field_names, unknown
from ctx.figures.spec import FigureSpec, parse
from ctx.report import text

TOP = {"figure", "inputs", "sha256"}


class ParamsError(ValueError):
    """A params file that cannot be used; the message lists every problem."""


@dataclass(frozen=True)
class Request:
    """A checked params file: the kind, its filled parameters, the inputs and the dataset file."""

    figure_id: str
    kind: Kind
    params: Any
    inputs: dict[str, Path]
    spec: Dataset


def _typed(value: Any, kind: Any) -> bool:
    """Whether a TOML value fits a parameter's type; a whole number is a float, never a bool."""
    if isinstance(value, bool):
        return kind is bool
    return isinstance(value, kind) or (kind is float and isinstance(value, int))


def _fields(raw: dict[str, Any], kind: Kind, problems: list[str]) -> dict[str, Any]:
    """The parameter values, each present unless it has a default, each of its field's type."""
    values = {}
    for f in dataclasses.fields(kind.params):
        if f.name not in raw:
            if f.default is dataclasses.MISSING:
                problems.append(f"{f.name}: missing")
            continue
        if not _typed(raw[f.name], f.type):
            problems.append(f"{f.name}: expected {getattr(f.type, '__name__', f.type)}, "
                            f"found {raw[f.name]!r}")  # fmt: skip
            continue
        values[f.name] = float(raw[f.name]) if f.type is float else raw[f.name]
    known = {*TOP, *field_names(kind.params)}
    problems += [f"{key}: not a parameter of {kind.figure}" for key in sorted(set(raw) - known)]
    return values


def bounds(kind: Kind, params: Any, spec: Dataset) -> list[str]:
    """Problems with parameter values the dataset file does not declare, or out of range."""
    names = [c.name for c in spec.contrasts]
    found = []
    # A kind the dataset is not offered has inputs its module never wrote (R1.3b).
    if not kinds.offered(kind, spec):
        found.append(f"{kind.figure}: {spec.slug}'s dataset file makes no such figure")
    if params.dataset != spec.slug:
        found.append(f"dataset: {params.dataset!r}, but the dataset file is {spec.slug!r}")
    for field in ("contrast", "x", "y"):
        value = getattr(params, field, None)
        if value is not None:
            found += unknown(field, value, names, "a contrast", spec.slug)
    return found + kind.check(params, spec)


def _files(raw: dict[str, Any], folder: Path, problems: list[str]) -> dict[str, Path]:
    """Each named input, present in the folder, and equal to its hash where one is given."""
    named, hashes = raw.get("inputs", {}), raw.get("sha256", {})
    if not isinstance(named, dict) or not isinstance(hashes, dict):
        problems.append("[inputs] and [sha256] are tables")
        return {}
    found = {}
    for name, file in named.items():
        path = folder / str(file)
        if not path.is_file():
            problems.append(f"inputs.{name}: no file {path}")
            continue
        found[name] = path
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        if name in hashes and hashes[name] != digest:
            problems.append(f"inputs.{name}: {path.name} has SHA-256 {digest}, "
                            f"the params file says {hashes[name]}")  # fmt: skip
    problems += [f"sha256.{n}: names no input" for n in sorted(set(hashes) - set(named))]
    return found


def load_params(path: Path, figure: str) -> Request:
    """Read and check a params file for the figure ``ctx figure`` was asked to draw.

    Args:
        path: The params file; input files are named relative to its folder.
        figure: The ID given on the command line, which the file must name too.

    Returns:
        The checked request.

    Raises:
        ParamsError: Listing every problem found, not only the first.
    """
    try:
        raw = tomllib.loads(path.read_text())
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise ParamsError(f"{path}: {error}") from error
    problems: list[str] = []
    if raw.get("figure") != figure:
        problems.append(f"figure: the file names {raw.get('figure')!r}, the command {figure!r}")
    parsed = parse(figure)
    kind = kinds.kind(parsed[1], parsed[2]) if parsed else None
    if kind is None:
        known = ", ".join(f"<cohort>_m*_{k.figure}" for k in kinds.KINDS)
        problems.append(f"{figure!r}: no figure kind has this ID; known: {known}")
        raise ParamsError(f"{path}:\n  " + "\n  ".join(problems))
    values = _fields(raw, kind, problems)
    inputs = _files(raw, path.parent, problems)
    spec = None
    if "dataset_file" in inputs:
        try:
            spec = load(inputs["dataset_file"])
        except DatasetError as error:
            problems.append(f"inputs.dataset_file: {error}")
    elif "dataset_file" not in raw.get("inputs", {}):
        problems.append("inputs.dataset_file: missing")
    params = None
    if spec is not None and not _missing(values, kind):
        params = kind.params(**values)
        problems += bounds(kind, params, spec)
        expected = set(kind.inputs(params, spec))
        given = set(raw.get("inputs", {}))
        problems += [f"inputs.{n}: missing" for n in sorted(expected - given)]
        problems += [f"inputs.{n}: not an input of {kind.figure}" for n in sorted(given - expected)]
    if problems or spec is None or params is None:
        raise ParamsError(f"{path}:\n  " + "\n  ".join(problems))
    return Request(figure, kind, params, inputs, spec)


def _missing(values: dict[str, Any], kind: Kind) -> bool:
    """Whether a required parameter has no usable value, so the dataclass cannot be filled."""
    return any(f.default is dataclasses.MISSING and f.name not in values
               for f in dataclasses.fields(kind.params))  # fmt: skip


def _comment(key: str) -> str:
    """One comment line of the params file, from the content file."""
    return f"# {text.word('figure', key)}"


def toml_value(value: Any) -> str:
    """A parameter as TOML writes it: a basic string, or a number."""
    return json.dumps(value) if isinstance(value, str) else repr(value)


def _input_key(name: str) -> str:
    """The content key of an input's comment; ``results_<contrast>`` takes ``results_each``."""
    known = set(text.raw("figure", "inputs"))
    if name in known:
        return name
    return f"{name.split('_', 1)[0]}_each"


def params_file(spec: FigureSpec) -> str:
    """The params file ``ctx figure`` reads: the figure, its parameters, inputs and hashes."""
    lines = [_comment("params.figure"), f"figure = {toml_value(spec.figure_id)}"]
    for name, value in dataclasses.asdict(spec.params).items():
        lines += [_comment(f"params.{name}"), f"{name} = {toml_value(value)}"]
    lines += ["", _comment("params.inputs"), "[inputs]"]
    for i in spec.inputs:
        lines += [_comment(f"inputs.{_input_key(i.name)}"), f"{i.name} = {toml_value(i.file)}"]
    lines += ["", _comment("params.sha256"), "[sha256]"]
    lines += [f"{i.name} = {toml_value(i.sha256)}" for i in spec.inputs]
    return "\n".join(lines) + "\n"
