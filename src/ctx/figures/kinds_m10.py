"""M10's figure kinds: drawn from the patient table and the tables ``ctx clinical`` wrote.

Every kind reads the patient table and the run's options; those that score a gene list also
read the stored values of the genes the page carries (``ctx.clinical.figure_tables``). The ER
curves are drawn only where the dataset file declares ER, and the box is the one it declares.
"""

from dataclasses import dataclass
from typing import Any

from ctx.dataset import ENDPOINTS, Dataset
from ctx.figures.kind import DATASET_FILE, Inputs, Kind, Plain, unknown


@dataclass(frozen=True)
class OneEndpoint:
    """ER-negative against ER-positive under one endpoint."""

    dataset: str
    endpoint: str


@dataclass(frozen=True)
class OneBox:
    """The build genes' expression by the two-level variable the dataset file declares."""

    dataset: str
    box: str


def _read(_: Any, __: Dataset) -> Inputs:
    """What every M10 figure reads: the dataset file, the patient table and the options."""
    return {"dataset_file": DATASET_FILE, "patients": ("clinical", "patients"),
            "run_options": ("clinical", "figure_options")}  # fmt: skip


def _scored(p: Any, spec: Dataset) -> Inputs:
    """``_read``'s inputs and the carried genes' stored values."""
    return {**_read(p, spec), "panel_values": ("clinical", "panel_values")}


def _clinical(spec: Dataset) -> list[Plain]:
    """The one figure of a kind every M10 dataset file makes."""
    return [Plain(spec.slug)] if spec.clinical else []


def _endpoints(spec: Dataset) -> list[OneEndpoint]:
    """One ER figure per endpoint, where the dataset file declares ER."""
    if not spec.clinical or "er" not in spec.clinical.variables:
        return []
    return [OneEndpoint(spec.slug, e) for e in ENDPOINTS]


def _boxes(spec: Dataset) -> list[OneBox]:
    """The dataset file's box."""
    return [OneBox(spec.slug, spec.clinical.box)] if spec.clinical else []


def _endpoint_check(p: OneEndpoint, spec: Dataset) -> list[str]:
    """The endpoint is one of the CDR's four, and the dataset file declares ER."""
    found = unknown("endpoint", p.endpoint, ENDPOINTS, "an endpoint of", spec.slug)
    if not _endpoints(spec):
        found.append(f"endpoint: {spec.slug} declares no er variable, so M10 draws no ER curves")
    return found


def _box_check(p: OneBox, spec: Dataset) -> list[str]:
    """The box is the one the dataset file declares, the only one M10 draws."""
    declared = [b.box for b in _boxes(spec)]
    return unknown("box", p.box, declared, "the clinical.box", spec.slug)


DRAW = "ctx.figures.draw_m10"
PAGE = "ctx.clinical.page"
KINDS = (
    Kind("clinical", "er", OneEndpoint, f"{PAGE}:_er_check", f"{DRAW}:er", _endpoints,
         lambda p: f"er_{p.endpoint.lower()}", _read, _endpoint_check),
    Kind("clinical", "expression", OneBox, f"{PAGE}:_expression", f"{DRAW}:expression", _boxes,
         lambda p: f"expression_{p.box}", _scored, _box_check),
    Kind("clinical", "km_score", Plain, f"{PAGE}:_scored", f"{DRAW}:km_score", _clinical,
         lambda _: "km_score", _scored),
    Kind("clinical", "scan_score", Plain, f"{PAGE}:_scored", f"{DRAW}:scan_score", _clinical,
         lambda _: "scan_score", _scored),
    Kind("clinical", "competing", Plain, f"{PAGE}:_competing", f"{DRAW}:competing", _clinical,
         lambda _: "competing", _read),
    Kind("clinical", "power", Plain, f"{PAGE}:_power", f"{DRAW}:power", _clinical,
         lambda _: "power", _scored),
    Kind("clinical", "confounding", Plain, f"{PAGE}:_context_groups", f"{DRAW}:confounding",
         _clinical, lambda _: "confounding", _read),
)  # fmt: skip
