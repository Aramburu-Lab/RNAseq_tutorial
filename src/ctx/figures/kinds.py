"""Which figures have a spec: every module's kinds, each with its parameters, inputs and drawer.

A kind names every figure of its sort the dataset file can make (``made``), so a module's
figure stem is matched against what the file declares, never parsed by guesswork: contrast
names and dataset slugs both hold underscores. Drawers are named by path and imported when
used, since the modules that draw import ``ctx.viz.save``, which imports this file.
"""

from typing import Any

from ctx.dataset import Dataset
from ctx.figures import (
    kinds_m1,
    kinds_m2,
    kinds_m3,
    kinds_m4,
    kinds_m5,
    kinds_m6,
    kinds_m7,
    kinds_m8,
    kinds_m9,
    kinds_m10,
)
from ctx.figures.kind import Kind
from ctx.pipeline import route_modules

KINDS: tuple[Kind, ...] = (
    *kinds_m1.KINDS,
    *kinds_m2.KINDS,
    *kinds_m3.KINDS,
    *kinds_m4.KINDS,
    *kinds_m5.KINDS,
    *kinds_m6.KINDS,
    *kinds_m7.KINDS,
    *kinds_m8.KINDS,
    *kinds_m9.KINDS,
    *kinds_m10.KINDS,
)


def kind(module: str, figure: str) -> Kind | None:
    """The kind of one module's figure name, or None when it has no spec yet."""
    return next((k for k in KINDS if (k.module, k.figure) == (module, figure)), None)


def offered(k: Kind, spec: Dataset) -> bool:
    """Whether the dataset is offered figures of ``k``: its route runs k's module, which makes one.

    The route is ``pipeline.route_modules``'s alone; a kind the dataset file lists none of
    (``made``) is offered none, unless its figures are named on the command line (``parsed``).
    """
    return k.module in route_modules(spec) and (k.parsed is not None or bool(k.made(spec)))


def matched(module: str, what: str, spec: Dataset) -> tuple[Kind, Any] | None:
    """The kind and parameters of a figure a module writes as ``<dataset>_<what>``, or None."""
    for k in KINDS:
        if k.module != module or not offered(k, spec):
            continue
        found = [k.parsed(what, spec)] if k.parsed else k.made(spec)
        for params in found:
            # A parsed description must read back as itself, or it names another kind's.
            if params is not None and k.what(params) == what:
                return k, params
    return None
