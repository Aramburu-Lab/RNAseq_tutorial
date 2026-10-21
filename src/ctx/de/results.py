"""Where M5's results are found: the one finder its readers and ``ctx run``'s plan share.

Light on purpose: ``ctx run --check`` asks it before any module is imported.
"""

import re
from pathlib import Path

from ctx.dataset import Dataset, DatasetError
from ctx.io.naming import MODULES


def newest(directory: Path, dataset: str, what: str) -> Path:
    """The newest ``<date>_<module>_<dataset>_<what>_v<N>.parquet`` by date, then version.

    Raises:
        DatasetError: If there is none.
    """
    pattern = re.compile(
        rf"(\d{{4}}-\d{{2}}-\d{{2}})_[a-z]+_{re.escape(dataset)}_{what}_v(\d+)\.parquet"
    )
    found = [
        (m.group(1), int(m.group(2)), p)
        for p in directory.glob("*.parquet")
        if (m := pattern.fullmatch(p.name))
    ]
    if not found:
        raise DatasetError(
            f"no {dataset} {what} table in {directory}; run the earlier module first"
        )
    return max(found, key=lambda x: (x[0], x[1]))[2]


def folder(root: Path) -> Path:
    """M5's output folder under a ctx root."""
    return root / "outputs" / f"{dict(MODULES)['de']:02d}_de"


def results(root: Path, slug: str, contrast: str) -> Path:
    """M5's newest results table of one contrast.

    Raises:
        DatasetError: If ``ctx de`` has not written it.
    """
    return newest(folder(root), slug, f"{contrast}_results")


def required(spec: Dataset, module: str) -> tuple[str, ...]:
    """The contrasts whose M5 results a module's command cannot start without.

    function reads every contrast's, ncrna its featured one's, a cross-study proteo its own
    contrast's; proteo's other designs note a missing one and draw the rest.
    """
    if module == "function":
        return tuple(c.name for c in spec.contrasts)
    if module == "ncrna" and spec.ncrna is not None:
        return (spec.ncrna.featured,)
    if module == "proteo" and spec.proteo and spec.proteo.design == "cross_dataset":
        return (spec.proteo.contrast,) if spec.proteo.contrast else ()
    return ()


def absent(root: Path, spec: Dataset, module: str) -> tuple[str, ...]:
    """The contrasts of ``required`` with no M5 results under the root."""
    missing = []
    for contrast in required(spec, module):
        try:
            results(root, spec.slug, contrast)
        except DatasetError:
            missing.append(contrast)
    return tuple(missing)
