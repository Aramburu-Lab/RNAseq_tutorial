"""Package wo13-de: where `ctx de` holds its memory, phase by phase (S.5).

A thread samples this process's RSS, and its children's, every 100 ms; wrappers around M5's
own steps name the phase each sample falls in. Prints, per phase, the peak RSS of the main
process and the largest child, and the peak traced by tracemalloc where asked.

Run: bash env.sh python docs/plans/logs/wo13_de/profile_de.py <root> <dataset> [--trace]
"""

import contextlib
import sys
import threading
import time
import tracemalloc
from collections import defaultdict
from pathlib import Path

import psutil

from ctx import runtime

settings = runtime.start(None, Path(sys.argv[1]))
from pydeseq2.dds import DeseqDataSet  # noqa: E402
from pydeseq2.ds import DeseqStats  # noqa: E402

from ctx.de import deseq  # noqa: E402 - after the runtime profile, as the CLI imports it
from ctx.de import run as module  # noqa: E402
from ctx.report.builder import ParameterCard  # noqa: E402

dataset, trace = sys.argv[2], "--trace" in sys.argv
phase = ["load"]
peaks: dict[str, list[int]] = defaultdict(lambda: [0, 0, 0])
order: list[str] = []
me = psutil.Process()
done = threading.Event()


def sample() -> None:
    """Record the main RSS, the largest child's, and their sum, under the current phase."""
    # Listing children walks every process on macOS and stalled the fits sevenfold when done
    # every sample; the list is refreshed every 2 s and only their RSS read each time.
    kids_of: list[psutil.Process] = []
    tick = 0
    while not done.is_set():
        if tick % 20 == 0:
            kids_of = me.children(recursive=True)
        tick += 1
        kids = []
        for child in kids_of:
            with contextlib.suppress(psutil.Error):
                kids.append(child.memory_info().rss)
        own = me.memory_info().rss
        p = peaks[phase[0]]
        p[0], p[1] = max(p[0], own), max(p[1], max(kids, default=0))
        p[2] = max(p[2], own + sum(kids))
        time.sleep(0.1)


def named(label: str, function):  # noqa: ANN001, ANN202 - a probe, not library code
    """Wrap ``function`` so its samples fall under ``label`` (numbered when repeated)."""

    def wrapper(*args, **kwargs):  # noqa: ANN002, ANN003, ANN202
        before = phase[0]
        name = f"{label}#{sum(o.startswith(label) for o in order) + 1}"
        order.append(name)
        phase[0] = name
        try:
            return function(*args, **kwargs)
        finally:
            if trace:
                # The run's peak so far, never reset: the phase where it rises holds the peak.
                now, peak = tracemalloc.get_traced_memory()
                print(f"{name}: traced now {now:,} B, run peak {peak:,} B", flush=True)
            phase[0] = before

    return wrapper


for method in ("fit_size_factors", "fit_genewise_dispersions", "fit_MAP_dispersions",
               "fit_LFC", "calculate_cooks", "refit", "vst_transform"):  # fmt: skip
    setattr(DeseqDataSet, method, named(method, getattr(DeseqDataSet, method)))
for method in ("summary", "lfc_shrink"):
    setattr(DeseqStats, method, named(method, getattr(DeseqStats, method)))
module.fit = named("fit", deseq.fit)
module.pca = named("pca", deseq.pca)
module._figures = named("figures", module._figures)
module.sex_check = named("sex_check", module.sex_check)
module.write = named("page", module.write)
if trace:
    tracemalloc.start()
threading.Thread(target=sample, daemon=True).start()
module.run(
    settings, (ParameterCard("--probe", "A probe run.", "none", "None."),), "profile", dataset,
    dataset_file=Path(__file__).parents[4] / "conf" / "datasets" / f"{dataset}.toml",
    min_count=10, lfc_shrink=True,
)  # fmt: skip
done.set()
time.sleep(0.05)
print(f"{'phase':<26}{'main RSS':>16}{'largest child':>16}{'main + children':>18}")
for name in ["load", *order]:
    own, kid, total = peaks.get(name, [0, 0, 0])
    print(f"{name:<26}{own:>16,}{kid:>16,}{total:>18,}")
usage = runtime.measure(settings.started)
print(f"getrusage: self {usage.peak_rss_self:,} B, children {usage.peak_rss_children:,} B, "
      f"wall {usage.wall_seconds:.1f} s")  # fmt: skip
