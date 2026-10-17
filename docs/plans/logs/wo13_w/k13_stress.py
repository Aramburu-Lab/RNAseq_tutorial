"""Package W, K.13: the heatmap-cycle test's browser plan, many drivers at once, traced.

Builds the test's page once, then starts ``drivers`` copies of tests/browser.mjs at the same
time with the test's own tasks and ``trace`` on, each given ``deadline`` seconds. A driver past
its deadline is reported with its last trace line and the CPU of every process of its browser,
then killed. Run with the CPU loaded (k13_repeat.sh's ``yes`` processes) to reproduce a loaded
suite.
Run: python k13_stress.py <out folder> <drivers> <deadline s> [rounds]
"""

import json
import shutil
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT / "tests"))
import polars as pl  # noqa: E402
from test_export_browser import CHROMIUM  # noqa: E402
from test_heatmap_cycle import FIXTURE, WATCH  # noqa: E402

from ctx.io.naming import stem  # noqa: E402
from ctx.quant.figures import correlation_heatmap  # noqa: E402
from ctx.report.builder import (  # noqa: E402
    FigureText,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    SelfCheck,
    write,
)
from ctx.viz.save import figure_block, write_figure  # noqa: E402


def build(out: Path) -> Path:
    """The test's page, as its fixture builds it."""
    table = pl.read_csv(FIXTURE)
    rho = table.select(pl.col("^rho_.*$")).to_numpy()
    data = correlation_heatmap(rho, table.select("sample", "cell_type", "batch"), "Made-up.")
    figure_stem = stem("quant", "cycle_correlation", 1, day=date(2026, 10, 14))
    text = FigureText("cycle", reading=("Each cell is one sample pair's rho.",))
    write_figure(data, out, figure_stem, title="cycle", data_releases={}, text=text)
    built = ModulePage(
        title="M4 Demo: cycle", purpose="Purpose.", input="Input.", method="Method.",
        parameters=(ParameterCard("--x", "A card.", "1", "Moves it."),),
        results=(Group("Figures", (figure_block(data, figure_stem, text),)),),
        interpretation="Interpretation.", questions=(SelfCheck("A question.", "An answer."),),
        pitfalls=("A pitfall.",),
        reproduce=Reproduction("ctx demo", 1.0, 1, 1, "laptop", "none", {"python": "3.12"}),
    )  # fmt: skip
    return write(built, out / f"{stem('quant', 'cycle_overview', 1, day=date(2026, 10, 14))}.html")


def plan(page: Path, downloads: Path) -> dict:
    """The test's tasks, traced."""
    tasks = [{"name": f"{w}_{toggle}", "url": page.as_uri(), "width": w, "height": 900,
              "act": WATCH % ("true" if toggle else "false"),
              "resize": [{"width": other, "act": WATCH % "false"},
                         {"width": w, "act": WATCH % "false"}]}
             for w, other, toggle in ((1280, 375, False), (375, 1280, True))]  # fmt: skip
    return {"chromium": CHROMIUM, "downloads": str(downloads), "tasks": tasks, "trace": True}


def tree(pid: int) -> tuple[str, set[int]]:
    """Every process under ``pid`` with its CPU and state, from ps, and their ids."""
    rows = subprocess.run(["ps", "-Ao", "pid,ppid,pcpu,stat,etime,comm"], capture_output=True,
                          text=True, check=True).stdout.splitlines()[1:]  # fmt: skip
    parent = {int(r.split()[0]): int(r.split()[1]) for r in rows}
    keep, out = {pid}, []
    for _ in range(4):
        keep |= {p for p, pp in parent.items() if pp in keep}
    for r in rows:
        if int(r.split()[0]) in keep:
            out.append(r[:160])
    return "\n".join(out), keep


def main() -> None:
    """Rounds of concurrent drivers; a summary line per driver."""
    out, drivers, deadline = Path(sys.argv[1]), int(sys.argv[2]), float(sys.argv[3])
    rounds = int(sys.argv[4]) if len(sys.argv) > 4 else 1  # noqa: PLR2004
    out.mkdir(parents=True, exist_ok=True)
    page = build(out)
    node = shutil.which("node")
    for r in range(rounds):
        runs = []
        for k in range(drivers):
            downloads = out / f"downloads_{r}_{k}"
            downloads.mkdir(exist_ok=True)
            err = (out / f"trace_{r}_{k}.txt").open("w")
            proc = subprocess.Popen([node, str(ROOT / "tests" / "browser.mjs")],
                                    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=err,
                                    text=True)  # fmt: skip
            assert proc.stdin is not None
            proc.stdin.write(json.dumps(plan(page, downloads)))
            proc.stdin.close()
            proc.stdin = None  # closed here, so communicate() writes nothing
            runs.append((k, proc, time.monotonic(), err))
        for k, proc, start, err in runs:
            try:
                left = max(1.0, deadline - (time.monotonic() - start))
                stdout, _ = proc.communicate(timeout=left)
                err.close()
                took = time.monotonic() - start
                results = json.loads(stdout) if proc.returncode == 0 else []
                drawn = [d["drawn"] for res in results for d in (res["acted"], *res["resized"])]
                print(f"round {r} driver {k}: exit {proc.returncode} in {took:.1f} s; "
                      f"draws after settling {drawn}", flush=True)  # fmt: skip
            except subprocess.TimeoutExpired:
                err.flush()
                last = (out / f"trace_{r}_{k}.txt").read_text().strip().splitlines()[-3:]
                listing, pids = tree(proc.pid)
                print(f"round {r} driver {k}: PAST {deadline:.0f} s; last trace {last}\n"
                      f"{listing}", flush=True)  # fmt: skip
                # The browser too: a killed driver leaves its Chromium running otherwise.
                subprocess.run(["kill", *map(str, sorted(pids))], check=False)
                err.close()


if __name__ == "__main__":
    main()
