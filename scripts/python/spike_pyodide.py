#!/usr/bin/env python3
# /// codebase-meta
# exempt: >-
#   One-off Phase 1 feasibility spike (TODO.md P1.5), not a pipeline step. It decides whether
#   the JupyterLite notebooks run pydeseq2 live; its record is a deliverable, read by no step.
# ///
"""Pyodide feasibility spike: can the JupyterLite notebooks run pydeseq2 live?

``decide()`` is the rule, committed before any run so the record cannot shape it.
``agreement()`` turns two fits into the two numbers the rule reads. The rest is stage A: Node
runs Pyodide from npm, a separate CPython process runs the same functions from the pixi
environment, and one long table records both, with the verdicts.

Usage, instructor-side, with the network; every cache stays outside the repository:
    npm install --prefix <scratch>/node pyodide@314.0.6
    pixi run python scripts/python/spike_pyodide.py --node-dir <scratch>/node
    --workdir <scratch>/run [--stage-b <browser timings json>]
"""

import argparse
import ast
import hashlib
import importlib.util
import inspect
import json
import math
import os
import platform
import subprocess
import sys
from collections.abc import Mapping, Sequence
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

from ctx.io.naming import stem
from ctx.io.readme_companion import write_table

Record = Mapping[tuple[str, str, str], str]
PYODIDE_VERSION = "314.0.6"
# plotly is pinned to the pixi environment's version, so the lite notebooks draw what ctx draws.
REQUIREMENTS = {"deseq": "pydeseq2==0.5.4", "fallback": "plotly==7.1.0"}
DESIGN = "~ cell_type + condition"
ALPHA = 0.05
SEED = 42
# (n_genes, donors_per_condition). Each donor gives both cell types: 6 donors per condition
# make 24 samples and 50 make 200. Taken from no cohort.
GRID = ((1000, 6), (20000, 6), (20000, 50))
TIMEOUT_S = 600
ROOT = Path(__file__).resolve().parents[2]
REPORTED = ["numpy", "scipy", "scikit-learn", "pandas", "pydeseq2", "anndata", "formulaic"]
REPORTED += ["typing-extensions", "matplotlib", "plotly"]


def agreement(
    reference_padj: Sequence[float],
    reference_lfc: Sequence[float],
    padj: Sequence[float],
    lfc: Sequence[float],
) -> tuple[int, float]:
    """Compare a fit with the CPython reference fit, gene by gene in the same order.

    A gene disagrees when its padj is NaN in one fit only, or when neither padj lies in the
    closed band [0.04, 0.06] and the two fall on different sides of 0.05. A gene with either
    padj inside the band is exempt, which is what the band is for.

    Args:
        reference_padj: BH-adjusted p-values of the CPython fit.
        reference_lfc: log2 fold changes of the CPython fit.
        padj: BH-adjusted p-values of the fit under test, same gene order.
        lfc: log2 fold changes of the fit under test, same gene order.

    Returns:
        The number of disagreeing genes, and max |dLFC| over genes with an LFC in both fits.
        An LFC present in one fit only makes max |dLFC| infinite.

    Raises:
        ValueError: If the four inputs are not the same length.
    """
    ref_p, new_p, ref_l, new_l = (
        np.asarray(x, dtype=float) for x in (reference_padj, padj, reference_lfc, lfc)
    )
    if not ref_p.shape == new_p.shape == ref_l.shape == new_l.shape:
        raise ValueError("the fits must cover the same genes in the same order")
    # Inside the function, so the source hash recorded with the spike covers the thresholds.
    low, alpha, high = 0.04, 0.05, 0.06
    one_nan = np.isnan(ref_p) != np.isnan(new_p)
    in_band = ((ref_p >= low) & (ref_p <= high)) | ((new_p >= low) & (new_p <= high))
    # NaN < alpha is False on both sides, so a gene NaN in both fits never lands here.
    other_side = ~one_nan & ~in_band & ((ref_p < alpha) != (new_p < alpha))
    if (np.isnan(ref_l) != np.isnan(new_l)).any():
        return int(one_nan.sum() + other_side.sum()), math.inf
    diff = np.abs(ref_l - new_l)
    diff = diff[~np.isnan(diff)]
    return int(one_nan.sum() + other_side.sum()), float(diff.max()) if diff.size else 0.0


def decide(record: Record) -> dict[str, tuple[str, str]]:
    """Apply the P1.5 rule to a spike record; the rule was fixed before the spike ran.

    Live pydeseq2 in JupyterLite at a grid size only if all of these hold: micropip installs
    pydeseq2; ``~ cell_type + condition`` fits; the fit agrees with a CPython fit of the same
    pydeseq2 version (no disagreeing gene as ``agreement()`` counts them, max |dLFC| <= 1e-3);
    wheels plus data are <= 50 MB, taken as 50,000,000 bytes; install plus fit take <= 120 s
    in the browser. Any failure gives the D4 fallback: precomputed grids plus live scipy tests.
    Stage A holding with no stage-B timing gives a provisional verdict. A missing or NaN
    metric fails its criterion, because every comparison is written so NaN reads as False.

    Args:
        record: (stage, component, field) to value, the spike's long table as a mapping.

    Returns:
        For each size in the record's grid, (verdict, why), with verdict ``live``,
        ``fallback`` or ``provisional``. Key ``plotting`` holds (``plotly`` or
        ``matplotlib``, why): the lite notebooks use matplotlib if plotly does not install.
    """

    def num(key: tuple[str, str, str]) -> float:
        value = record.get(key)
        return math.nan if value is None else float(value)

    # Inside the function, so the source hash recorded with the spike covers the thresholds.
    max_dlfc, max_bytes, max_seconds = 1e-3, 50_000_000, 120.0
    a = "stage_a"
    plotly = record.get((a, "plotly", "status")) == "installed"
    verdicts = {
        "plotting": ("plotly", "plotly installed")
        if plotly
        else ("matplotlib", "plotly did not install")
    }
    for size in record[("meta", "grid", "sizes")].split(";"):
        fit = (a, f"fit_{size}")
        wheels_data = num((a, "pydeseq2_wheels", "total_bytes")) + num(
            ("meta", f"input_{size}", "total_bytes")
        )
        browser = ("stage_b", f"browser_{size}", "install_plus_fit_s")
        criteria = (
            (
                record.get((a, "pydeseq2", "status")) == "installed",
                "micropip did not install pydeseq2",
            ),
            (
                record.get((a, "pydeseq2", "version"))
                == record.get(("cpython", "pydeseq2", "version")),
                "the CPython reference is not the same pydeseq2 version",
            ),
            (
                record.get((*fit, "status")) == "ok",
                f"~ cell_type + condition fit: {record.get((*fit, 'status'), 'missing')}",
            ),
            (
                num((*fit, "padj_side_disagreements_n")) == 0,
                f"{num((*fit, 'padj_side_disagreements_n')):g} genes on the other side of 0.05",
            ),
            (
                num((*fit, "max_abs_dlfc_log2")) <= max_dlfc,
                f"max |dLFC| {num((*fit, 'max_abs_dlfc_log2')):g} > {max_dlfc:g}",
            ),
            (wheels_data <= max_bytes, f"wheels plus data {wheels_data:.0f} > {max_bytes} bytes"),
            (
                browser not in record or num(browser) <= max_seconds,
                f"install plus fit {num(browser):g} s in the browser > {max_seconds:g} s",
            ),
        )
        failed = [why for held, why in criteria if not held]
        if failed:
            verdicts[size] = ("fallback", "; ".join(failed))
        elif browser not in record:
            verdicts[size] = ("provisional", "stage A holds; the stage-B browser timing is missing")
        else:
            verdicts[size] = ("live", "every criterion holds")
    return verdicts


def size_label(n_genes: int, donors_per_condition: int) -> str:
    """Name a grid size genes x samples, e.g. ``1000x24``; every donor gives two samples."""
    return f"{n_genes}x{4 * donors_per_condition}"


def counts_call(n_genes: int, donors_per_condition: int) -> str:
    """The ``counts()`` call behind one grid size, with every argument explicit."""
    return f"counts(seed={SEED}, n_genes={n_genes}, donors_per_condition={donors_per_condition})"


def input_csv(call: str) -> tuple[bytes, bytes]:
    """Run a recorded ``counts()`` call and serialise what it returns, as both arms read it.

    Args:
        call: A string from ``counts_call()``.

    Returns:
        The counts CSV (``gene_id``, then one column per sample) and the sample-sheet CSV.

    Raises:
        ValueError: If the call is not ``counts(...)`` with every argument given by keyword.
    """
    # tests/synth.py is test code, outside the wheel, so it is loaded by path (P1.3).
    spec = importlib.util.spec_from_file_location("spike_synth", ROOT / "tests" / "synth.py")
    if spec is None or spec.loader is None:
        raise ValueError("tests/synth.py cannot be loaded")
    synth = importlib.util.module_from_spec(spec)
    sys.modules["spike_synth"] = synth
    spec.loader.exec_module(synth)
    node = ast.parse(call, mode="eval").body
    if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and not node.args):
        raise ValueError(f"not a keyword-only call: {call!r}")
    kwargs = {str(k.arg): ast.literal_eval(k.value) for k in node.keywords}
    declared = set(inspect.signature(synth.counts).parameters)
    if node.func.id != "counts" or set(kwargs) != declared:
        raise ValueError(f"{call!r} must set every argument of counts(): {sorted(declared)}")
    data = synth.counts(**kwargs)
    return data.counts.write_csv().encode(), data.samples.write_csv().encode()


def source_sha256(function: Any) -> str:
    """SHA-256 of a function's source text, as ``inspect.getsource`` returns it."""
    return hashlib.sha256(inspect.getsource(function).encode()).hexdigest()


def deseq_fit(  # noqa: PLR0913 - one call signature for both runtimes
    counts_csv: str,
    samples_csv: str,
    results_csv: str,
    design: str,
    contrast: list[str],
    backend: str,
) -> None:
    """Fit a pydeseq2 design and write one contrast's results.

    Pyodide executes this function from its source text and CPython imports it, so both arms
    run the same code; the imports sit inside for that reason.

    Args:
        counts_csv: Genes by samples, first column ``gene_id``.
        samples_csv: The sample sheet, first column ``sample``.
        results_csv: Where baseMean, log2FoldChange, pvalue and padj per gene go.
        design: The formula, e.g. ``~ cell_type + condition``.
        contrast: [factor, tested level, reference level].
        backend: The joblib backend pydeseq2 runs under. ``loky`` is pydeseq2's default.
            ``in_process`` is registered here: a sequential backend that accepts the
            ``inner_max_num_threads`` pydeseq2 always sets, which Pyodide needs because it
            has no loky and joblib's threading fallback rejects that argument.
    """
    import pandas as pd  # noqa: PLC0415 - runs from source inside Pyodide
    from joblib import register_parallel_backend  # noqa: PLC0415
    from joblib.parallel import SequentialBackend  # noqa: PLC0415
    from pydeseq2.dds import DeseqDataSet  # noqa: PLC0415
    from pydeseq2.default_inference import DefaultInference  # noqa: PLC0415
    from pydeseq2.ds import DeseqStats  # noqa: PLC0415

    class InProcess(SequentialBackend):
        supports_inner_max_num_threads = True

    register_parallel_backend("in_process", InProcess)
    counts = pd.read_csv(counts_csv, index_col=0).T
    samples = pd.read_csv(samples_csv, index_col=0).loc[counts.index]
    # One CPU in both arms: Pyodide has no worker processes, and the arms must not differ.
    inference = DefaultInference(n_cpus=1, backend=backend)
    dds = DeseqDataSet(
        counts=counts, metadata=samples, design=design, inference=inference, quiet=True
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=contrast, inference=inference, quiet=True)
    stats.summary()
    stats.results_df[["baseMean", "log2FoldChange", "pvalue", "padj"]].to_csv(results_csv)


def scipy_tests(counts_csv: str, samples_csv: str, results_csv: str, cell_type: str) -> None:
    """The fallback's live test: Mann-Whitney U on CPM, case vs control in one cell type, BH.

    Runs in both arms from the same source, like ``deseq_fit``. A gene whose values are all
    tied has no p-value and stays out of the BH set.

    Args:
        counts_csv: Genes by samples, first column ``gene_id``.
        samples_csv: The sample sheet, first column ``sample``.
        results_csv: Where pvalue and padj per gene go.
        cell_type: The cell type whose case and control samples are compared.
    """
    import pandas as pd  # noqa: PLC0415 - runs from source inside Pyodide
    from scipy.stats import false_discovery_control, mannwhitneyu  # noqa: PLC0415

    counts = pd.read_csv(counts_csv, index_col=0)
    samples = pd.read_csv(samples_csv, index_col=0).loc[counts.columns]
    cpm = counts / counts.sum(axis=0) * 1e6
    keep = samples["cell_type"] == cell_type
    case = cpm.loc[:, keep & (samples["condition"] == "case")].to_numpy()
    control = cpm.loc[:, keep & (samples["condition"] == "control")].to_numpy()
    out = pd.DataFrame(
        {"pvalue": mannwhitneyu(case, control, axis=1, alternative="two-sided").pvalue},
        index=counts.index,
    )
    tested = out["pvalue"].notna()
    out.loc[tested, "padj"] = false_discovery_control(out.loc[tested, "pvalue"], method="bh")
    out.to_csv(results_csv)


# arm: (function, extra arguments, micropip requirement, Pyodide packages loaded, must not import)
# deseq is the rule's arm: pydeseq2 exactly as the CLI calls it. deseq_shim is recorded beside
# it and decide() never reads it; it measures what the in-process backend would buy.
CONTRAST = ["condition", "case", "control"]
ARMS: dict[str, tuple[Any, tuple[Any, ...], str, list[str], str]] = {
    "deseq": (deseq_fit, (DESIGN, CONTRAST, "loky"), REQUIREMENTS["deseq"], [], ""),
    "deseq_shim": (deseq_fit, (DESIGN, CONTRAST, "in_process"), REQUIREMENTS["deseq"], [], ""),
    "fallback": (
        scipy_tests,
        ("type_a",),
        REQUIREMENTS["fallback"],
        ["scipy", "pandas"],
        "pydeseq2",
    ),
}


def _results(path: Path) -> pl.DataFrame:
    """A results CSV written by pandas, gene id first and empty cells as NaN."""
    frame = pl.read_csv(path, infer_schema_length=None)
    frame = frame.rename({frame.columns[0]: "gene_id"})
    return frame.with_columns(pl.exclude("gene_id").cast(pl.Float64).fill_null(math.nan))


def compare(arm: str, reference_csv: Path, results_csv: Path) -> dict[str, Any]:
    """Compare one Pyodide run with its CPython reference, gene by gene.

    Args:
        arm: ``deseq``, compared with ``agreement()``, or ``fallback``.
        reference_csv: The CPython results.
        results_csv: The Pyodide results.

    Returns:
        Named numbers; each name's suffix gives its unit.

    Raises:
        ValueError: If the two runs did not test the same genes in the same order.
    """
    ref, new = _results(reference_csv), _results(results_csv)
    if not ref["gene_id"].equals(new["gene_id"]):
        raise ValueError(f"{reference_csv} and {results_csv} cover different genes")
    out: dict[str, Any] = {"genes_n": ref.height}
    for name, frame in (("cpython", ref), ("pyodide", new)):
        out[f"padj_below_alpha_{name}_n"] = int((frame["padj"] < ALPHA).sum())
    # A padj NaN in one fit only is how a shifted independent-filtering threshold shows.
    for column in ("pvalue", "padj"):
        a, b = ref[column].to_numpy(), new[column].to_numpy()
        out[f"{column}_nan_mismatch_n"] = int((np.isnan(a) != np.isnan(b)).sum())
        out[f"max_abs_d{column}"] = float(np.nanmax(np.abs(a - b), initial=0.0))
    if arm.startswith("deseq"):
        n, dlfc = agreement(
            ref["padj"].to_numpy(),
            ref["log2FoldChange"].to_numpy(),
            new["padj"].to_numpy(),
            new["log2FoldChange"].to_numpy(),
        )
        return out | {"padj_side_disagreements_n": n, "max_abs_dlfc_log2": dlfc}
    out["padj_below_alpha_mismatch_n"] = int(((ref["padj"] < ALPHA) != (new["padj"] < ALPHA)).sum())
    return out


Row = tuple[str, str, str, str]


def _text(value: Any) -> str:
    """A value as the long table holds it; floats keep every digit."""
    return repr(value) if isinstance(value, float) else str(value)


def _run_rows(arm: str, size: str, ref: Mapping[str, Any], new: Mapping[str, Any]) -> list[Row]:
    """Rows for one arm at one size: both runtimes' status and cost, and their comparison."""
    from pyodide_harness import wheels  # noqa: PLC0415 - sibling script, see spike()

    component = {"deseq": "fit", "deseq_shim": "fit_shim", "fallback": "scipy_tests"}[arm]
    component = f"{component}_{size}"
    process = new["process"]
    status = new.get("run", {}).get("status") or (
        "timeout" if process["status"] == "timeout" else "not_run: install did not complete"
    )
    rows = [("cpython", component, "status", ref.get("status") or ref["process"]["status"])]
    rows += [("cpython", component, "error", ref["error"])] if "error" in ref else []
    rows += [("cpython", component, f"process_{k}", _text(v)) for k, v in ref["process"].items()]
    rows += [("cpython", component, "python_s", _text(ref.get("python_s")))]
    rows += [("stage_a", component, "status", status)]
    rows += [("stage_a", component, "error", new["run"]["error"])] if status == "error" else []
    phases = new.get("phases_s", {})
    rows += [("stage_a", component, f"node_{k}_s", _text(v)) for k, v in phases.items()]
    if {"install", "run"} <= phases.keys():
        total = phases["install"] + phases["run"]
        rows.append(("stage_a", component, "node_install_plus_run_s", _text(total)))
    rows += [("stage_a", component, f"process_{k}", _text(v)) for k, v in process.items()]
    downloaded = sum(wheels(new.get("fetches", [])).values())
    rows.append(("stage_a", component, "wheels_downloaded_bytes", str(downloaded)))
    if ref.get("status") == status == "ok":
        measured = compare(arm, ref["results"], new["results"])
        rows += [("stage_a", component, k, _text(v)) for k, v in measured.items()]
    return rows


def _install_rows(arm: str, run: Mapping[str, Any]) -> list[Row]:
    """Package status, versions and every downloaded file, from one arm's first run."""
    from pyodide_harness import wheels  # noqa: PLC0415 - sibling script, see spike()

    group = "pydeseq2_wheels" if arm == "deseq" else "fallback_wheels"
    files = wheels(run.get("fetches", []))
    rows = [("stage_a", f"{group}:{name}", "bytes", str(b)) for name, b in sorted(files.items())]
    rows += [("stage_a", group, "total_bytes", str(sum(files.values())))]
    rows += [("stage_a", group, "files_n", str(len(files)))]
    install = run.get("install", {})
    for name, step in install.items():
        if isinstance(step, dict) and "status" in step:
            package = "pydeseq2_absent" if name == "absent" else name
            rows += [("stage_a", package, k, str(v)) for k, v in step.items()]
    for name, (version, source) in sorted(install.get("distributions", {}).items()):
        rows.append(("stage_a", f"{arm}_dist:{name}", "version", f"{version} ({source})"))
    return rows


def _meta_rows(node_dir: Path) -> list[Row]:
    """Tools, rule hashes and the design: what a reader needs to rerun the spike."""
    lock = json.loads((node_dir / "node_modules" / "pyodide" / "pyodide-lock.json").read_text())
    npm = json.loads((node_dir / "package-lock.json").read_text())["packages"]
    node = subprocess.run(["node", "--version"], capture_output=True, text=True, check=True)
    return [
        ("meta", "pyodide", "version", PYODIDE_VERSION),
        ("meta", "pyodide", "npm_version", npm["node_modules/pyodide"]["version"]),
        ("meta", "pyodide", "npm_integrity", npm["node_modules/pyodide"]["integrity"]),
        ("meta", "pyodide", "python_version", lock["info"]["python"]),
        ("meta", "node", "version", node.stdout.strip()),
        ("meta", "host", "platform", platform.platform()),
        ("meta", "host", "cpu_count_n", str(os.cpu_count())),
        ("meta", "rule", "decide_source_sha256", source_sha256(decide)),
        ("meta", "rule", "agreement_source_sha256", source_sha256(agreement)),
        ("meta", "design", "formula", DESIGN),
        ("meta", "design", "contrast", "condition: case vs control"),
        ("meta", "design", "fallback_test", "Mann-Whitney U on CPM, type_a, case vs control; BH"),
        ("meta", "timeout", "limit_s", str(TIMEOUT_S)),
        ("meta", "grid", "sizes", ";".join(size_label(g, d) for g, d in GRID)),
    ]


def spike(work: Path, node_dir: Path, stage_b: Mapping[str, Mapping[str, Any]]) -> list[Row]:
    """Write the inputs, run both arms in both runtimes at every size, and collect the rows.

    Args:
        work: An empty scratch directory.
        node_dir: The npm prefix holding ``node_modules/pyodide``.
        stage_b: Browser timings as {component: {field: value}}, or empty when not run.

    Returns:
        The record's rows, the decision last.

    Raises:
        ValueError: If two rows share a (stage, component, field) key.
    """
    # A sibling script: importable when this file runs as a script, which is the only time
    # spike() runs. The tests load this file by path and never reach this line.
    import pyodide_harness as harness  # noqa: PLC0415

    rows = _meta_rows(node_dir)
    (work / "inputs").mkdir()
    for n_genes, donors in GRID:
        size, call = size_label(n_genes, donors), counts_call(n_genes, donors)
        blobs = dict(zip(("counts", "samples"), input_csv(call), strict=True))
        rows.append(("meta", f"input_{size}", "counts_call", call))
        for kind, blob in blobs.items():
            (work / "inputs" / f"{size}_{kind}.csv").write_bytes(blob)
            rows.append(
                ("meta", f"input_{size}", f"{kind}_sha256", hashlib.sha256(blob).hexdigest())
            )
            rows.append(("meta", f"input_{size}", f"{kind}_bytes", str(len(blob))))
        rows.append(("meta", f"input_{size}", "total_bytes", str(sum(map(len, blobs.values())))))
    versions: dict[str, str] = {}
    for arm, (function, extra, requirement, loaded, absent) in ARMS.items():
        for i, (n_genes, donors) in enumerate(GRID):
            size = size_label(n_genes, donors)
            host = [str(work / "inputs" / f"{size}_{k}.csv") for k in ("counts", "samples")]
            ref_dir = work / f"cpython_{arm}_{size}"
            args = [*host, str(ref_dir / "results.csv"), *extra]
            ref = harness.in_cpython(Path(__file__), function, args, ref_dir, TIMEOUT_S)
            ref["results"] = ref_dir / "results.csv"
            versions |= ref.get("versions", {})
            job = {
                "function": function,
                "args": ["/counts.csv", "/samples.csv", "/results.csv", *extra],
                "inputs": {"/counts.csv": host[0], "/samples.csv": host[1]},
                "output": "/results.csv",
                "requirements": [requirement],
                "loaded": loaded,
                "absent": absent,
            }
            new = harness.in_pyodide(job, work / f"pyodide_{arm}_{size}", node_dir, TIMEOUT_S)
            rows += _run_rows(arm, size, ref, new)
            # The shim arm installs the same wheels as deseq; its install rows would repeat keys.
            rows += _install_rows(arm, new) if i == 0 and arm != "deseq_shim" else []
    rows += [("cpython", name, "version", v) for name, v in sorted(versions.items())]
    rows += [
        ("stage_b", c, f, _text(v)) for c, fields in stage_b.items() for f, v in fields.items()
    ]
    verdicts = decide({(s, c, f): v for s, c, f, v in rows})
    if not stage_b:
        # Stage B runs only where stage A holds; decide() reads no row of this kind.
        pending = [k for k, (verdict, _) in verdicts.items() if verdict == "provisional"]
        why = f"sizes awaiting it: {', '.join(pending)}" if pending else "no size passed stage A"
        rows.append(("stage_b", "run", "status", f"not run: {why}"))
    for key, (verdict, why) in verdicts.items():
        rows += [("decision", key, "verdict", verdict), ("decision", key, "why", why)]
    if len({row[:3] for row in rows}) != len(rows):
        raise ValueError("two rows share a (stage, component, field) key")
    return rows


COLUMNS = {
    "stage": (
        "none",
        "Source of the value: meta (inputs, tools, rule hashes), cpython (the reference in the "
        "pixi environment), stage_a (Pyodide 314.0.6 in Node), stage_b (Pyodide in a browser), "
        "decision (decide() applied to every row above it).",
    ),
    "component": (
        "none",
        "What was measured: a package, a grid size (fit_<genes>x<samples> for pydeseq2 as the "
        "CLI calls it, fit_shim_<size> for pydeseq2 under the in-process joblib backend, which "
        "decide() does not read, scipy_tests_<size> for the fallback), an input "
        "(input_<size>) or a downloaded file (<group>:<file name>).",
    ),
    "field": (
        "none",
        "The quantity. Its suffix is the unit of value: _s seconds, _bytes bytes, _n count, "
        "_log2 log2 fold change; any other field is text or a dimensionless number.",
    ),
    "value": ("given by the field suffix", "The recorded value, as text."),
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Read the command line.

    Args:
        argv: Arguments without the program name; ``sys.argv[1:]`` when None.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--node-dir", type=Path, help="npm prefix with pyodide installed.")
    parser.add_argument("--workdir", type=Path, help="Empty scratch directory for the runs.")
    parser.add_argument("--outdir", type=Path, default=ROOT / "deliverables")
    parser.add_argument("--stage-b", type=Path, help="JSON of browser timings, if measured.")
    parser.add_argument("--day", type=date.fromisoformat, default=date.today())
    parser.add_argument("--reference", type=Path, help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    """Run the spike and write the record with its README, or serve as the CPython child.

    Args:
        argv: Arguments without the program name; ``sys.argv[1:]`` when None.

    Returns:
        0 once the record, or the child's ``out.json``, is written.
    """
    args = parse_args(argv)
    if args.reference:
        from pyodide_harness import reference  # noqa: PLC0415 - sibling script, see spike()

        reference(args.reference, [f for f, *_ in ARMS.values()], REPORTED)
        return 0
    if args.node_dir is None or args.workdir is None:
        raise SystemExit("--node-dir and --workdir are required")
    args.workdir.mkdir(parents=True, exist_ok=True)
    if any(args.workdir.iterdir()):
        raise SystemExit(f"{args.workdir} is not empty; every run starts from fresh caches")
    stage_b = json.loads(args.stage_b.read_text()) if args.stage_b else {}
    frame = pl.DataFrame(
        spike(args.workdir, args.node_dir, stage_b), schema=list(COLUMNS), orient="row"
    )
    written = write_table(
        frame,
        args.outdir,
        stem("lite_build", "pyodide_spike", 1, day=args.day),
        columns=COLUMNS,
        data_releases={
            "Pyodide": f"{PYODIDE_VERSION}, cdn.jsdelivr.net/pyodide/v{PYODIDE_VERSION}/full/",
            "PyPI": f"the wheels listed in the table, downloaded {args.day.isoformat()}",
        },
        formats=("csv",),
    )
    print(*written, sep="\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
