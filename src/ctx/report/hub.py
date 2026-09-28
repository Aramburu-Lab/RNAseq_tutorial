"""``ctx report``: the hub page linking the newest page of each module for one dataset.

Links are relative, so the hub works after ``outputs/`` is copied elsewhere or opened from
a zip; nothing points to the web.
"""

import hashlib
import math
import os
import re
from datetime import date
from pathlib import Path
from typing import Any

import polars as pl

from ctx import runtime
from ctx.dataset import load
from ctx.de.sanity import VERDICT_COLUMNS, decide, parse
from ctx.io.naming import MODULES, stem
from ctx.io.readme_companion import write_table
from ctx.report.builder import HubLink, HubPage, HubRule, KeyNumber, ParameterCard, write

MODULE = "report"
TITLES = {
    "setup": "M0 Setup",
    "gather": "M1 Gather",
    "qc": "M2 Read QC",
    "process": "M3 Trim, align, count",
    "quant": "M4 Quantification",
    "de": "M5 Differential expression",
    "function": "M6 Function",
    "proteo": "M7 RNA versus protein",
    "ncrna": "M8 Non-coding RNA",
    "lite_build": "JupyterLite build",
}
# The module is one of the MODULES tokens, longest first, so "lite_build" is not read as "lite".
_TOKENS = "|".join(sorted((t for t, _ in MODULES), key=len, reverse=True))
PAGE = re.compile(
    rf"(\d{{4}}-\d{{2}}-\d{{2}})_(?P<module>{_TOKENS})_(?P<what>.+)_v(?P<version>\d+)\.html"
)


def _newest(directory: Path, what: str) -> Path | None:
    """The newest page named ``..._<what>_v<N>.html`` in a module directory."""
    found = [
        (m.group(1), int(m.group("version")), p)
        for p in directory.glob("*.html")
        if (m := PAGE.fullmatch(p.name)) and m.group("what") == what
    ]
    return max(found, key=lambda x: (x[0], x[1]))[2] if found else None


def _record(outputs: Path, dataset: str) -> tuple[dict[str, float], dict[str, str]]:
    """Every metric the newest M5 and M6 metrics tables hold, and their SHA-256."""
    record: dict[str, float] = {}
    read: dict[str, str] = {}
    for module in ("de", "function"):
        folder = outputs / f"{dict(MODULES)[module]:02d}_{module}"
        pattern = re.compile(
            rf"(\d{{4}}-\d{{2}}-\d{{2}})_{module}_{re.escape(dataset)}_metrics_v(\d+)\.parquet"
        )
        found = [
            (m.group(1), int(m.group(2)), p)
            for p in folder.glob("*.parquet")
            if (m := pattern.fullmatch(p.name))
        ]
        if found:
            path = max(found, key=lambda x: (x[0], x[1]))[2]
            record |= dict(pl.read_parquet(path).select("metric", "value").rows())
            read[f"sha256:{path.name}"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return record, read


def sanity(
    settings: runtime.Settings, dataset: str, dataset_file: Path
) -> tuple[list[dict[str, Any]], list[Path]]:
    """Verdicts of the dataset file's pre-registered rules on the recorded metrics."""
    spec = load(dataset_file)
    if not spec.sanity:
        return [], []
    record, read = _record(settings.root / "outputs", dataset)
    rows = decide(record, parse(spec.sanity))
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    releases = {
        "dataset file": str(spec.path),
        f"sha256:{spec.path.name}": hashlib.sha256(spec.path.read_bytes()).hexdigest(),
        **read,
    }
    frame = pl.DataFrame(
        rows,
        schema={
            "rule": pl.String,
            "text": pl.String,
            "condition": pl.String,
            "left": pl.Float64,
            "right": pl.Float64,
            "verdict": pl.String,
            "source": pl.String,
            "rule_verdict": pl.String,
        },
    )
    table_stem = stem(MODULE, f"{dataset}_sanity", 1, day=date.today())
    files = write_table(
        frame.select(list(VERDICT_COLUMNS)),
        out,
        table_stem,
        columns=VERDICT_COLUMNS,
        data_releases=releases,
    )
    return rows, files


def _measured(value: Any) -> str:
    """A recorded value as the hub prints it; a missing one says so."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "not measured"
    # A count prints whole: .4g would show 10,523 genes as 1.052e+04.
    if float(value).is_integer():
        return f"{int(value):,}"
    return f"{value:.4g}"


def _rules(rows: list[dict[str, Any]]) -> tuple[HubRule, ...]:
    """One entry per rule, in file order, with each condition's measured value and verdict."""
    rules: dict[str, HubRule] = {}
    for r in rows:
        condition = (r["condition"], _measured(r["left"]), r["verdict"])
        if r["rule"] in rules:
            known = rules[r["rule"]]
            rules[r["rule"]] = HubRule(known.rule, known.text, known.verdict,
                                       (*known.conditions, condition))  # fmt: skip
        else:
            rules[r["rule"]] = HubRule(r["rule"], r["text"], r["rule_verdict"], (condition,))
    return tuple(rules.values())


def run(
    settings: runtime.Settings,
    cards: tuple[ParameterCard, ...],
    dataset: str,
    dataset_file: Path | None = None,
) -> list[Path]:
    """Write the hub for one dataset under ``outputs/09_report``, with the sanity verdicts.

    Args:
        settings: The resolved runtime.
        cards: Parameter cards for every option of this command.
        dataset: Slug naming the dataset.
        dataset_file: The dataset file; its ``[sanity]`` rules are decided when given.

    Returns:
        The hub page, then the verdict table's files.
    """
    rows, files = sanity(settings, dataset, dataset_file) if dataset_file else ([], [])
    outputs = settings.root / "outputs"
    out = outputs / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    out.mkdir(parents=True, exist_ok=True)
    links, missing = [], []
    for module, number in MODULES:
        if module == MODULE:
            continue
        what = "environment" if module == "setup" else f"{dataset}_overview"
        page = _newest(outputs / f"{number:02d}_{module}", what)
        if page is None:
            missing.append(TITLES.get(module, module))
            continue
        href = Path(os.path.relpath(page, out)).as_posix()
        links.append(HubLink(TITLES.get(module, module), page.name, href, page.name[:10]))
    rules = _rules(rows)
    gated = [r for r in rules if r.verdict != "RECORDED"]
    keys: tuple[KeyNumber, ...] = (
        KeyNumber("Modules run", f"{len(links)} of {len(links) + len(missing)}"),
    )
    if gated:
        passed = sum(r.verdict == "PASS" for r in gated)
        keys += (
            KeyNumber("Rules passed", f"{passed} of {len(gated)}", "recorded rules not counted"),
        )
    hub = HubPage(
        title=f"Course hub: {dataset}",
        intro=f"The newest page of every module run on {dataset}, in course order.",
        parameters=cards,
        links=tuple(links),
        missing=tuple(missing),
        rules=rules,
        key_numbers=keys,
    )
    return [write(hub, out / f"{stem(MODULE, f'{dataset}_hub', 1, day=date.today())}.html"), *files]
