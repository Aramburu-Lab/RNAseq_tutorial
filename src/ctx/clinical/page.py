"""The M10 page.

Patients and endpoints, the ER check, expression by clinical variable, a gene score's survival
(live), the scan (live), power, the per-gene screen, confounding and treatment. Words come
from ``content/clinical.toml``; numbers from ``analyses`` and ``survival``.
"""

from collections.abc import Mapping
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl
from scipy.stats import norm

from ctx import runtime
from ctx.clinical import (
    ClinicalError,
    analyses,
    cutpoint,
    figure_tables,
    figures,
    harmonise,
    panel,
    score,
    survival,
)
from ctx.dataset import ENDPOINTS, Dataset
from ctx.de.run import Writer
from ctx.gather import gdc
from ctx.io.naming import stem
from ctx.report import genes, text
from ctx.report.builder import (
    Group,
    KeyNumber,
    ModulePage,
    ParameterCard,
    Reproduction,
    Result,
    ResultTable,
    write,
)
from ctx.report.entry import write_entry
from ctx.viz.style import KM

PAGE = "clinical"
# Patient-table columns a model can stratify on.
CATEGORICAL = ("er", "pr", "her2", "histology", "menopause", "sex", "plate", "site")
ER = ("positive", "negative")
SCREEN_ROWS = 20
# Confounding-map columns every cohort has, after its own categorical variables.
DESIGN = ("age_group", "plate", "site")
# GDC's sex_at_birth values the sex sentence counts.
SEXES = ("female", "male")
# Liu et al. 2018 Table 3's use marks as a dataset file's [clinical.published.table3_use] copies
# them: a tick, a tick with an asterisk, a cross.
USES = ("recommended", "caution", "not_recommended")


@dataclass(frozen=True)
class Context:
    """What the page is built from: the run's options and the tables it read."""

    spec: Dataset
    options: Any
    endpoint: str
    clinical: harmonise.Clinical
    values: np.ndarray
    gene_table: pl.DataFrame
    carried: panel.Panel
    chosen: list[int]
    source: str
    filtered: tuple[int, int]
    settings: runtime.Settings
    command: str
    cards: tuple[ParameterCard, ...]


@dataclass(frozen=True)
class Scored:
    """The build list's score on the endpoint's patients, and everything the page draws of it."""

    out: analyses.Outcome
    score: np.ndarray
    labels: np.ndarray
    names: tuple[str, ...]
    km: pl.DataFrame
    groups: pl.DataFrame
    logrank: tuple[float, int, float]
    rmst: float
    models: tuple[analyses.Model, analyses.Model]
    scan: cutpoint.Scan


def fill(key: str, **values: Any) -> str:
    """A template of this page's content file, filled."""
    return text.fill(PAGE, key, **values)


def word(key: str) -> str:
    """A fixed text of this page's content file."""
    return text.word(PAGE, key)


def number(value: float | None, spec: str = ".2g") -> str:
    """A number as the page prints it; a missing one as the page's word for none."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return word("none")
    return format(value, spec)


def listed(names: list[str]) -> str:
    """Names as a sentence lists them: "age and stage"."""
    if len(names) < 2:  # noqa: PLR2004 - one name, or none
        return "".join(names)
    return fill("results.and", rest=", ".join(names[:-1]), last=names[-1])


def split_names(rule: str) -> tuple[str, ...]:
    """The groups of a split, lowest first, in page words."""
    keys = (
        ("low", "high") if rule == "median" else ("tertile_low", "tertile_middle", "tertile_high")
    )
    return tuple(word(f"split.{k}") for k in keys)


def slots(rule: str) -> tuple[str, str, str]:
    """The live figure's three groups: the split's own, the median's middle slot left empty."""
    names = split_names(rule)
    middle = word("split.tertile_middle") if rule == "median" else names[1]
    return names[0], middle, names[-1]


def strata(rows: pl.DataFrame, names: tuple[str, ...]) -> np.ndarray | None:
    """One stratum label per row from the strata columns, or None without any."""
    if not names:
        return None
    return np.array(["|".join(str(v) for v in r) for r in rows.select(names).rows()])


def scored(ctx: Context) -> Scored:
    """The build list's score and its analyses, as the page computes them by default."""
    p, options = ctx.clinical.patients, ctx.options
    out = analyses.outcome(p, ctx.endpoint)
    analysed = p[ctx.endpoint].is_not_null().to_numpy()
    z = score.zscores(ctx.values[ctx.chosen][:, analysed])
    value = score.score(z, np.ones(len(ctx.chosen)))
    names = split_names(options.split)
    labels = np.array(names)[score.split(value, options.split)]
    months, event = out.months, out.event
    summary = analyses.groups(months, event, labels, names, options.horizon)
    rmst = summary["rmst"]
    extra, complete = analyses.covariates(out.patients, options.covariates)
    layers = strata(out.patients, options.strata)
    alone = analyses.cox_model(months, event, value[:, None], layers)
    both = np.concatenate([value[:, None], extra], axis=1)[complete]
    adjusted = analyses.cox_model(months[complete], event[complete], both,
                                  None if layers is None else layers[complete])  # fmt: skip
    return Scored(out, value, labels, names, analyses.km_frame(months, event, labels, names),
                  summary, survival.logrank(months, event, labels),
                  float(rmst[-1] - rmst[0]), (alone, adjusted),
                  cutpoint.scan(months, event, value))  # fmt: skip


def variables(ctx: Context) -> tuple[str, ...]:
    """The patient-table variables the dataset file declares for this cohort."""
    return ctx.spec.clinical.variables if ctx.spec.clinical else ()


def missing_columns(ctx: Context) -> tuple[str, ...]:
    """The columns the missing-values table counts.

    The declared variables, PR and HER2 where ER is one, and the treatment.
    """
    receptors = ("pr", "her2") if "er" in variables(ctx) else ()
    return (*variables(ctx), *receptors, "therapy")


def confounders(ctx: Context) -> tuple[str, ...]:
    """The confounding map's variables: the declared categorical ones, then ``DESIGN``."""
    return (*(v for v in variables(ctx) if v != "age"), *DESIGN)


def endpoint_note(ctx: Context) -> tuple[str, ...]:
    """Why the models fit their endpoint: who chose it, then Table 3's marks for the cohort.

    The first paragraph names the choice and Table 3's mark and counts for it; the second the
    endpoints Table 3 marks for caution or against, with its reason, or that it recommends all.

    Raises:
        ClinicalError: If the dataset file's copy of Table 3 holds a mark not in ``USES``.
    """
    assert ctx.spec.clinical and ctx.spec.gdc  # noqa: S101 - run() refuses a dataset without them
    clinical, kind = ctx.spec.clinical, ctx.spec.gdc.cdr_type or ctx.spec.gdc.project
    chosen = (fill("patients.endpoint.file", endpoint=ctx.endpoint)
              if ctx.endpoint == clinical.endpoint
              else fill("patients.endpoint.option", endpoint=ctx.endpoint,
                        default=clinical.endpoint))  # fmt: skip
    marks = clinical.published.get("table3_use", {})
    counts = clinical.published.get("table3", {}).get(ctx.endpoint)
    if not all(marks.get(e) for e in ENDPOINTS) or not counts:
        return (f"{chosen} {fill('patients.endpoint.uncopied', type=kind)}",)
    if bad := {e: marks[e] for e in ENDPOINTS if marks[e] not in USES}:
        raise ClinicalError(f"{ctx.spec.path.name} [clinical.published.table3_use]: {bad}, "
                            f"each one of {USES}")  # fmt: skip
    use = fill("patients.endpoint.use", endpoint=ctx.endpoint, type=kind, events=counts[0],
               censored=counts[1], use=word(f"patients.use.{marks[ctx.endpoint]}"))  # fmt: skip
    caution = [e for e in ENDPOINTS if marks[e] == "caution"]
    against = [e for e in ENDPOINTS if marks[e] == "not_recommended"]
    if caution and not marks.get("caution"):
        raise ClinicalError(f"{ctx.spec.path.name} [clinical.published.table3_use]: copy "
                            f"Table 3's Explanation/Caution for {caution} as caution")  # fmt: skip
    rest = [fill("patients.endpoint.caution", endpoints=listed(caution), type=kind,
                 why=marks["caution"])] if caution else []  # fmt: skip
    if against:
        rest.append(fill("patients.endpoint.against", endpoints=listed(against), type=kind))
    return (f"{chosen} {use}", " ".join(rest or [fill("patients.endpoint.all", type=kind)]))


def sex_flag(ctx: Context) -> tuple[str, ...]:
    """The less frequent sex, kept and flagged, where the dataset file does not declare sex.

    A cohort that declares sex tests it as a variable, so its men and women are no flag; one
    that does not, a cancer nearly all of one sex, names the few of the other. Empty without.
    """
    if "sex" in variables(ctx):
        return ()
    p = ctx.clinical.patients
    counts = {s: p.filter(pl.col("sex") == s).height for s in SEXES}
    fewer = min(SEXES, key=counts.__getitem__)
    if counts[fewer] == 0 or counts[fewer] == max(counts.values()):
        return ()
    people = word(f"patients.people.{fewer}")
    return (fill("patients.sex_flag", count=counts[fewer], people=people),)


def _patients(writer: Writer, ctx: Context) -> Group:
    """The flow to the analysed patients, each endpoint beside Table 3, the missing values."""
    clinical, p = ctx.clinical, ctx.clinical.patients
    published = ctx.spec.clinical.published if ctx.spec.clinical else {}
    writer.table(p, "patients", figures.columns("columns.patients"))
    flow = ResultTable(word("patients.flow_title"), text.words(PAGE, "patients.flow_columns"),
                       tuple((word(f"patients.flow.{k}"), f"{n:,}") for k, n in clinical.flow),
                       numeric=(1,))  # fmt: skip
    table3 = published.get("table3", {})
    rows = []
    for e in ENDPOINTS:
        analysed = harmonise.analysed(p, e)
        liu = table3.get(e)
        rows.append((e, f"{clinical.flags[e][0]:,}", f"{clinical.flags[e][1]:,}",
                     f"{liu[0]:,}; {liu[1]:,}" if liu else word("none"),
                     f"{analysed.height:,}", f"{int(analysed[e].sum()):,}"))  # fmt: skip
    endpoints = ResultTable(fill("patients.endpoints_title", source=table3.get("source", "")),
                            text.words(PAGE, "patients.endpoints_columns"), tuple(rows),
                            numeric=(1, 2, 4, 5))  # fmt: skip
    gaps = tuple((word(f"patients.variables.{c}"), f"{p[c].null_count():,}",
                  "; ".join(f"{r} {n:,}" for r, n in harmonise.reasons(p, c)) or word("none"))
                 for c in missing_columns(ctx))  # fmt: skip
    heads = text.words(PAGE, "patients.missing_columns")
    missing = ResultTable(word("patients.missing_title"), heads, gaps, numeric=(1,))
    noted = {k: int(p[k].sum()) for k in harmonise.NOTIFICATIONS}
    notified = p.filter(pl.any_horizontal(*harmonise.NOTIFICATIONS)).height
    note = " ".join((fill("patients.notified", any=notified, n=p.height, **noted),
                     *sex_flag(ctx)))  # fmt: skip
    return Group(word("patients.title"), (flow, endpoints, *endpoint_note(ctx), missing, note),
                 lead=fill("patients.lead", n=p.height))  # fmt: skip


def _er_check(writer: Writer, ctx: Context) -> list[Result]:
    """ER-negative against ER-positive under all four endpoints, beside Liu et al.'s p."""
    published = ctx.spec.clinical.published.get("er_logrank", {}) if ctx.spec.clinical else {}
    named = {v: word(f"er.{v}") for v in ER}
    blocks, rows = [], []
    for e in ENDPOINTS:
        out = analyses.outcome(ctx.clinical.patients, e)
        called = out.patients["er"].to_numpy()
        keep = np.isin(called, ER)
        months, event = out.months[keep], out.event[keep]
        labels = np.array([named[v] for v in called[keep]])
        chi, _, p = survival.logrank(months, event, labels)
        model = analyses.cox_model(
            months, event, (called[keep] == "negative").astype(float)[:, None]
        )
        liu = published.get(e)
        values = {"endpoint": e, "n": len(months), "events": int(event.sum()), "chi2": chi,
                  "p": number(p), "liu": number(liu), "hr": number(model.hr, ".2f")}  # fmt: skip
        data = figures.km(analyses.km_frame(months, event, labels, list(named.values())),
                          list(named.values()), fill("figures.er.caption", **values),
                          ctx.options.horizon)  # fmt: skip
        blocks.append(writer.figure(data, f"er_{e.lower()}", figures.words("er", **values)))
        per = [(int((called[keep] == v).sum()), int(event[called[keep] == v].sum())) for v in ER]
        rows.append((e, f"{per[0][0]:,}; {per[1][0]:,}", f"{per[0][1]}; {per[1][1]}", number(p),
                     number(liu), fill("er.hr_value", hr=number(model.hr, ".2f"),
                                       low=number(model.low, ".2f"),
                                       high=number(model.high, ".2f"))))  # fmt: skip
    table = ResultTable(
        fill("er.table_title", source=published.get("source", "")),
        text.words(PAGE, "er.columns"),
        tuple(rows),
        numeric=(1, 2, 3, 4),
    )
    return [
        Group(word("er.title"), tuple(blocks), lead=word("er.lead"), tabs=True),
        Group(word("er.table_group"), (table,)),
    ]


def _dunn(pairs: pl.DataFrame) -> ResultTable:
    """Dunn's pairs of every Kruskal-Wallis row, BH within each test's pairs."""
    rows = tuple((r["gene"], word(f"patients.variables.{r['variable']}"),
                  fill("expression.dunn_pair", a=r["level_a"], b=r["level_b"]), f"{r['z']:.2f}",
                  number(r["p"]), number(r["p_adj"]))
                 for r in pairs.iter_rows(named=True))  # fmt: skip
    return ResultTable(word("expression.dunn_title"), text.words(PAGE, "expression.dunn_columns"),
                       rows, numeric=(3, 4, 5))  # fmt: skip


def _expression(writer: Writer, ctx: Context) -> Group:
    """The build genes against each declared variable, Dunn's pairs, and the box."""
    p = ctx.clinical.patients
    symbols = ctx.gene_table["gene_name"].to_list()
    values = {symbols[r]: ctx.values[r] / score.SCALE for r in ctx.chosen}
    tests, pairs = analyses.expression_tests(p, values, analyses.plans(variables(ctx)))
    writer.table(tests, "expression_tests", figures.columns("columns.tests"))
    writer.table(pairs, "dunn_pairs", figures.columns("columns.dunn"))
    rows = tuple((r["gene"], word(f"patients.variables.{r['variable']}"),
                  word(f"expression.tests.{r['test']}"), f"{r['n']:,}", f"{r['effect']:.3f}",
                  number(r["p"]), number(r["q"])) for r in tests.iter_rows(named=True))  # fmt: skip
    table = ResultTable(word("expression.table_title"), text.words(PAGE, "expression.columns"),
                        rows, numeric=(3, 4, 5, 6))  # fmt: skip
    split = ctx.spec.clinical.box if ctx.spec.clinical else "er"
    data = figures.box(p, values, split)
    first = data.caption.tests[0]
    auc = (first.r + 1) / 2
    box = writer.figure(
        data,
        f"expression_{split}",
        figures.words(f"expression.{split}", gene=first.feature, auc=auc),
    )
    shown = (box, table, _dunn(pairs)) if pairs.height else (box, table)
    return Group(word("expression.title"), shown)


def _models(result: Scored, ctx: Context) -> ResultTable:
    """The Cox models of the score, alone and adjusted: one row each."""
    covariates = listed([word(f"patients.variables.{c}").lower() for c in ctx.options.covariates])
    rows = []
    for name, m in zip((word("models.alone"), fill("models.adjusted", covariates=covariates)),
                       result.models, strict=True):  # fmt: skip
        rows.append((name, f"{m.patients:,}", f"{m.events:,}", number(m.hr, ".2f"),
                     fill("models.ci", low=number(m.low, ".2f"), high=number(m.high, ".2f")),
                     number(m.p), number(m.ph_p), number(m.c_index, ".3f")))  # fmt: skip
    return ResultTable(fill("models.title", endpoint=ctx.endpoint),
                       text.words(PAGE, "models.columns"), tuple(rows),
                       numeric=(1, 2, 3, 5, 6, 7))  # fmt: skip


def _groups_table(result: Scored, ctx: Context) -> ResultTable:
    """Per group: patients, events, median with its interval, RMST."""
    rows = tuple((r["group"], f"{r['patients']:,}", f"{r['events']:,}",
                  number(r["median"], ".1f"),
                  fill("models.ci", low=number(r["low"], ".1f"), high=number(r["high"], ".1f")),
                  f"{r['rmst']:.1f}") for r in result.groups.iter_rows(named=True))  # fmt: skip
    return ResultTable(fill("groups.title", endpoint=ctx.endpoint, horizon=ctx.options.horizon),
                       text.words(PAGE, "groups.columns"), rows, numeric=(1, 2, 3, 5))  # fmt: skip


def values(result: Scored, ctx: Context) -> dict[str, Any]:
    """The numbers the live figures' words print, as the page script fills them too."""
    chi, df, p = result.logrank
    alone, adjusted = result.models
    return {
        "endpoint": ctx.endpoint, "genes": len(ctx.chosen), "source": ctx.source,
        "n": len(result.out.months), "events": int(result.out.event.sum()),
        "split": word(f"splits.{ctx.options.split}"), "chi2": f"{chi:.2f}", "df": df,
        "p": number(p), "horizon": f"{ctx.options.horizon:g}", "rmst": f"{result.rmst:.1f}",
        "hr": number(alone.hr, ".2f"), "hr_adj": number(adjusted.hr, ".2f"),
        "corrected": number(result.scan.corrected), "min_p": number(float(result.scan.p.min())),
        "best": result.scan.best, "perms": cutpoint.PERMUTATIONS, "low": cutpoint.LOW,
        "high": cutpoint.HIGH, "cuts": len(result.scan.cuts), "alpha": analyses.ALPHA,
        "floor": number(1 / (cutpoint.PERMUTATIONS + 1)), "seed": cutpoint.SEED,
    }  # fmt: skip


def live(result: Scored, ctx: Context, screen: pl.DataFrame) -> dict[str, Any]:
    """The survival panel's settings and data.

    The panel's genes as the encoded block, every patient's endpoints in days and covariates,
    the per-gene screen's rows, which the screen's table follows for a scored list, the
    captions the page refills for it, and the words the page script fills.
    """
    p, carried, options = ctx.clinical.patients, ctx.carried, ctx.options
    layers = strata(p, options.strata)
    ids = ctx.gene_table["gene_id"].to_list()
    symbols = ctx.gene_table["gene_name"].to_list()
    words = {k: text.raw(PAGE, f"live.{k}") for k in text.raw(PAGE, "live")}
    return {
        "block": carried.block, "genes": carried.genes, "patients": carried.patients,
        "numbers": [genes.number(ids[r]) for r in carried.rows],
        "symbols": [symbols[r] for r in carried.rows], "build": list(range(len(ctx.chosen))),
        "endpoints": {e: {"days": p[f"{e}_days"].to_list(), "event": p[e].to_list()}
                      for e in ENDPOINTS},
        "covariates": [p[c].to_list() for c in options.covariates],
        "strata": None if layers is None else layers.tolist(),
        "endpoint": ctx.endpoint, "split": options.split, "horizon": options.horizon,
        "days_per_month": harmonise.DAYS_PER_MONTH, "z": float(norm.ppf(0.975)),
        "scan": {"low": cutpoint.LOW, "high": cutpoint.HIGH,
                 "permutations": cutpoint.PERMUTATIONS, "seed": cutpoint.SEED},
        "names": {rule: list(split_names(rule)) for rule in score.SPLITS},
        "screen": screen_rows(screen), "screen_rows": SCREEN_ROWS, "alpha": analyses.ALPHA,
        "captions": {"km": text.raw(PAGE, "figures.km_score.caption"),
                     "scan": text.raw(PAGE, "figures.scan_score.caption")},
        "risk_columns": KM["risk_columns"], "words": words, "source": ctx.source,
        "none": word("none"),
        "leads": {"km": text.raw(PAGE, "figures.km_score.lead"),
                  "scan": text.raw(PAGE, "figures.scan_score.lead"),
                  "scan_floor": text.raw(PAGE, "figures.scan_score.lead_floor")},
        "defaults": {
            "endpoints": [[e, e] for e in ENDPOINTS], "endpoint": ctx.endpoint,
            "splits": [[s, word(f"splits.{s}")] for s in score.SPLITS], "split": options.split,
            "genes": "\n".join(symbols[r] for r in ctx.chosen),
            "summary": fill("live.summary", **values(result, ctx)),
            "notice": fill("live.notice", genes=len(ctx.chosen), source=ctx.source),
            "panel": fill("live.panel", genes=carried.genes, patients=carried.patients),
            "screen_rows": SCREEN_ROWS,
        },
    }  # fmt: skip


def _scored(writer: Writer, ctx: Context, result: Scored, screen: pl.DataFrame
            ) -> list[Result]:  # fmt: skip
    """The score's curves, groups and models (live), then its cut-point scan (live)."""
    shown = values(result, ctx)
    km = figures.km(result.km, list(slots(ctx.options.split)),
                    fill("figures.km_score.caption", **shown), ctx.options.horizon,
                    live=live(result, ctx, screen))  # fmt: skip
    block = writer.figure(km, "km_score", figures.words("km_score", **shown))
    models = pl.DataFrame([vars(m) for m in result.models])
    models = models.with_columns(model=pl.Series(["alone", "adjusted"]))
    writer.table(models, "cox_score", figures.columns("columns.models"))
    writer.table(result.groups, "groups_score", figures.columns("columns.groups"))
    scan = figures.scan(result.scan, fill("figures.scan_score.caption", **shown),
                        live={"kind": "follow", "part": "scan"})  # fmt: skip
    lead = "figures.scan_score.lead" if result.scan.beaten else "figures.scan_score.lead_floor"
    scan_words = replace(figures.words("scan_score", **shown), lead=fill(lead, **shown))
    scan_block = writer.figure(scan, "scan_score", scan_words)
    return [Group(word("score.title"), (block, _groups_table(result, ctx), _models(result, ctx)),
                  lead=fill("score.lead", **shown)),
            Group(word("scan.title"), (scan_block,))]  # fmt: skip


def _competing(writer: Writer, ctx: Context) -> Group:
    """Death from the cancer and from other causes as competing risks, beside 1 - KM of DSS."""
    competed = figures.Competing.of(ctx.clinical.patients)
    at = competed.at(ctx.options.horizon)
    shown = {"n": len(competed.months), "cancer": int((competed.cause == figures.CANCER).sum()),
             "other": int((competed.cause == figures.OTHER).sum()), "unknown": competed.unknown,
             "horizon": f"{ctx.options.horizon:g}", "cif": at["cancer"],
             "naive": at["naive"]}  # fmt: skip
    data = figures.competing(competed, fill("figures.competing.caption", **shown))
    block = writer.figure(data, "competing", figures.words("competing", **shown))
    return Group(word("competing.title"), (block,))


def _power(writer: Writer, ctx: Context, result: Scored) -> Group:
    """Power over events and hazard ratio, this endpoint's events marked."""
    events = int(result.out.event.sum())
    shown = {"events": events, "endpoint": ctx.endpoint, "alpha": analyses.ALPHA,
             "target": analyses.TARGET_POWER, "hr": analyses.detectable(events)}  # fmt: skip
    data = figures.power(events, fill("figures.power_map.caption", **shown))
    block = writer.figure(data, "power", figures.words("power_map", **shown))
    return Group(word("power.title"), (block,))


def screen_table(ctx: Context, result: Scored) -> pl.DataFrame:
    """Cox per panel gene, alone and adjusted, BH over the panel, in the panel's order."""
    p = ctx.clinical.patients
    analysed = p[ctx.endpoint].is_not_null().to_numpy()
    rows = ctx.carried.rows
    values_ = ctx.values[rows][:, analysed].astype(float) / score.SCALE
    extra, complete = analyses.covariates(result.out.patients, ctx.options.covariates)
    found = analyses.screen(result.out, values_, extra, complete)
    ids, symbols = ctx.gene_table["gene_id"].to_list(), ctx.gene_table["gene_name"].to_list()
    names = pl.DataFrame(
        {"gene_id": [ids[r] for r in rows], "gene_name": [symbols[r] for r in rows]}
    )
    return names.hstack(found)


def screen_cells(row: Mapping[str, Any]) -> tuple[str, ...]:
    """One screen row as the table prints it, here and in the live panel's rows."""
    return (row["gene_name"], number(row["hr_alone"], ".2f"), number(row["q_alone"]),
            number(row["hr_adjusted"], ".2f"),
            fill("models.ci", low=number(row["low_adjusted"], ".2f"),
                 high=number(row["high_adjusted"], ".2f")), number(row["q_adjusted"]))  # fmt: skip


def screen_rows(table: pl.DataFrame) -> list[list[Any]]:
    """The live panel's screen: panel row, adjusted q (None where refused), printed cells, sign.

    In the panel's order, so a scored list's rows are read by its panel rows; the page sorts
    them by q as the table does, the refused last. The sign is the side of 1 of the adjusted
    hazard ratio, the one the screen ranks by: a gene taken from the screen into a list adds
    its z where higher expression goes with more risk and subtracts it where with less, so a
    protective gene does not cancel a harmful one in the signed mean (E.3). None where refused.
    """
    out = []
    for i, r in enumerate(table.iter_rows(named=True)):
        q, hr = r["q_adjusted"], r["hr_adjusted"]
        kept = q is not None and np.isfinite(q) and hr is not None and np.isfinite(hr)
        out.append([i, q if kept else None, list(screen_cells(r)),
                    (1 if hr >= 1 else -1) if kept else None])  # fmt: skip
    return out


def _screen(writer: Writer, ctx: Context, table: pl.DataFrame) -> Group:
    """The screen's table file, its top rows, then the build list's."""
    writer.table(table, "cox_screen", figures.columns("columns.screen"))
    top = table.sort("q_adjusted", nulls_last=True).head(SCREEN_ROWS)
    shown = pl.concat([top, table[: len(ctx.chosen)]]).unique("gene_id", maintain_order=True)
    cells = tuple(screen_cells(r) for r in shown.iter_rows(named=True))
    passing = int((table["q_adjusted"] < analyses.ALPHA).sum())
    lead = fill("screen.lead", genes=table.height, n=passing, alpha=analyses.ALPHA,
                endpoint=ctx.endpoint)  # fmt: skip
    return Group(word("screen.title"), (ResultTable(fill("screen.table_title", rows=SCREEN_ROWS),
                                                    text.words(PAGE, "screen.columns"), cells,
                                                    numeric=(1, 2, 3, 5), css="screen"),),
                 lead=lead)  # fmt: skip


def _context_groups(writer: Writer, ctx: Context) -> list[Result]:
    """Confounding between patient variables, and who was treated, both descriptive."""
    p = ctx.clinical.patients
    data, top = figures.confounding(p, confounders(ctx))
    shown = {"a": top["row_variable"], "b": top["column_variable"], "v": top["v"], "n": p.height}
    heat = writer.figure(data, "confounding", figures.words("confounding_map", **shown))
    split = ctx.spec.clinical.box if ctx.spec.clinical else "er"
    levels, names = figures.two_levels(p, split)
    treated = (p["therapy"] == "yes").fill_null(False).to_numpy()
    rows = []
    for name in (*reversed(names), None):
        member = np.array([v == name for v in levels])
        yes = int((member & treated).sum())
        label = name or word(f"treatment.unknown.{split}")
        rows.append((label, f"{member.sum():,}", f"{yes:,}", f"{member.sum() - yes:,}"))
    named = {"therapy": ctx.spec.clinical.therapy if ctx.spec.clinical else "",
             "variable": word(f"patients.variables.{split}")}  # fmt: skip
    heads = (named["variable"], *text.words(PAGE, "treatment.columns"))
    table = ResultTable(fill("treatment.table_title", **named), heads, tuple(rows),
                        numeric=(1, 2, 3))  # fmt: skip
    return [
        Group(word("confounding.title"), (heat,)),
        Group(fill("treatment.title", **named), (table, word("treatment.note"))),
    ]


def _keys(ctx: Context, result: Scored) -> tuple[KeyNumber, ...]:
    """Patients, the endpoint's events, median follow-up and the adjusted HR of the score."""
    out, adjusted = result.out, result.models[1]
    events = int(out.event.sum())
    follow = survival.follow_up(out.months, out.event)
    return (
        KeyNumber(word("keys.patients"), f"{ctx.clinical.patients.height:,}"),
        KeyNumber(fill("keys.events", endpoint=ctx.endpoint), f"{events:,}",
                  fill("keys.events_detail", pct=100 * events / len(out.months))),
        KeyNumber(word("keys.follow_up"), fill("keys.follow_up_value", months=follow),
                  word("keys.follow_up_detail")),
        KeyNumber(word("keys.score_hr"), number(adjusted.hr, ".2f"),
                  fill("keys.score_hr_detail", low=number(adjusted.low, ".2f"),
                       high=number(adjusted.high, ".2f"))),
    )  # fmt: skip


def notes(ctx: Context) -> list[Result]:
    """The data's sources: the GDC acknowledgement, the CDR's licence, the filter's drops."""
    below, flat = ctx.filtered
    assert ctx.spec.gdc is not None  # noqa: S101 - run() refuses a dataset without [gdc]
    return [*gdc.acknowledgement(ctx.spec.gdc), word("cdr.licence"),
            fill("panel.note", kept=ctx.values.shape[0], below=below, flat=flat,
                 carried=ctx.carried.genes, minimum=ctx.options.min_log2cpm)]  # fmt: skip


def build(writer: Writer, ctx: Context) -> list[Path]:
    """Write every table and figure, then the page; the page first in the returned list."""
    result = scored(ctx)
    screen = screen_table(ctx, result)
    figure_tables.write(writer, ctx)
    results: list[Result] = [
        *notes(ctx), _patients(writer, ctx),
        *(_er_check(writer, ctx) if "er" in variables(ctx) else []), _expression(writer, ctx),
        *_scored(writer, ctx, result, screen), _competing(writer, ctx),
        _power(writer, ctx, result), _screen(writer, ctx, screen),
        *_context_groups(writer, ctx),
    ]  # fmt: skip
    usage, budget = runtime.measure(ctx.settings.started), runtime.budget(PAGE)
    built = ModulePage(
        title=fill("title", dataset=ctx.spec.name), **text.skeleton(PAGE),
        parameters=ctx.cards, results=tuple(results), key_numbers=_keys(ctx, result),
        reproduce=Reproduction(
            command=ctx.command, wall_seconds=usage.wall_seconds,
            peak_rss_self=usage.peak_rss_self, peak_rss_children=usage.peak_rss_children,
            profile=f"{ctx.settings.profile}, {ctx.settings.threads} threads",
            budget=fill("budget", minutes=budget.minutes, ram_gb=budget.ram_gb),
            versions=runtime.versions(),
        ),
    )  # fmt: skip
    page = writer.out / f"{stem(PAGE, f'{writer.dataset}_overview', 1, day=writer.day)}.html"
    return [write(built, page), write_entry(writer.out, "clinical"), *writer.files]
