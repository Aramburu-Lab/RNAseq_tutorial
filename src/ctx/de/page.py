"""M5 page copy: the words around the numbers, each from ``content/de.toml``, filled by the run.

Nothing here names a cohort or writes a sentence. A self-check tied to a feature of the design
(a batch term, a paired contrast) is asked only when the dataset file declares that feature. A
contrast reads by its plain label, never its identifier (Rytis 2026-10-01).
"""

import math
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from typing import Any

from ctx.dataset import Column, Contrast, Dataset
from ctx.de.deseq import TOP_GENES, Fit
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.report import text
from ctx.report.builder import FigureText, SelfCheck
from ctx.viz.annotate import Permanova

PAGE = "de"
# D-018, measured on tests/synth.py design_counts, seeds 0 to 199: the share of genes at
# padj < 0.05 that were false, against what BH promised, and after the call's fold-change
# filter. Not this run's numbers; the page says synthetic.
SYNTHETIC = {"per_group": 6, "false_pct": 8.4, "bh_pct": 4.5, "filtered_pct": 3.2}


@dataclass(frozen=True)
class Words:
    """How one contrast reads on the page: its label and its two levels in plain words."""

    label: str
    tested: str
    reference: str


def level_words(columns: Sequence[Column], variable: str, level: str) -> str:
    """A level as the dataset's own table spells it, sentence case kept for acronyms.

    The column's value map sends each SDRF value to its curated level; read backwards, a level
    is the SDRF's words: ``steady_state`` is "steady state", ``gcsf`` "G-CSF". A level with no
    map reads with spaces for underscores.
    """
    found = next((c.values for c in columns if c.name == variable and c.values), None)
    words = next((k for k, v in (found or {}).items() if v == level and k), None)
    if words is None:
        return level.replace("_", " ")
    # "Steady state" mid-sentence is "steady state"; "G-CSF" and "HSCT" stay as they are.
    return words if any(c.isupper() for c in words[1:]) else words[:1].lower() + words[1:]


def words(contrast: Contrast, columns: Sequence[Column] = ()) -> Words:
    """The contrast's plain label, its own or its two levels', and the levels in words."""
    variable = contrast.variable if contrast.split is None else ""
    tested = level_words(columns, variable, contrast.tested)
    reference = level_words(columns, variable, contrast.reference)
    label = contrast.label or text.fill(PAGE, "contrast.label", tested=tested, reference=reference)
    return Words(label, tested, reference)


def every(contrasts: Sequence[Contrast], columns: Sequence[Column]) -> dict[str, Words]:
    """Each contrast in words, by name; two with no label that read alike are told apart.

    A file that labels none of its contrasts can declare two of the same levels (Cassetta's
    breast files, the synthetic fixture), and two bars of one name would leave one unseen. Such
    a pair reads as its levels with its identifier after them, until the file labels them.
    """
    found = {c.name: words(c, columns) for c in contrasts}
    seen = Counter(w.label for w in found.values())
    return {
        c.name: replace(found[c.name], label=text.fill(PAGE, "contrast.alike",
                                                       label=found[c.name].label, name=c.name))
        if seen[found[c.name].label] > 1 and not c.label else found[c.name]
        for c in contrasts
    }  # fmt: skip


def skeleton() -> dict[str, Any]:
    """The page's fixed sections, their numbers filled, under ``ModulePage``'s names."""
    return text.skeleton(PAGE, alpha=ALPHA, lfc_min=f"{LFC_MIN:g}", **SYNTHETIC)


def _terms(design: str) -> set[str]:
    """The terms of a design formula, e.g. ``{"batch", "group"}`` for ``~batch + group``."""
    return set(re.findall(r"[A-Za-z_]\w*", design))


def _check(key: str, **values: Any) -> SelfCheck:
    """One self-check from ``[questions.<key>]``."""
    return SelfCheck(text.fill(PAGE, f"questions.{key}.prompt", **values),
                     text.word(PAGE, f"questions.{key}.answer"))  # fmt: skip


def questions(spec: Dataset, fits: Mapping[str, Fit]) -> tuple[SelfCheck, ...]:
    """The self-checks this dataset's design can answer, with its own labels and counts."""
    checks = []
    assert spec.de is not None  # noqa: S101 - M5 runs only with [de]
    for names in spec.de.lessons.values():
        batch = [
            (a, b)
            for i, a in enumerate(names)
            for b in names[i + 1 :]
            if set(fits[a].samples["sample"]) == set(fits[b].samples["sample"])
            and _terms(fits[a].contrast.design) ^ _terms(fits[b].contrast.design) == {"batch"}
        ]
        if batch:
            a, b = (words(fits[n].contrast, spec.columns).label for n in batch[0])
            checks.append(_check("batch", a=a, b=b, n=fits[batch[0][0]].samples.height))
            break
    if any(c.paired for c in spec.contrasts):
        checks.append(_check("paired"))
    return (*checks, _check("na"))


def figure(name: str, key: str, **values: Any) -> str:
    """One figure's word at ``figures.<name>.<key>``, filled."""
    return text.fill(PAGE, f"figures.{name}.{key}", **values)


def reading(name: str, key: str = "reading", **values: Any) -> tuple[str, ...]:
    """A figure's reading guide, each line filled."""
    lines = text.raw(PAGE, f"figures.{name}.{key}")
    return tuple(line.format_map(values) for line in lines)


def subset_words(filters: Mapping[str, Sequence[str]], columns: Sequence[Column]) -> str:
    """The samples a figure uses, in words: each level as the SDRF spells it."""
    return "; ".join(
        figure("pca", "filter", column=k.replace("_", " "),
               values=", ".join(level_words(columns, k, level) for level in v))
        for k, v in filters.items()
    )  # fmt: skip


def pca_text(name: str, subset: str, factor: str, symbol: str | None,
             share: Sequence[float], test: tuple[str, Permanova] | None) -> FigureText:  # fmt: skip
    """Title, lead and reading guide of one PCA; ``test`` is the tested factor and its test.

    The colour and the test can name different factors: a factor with one level colours the
    points but cannot be tested, and the next factor is.
    """
    shown = factor.replace("_", " ")
    shares = {"pc1": share[0], "pc2": share[1]}
    lines = [figure("pca", "reading_points", top=TOP_GENES)]
    lines.append(figure("pca", "reading_shape", colour=shown, symbol=symbol.replace("_", " "),
                        **shares) if symbol else figure("pca", "reading_colour", colour=shown,
                                                        **shares))  # fmt: skip
    lead = ""
    if test is not None:
        tested, result = test[0].replace("_", " "), test[1]
        p = figure("pca", "no_p") if math.isnan(result.p) else figure("pca", "p", p=result.p)
        lines.append(figure("pca", "reading_test", factor=tested, r2=result.r2, p=p))
        lead = figure("pca", "lead", factor=tested, r2=result.r2)
    title = (figure("pca", "title_subset", subset=subset, colour=shown) if subset
             else figure("pca", "title", colour=shown))  # fmt: skip
    return FigureText(title=title, lead=lead, reading=tuple(lines), tab=name)


def volcano_text(said: Words) -> FigureText:
    """Title of the volcano; its lead and reading guide come live from the figure."""
    return FigureText(title=figure("volcano", "title", label=said.label),
                      tab=figure("volcano", "tab"))  # fmt: skip


def ma_text(said: Words, shrunk: bool = True) -> FigureText:
    """Title, lead and reading guide of the MA plot; the raw estimate alone without shrinkage."""
    suffix = "" if shrunk else "_raw"
    return FigureText(
        title=figure("ma", "title", label=said.label),
        lead=figure("ma", f"lead{suffix}"),
        reading=reading("ma", f"reading{suffix}"),
        tab=figure("ma", "tab"),
    )


def phist_text(
    said: Words, tested: int, low: int, flat: float, first: int, width: float
) -> FigureText:
    """Title, lead and reading guide of the p-value histogram.

    ``low`` counts p below ALPHA and ``first`` the first bar; they agree only when the bin
    width equals ALPHA, so each sentence names its own.
    """
    ratio = figure("phist", "times", ratio=first / flat) if flat else figure("phist", "far")
    return FigureText(
        title=figure("phist", "title", label=said.label),
        lead=figure("phist", "lead", low=low, alpha=ALPHA, flat=flat),
        reading=reading("phist", width=width, tested=tested, flat=flat, first=first, ratio=ratio),
        tab=figure("phist", "tab"),
    )


def dispersion_text(said: Words) -> FigureText:
    """Title, lead and reading guide of the dispersion plot."""
    return FigureText(
        title=figure("dispersion", "title", label=said.label),
        lead=figure("dispersion", "lead"),
        reading=reading("dispersion"),
        tab=figure("dispersion", "tab"),
    )


def calls_text(lesson: str, labels: Mapping[str, str], fits: Mapping[str, Fit]) -> FigureText:
    """Title, lead and reading guide of one lesson's call counts; ``labels`` by contrast."""
    counts = {n: int((fits[n].results["call"] != "ns").sum()) for n in labels}
    most = max(counts, key=lambda n: counts[n])
    return FigureText(
        title=figure("calls", "title", lesson=lesson.lower()),
        lead=figure("calls", "lead", label=labels[most], n=counts[most]),
        reading=reading("calls", alpha=ALPHA, lfc_min=f"{LFC_MIN:g}"),
    )


def lfc_text(x: str, y: str, rho: float) -> FigureText:
    """Title, lead and reading guide of one contrast's fold changes against another's."""
    return FigureText(
        title=figure("lfc", "title", x=x, y=y),
        lead=figure("lfc", "lead", rho=rho),
        reading=reading("lfc", x=x, y=y),
    )


def box_text(said: Words, n_genes: int) -> FigureText:
    """Title, lead and reading guide of the named-gene boxes."""
    return FigureText(
        title=figure("box", "title", label=said.label),
        reading=reading("box", n_genes=n_genes, tested=said.tested),
    )
