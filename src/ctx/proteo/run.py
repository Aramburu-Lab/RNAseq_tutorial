"""``ctx proteo``: read both layers of a dataset, measure, decide the rules, write the page.

Across samples, RNA is gather's counts of the chosen tumours (``gdc.chosen``) as log2(CPM + 1)
and protein the freeze's GCT; across studies, RNA is M5's results for the declared contrast and
protein the copy-number table. Either way each layer is one row per gene before it is joined.
"""

from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Any

import numpy as np
import polars as pl

from ctx import concordance, runtime
from ctx.dataset import Dataset, DatasetError, ProteoSpec, load
from ctx.de.filters import log2cpm
from ctx.de.run import Writer, newest
from ctx.de.sanity import decide, parse
from ctx.function.lite import read_gmt
from ctx.gather import gdc
from ctx.io.naming import MODULES
from ctx.proteo import ProteoError, folds, measure, tables
from ctx.proteo.annotate import pair_samples
from ctx.report.builder import ParameterCard

MODULE = "proteo"
GATHER = ("samplesheet", "counts", "genes")


@dataclass(frozen=True)
class Options:
    """Every option of ``ctx proteo`` but the dataset.

    Attributes:
        dataset_file: The dataset file, with ``[proteo]``.
        method: The per-gene correlation the class figure and tests use, one of ``METHODS``.
        min_quantified: Share of samples a protein needs a value in to be tested.
        lfc_min: Absolute log2 fold change a protein needs to be called changed.
        gene_class: A declared class drawn alone against the others, or ``all``.
    """

    dataset_file: Path
    method: str = "spearman"
    min_quantified: float = measure.MIN_QUANTIFIED
    lfc_min: float = 1.0
    gene_class: str = "all"


@dataclass
class Measured:
    """What either design measured, for the figures, the page and the rules.

    Attributes:
        protein: The protein layer, one row per gene, before the join.
        groups: Each protein sample's group: its plex, or its cell type.
        genes: One row per gene in both layers, with every per-gene number and its class.
        abundance: The RNA abundance each gene's tertile is cut on, in ``genes`` order.
        metrics: (metric, value, detail) rows the sanity rules read.
        notes: Lines the page prints before the figures.
        extra: Design-specific values: the matrices, the fold-change interval, the n.
    """

    protein: tables.Layer
    groups: list[str]
    genes: pl.DataFrame
    abundance: np.ndarray
    metrics: list[tuple[str, float, str]] = field(default_factory=list)
    notes: list[tuple[str, dict[str, object]]] = field(default_factory=list)
    extra: dict[str, Any] = field(default_factory=dict)


def resource(settings: runtime.Settings, spec: Dataset, name: str) -> Path:
    """Where ``ctx setup`` puts a resource of the dataset file."""
    return settings.root / spec.resource(name).dest


def classes(settings: runtime.Settings, spec: Dataset) -> dict[str, set[str]]:
    """Each declared gene class as the union of its sets, read from the GMT resources.

    Raises:
        ProteoError: If a GMT is missing or lacks a named set.
    """
    assert spec.proteo is not None  # noqa: S101 - run() refuses a dataset without [proteo]
    found = {}
    for name, gene_class in spec.proteo.classes.items():
        path = resource(settings, spec, gene_class.resource)
        if not path.is_file():
            raise ProteoError(f"{path} is missing; run ctx setup --dataset-file {spec.path}")
        sets = read_gmt(path)
        if missing := [s for s in gene_class.sets if s not in sets]:
            raise ProteoError(f"{path.name}: no set {missing}")
        found[name] = {g for s in gene_class.sets for g in sets[s]}
    return found


def unique_symbols(names: list[str]) -> tuple[np.ndarray, int]:
    """Rows whose gene symbol no other row shares, and how many symbols were shared."""
    counts = pl.Series("name", names).value_counts()
    twice = set(counts.filter(pl.col("count") > 1)["name"].to_list())
    return np.array([n not in twice and bool(n) for n in names]), len(twice)


def _labelled(measured: Measured, known: dict[str, set[str]]) -> int:
    """Add each gene's class; return how many genes sat in two classes."""
    labels, twice = measure.membership(measured.genes["gene"].to_list(), known)
    measured.genes = measured.genes.with_columns(gene_class=pl.Series(labels))
    return twice


def across(settings: runtime.Settings, spec: Dataset, options: Options) -> Measured:
    """One cohort at both layers: per-gene correlation across the samples in both.

    Raises:
        ProteoError: If no sample is in both layers.
    """
    p, gather = spec.proteo, settings.root / "outputs" / f"{dict(MODULES)['gather']:02d}_gather"
    assert p is not None and p.member and p.plex and spec.gdc  # noqa: S101 - load() checks
    found = {w: pl.read_parquet(newest(gather, spec.slug, w)) for w in GATHER}
    protein, annotation = tables.read_gct(resource(settings, spec, p.protein), p.member)
    sheet = gdc.chosen(found["samplesheet"], spec.gdc.tumour)
    pairs, by_sample_id = pair_samples(protein.samples, sheet, p.sample_prefix)
    both = [j for j, _ in pairs]
    samples = [s for _, s in pairs]
    sizes = dict(sheet.select("sample", "library_size").rows())
    rna_all = log2cpm(found["counts"].select(samples).to_numpy(),
                      np.array([sizes[s] for s in samples]))  # fmt: skip
    names = found["genes"]["gene_name"].to_list()
    single, twice = unique_symbols(names)
    rna_row = {n: i for i, n in enumerate(names) if single[i]}
    genes = [g for g in protein.genes if g in rna_row]
    rows = protein.row()
    rna = rna_all[[rna_row[g] for g in genes]]
    values = protein.values[[rows[g] for g in genes]][:, both]
    frame = measure.per_gene(tuple(genes), rna, values, options.min_quantified)
    out = Measured(protein, annotation[p.plex].to_list(), frame, rna.mean(axis=1))
    out.notes.append(("rna_twice", {"twice": twice}))
    out.notes.append(("paired", {"n": len(both), "n_columns": len(protein.samples),
                                 "by_sample_id": by_sample_id}))  # fmt: skip
    out.extra |= {"rna": rna, "values": values, "n_samples": len(both),
                  "annotation": annotation[both], "pairs": pairs,
                  "sheet": gdc.chosen(found["samplesheet"])}  # fmt: skip
    if p.tracks:
        out.extra["colour"] = [str(v) for v in annotation[p.tracks[0]]]
    return out


def contrast_metrics(settings: runtime.Settings, spec: Dataset, measured: Measured,
                     options: Options) -> None:  # fmt: skip
    """F4 for every declared contrast whose M5 results exist; a missing one is noted as such."""
    de_dir = settings.root / "outputs" / f"{dict(MODULES)['de']:02d}_de"
    found = []
    for contrast in spec.contrasts:
        path = folds.results_file(de_dir, spec.slug, contrast)
        if path is None:
            measured.notes.append(("f4_missing", {"contrast": contrast.name}))
            continue
        measured.extra.setdefault("inputs", {})[f"M5 results {contrast.name}"] = path.name
        made = folds.contrast_folds(
            pl.read_parquet(path), measured.extra["sheet"], measured.extra["pairs"],
            (measured.protein.genes, measured.protein.values), contrast, order=dict(spec.order),
            limits=(options.min_quantified, options.lfc_min),
        )  # fmt: skip
        measured.metrics += made.metrics()
        found.append(made)
    measured.extra["folds"] = found


def cross(settings: runtime.Settings, spec: Dataset, options: Options) -> Measured:
    """RNA and protein from two studies: fold changes of one contrast, gene by gene."""
    p = spec.proteo
    assert p is not None and p.tested and p.reference and p.contrast  # noqa: S101 - load() checks
    protein, n_tested, n_reference = tables.read_copies(resource(settings, spec, p.protein),
                                                        p.tested, p.reference)  # fmt: skip
    group = np.arange(protein.values.shape[1]) < n_tested
    contrast = next(c for c in spec.contrasts if c.name == p.contrast)
    de_dir = settings.root / "outputs" / f"{dict(MODULES)['de']:02d}_de"
    results = pl.read_parquet(newest(de_dir, spec.slug, f"{p.contrast}_results")).filter(
        pl.col("log2FoldChange").is_not_null() & pl.col("gene_name").is_not_null()
    )
    single, twice = unique_symbols(results["gene_name"].to_list())
    rna = results.filter(pl.Series(single)).select(
        gene="gene_name", rna_lfc="log2FoldChange", rna_padj="padj", rna_mean="baseMean"
    )
    layer = folds.protein_calls(protein.genes, protein.values, group, options.min_quantified,
                                options.lfc_min)  # fmt: skip
    joined = rna.join(layer, on="gene", how="inner", maintain_order="right")
    joined = joined.with_columns(difference=pl.col("protein_lfc") - pl.col("rna_lfc"))
    groups = [contrast.tested] * n_tested + [contrast.reference] * n_reference
    out = Measured(protein, groups, joined, np.log2(joined["rna_mean"].to_numpy() + 1))
    out.notes.append(("rna_twice", {"twice": twice}))
    out.extra |= {"contrast": contrast, "n_tested": n_tested, "n_reference": n_reference}
    return out


def fold_metrics(measured: Measured) -> None:
    """The fold-change concordance and the granule shift, over the tested genes in both."""
    tested = measured.genes.filter(pl.col("tested") & pl.col("protein_lfc").is_not_nan())
    rho, low, high, _ = concordance.spearman(tested["rna_lfc"], tested["protein_lfc"])
    measured.extra |= {"fold": (rho, low, high), "tested": tested}
    measured.metrics += [(f"fold.spearman{k}", v if v is not None else np.nan,
                          f"{tested.height} genes, bootstrap 1000, seed 42")
                         for k, v in (("", rho), (".low", low), (".high", high))]  # fmt: skip
    for name in sorted(set(tested["gene_class"]) - {measure.OTHER, ""}):
        inside = tested.filter(pl.col("gene_class") == name)["difference"].to_numpy()
        other = tested.filter(pl.col("gene_class") == measure.OTHER)["difference"].to_numpy()
        p, r = measure.shifted(inside, other)
        detail = f"{inside.size} against {other.size} genes, one-sided Mann-Whitney"
        measured.metrics += [(f"class.{name}.p", p, detail),
                             (f"class.{name}.rank_biserial", r, detail)]  # fmt: skip
        measured.extra[f"shift.{name}"] = (p, r, inside.size, other.size)


def correlation_metrics(measured: Measured, spec: ProteoSpec, method: str) -> None:
    """The per-gene medians, the four genes' r, and the class tests the rules read."""
    tested = measured.genes.filter(pl.col("tested"))
    short = "rho" if method == "spearman" else "r"
    measured.metrics.append(("pergene.pearson.median", measure.median(tested["r"]),
                             f"{tested['r'].drop_nulls().len()} genes"))  # fmt: skip
    measured.metrics.append(("pergene.spearman.median", measure.median(tested["rho"]),
                             f"{tested['rho'].drop_nulls().len()} genes"))  # fmt: skip
    published = spec.published.get("mertins", {})
    # ERBB2 is both a published gene and a worked one: one metric per gene.
    for gene in dict.fromkeys([*(k for k in published if k.isupper()), *spec.genes]):
        row = measured.genes.filter((pl.col("gene") == gene) & pl.col("tested"))
        value = float(row["r"][0]) if row.height and row["r"][0] is not None else np.nan
        measured.metrics.append((f"gene.{gene}.pearson", value, "per-gene Pearson r"))
    groups = {name: tested.filter(pl.col("gene_class") == name)[short].drop_nulls().to_numpy()
              for name in [*spec.classes, measure.OTHER]}  # fmt: skip
    measured.extra["class_test"] = measure.class_test(groups)
    for name, values in groups.items():
        measured.metrics.append((f"class.{name}.median", float(np.median(values)) if values.size
                                 else np.nan, f"{values.size} genes, {method}"))  # fmt: skip
    if {"complex", "metabolic"} <= set(groups):
        p, r = measure.shifted(groups["metabolic"], groups["complex"])
        measured.metrics.append(("class.complex_below_metabolic.p", p,
                                 f"one-sided Mann-Whitney, rank-biserial {r:.3f}"))  # fmt: skip
        measured.extra["complex_below"] = (p, r)


def verdicts(spec: ProteoSpec, metrics: list[tuple[str, float, str]]) -> list[dict[str, object]]:
    """The registered rules decided on this run's metrics."""
    return decide({m: v for m, v, _ in metrics}, parse(spec.sanity))


def run(settings: runtime.Settings, cards: tuple[ParameterCard, ...], command: str,
        dataset: str, options: Options) -> list[Path]:  # fmt: skip
    """Build the M7 tables, figures and page.

    Returns:
        The page, then every table and figure file.

    Raises:
        DatasetError: If the dataset file names another slug or lacks [proteo].
        ProteoError: If an option is out of range or an input cannot make the analysis.
    """
    from ctx.proteo import page  # noqa: PLC0415 - figures pull in plotting; after the profile

    spec = load(options.dataset_file)
    if spec.slug != dataset or spec.proteo is None:
        raise DatasetError(f"{options.dataset_file}: needs slug {dataset!r} and [proteo]")
    if options.method not in measure.METHODS:
        raise ProteoError(f"--method {options.method!r}: one of {measure.METHODS}")
    if not 0 < options.min_quantified <= 1:
        raise ProteoError(f"--min-quantified {options.min_quantified}: a share in (0, 1]")
    if options.gene_class not in {"all", *spec.proteo.classes}:
        raise ProteoError(
            f"--gene-class {options.gene_class!r}: all or {list(spec.proteo.classes)}"
        )
    known = classes(settings, spec)
    design = across if spec.proteo.design == "across_samples" else cross
    measured = design(settings, spec, options)
    measured.notes.append(("classes_note", {"twice": _labelled(measured, known)}))
    if spec.proteo.design == "across_samples":
        correlation_metrics(measured, spec.proteo, options.method)
        contrast_metrics(settings, spec, measured, options)
    else:
        fold_metrics(measured)
    out = settings.root / "outputs" / f"{dict(MODULES)[MODULE]:02d}_{MODULE}"
    releases = {"dataset file": str(spec.path), "protein": spec.resource(spec.proteo.protein).dest,
                **{f"sha256:{r.name}": r.sha256 for r in spec.resources},
                **measured.extra.get("inputs", {})}  # fmt: skip
    writer = Writer(out, dataset, date.today(), releases)
    return page.build(writer, page.Context(spec, options, measured, settings, command, cards))
