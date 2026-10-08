"""M8's miRNA targets: TargetScan 8.0 predictions and miRTarBase strong evidence, kept apart.

Validated means a row of miRTarBase's strong-evidence file, a reporter assay or a Western blot
(plan decision 2), for a human miRNA and a human gene. Predicted means a gene of TargetScan's
default predictions for the miRNA; its score is TargetScan's cumulative weighted context++
score, the sum over the gene's sites, lower being stronger. A gene's class for one miRNA is
decided in one order: validated first, then predicted at the preset, then a non-target only if
TargetScan predicts no site of that miRNA in it at all, so a weak prediction is never counted
as a non-target. Every class is cut to the genes M5 tested and TargetScan's universe.
"""

import io
import zipfile
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import polars as pl

from ctx.dataset import Dataset
from ctx.ncrna import NcrnaError
from ctx.process.reference import sha256

HUMAN_TAX = 9606
HUMAN = "hsa"
RESOURCES = ("targetscan_predictions", "targetscan_families", "targetscan_genes",
             "mirtarbase_strong")  # fmt: skip
VALIDATED, PREDICTED, NON_TARGET = "validated", "predicted only", "non-target"
CLASSES = (VALIDATED, PREDICTED, NON_TARGET)
# Cumulative weighted context++ cutoffs the page offers; the default is the second.
PRESETS = (-0.1, -0.2, -0.4)
DEFAULT_PRESET = PRESETS[1]


@dataclass(frozen=True)
class Targets:
    """Both sources read, human rows only, each gene by its unversioned Ensembl ID or symbol.

    Attributes:
        predicted: mirna, gene_id, gene_name, score (cumulative weighted context++), sites.
        validated: mirna, gene_name, experiments, references (distinct PMIDs).
        names: mimat to miRBase name, from TargetScan's family file.
        universe: TargetScan's human genes, unversioned Ensembl IDs.
    """

    predicted: pl.DataFrame
    validated: pl.DataFrame
    names: Mapping[str, str]
    universe: frozenset[str]


def located(root: Path, spec: Dataset) -> tuple[dict[str, Path], list[str]]:
    """Each target resource on disk with the SHA-256 the dataset file declares, or the problems."""
    paths, problems = {}, []
    declared = {r.name: r for r in spec.resources}
    for name in RESOURCES:
        resource = declared.get(name)
        if resource is None:
            problems.append(f"{spec.path.name} declares no {name} resource")
            continue
        path = root / resource.dest
        if not path.is_file():
            problems.append(f"{name} is not at {path}; run ctx setup with this dataset file")
        elif (digest := sha256(path)) != resource.sha256:
            problems.append(f"{path} has SHA-256 {digest}, not the declared {resource.sha256}")
        else:
            paths[name] = path
    return paths, problems


def _zipped(path: Path) -> bytes:
    """The one member of a zip archive."""
    with zipfile.ZipFile(path) as archive:
        members = archive.namelist()
        if len(members) != 1:
            raise NcrnaError(f"{path.name}: {len(members)} members, not one")
        return archive.read(members[0])


def _unversioned(column: str) -> pl.Expr:
    """An Ensembl ID without its version."""
    return pl.col(column).str.replace(r"\.\d+$", "")


def read_predictions(path: Path) -> pl.DataFrame:
    """TargetScan's default predictions: one row per human miRNA and gene, sites summed."""
    columns = {"Gene ID": pl.String, "Gene Symbol": pl.String, "Gene Tax ID": pl.Int64,
               "miRNA": pl.String, "weighted context++ score": pl.Float64}  # fmt: skip
    frame = pl.read_csv(io.BytesIO(_zipped(path)), separator="\t", columns=list(columns),
                        schema_overrides=columns, null_values=["NULL"])  # fmt: skip
    human = frame.filter((pl.col("Gene Tax ID") == HUMAN_TAX)
                         & pl.col("miRNA").str.starts_with(f"{HUMAN}-"))  # fmt: skip
    return (human.group_by("miRNA", "Gene ID")
            .agg(pl.col("Gene Symbol").first(),
                 pl.col("weighted context++ score").sum().alias("score"),
                 pl.len().alias("sites"))
            .select(mirna="miRNA", gene_id=_unversioned("Gene ID"), gene_name="Gene Symbol",
                    score="score", sites="sites")
            .sort("mirna", "gene_id"))  # fmt: skip


def read_names(path: Path) -> dict[str, str]:
    """Human MIMAT accession to miRBase name, from TargetScan's family file."""
    frame = pl.read_csv(io.BytesIO(_zipped(path)), separator="\t", infer_schema=False)
    human = frame.filter(pl.col("Species ID") == str(HUMAN_TAX))
    return dict(human.select("MiRBase Accession", "MiRBase ID").unique().rows())


def read_universe(path: Path) -> frozenset[str]:
    """TargetScan's human genes, unversioned: the genes any prediction could name."""
    frame = pl.read_csv(io.BytesIO(_zipped(path)), separator="\t", infer_schema=False)
    human = frame.filter(pl.col("Species ID") == str(HUMAN_TAX))
    return frozenset(human.select(_unversioned("Gene ID"))["Gene ID"].to_list())


def read_validated(path: Path) -> pl.DataFrame:
    """The miRTarBase strong-evidence rows of a human miRNA and gene, one row per pair."""
    text = path.read_text(encoding="utf-8-sig")
    frame = pl.read_csv(io.StringIO(text), infer_schema=False)
    human = frame.filter((pl.col("Species (miRNA)") == HUMAN)
                         & (pl.col("Species (Target Gene)") == HUMAN))  # fmt: skip
    return (human.group_by("miRNA", "Target Gene")
            .agg(pl.col("Experiments").str.split("//").explode(empty_as_null=False).unique()
                 .sort().str.join("; "),
                 pl.col("References (PMID)").n_unique().alias("references"))
            .select(mirna="miRNA", gene_name="Target Gene", experiments="Experiments",
                    references="references")
            .sort("mirna", "gene_name"))  # fmt: skip


def read(paths: Mapping[str, Path]) -> Targets:
    """All four target files."""
    return Targets(read_predictions(paths["targetscan_predictions"]),
                   read_validated(paths["mirtarbase_strong"]),
                   read_names(paths["targetscan_families"]),
                   read_universe(paths["targetscan_genes"]))  # fmt: skip


def classes(found: Targets, mirna: str, tested: pl.DataFrame, preset: float) -> pl.DataFrame:
    """Each tested gene in TargetScan's universe, classed for one miRNA.

    Args:
        found: The target files.
        mirna: The miRBase name, e.g. hsa-miR-210-3p.
        tested: gene_id (unversioned) and gene_name of the genes M5 tested.
        preset: The cumulative weighted context++ score a predicted target must reach.

    Returns:
        gene_id, gene_name, score (null where TargetScan names no site) and ``target``, one of
        ``CLASSES`` or null for a weaker prediction, which is neither.
    """
    genes = tested.filter(pl.col("gene_id").is_in(list(found.universe)))
    predicted = found.predicted.filter(pl.col("mirna") == mirna).select("gene_id", "score")
    validated = set(found.validated.filter(pl.col("mirna") == mirna)["gene_name"])
    frame = genes.join(predicted, on="gene_id", how="left", maintain_order="left")
    target = (
        pl.when(pl.col("gene_name").is_in(list(validated))).then(pl.lit(VALIDATED))
        .when(pl.col("score") <= preset).then(pl.lit(PREDICTED))
        .when(pl.col("score").is_null()).then(pl.lit(NON_TARGET))
        .otherwise(pl.lit(None, dtype=pl.String))
    )  # fmt: skip
    return frame.with_columns(target=target)


def shift(values: np.ndarray, background: np.ndarray) -> dict[str, float]:
    """A class's fold changes against the non-targets': KS D, Mann-Whitney p and median shift."""
    from scipy.stats import ks_2samp, mannwhitneyu  # noqa: PLC0415 - only the target tests

    if values.size == 0 or background.size == 0:
        return {"n": float(values.size), "shift": np.nan, "d": np.nan, "p_ks": np.nan,
                "p_mw": np.nan}  # fmt: skip
    ks = ks_2samp(values, background)
    return {"n": float(values.size), "shift": float(np.median(values) - np.median(background)),
            "d": float(ks.statistic), "p_ks": float(ks.pvalue),
            "p_mw": float(mannwhitneyu(values, background).pvalue)}  # fmt: skip


def derived(found: Targets, mirnas: Collection[str], genes: pl.DataFrame) -> pl.DataFrame:
    """The target rows the kit ships: the featured miRNAs' targets, each with its source.

    One row per miRNA, gene and source, validated and predicted never merged; the raw files
    stay linked, not shipped (plan decision 3).
    """
    names = genes.select("gene_id", "gene_name")
    predicted = (found.predicted.filter(pl.col("mirna").is_in(list(mirnas)))
                 .join(names, on="gene_id", how="inner", suffix="_v36")
                 .select("mirna", "gene_id", gene_name="gene_name_v36", source=pl.lit("TargetScan"),
                         score="score", evidence=pl.format("{} sites", "sites")))  # fmt: skip
    validated = (found.validated.filter(pl.col("mirna").is_in(list(mirnas)))
                 .join(names, on="gene_name", how="inner")
                 .select("mirna", "gene_id", "gene_name", source=pl.lit("miRTarBase"),
                         score=pl.lit(None, dtype=pl.Float64), evidence="experiments"))  # fmt: skip
    return pl.concat([validated, predicted]).sort("mirna", "source", "gene_name")
