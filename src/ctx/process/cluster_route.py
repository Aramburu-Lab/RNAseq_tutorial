"""The cluster's alignment as this run's stages: its indices, alignments and counter (ED.6e-4).

``ctx.process.run`` calls these in place of building an index, aligning and running
featureCounts; MultiQC, the checks, the tables and the page then run as on the local route.
"""

from collections.abc import Mapping
from pathlib import Path

from ctx.process import ProcessError, align, cluster, count, fc_compact, reference
from ctx.process.stages import Built, Counter, Inputs, Options

# The STAR index each featured setting was aligned to on the cluster.
SETTING_INDEX = {"untrimmed": "main", "cut": "cut"}


def _check(problems: list[str]) -> None:
    """Stop on any problem with the copied files, every one listed.

    Raises:
        ProcessError: If any problem was found.
    """
    if problems:
        raise ProcessError("the cluster's alignment cannot be read:\n  " + "\n  ".join(problems))


def _source(inputs: Inputs) -> cluster.Source:
    """The run's cluster source; ``resolve`` sets it on this route."""
    if inputs.source is None:
        raise ProcessError("no --align-dir: this run aligns here")
    return inputs.source


def indices(inputs: Inputs, options: Options, parent: Path) -> dict[str, Built]:
    """The cluster's STAR indices, named by their folders, and the references they index.

    The whole reference is ``inputs.reference``; the half is rebuilt from block E's facts.
    Neither index is on this machine, so neither was built or timed here.
    """
    source = _source(inputs)
    found = {"main": Built(inputs.reference, Path(source.index), built=False, seconds=0.0)}
    if "cut" in options.settings:
        half = reference.from_facts(source.reference("cut"), inputs.gtf, parent)
        _check([] if source.cut_index else [f"{source.folder} names no cut index"])
        found["cut"] = Built(half, Path(source.cut_index), built=False, seconds=0.0)
    return found


def aligned(inputs: Inputs, options: Options, parent: Path
            ) -> tuple[align.Trimmed, align.Alignment, dict[str, align.Alignment]]:  # fmt: skip
    """Fastp's reports, the main alignment and the featured sample's, as keyed folders.

    Raises:
        ProcessError: Listing every absent or mismatched file.
    """
    source = _source(inputs)
    paired = inputs.files["mate"].max() == 2  # noqa: PLR2004 - a pair's second mate
    samples = inputs.files["sample"].unique(maintain_order=True).to_list()
    trimmed, problems = cluster.trimmed(source, inputs.files, parent, options.fastp)
    main, found = cluster.alignment({s: source.sample_dir(s) for s in samples}, parent,
                                    bool(paired), source.index)  # fmt: skip
    problems += found
    featured = {}
    for setting in options.settings:
        index = source.index if SETTING_INDEX[setting] == "main" else source.cut_index
        one, found = cluster.alignment({inputs.featured: source.featured_dir(inputs.featured,
                                                                             setting)},
                                       parent, bool(paired), index)  # fmt: skip
        problems += found
        if one is not None:
            featured[setting] = one
    _check(problems)
    assert trimmed is not None and main is not None  # noqa: S101 - _check stops on either missing
    return trimmed, main, featured


def counter(inputs: Inputs, source: cluster.Source, parent: Path) -> Counter:
    """The cluster's per-sample featureCounts files, joined for each setting asked."""

    def counts(which: str, alignment: align.Alignment,
               settings: Mapping[str, count.FcOptions]) -> count.Counts:  # fmt: skip
        """One setting group's joined counts; ``cut`` reads the featured sample's half."""
        if which == "cut":
            cut = source.featured_dir(inputs.featured, "cut")
            dirs, gtf, kept = {inputs.featured: cut}, reference.ANNOTATION, cut
        else:
            dirs = {s: source.sample_dir(s) for s in alignment.samples}
            gtf, kept = source.gtf, source.folder
        return cluster.Counter(dirs, alignment, parent, gtf, kept / fc_compact.ANNOTATION)(settings)

    return counts
