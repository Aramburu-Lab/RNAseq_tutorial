"""``[counts]``: a count matrix from a repository, for ``ctx gather --counts``."""

import re
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from ctx.dataset.checks import _need

# What a count matrix's values are, as [counts] kind declares it: raw counts, or RSEM's
# expected_count, which ctx gather rounds (ctx.gather.counts_rsem).
COUNT_KINDS = ("raw", "rsem_expected")


@dataclass(frozen=True)
class CountsSpec:
    """Where a count matrix from a repository comes from, for ``ctx gather --counts``.

    Attributes:
        source: The repository, accession and file, as the outputs' READMEs name them.
        url: Where the file was downloaded from; empty when it is the student's own.
        annotation: The annotation release the source says its gene IDs come from, quoted
            from the source; empty where it names none.
        kind: One of ``COUNT_KINDS``: ``raw``, or ``rsem_expected`` for RSEM's expected_count.
        sample_pattern: A regular expression with one group, matched against each RSEM file's
            name, whose group is the sample; empty takes the name less its RSEM ending.
    """

    source: str
    url: str = ""
    annotation: str = ""
    kind: str = COUNT_KINDS[0]
    sample_pattern: str = ""


def _counts(raw: Mapping[str, Any], problems: list[str]) -> CountsSpec | None:
    """Parse [counts]: the source is required; url, annotation, kind and sample_pattern not."""
    source = _need(raw, "source", str, "counts", problems)
    optional = {k: raw.get(k, "") for k in ("url", "annotation", "kind", "sample_pattern")}
    problems += [f"counts.{k}: expected str, found {v!r}" for k, v in optional.items()
                 if not isinstance(v, str)]  # fmt: skip
    problems += [f"counts.{k}: unknown key" for k in sorted(set(raw) - {"source", *optional})]
    kind = optional["kind"] or COUNT_KINDS[0]
    if kind not in COUNT_KINDS:
        problems.append(f"counts.kind: {kind!r} is not one of {list(COUNT_KINDS)}")
    if isinstance(pattern := optional["sample_pattern"], str) and pattern:
        try:
            if re.compile(pattern).groups < 1:
                problems.append(f"counts.sample_pattern: {pattern!r} has no group")
        except re.error as error:
            problems.append(f"counts.sample_pattern: {pattern!r}: {error}")
    if source is None or not all(isinstance(v, str) for v in optional.values()):
        return None
    return CountsSpec(source, str(optional["url"]), str(optional["annotation"]), str(kind),
                      str(pattern))  # fmt: skip
