"""M9's figure kinds: the evidence view and each card's panels, drawn from M9's tables.

A card's gene comes from the dataset file's box genes, the featured contrast's top calls, or
``--genes`` and ``--genes-from``, so no kind can list the cards; each reads its gene from the
figure's description instead (``Kind.parsed``), as M9 names it (``ctx.genecard.page._slug``).
"""

import re
from collections import Counter
from collections.abc import Callable
from dataclasses import dataclass

from ctx.dataset import Dataset
from ctx.figures.kind import DATASET_FILE, Inputs, Kind, Plain

GC = "genecard"
# A card's name in an output name: its symbol in lower case, then its Ensembl number.
CARD = re.compile(r"[a-z0-9_]*_[0-9]+")
# The own-sample panels, whose values table names each contrast by its words.
OWN = ("box", "contrasts")


@dataclass(frozen=True)
class OneCard:
    """One panel of one gene's card."""

    dataset: str
    gene: str


def _evidence_view(spec: Dataset) -> list[Plain]:
    """The evidence view, which M9 draws for two contrasts or more."""
    return [Plain(spec.slug)] if len(spec.contrasts) >= 2 else []  # noqa: PLR2004 - see page


def words_unique(spec: Dataset) -> bool:
    """Whether every contrast has words of its own, so the values table names it once."""
    from ctx.genecard import boxes  # noqa: PLC0415 - genecard imports ctx.viz.save, then this

    return max(Counter(boxes.names(spec).values()).values(), default=1) == 1


def _card_check(kind: str) -> Callable[[OneCard, Dataset], list[str]]:
    """The problems of a card's parameters: a gene named as M9 names it, contrasts told apart."""

    def check(p: OneCard, spec: Dataset) -> list[str]:
        """See ``_card_check``."""
        found = [] if CARD.fullmatch(p.gene) else [
            f"gene: {p.gene!r} is not a card's name, its symbol in lower case, '_', its number"
        ]  # fmt: skip
        if kind in OWN and not words_unique(spec):
            found.append(f"{spec.slug}: two contrasts share their words, which M9's values name")
        return found

    return check


def _card_parsed(kind: str) -> Callable[[str, Dataset], OneCard | None]:
    """The card a description names for ``kind``, or None for another kind's."""
    check = _card_check(kind)

    def parsed(what: str, spec: Dataset) -> OneCard | None:
        """See ``_card_parsed``."""
        if not what.startswith(f"{kind}_"):
            return None
        params = OneCard(spec.slug, what[len(kind) + 1 :])
        return None if check(params, spec) else params

    return parsed


def _tables(*names: str) -> Callable[[OneCard | Plain, Dataset], Inputs]:
    """M9's tables by name, and the dataset file."""
    return lambda _, __: {**{n: (GC, n) for n in names}, "dataset_file": DATASET_FILE}


def _options(_: OneCard, __: Dataset) -> Inputs:
    """The immune panel's table and the run's HPA release."""
    return {"immune": (GC, "immune"), "genecard_options": (GC, "figure_options"),
            "dataset_file": DATASET_FILE}  # fmt: skip


def _card_kind(kind: str, plotting: str, inputs: Callable[[OneCard, Dataset], Inputs]) -> Kind:
    """One panel of a card, its gene read from the description."""
    return Kind(GC, kind, OneCard, plotting, f"ctx.figures.draw_m9:{kind}", lambda _: [],
                lambda p: f"{kind}_{p.gene}", inputs, _card_check(kind),
                parsed=_card_parsed(kind))  # fmt: skip


KINDS = (
    Kind(GC, "evidence_view", Plain, "ctx.genecard.figures:evidence_view",
         "ctx.figures.draw_m9:evidence_view", _evidence_view, lambda _: "evidence_view",
         _tables("evidence_listed")),
    _card_kind("box", "ctx.genecard.boxes:box", _tables("own_values", "evidence", "own_family")),
    _card_kind("contrasts", "ctx.genecard.boxes:contrasts", _tables("own_values")),
    _card_kind("tissues", "ctx.genecard.figures:tissues", _tables("tissues")),
    _card_kind("exons", "ctx.genecard.figures:exons",
               _tables("exons", "tissues", "exons_dropped")),
    _card_kind("junctions", "ctx.genecard.figures:junctions", _tables("junctions", "mane")),
    _card_kind("immune", "ctx.genecard.figures:immune", _options),
)  # fmt: skip
