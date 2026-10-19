"""M5's figure kinds: drawn from the tables ``ctx de`` wrote, and gather's for the samples."""

from dataclasses import dataclass

from ctx.dataset import Dataset
from ctx.de.figures import BASES
from ctx.de.stats import ALPHA, LFC_MIN
from ctx.figures.kind import DATASET_FILE, Inputs, Kind, OneContrast, featured, slug, unknown

SIZE_FACTORS = ("de", "size_factors")


@dataclass(frozen=True)
class TwoContrasts:
    """One contrast's fold changes against another's."""

    dataset: str
    x: str
    y: str


@dataclass(frozen=True)
class Lesson:
    """One lesson's call counts at one setting of the call rule (``ctx.de.figures.Setting``)."""

    dataset: str
    lesson: str
    basis: str = "padj"
    alpha: float = ALPHA
    lfc_min: float = LFC_MIN


@dataclass(frozen=True)
class OnePca:
    """One PCA of ``[[de.pca]]``, by its name."""

    dataset: str
    pca: str


def _one(p: OneContrast, _: Dataset) -> Inputs:
    """One contrast's results, the size factors that name its samples, the dataset file."""
    return {"results": ("de", f"{p.contrast}_results"), "size_factors": SIZE_FACTORS,
            "dataset_file": DATASET_FILE}  # fmt: skip


def _lesson_names(p: Lesson, spec: Dataset) -> tuple[str, ...]:
    """The contrasts of one lesson, as the dataset file lists them; none for an unknown one."""
    lessons = spec.de.lessons if spec.de else {}
    return next((names for title, names in lessons.items() if title == p.lesson), ())


def _lesson_inputs(p: Lesson, spec: Dataset) -> Inputs:
    """Every contrast of the lesson's results, the size factors and the dataset file."""
    found = {f"results_{n}": ("de", f"{n}_results") for n in _lesson_names(p, spec)}
    return {**found, "size_factors": SIZE_FACTORS, "dataset_file": DATASET_FILE}


def _pair_inputs(p: TwoContrasts, _: Dataset) -> Inputs:
    """Both contrasts' results, the size factors and the dataset file."""
    return {"results_x": ("de", f"{p.x}_results"), "results_y": ("de", f"{p.y}_results"),
            "size_factors": SIZE_FACTORS, "dataset_file": DATASET_FILE}  # fmt: skip


def _pca_inputs(_: OnePca, __: Dataset) -> Inputs:
    """Gather's sheet and counts, from which M5's one VST is computed again, and the file."""
    return {"samplesheet": ("gather", "samplesheet"), "counts": ("gather", "counts"),
            "dataset_file": DATASET_FILE}  # fmt: skip


def _box_inputs(p: OneContrast, spec: Dataset) -> Inputs:
    """The contrast's fit as M5 wrote it, and gather's counts and genes the boxes read."""
    return {**_one(p, spec), "counts": ("gather", "counts"), "genes": ("gather", "genes")}


def _dispersion_inputs(p: OneContrast, spec: Dataset) -> Inputs:
    """The contrast's fit and the dispersions M5 wrote beside it."""
    return {**_one(p, spec), "dispersions": ("de", f"{p.contrast}_dispersions")}


def _lessons(spec: Dataset) -> list[Lesson]:
    """Every lesson of ``[de.lessons]`` at the build's setting."""
    return [Lesson(spec.slug, title) for title in (spec.de.lessons if spec.de else {})]


def _pairs(spec: Dataset) -> list[TwoContrasts]:
    """Every pair of ``[de] scatters``."""
    return [TwoContrasts(spec.slug, x, y) for x, y in (spec.de.scatters if spec.de else ())]


def _pcas(spec: Dataset) -> list[OnePca]:
    """Every PCA of ``[[de.pca]]``."""
    return [OnePca(spec.slug, p.name) for p in (spec.de.pca if spec.de else ())]


def _boxed(spec: Dataset) -> list[OneContrast]:
    """The box contrast, where the dataset file names one and genes for it."""
    de = spec.de
    return [OneContrast(spec.slug, de.box_contrast)] if de and de.box_contrast else []


def _lesson_check(p: Lesson, spec: Dataset) -> list[str]:
    """The lesson is the dataset file's, and the setting one the call rule can take."""
    found = unknown("lesson", p.lesson, spec.de.lessons if spec.de else {}, "a lesson",
                    spec.slug)  # fmt: skip
    if p.basis not in BASES:
        found.append(f"basis: {p.basis!r} is not one of {', '.join(BASES)}")
    if not 0 < p.alpha <= 1:
        found.append(f"alpha: {p.alpha}, expected above 0 and at most 1")
    if p.lfc_min < 0:
        found.append(f"lfc_min: {p.lfc_min}, expected 0 or more")
    return found


def _pca_check(p: OnePca, spec: Dataset) -> list[str]:
    """The PCA is one of the dataset file's."""
    return unknown("pca", p.pca, [q.name for q in spec.de.pca] if spec.de else [],
                   "a [[de.pca]] name", spec.slug)  # fmt: skip


def lesson_slug(title: str) -> str:
    """A lesson's title as M5's output names carry it (``ctx.de.run._lessons``)."""
    return slug(title)


KINDS = (
    Kind("de", "volcano", OneContrast, "ctx.de.run:volcano", "ctx.figures.draw_m5:volcano",
         featured, lambda p: f"volcano_{p.contrast}", _one),
    Kind("de", "ma", OneContrast, "ctx.de.run:_ma", "ctx.figures.draw_m5:ma", featured,
         lambda p: f"ma_{p.contrast}", _one),
    Kind("de", "phist", OneContrast, "ctx.de.run:p_histogram", "ctx.figures.draw_m5:phist",
         featured, lambda p: f"phist_{p.contrast}", _one),
    Kind("de", "lfc", TwoContrasts, "ctx.de.figures:lfc_scatter", "ctx.figures.draw_m5:lfc",
         _pairs, lambda p: f"lfc_{p.x}_vs_{p.y}", _pair_inputs),
    Kind("de", "calls", Lesson, "ctx.de.figures:calls_bars", "ctx.figures.draw_m5:calls",
         _lessons, lambda p: f"calls_{lesson_slug(p.lesson)}", _lesson_inputs, _lesson_check),
    Kind("de", "pca", OnePca, "ctx.de.run:_pca_figure", "ctx.figures.draw_m5:pca", _pcas,
         lambda p: f"pca_{p.pca}", _pca_inputs, _pca_check),
    Kind("de", "box", OneContrast, "ctx.de.run:_box", "ctx.figures.draw_m5:box", _boxed,
         lambda p: f"box_{p.contrast}", _box_inputs),
    Kind("de", "dispersion", OneContrast, "ctx.de.run:_dispersion",
         "ctx.figures.draw_m5:dispersion", featured, lambda p: f"dispersion_{p.contrast}",
         _dispersion_inputs),
)  # fmt: skip
