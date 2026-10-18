"""Read each scratch root's M10 outputs and check them against the cohort (wo13-v, R1.2).

For every cohort under ``.wo13v_roots/<label>/``: every text file of ``outputs/12_clinical`` (the
page with the plotly bundle cut out, figure pages, SVGs, CSVs, READMEs) is searched for another
cohort's cancer words; the page's visible text is printed for the sentences R1.2 changes; the
Table 3 counts the page prints are compared with the CDR's own, and the sex and Table 3
sentences with the patient table and the dataset file. Usage: check_pages.py <label>
"""

import re
import sys
import tomllib
from html.parser import HTMLParser
from pathlib import Path

import polars as pl

from ctx.report.builder import plotly_js

WT = Path(__file__).resolve().parents[4]
# Words naming each cohort's cancer; a page may print its own and no other's.
WORDS = {
    "tcga_brca": ("breast", "BRCA", "mammary"),
    "tcga_kirc": ("kidney", "renal", "KIRC", "clear cell", "clear-cell", "ccRCC"),
    "tcga_luad": ("lung", "LUAD", "adenocarcinoma"),
    "tcga_lusc": ("lung", "LUSC", "squamous"),
    "tcga_hnsc": ("head and neck", "HNSC", "squamous", "larynx", "laryngeal"),
}
# Concepts of one cancer type the page should print only where the dataset file declares them.
BREAST_TERMS = ("endocrine", "oestrogen", "estrogen", "ER-positive", "ER-negative")
TEXT = (".html", ".svg", ".csv", ".md")


class Visible(HTMLParser):
    """The text a reader sees: everything outside script and style."""

    def __init__(self) -> None:  # noqa: D107
        super().__init__()
        self.skip, self.parts = 0, []

    def handle_starttag(self, tag, attrs):  # noqa: D102, ARG002
        if tag in ("script", "style"):
            self.skip += 1

    def handle_endtag(self, tag):  # noqa: D102
        if tag in ("script", "style"):
            self.skip -= 1
        if tag in ("p", "li", "td", "th", "h2", "h3", "h4", "tr", "figcaption", "div"):
            self.parts.append("\n")

    def handle_data(self, data):  # noqa: D102
        if not self.skip:
            self.parts.append(data)


def visible(html: str) -> str:
    """The page's visible text, one block per line."""
    parser = Visible()
    parser.feed(html)
    return re.sub(r"\n\s*\n+", "\n", "".join(parser.parts))


def foreign(slug: str) -> list[str]:
    """The cancer words of every other cohort that are not also this cohort's own."""
    own = {w.lower() for w in WORDS[slug]}
    return sorted({w for s, ws in WORDS.items() if s != slug for w in ws if w.lower() not in own})


def hits(text: str, words: list[str] | tuple[str, ...]) -> dict[str, int]:
    """Whole-word counts of each word that occurs.

    A code such as BRCA is case-sensitive, so a gene symbol (BRCA1) or a base64 block does not
    count; a word is case-blind.
    """
    edge = "A-Za-z0-9+/"
    found = {w: len(re.findall(rf"(?<![{edge}]){re.escape(w)}(?![{edge}])", text,
                               0 if w.isupper() else re.I))
             for w in words}  # fmt: skip
    return {w: n for w, n in found.items() if n}


def check(root: Path, slug: str) -> bool:
    """Print one cohort's findings; True when every check passes."""
    out = root / "outputs" / "12_clinical"
    spec = tomllib.loads((WT / "conf" / "datasets" / f"{slug}.toml").read_text())
    (page,) = out.glob(f"*_clinical_{slug}_overview_v1.html")
    html = page.read_text(encoding="utf-8")
    bundle = plotly_js()
    ok = True
    named: dict[str, dict[str, int]] = {}
    for path in sorted(out.iterdir()):
        if path.suffix not in TEXT:
            continue
        body = path.read_text(encoding="utf-8", errors="replace").replace(bundle, "")
        if found := hits(body, foreign(slug)):
            named[path.name] = found
    print(f"=== {slug}: {len([p for p in out.iterdir() if p.suffix in TEXT])} text files read")
    print(f"other cohorts' cancer words: {named or 'none'}")
    ok &= not named
    shown = visible(html)
    print(f"breast terms on the page's visible text: {hits(shown, BREAST_TERMS) or 'none'}")
    for line in shown.splitlines():
        if re.search(r"Table 3|flagged by sex|GDC notifications|Liu et al|follow-up ends|"
                     r"systemic therapy|ER-negative", line):  # fmt: skip
            print("  |", line.strip())
    ok &= table3(shown, spec)
    ok &= sexes(out, slug, shown, spec)
    return ok


def table3(shown: str, spec: dict) -> bool:
    """The endpoints table's CDR counts equal Liu's on every endpoint.

    The Table 3 sentence's counts are the dataset file's for the page's endpoint.
    """
    liu = spec["clinical"]["published"]["table3"]
    ok = True
    for e in ("OS", "DSS", "DFI", "PFI"):
        events, censored = liu[e]
        row = f"{e}\n{events:,}\n{censored:,}\n{events:,}; {censored:,}\n"
        equal = row in shown
        print(f"  Table 3 {e}: CDR {events:,} events, {censored:,} censored, equal to Liu: {equal}")
        ok &= equal
    endpoint = spec["clinical"]["endpoint"]
    events, censored = liu[endpoint]
    said = f"with {events:,} events and {censored:,} censored"
    print(f"  endpoint sentence carries {endpoint}'s Table 3 counts ({said}): {said in shown}")
    return ok and said in shown


def sexes(out: Path, slug: str, shown: str, spec: dict) -> bool:
    """The sex sentence is printed only where sex is not a declared variable.

    It carries the patient table's own count of the less frequent sex.
    """
    (table,) = out.glob(f"*_clinical_{slug}_patients_v1.parquet")
    sex = pl.read_parquet(table)["sex"].value_counts().rows()
    counts = {k: n for k, n in sex if k in ("female", "male")}
    print(f"  patient table sex counts: {dict(sorted((str(k), n) for k, n in sex))}")
    said = re.findall(r"([\d,]+) (men|women) are kept and flagged by sex", shown)
    declared = "sex" in spec["clinical"]["variables"]
    if declared:
        print(f"  sex is a declared variable; sex sentences printed: {said or 'none'}")
        return not said
    fewer = min(counts, key=counts.__getitem__)
    want = [(f"{counts[fewer]:,}", {"male": "men", "female": "women"}[fewer])]
    print(f"  sex is not declared; printed {said}, the patient table gives {want}")
    return said == want


def main(label: str) -> None:
    """Check every cohort built under the label; exit 1 on any failure."""
    roots = WT / ".wo13v_roots" / label
    results = {slug: check(roots / slug, slug) for slug in WORDS if (roots / slug).exists()}
    print("\nPASS" if all(results.values()) else f"\nFAIL: {results}")
    sys.exit(0 if all(results.values()) else 1)


if __name__ == "__main__":
    main(sys.argv[1])
