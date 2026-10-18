"""Add Liu et al. 2018 Table 3's use marks to each M10 dataset file and the GDC template (wo13-v).

The marks are copied from liu_table3.log, this folder's parse of PubMed Central's XML; each
block goes after its file's [clinical.published.table3] DSS line, exactly one match asserted.
"""

from pathlib import Path

WT = Path(__file__).resolve().parents[4]
# Type: (DSS line of the file, OS, PFI, DFI, DSS marks, Explanation/Caution), from the log.
MARKS = {
    "tcga_brca": (
        "BRCA",
        "DSS = [83, 995]",
        "caution",
        "recommended",
        "recommended",
        "caution",
        "need a longer follow-up for OS and DSS",
    ),
    "tcga_kirc": (
        "KIRC",
        "DSS = [110, 415]",
        "recommended",
        "recommended",
        "caution",
        "recommended",
        "number of events is small",
    ),
    "tcga_luad": (
        "LUAD",
        "DSS = [116, 370]",
        "recommended",
        "recommended",
        "recommended",
        "recommended",
        "",
    ),
    "tcga_lusc": (
        "LUSC",
        "DSS = [91, 361]",
        "recommended",
        "recommended",
        "recommended",
        "recommended",
        "",
    ),
    "tcga_hnsc": (
        "HNSC",
        "DSS = [130, 372]",
        "recommended",
        "recommended",
        "recommended",
        "recommended",
        "",
    ),
}
COMMENT = (
    "\n# Liu et al. 2018 Table 3's Use column for {type}, read through PubMed Central E-utilities\n"
    "# (PMC6066282) on 2026-10-18 (docs/plans/logs/wo13_v/liu_table3.log): a tick is recommended,\n"
    "# a tick with an asterisk caution, a cross not_recommended; caution is the table's\n"
    "# Explanation/Caution, verbatim. M10 prints them beside its endpoint.\n"
)
TEMPLATE = (
    "\n# EDIT: Liu 2018 Table 3's Use column for this type, from the same row: a tick is\n"
    "# recommended, a tick with an asterisk caution, a cross not_recommended; caution is the\n"
    "# row's Explanation/Caution, verbatim, where any mark is caution. Left empty, M10 says the\n"
    "# dataset file holds no copy.\n"
    "[clinical.published.table3_use]\n"
    'source = "Liu et al. 2018 Cell, Table 3, XXXX: recommended use of each endpoint"\n'
    'OS = ""\nPFI = ""\nDFI = ""\nDSS = ""\ncaution = ""\n'
)


def block(kind: str, os_: str, pfi: str, dfi: str, dss: str, caution: str) -> str:
    """One dataset file's table, its source naming the type."""
    lines = [COMMENT.format(type=kind) + "[clinical.published.table3_use]",
             f'source = "Liu et al. 2018 Cell, Table 3, {kind}: recommended use of each endpoint"',
             f'OS = "{os_}"', f'PFI = "{pfi}"', f'DFI = "{dfi}"', f'DSS = "{dss}"']  # fmt: skip
    if caution:
        lines.append(f'caution = "{caution}"')
    return "\n".join(lines) + "\n"


def insert(path: Path, anchor: str, text: str) -> None:
    """Put ``text`` after the one line equal to ``anchor``."""
    body = path.read_text()
    assert body.count(anchor + "\n") == 1, (path, anchor)
    assert "table3_use" not in body, f"{path} has the table already; this script runs once"
    path.write_text(body.replace(anchor + "\n", anchor + "\n" + text, 1))


for slug, (kind, anchor, *marks, caution) in MARKS.items():
    insert(WT / "conf" / "datasets" / f"{slug}.toml", anchor, block(kind, *marks, caution))
insert(WT / "conf" / "dataset_template_gdc.toml", "DSS = [0, 0]", TEMPLATE)
