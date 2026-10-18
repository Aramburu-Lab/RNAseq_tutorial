"""wo13-r1 proof: the edits a student makes to the copied GDC template, each an EDIT line.

Run once on ~/ctx_r1_clone/conf/datasets/tcga_hnsc.toml, before any ctx command reads it. The
contrast n values are left at the template's zeros: `ctx run --check` gives them.
"""

from datetime import datetime
from pathlib import Path

PATH = Path.home() / "ctx_r1_clone" / "conf" / "datasets" / "tcga_hnsc.toml"
REGISTERED = datetime.now().strftime("%Y-%m-%d %H:%M")
YE = (
    "Ye et al. 2008 BMC Genomics, doi:10.1186/1471-2164-9-69 (PMC2262071, read in full "
    "2026-10-18): 53 oral tongue squamous carcinomas against 22 matching normal tissues"
)
AVISSAR = (
    "Avissar et al. 2009 Clin Cancer Res, doi:10.1158/1078-0432.CCR-08-3131 (PMC2669849, "
    "read in full 2026-10-18): qRT-PCR on 99 HNSCC tumours against 14 non-diseased head "
    "and neck epithelia, normals unpaired"
)
SITES = (
    "; that cohort is the tongue alone, while TCGA-HNSC spans oral cavity, pharynx and "
    "larynx. Registered blind: no ctx module had run on these files"
)

EDITS = [
    ('slug = "tcga_xxxx"', 'slug = "tcga_hnsc"'),
    ('label = "TCGA project"', 'label = "TCGA head and neck squamous"'),
    ("project TCGA-XXXX, through", "project TCGA-HNSC, through"),
    ('project = "TCGA-XXXX"', 'project = "TCGA-HNSC"'),
    ('cdr_type = "XXXX"', 'cdr_type = "HNSC"'),
    ('ids = ["TCGA-XXXX"]', 'ids = ["TCGA-HNSC"]'),
    ('genes = ["MKI67"]', 'genes = ["MMP1", "CRNN", "CDKN2A", "MKI67"]'),
    ('variables = ["stage", "age", "sex"]', 'variables = ["stage", "grade", "age", "sex"]'),
    ("Table 3, XXXX: events", "Table 3, HNSC: events"),
    ("OS = [0, 0]", "OS = [223, 305]"),
    ("PFI = [0, 0]", "PFI = [198, 330]"),
    ("DFI = [0, 0]", "DFI = [28, 106]"),
    ("DSS = [0, 0]", "DSS = [130, 372]"),
    ("mirnas = []", 'mirnas = ["hsa-miR-221-3p", "hsa-miR-375"]'),
]


def rule(rid: str, gene: str, up: bool, what: str) -> str:
    """One [[sanity.rules]] entry on M5's metrics for a gene, at M5's own thresholds."""
    op, value = (">=", 1.0) if up else ("<=", -1.0)
    verb = "rises" if up else "falls"
    return (
        f'\n[[sanity.rules]]\nid = "{rid}"\n'
        f'text = "{gene} {verb} in tumours against their paired normals."\n'
        f'source = "{YE}; {what}{SITES}"\n'
        f'conditions = [{{ metric = "gene.tumour_normal.{gene}.lfc", op = "{op}", '
        f'value = {value} }}, {{ metric = "gene.tumour_normal.{gene}.padj", op = "<", '
        f"value = 0.05 }}]\n"
    )


def mirna_rule(rid: str, name: str, up: bool, what: str) -> str:
    """One [[ncrna.sanity.rules]] entry on M8's miRNA metrics, at M5's thresholds."""
    op, value = (">=", 1.0) if up else ("<=", -1.0)
    verb = "rises" if up else "falls"
    return (
        f'\n[[ncrna.sanity.rules]]\nid = "{rid}"\n'
        f'text = "{name} {verb} in tumours against their paired normals."\n'
        f'source = "{AVISSAR}; {what}. Registered blind: no ctx module had run on these files"\n'
        f'conditions = [{{ metric = "mirna.{name}.lfc", op = "{op}", value = {value} }}, '
        f'{{ metric = "mirna.{name}.padj", op = "<", value = 0.05 }}]\n'
    )


def main() -> None:
    """Apply every edit once, then register the rules with the time they were written."""
    text = PATH.read_text()
    for old, new in EDITS:
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    old = '[sanity]\nregistered_on = ""\n'
    assert text.count(old) == 1
    rules = (
        rule(
            "HNSC.mmp1",
            "MMP1",
            True,
            "MMP1 up 57.6-fold, FDR 0 (Table 3), in the 35-gene "
            "signature at fold change over 4 and FDR under 0.0001 (Table 1)",
        )
        + rule(
            "HNSC.mmp9",
            "MMP9",
            True,
            "MMP9 up 4.08-fold, FDR 5.85e-05 (Table 3), "
            "confirmed by qRT-PCR on 25 tumours and 12 normals (Figure 4)",
        )
        + rule(
            "HNSC.krt4",
            "KRT4",
            False,
            "KRT4 among the signature's down-regulated "
            "genes at fold change over 4 and FDR under 0.0001 (Table 1)",
        )
        + rule(
            "HNSC.crnn",
            "CRNN",
            False,
            "CRNN among the signature's down-regulated "
            "genes at fold change over 4 and FDR under 0.0001 (Table 1)",
        )
    )
    mirna = (
        "\n[ncrna.sanity]\n"
        "# Registered before the first ctx run on this project; the miRNA names are "
        "GDC's miRBase 21 spellings.\n"
        f'registered_on = "{REGISTERED}"\n'
        + mirna_rule(
            "HNSC.mir375",
            "hsa-miR-375",
            False,
            "miR-375 the one miRNA down in "
            "tumours, about 22-fold on the array, and down in the qRT-PCR series",
        )
        + mirna_rule(
            "HNSC.mir221",
            "hsa-miR-221-3p",
            True,
            "miR-221 up in tumours, the numerator of the miR-221:miR-375 ratio",
        )
    )
    text = text.replace(old, f'[sanity]\nregistered_on = "{REGISTERED}"\n{rules}')
    # [ncrna.sanity] after [ncrna.mirna], whose table ends before the [sanity] comments.
    anchor = "# The hub's rules, decided once M5 and M6 have run."
    assert text.count(anchor) == 1
    text = text.replace(anchor, mirna.lstrip("\n") + "\n" + anchor)
    PATH.write_text(text)
    print(f"registered_on {REGISTERED}")


if __name__ == "__main__":
    main()
