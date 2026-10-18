"""Print Liu et al. 2018 Table 3's rows for the cohorts M10 builds (wo13-v, R1.2).

Reads the PubMed Central XML given as the one argument (efetch db=pmc id=PMC6066282, no email or
tool parameter sent) and prints each row of the table and the table's footnotes as the XML gives
them, so the dataset files' recommendation tables are copied from the paper, not from memory.
"""

import sys
import xml.etree.ElementTree as ET

TYPES = {"BRCA", "KIRC", "LUAD", "LUSC", "HNSC"}


def cells(row: ET.Element) -> list[str]:
    """The text of every cell of one table row."""
    return ["".join(c.itertext()).strip() for c in row if c.tag in ("td", "th")]


def main(path: str) -> None:
    """Print the header rows, the five cohorts' rows and the footnotes of Table 3."""
    root = ET.parse(path).getroot()  # noqa: S314 - NCBI's own XML, read once
    for wrap in root.iter("table-wrap"):
        label = "".join(wrap.findtext("label") or "")
        if label.strip() != "Table 3":
            continue
        print("# caption:", " ".join("".join(wrap.find("caption").itertext()).split()))
        for i, row in enumerate(wrap.iter("tr")):
            got = cells(row)
            if i < 2 or (got and got[0] in TYPES):  # noqa: PLR2004 - the two header rows
                print(got)
        for foot in wrap.iter("table-wrap-foot"):
            print("# foot:", " ".join("".join(foot.itertext()).split()))


if __name__ == "__main__":
    main(sys.argv[1])
