# ruff: noqa: SIM115, PLR2004 - a proof record, run once and kept as run
"""Package H: copy tcga_kirc.toml's GTEx V11 block into the two lung files, once each."""

from pathlib import Path

D = Path(__file__).resolve().parents[3] / "conf" / "datasets"
kirc = (D / "tcga_kirc.toml").read_text()
start = kirc.index("# GTEx V11 for M9's gene cards")
block = kirc[start : kirc.index("# MSigDB Hallmark for M6", start)]
assert block.count("[[resources]]") == 4, block
anchor = "# TargetScan 8.0 (McGeary"
for name in ("tcga_luad.toml", "tcga_lusc.toml"):
    path = D / name
    text = path.read_text()
    assert text.count(anchor) == 1 and "gtex_gene_tpm" not in text, name
    path.write_text(text.replace(anchor, block + anchor))
    print(name, "inserted", block.count("\n"), "lines")
