"""Track B at Step-1 size: 60,000 genes by 171 samples, the Step-1 contrasts at their real n.

Marked ``budget`` and deselected by default; run once with ``-m budget -s`` and paste the
printed measurements into TODO.md. The synthetic design mirrors Step-1's layout (cell types,
batches, individuals shared across cell types and conditions, repeated samples), not its
biology.
"""

import os
import re
import subprocess
import sys
import time
from datetime import date
from pathlib import Path

import pytest
from synth import design_counts, write_gather_outputs, write_gene_sets

from ctx.dataset import Contrast, load
from ctx.de.deseq import fit
from ctx.runtime import BUDGETS

CTX = Path(sys.executable).parent / "ctx"
N_GENES = 60_000
# (condition, batch, n, cell type, individuals): a shared prefix is a shared individual.
GROUPS = (
    ("steady_state", "short", 13, "NDN", "A"), ("steady_state", "short", 1, "NDN", "A"),
    ("gcsf", "short", 13, "NDN", "A"), ("steady_state", "long", 4, "NDN", "B"),
    ("gcsf", "long", 4, "NDN", "B"), ("hsct", "short", 5, "NDN", "C"),
    ("hsct", "short", 5, "NDN", "C"), ("hsct", "short", 1, "NDN", "C"),
    ("hsct", "long", 2, "NDN", "F"), ("pdac", "long", 14, "NDN", "D"),
    ("ipmn", "long", 11, "NDN", "E"),
    ("steady_state", "short", 8, "Monocyte", "A"), ("steady_state", "short", 1, "Monocyte", "A"),
    ("gcsf", "short", 8, "Monocyte", "A"), ("hsct", "short", 4, "Monocyte", "C"),
    ("hsct", "short", 4, "Monocyte", "C"), ("hsct", "short", 2, "Monocyte", "C"),
    ("pdac", "long", 12, "Monocyte", "D"), ("ipmn", "long", 10, "Monocyte", "E"),
    ("steady_state", "long", 4, "Monocyte", "B"), ("gcsf", "long", 4, "Monocyte", "B"),
    ("hsct", "long", 3, "Monocyte", "F"),
    ("gcsf", "short", 13, "LDN", "A"), ("gcsf", "long", 4, "LDN", "B"),
    ("hsct", "short", 10, "LDN", "C"), ("pdac", "long", 11, "LDN", "D"),
)  # fmt: skip
NINE = ", ".join(f'"A{i:03d}"' for i in range(9))
NDN_SHORT = 'cell_type = ["NDN"], batch = ["short"]'
MONO_SHORT = 'cell_type = ["Monocyte"], batch = ["short"]'
NINE_NDN = f"{NDN_SHORT}, donor = [{NINE}]"
NULL = f'{NDN_SHORT}, condition = ["steady_state"]'
SPLIT = 'split = { seed = 42, strata = "sex" }'
SS, PAIR, ADJ = "steady_state", "~donor + group", "~batch + group"
CONTRASTS = (
    ("ndn_gcsf", SS, "gcsf", 13, 13, "~group", False, NDN_SHORT, ""),
    ("mono_gcsf", SS, "gcsf", 8, 8, "~group", False, MONO_SHORT, ""),
    ("ndn_hsct_early", SS, "hsct", 13, 5, "~group", False, NDN_SHORT, ""),
    ("mono_hsct_early", SS, "hsct", 8, 4, "~group", False, MONO_SHORT, ""),
    ("ndn_gcsf_paired", SS, "gcsf", 9, 9, PAIR, True, NINE_NDN, ""),
    ("ndn_gcsf_unpaired9", SS, "gcsf", 9, 9, "~group", False, NINE_NDN, ""),
    ("ndn_pdac_naive", SS, "pdac", 17, 14, "~group", False, 'cell_type = ["NDN"]', ""),
    ("ndn_pdac_adjusted", SS, "pdac", 17, 14, ADJ, False, 'cell_type = ["NDN"]', ""),
    (
        "ndn_pdac_within",
        SS,
        "pdac",
        4,
        14,
        "~group",
        False,
        'cell_type = ["NDN"], batch = ["long"]',
        "",
    ),
    ("ndn_null_split", SS, SS, 7, 6, "~group", False, NULL, SPLIT),
)
SPEC = f"""
[dataset]
slug = "budget"

[[resources]]
name = "sets"
urls = ["http://127.0.0.1/sets.gmt"]
bytes = 0
sha256 = "{"0" * 64}"
dest = "resources/sets/budget.gmt"

[guard]
ids = []
patterns = []
""" + "".join(
    f"""
[[contrasts]]
name = "{name}"
variable = "condition"
reference = "{ref}"
tested = "{tested}"
n = [{a}, {b}]
design = "{design}"
paired = {str(paired).lower()}
filters = {{ {filters} }}
{split}"""
    for name, ref, tested, a, b, design, paired, filters, split in CONTRASTS
) + """
[de]
featured = "ndn_gcsf"
box_contrast = "ndn_hsct_early"
box_genes = ["G00001", "G00002", "G00003"]
scatters = [["ndn_pdac_naive", "ndn_pdac_adjusted"], ["ndn_hsct_early", "mono_hsct_early"]]

[de.lessons]
"Batch" = ["ndn_pdac_naive", "ndn_pdac_adjusted", "ndn_pdac_within"]
"Pairing" = ["ndn_gcsf_unpaired9", "ndn_gcsf_paired"]

[[de.pca]]
name = "pb"
factors = ["cell_type", "condition", "batch"]
strata = "donor"

[[de.pca]]
name = "ss_ndn"
filters = { cell_type = ["NDN"], condition = ["steady_state"] }
factors = ["batch"]
one_per = "donor"

[function]
gene_sets = "sets"
licence = "Synthetic gene sets; no licence."
featured_set = "UP_SET"
running = ["ndn_hsct_early", "UP_SET"]
"""  # fmt: skip


def _measured(html: str) -> tuple[float, int, int]:
    wall = float(re.search(r'<td class="wall">([\d.]+) s</td>', html).group(1))
    self_rss = int(
        re.search(r'<td class="rss-self">[\d.]+ MiB \((\d+) bytes\)</td>', html).group(1)
    )
    children = int(
        re.search(r'<td class="rss-children">[\d.]+ MiB \((\d+) bytes\)</td>', html).group(1)
    )
    return wall, self_rss, children


@pytest.mark.budget
def test_track_b_at_tier1_size(tmp_path):
    spec_path = tmp_path / "budget.toml"
    spec_path.write_text(SPEC)
    spec = load(spec_path)
    synthetic = design_counts(GROUPS, n_genes=N_GENES, batch_lfc=1.0)
    assert synthetic.samples.height == 171
    write_gather_outputs(tmp_path, "budget", synthetic, date.today())
    write_gene_sets(tmp_path / spec.resource("sets").dest, n_genes=N_GENES, n_random=48)
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    threads = min(4, os.cpu_count() or 1)
    for command, folder in (("de", "05_de"), ("function", "06_function")):
        done = subprocess.run(
            [str(CTX), "--root", str(tmp_path), command, "--dataset", "budget",
             "--dataset-file", str(spec_path)],
            cwd=tmp_path, env=env, capture_output=True, text=True, check=False,
        )  # fmt: skip
        assert done.returncode == 0, done.stderr[-3000:]
        (page,) = (tmp_path / "outputs" / folder).glob("*_budget_overview_v1.html")
        wall, self_rss, children = _measured(page.read_text())
        bound = self_rss + threads * children
        size = page.stat().st_size
        print(f"\n{command}: wall {wall:.1f} s, self {self_rss:,} B, largest child {children:,} B, "
              f"bound {bound:,} B ({threads} threads), page {size:,} B")  # fmt: skip
        assert wall < BUDGETS["B"].minutes * 60
        assert bound < BUDGETS["B"].ram_gb * 10**9
        assert size < 15 * 10**6
    contrast = next(c for c in spec.contrasts if c.name == "ndn_pdac_naive")
    assert isinstance(contrast, Contrast)
    for n in (1, threads):
        started = time.perf_counter()
        fit(synthetic.counts, synthetic.samples, contrast, spec.order, threads=n)
        seconds = time.perf_counter() - started
        print(f"ndn_pdac_naive at {N_GENES:,} genes, n_cpus {n}: {seconds:.1f} s")
