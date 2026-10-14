"""Package H: each ccRCC STAR-Counts file's first-strand share, by gather's own reader and rule."""

import sys
from collections import Counter
from pathlib import Path

from ctx.dataset import load
from ctx.gather import gdc
from ctx.qc import registry

spec = load(Path(sys.argv[1]))
assert spec.gdc is not None
folder = gdc.read_folder(Path("/Users/rytis/ctx_kidney/data/raw/cptac_ccrcc_star"), spec.gdc)
print("files read", folder.files.height, "dropped by site", folder.dropped_site)
metric = registry.load(registry.resolve(None)).metrics[gdc.STRAND_METRIC]
stats = gdc.strand(gdc.read_counts(folder, spec.gdc).stats, metric)
print("calls", dict(Counter(stats["strand_call"])))
share = stats["first_strand_share"]
print(f"first-strand share min {share.min():.4f} max {share.max():.4f} median {share.median():.4f}")
