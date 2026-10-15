"""Package M9, H.2: the one HPA request that names the lines of HPA's lung cancer group.

It is the request ``ctx genecard`` makes for NKX2-1's cell-line row (``run._lines``), sent
through M9's own cache under M9's own key, so the online run reads it from the cache and
never sends it twice. It prints every line of the group, which the dataset files then name.
"""

from pathlib import Path

from ctx.genecard import hpa
from ctx.genecard.cache import Cache

ROOT = Path.home() / "ctx_lung"
GENE, GROUP = "ENSG00000136352", "cell_RNA_lung_cancer"  # NKX2-1

cache = Cache(ROOT)
answer = cache.get("hpa", hpa.RELEASE, GENE, f"lines_{GROUP}",
                   hpa.search_url(hpa.RELEASE, GENE, ("g", "eg", GROUP)))  # fmt: skip
print(answer.url, answer.retrieved, answer.error or "ok")
if answer.body is not None:
    lines = hpa.parse_lines(answer.body, GENE, GROUP)
    print(f"{len(lines.values)} lines in {GROUP}, {len(answer.body)} B")
    for name, value in lines.values:
        print(f"{name}\t{value}")
