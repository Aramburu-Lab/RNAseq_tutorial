"""wo13-r1 proof: one block of the guide's TCGA part, as the default guide page prints it.

Only the placeholders a student fills are replaced: the project, the slug, the clone and the
course repository. Usage: make_paste.py <block key> [<block key> ...] > paste.sh
"""

import sys
from pathlib import Path

from ctx.dataset import load
from ctx.report import guide_gdc, text

WT = Path("/Users/rytis/Scripts/Pipeline/RNAseq_tutorial/.claude/worktrees/agent-a00f69bc2d903ff45")
FILLS = {
    "<project>": "TCGA-HNSC",
    "<slug>": "tcga_hnsc",
    "<your clone of the course repository>": str(Path.home() / "ctx_r1_clone"),
    "<course repository>": str(WT),
}


def main() -> None:
    """Print the named blocks, filled, in the order given."""
    group = guide_gdc.group(load(WT / "conf" / "datasets" / "four_runs.toml"), 4)
    blocks = {i.title: i.text for i in group.items if hasattr(i, "text")}
    for key in sys.argv[1:]:
        body = blocks[text.word("guide", f"gdc.{key}")]
        for placeholder, value in FILLS.items():
            body = body.replace(placeholder, value)
        sys.stdout.write(body)


if __name__ == "__main__":
    main()
