# /// codebase-meta
# exempt: >-
#   build tooling; the figures image's check, run by containers.yml, in a scratch folder
# ///
"""Prove the figures image: run one Reproduce panel's Docker tab and compare the redraw.

    python3 containers/prove_figures.py --image <ref> --work <empty folder> [--root <ctx root>]

Standard library only, so it runs on a bare CI runner. Without ``--root`` it first writes the
fixture root inside the image (containers/figures_root.py), so the module's own figure and the
redraw come from one environment. It then takes the volcano's Docker tab from M5's page, as a
reader downloads it, runs it in a clean folder with ``CTX_FIGURES_IMAGE`` naming the image, and
compares the redrawn PDF, SVG and CSV files with the module's byte for byte. It prints each
file's result and the font names each PDF embeds and each SVG names; it exits 1 on any
difference, or on a file missing from either side.
"""

import argparse
import base64
import html
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PLATFORM = "linux/amd64"
PANEL = re.compile(
    r'<code id="ctx-repro-\d+-docker">(.*?)</code>'
    r'.*?href="data:text/x-shellscript;base64,([^"]+)"',
    re.S,
)
CHECK = re.compile(r"^printf '%s  %s\\n' [0-9a-f]{64} \?{4}-\?{2}-\?{2}_(\S+)\.csv", re.M)
PDF_FONT = re.compile(rb"/BaseFont\s*/([^\s/<>\[\]()]+)")
SVG_FONT = re.compile(r"font-family[:=]\s*([^;\"]+)|font:\s*[^;\"]*?px\s+([^;\"]+)")
COMPARED = (".pdf", ".svg", ".csv")


def _run(command: list[str], **kwargs: object) -> None:
    """Run a command, echoed; exit with its code when it fails."""
    print("+", " ".join(command), flush=True)
    done = subprocess.run(command, check=False, **kwargs)  # type: ignore[call-overload]
    if done.returncode != 0:
        sys.exit(done.returncode)


def build_root(image: str, root: Path) -> None:
    """Write the fixture root inside the image, every path mounted where the host has it."""
    root.parent.mkdir(parents=True, exist_ok=True)
    user = f"{os.getuid()}:{os.getgid()}"
    _run(["docker", "run", "--rm", "--platform", PLATFORM, "-u", user, "-e", "HOME=/tmp",
          "-v", f"{REPO}:{REPO}:ro", "-v", f"{root.parent}:{root.parent}", "-w", str(root.parent),
          image, "python", str(REPO / "containers" / "figures_root.py"), str(root)])  # fmt: skip


def volcano_tab(root: Path) -> tuple[str, str]:
    """The volcano's Docker tab on M5's page: the downloaded script and the figure's stem."""
    (page,) = (root / "outputs" / "05_de").glob("*_overview_v1.html")
    for shown, encoded in PANEL.findall(page.read_text()):
        pasted = html.unescape(shown)
        if "_m5_volcano --params" in pasted:
            match = CHECK.search(pasted)
            if match is None:
                sys.exit("the volcano's Docker tab has no data check line")
            return base64.b64decode(encoded).decode(), match.group(1)
    sys.exit(f"{page}: no Docker tab draws a volcano")


def fonts(path: Path) -> list[str]:
    """The font names a PDF embeds, or the font families an SVG names, sorted and distinct."""
    if path.suffix == ".pdf":
        return sorted({m.decode() for m in PDF_FONT.findall(path.read_bytes())})
    return sorted({(a or b).strip() for a, b in SVG_FONT.findall(path.read_text())})


def compare(module: Path, redrawn: Path, stem: str) -> bool:
    """Print each compared file's result and fonts; True when every pair is byte-equal."""
    # The day each was drawn leads both names; the rest must match.
    own = {p.name[10:]: p for p in module.glob(f"*_{stem}*") if p.suffix in COMPARED}
    new = {p.name[10:]: p for p in redrawn.glob(f"*_{stem}*") if p.suffix in COMPARED}
    equal = bool(own) and own.keys() == new.keys()
    for name in sorted(own.keys() | new.keys()):
        a, b = own.get(name), new.get(name)
        if a is None or b is None:
            print(f"MISSING  {name}: module {a is not None}, redraw {b is not None}")
            continue
        same = a.read_bytes() == b.read_bytes()
        equal &= same
        sizes = f"{a.stat().st_size} B and {b.stat().st_size} B"
        print(f"{'EQUAL' if same else 'DIFFER'}  {name}: {sizes}")
        if name.endswith((".pdf", ".svg")):
            print(f"         fonts, module: {fonts(a)}; redraw: {fonts(b)}")
    return equal


def main() -> None:
    """Build or take the root, run the Docker tab, compare; exit 1 on any difference."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--image", required=True, help="the figures image to prove")
    parser.add_argument("--work", required=True, type=Path, help="an empty scratch folder")
    parser.add_argument("--root", type=Path, help="an existing ctx root; else one is built")
    args = parser.parse_args()
    work = args.work.resolve()
    if work.exists() and any(work.iterdir()):
        sys.exit(f"{work}: not empty")
    root = args.root.resolve() if args.root else work / "root"
    if args.root is None:
        build_root(args.image, root)
    script, stem = volcano_tab(root)
    reader = work / "reader"
    reader.mkdir(parents=True)
    (reader / "reproduce.sh").write_text(script)
    env = {**os.environ, "CTX_FIGURES_IMAGE": args.image}
    _run(["bash", "reproduce.sh"], cwd=reader, env=env)
    (redrawn,) = reader.glob("repro_*")
    if not compare(root / "outputs" / "05_de", redrawn, stem):
        sys.exit(1)
    print("every compared file is byte-equal")


if __name__ == "__main__":
    main()
