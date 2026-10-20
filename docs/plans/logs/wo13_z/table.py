"""wo13-z: every rebuilt page's bytes on main's side and on this package's.

The parts of the pages within 1 MB of the cap on either side follow the table.

Run: python table.py. Reads build.py's scratch roots; a page is a module's overview page that
the run wrote (gather tables were copied in, so a group's gather page is the real root's own,
except Step-1's, which gather wrote).
"""

from pathlib import Path

from parts import parts

SCRATCH = Path(__file__).resolve().parents[4] / "outputs" / "scratch_z"
CAP = 15 * 10**6
NEAR = 10**6


def pages(side: str) -> dict[str, Path]:
    """Each written overview page of one side, keyed by its name without the day."""
    found = {}
    for path in sorted((SCRATCH / side).glob("*/outputs/[0-9][0-9]_*/*_overview_v1.html")):
        if path.parent.name == "01_gather" and path.name[:10] != "2026-10-20":
            continue
        found[path.name[11:].removesuffix("_overview_v1.html")] = path
    return found


def main() -> None:
    """Print the size table, then the parts of the pages near the cap."""
    main_pages, z_pages = pages("main"), pages("z")
    print(f"{'page':<34} {'main':>12} {'Z':>12} {'Z minus main':>13}  status on Z")
    for name in sorted(main_pages.keys() | z_pages.keys()):
        before = main_pages[name].stat().st_size if name in main_pages else None
        after = z_pages[name].stat().st_size if name in z_pages else None
        delta = f"{after - before:>+13,}" if before is not None and after is not None else ""
        state = "" if after is None else ("written" if after <= CAP else "OVER")
        print(f"{name:<34} {before or 0:>12,} {after or 0:>12,} {delta}  {state}")
    near = sorted(n for n in main_pages.keys() | z_pages.keys()
                  if any(p[n].stat().st_size >= CAP - NEAR for p in (main_pages, z_pages)
                         if n in p))  # fmt: skip
    for name in near:
        sides = [(s, p[name]) for s, p in (("main", main_pages), ("Z", z_pages)) if name in p]
        found = {s: parts(path.read_text()) for s, path in sides}
        print(f"\n{name}")
        for key in found[sides[0][0]]:
            print(f"  {key:<42}" + "".join(f" {found[s][key]:>12,}" for s, _ in sides))


if __name__ == "__main__":
    main()
