"""Package P, E.4: render the page script from two template folders and compare the bytes.

Run: python render_page_js.py <old templates dir> <new templates dir>
The environment's options are builder._ENV's (trim_blocks, lstrip_blocks).
"""

import difflib
import hashlib
import sys

from jinja2 import Environment, FileSystemLoader


def render(folder: str, live: bool) -> str:
    """The page script as base.html.j2 includes it."""
    env = Environment(loader=FileSystemLoader(folder), trim_blocks=True, lstrip_blocks=True)
    return env.get_template("_page.js.j2").render(live=live)


for live in (False, True):
    old, new = render(sys.argv[1], live), render(sys.argv[2], live)
    same = old == new
    digest = hashlib.sha256(new.encode()).hexdigest()[:16]
    print(f"live={live}: old {len(old):,} B, new {len(new):,} B, identical {same}, sha {digest}")
    if not same:
        diff = difflib.unified_diff(old.splitlines(), new.splitlines(), lineterm="", n=1)
        print("\n".join(list(diff)[:60]))
