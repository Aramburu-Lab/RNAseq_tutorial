"""Package P, E.4: cut _page.js.j2 into partials by line range, the text moved unchanged.

Each partial is its lines inside raw tags; the page script includes it where the lines stood.
Under the builder's trim_blocks the newline after a comment, an endraw or an include is
dropped and the one after a raw tag is kept, so each partial opens with its raw tag on its
first code line, and the page script reopens raw on a line whose newline is the blank line
the cut removed. The rendered script then keeps every byte in order (render_page_js.py).
Run once, from the worktree.
"""

from pathlib import Path

T = Path("src/ctx/report/templates")
ORIGINAL_LINES = 1100
lines = (T / "_page.js.j2").read_text().splitlines(keepends=True)

# name: (first line, last line, 1-based inclusive, what it holds), in file order.
PARTS = {
    "_page_theme.js.j2": (22, 62, "The page's templates, the theme toggle and the tooltips."),
    "_page_fit.js.j2": (322, 473, "Figure settling: each plot follows its box and swaps layouts."),
    "_page_export.js.j2": (474, 528, "A figure's PDF and SVG buttons for a moved view; its rows."),
    "_page_brush.js.j2": (529, 601, "Controls: the linked brush and M8's pair select."),
    "_page_pick.js.j2": (602, 823, "Controls: M5's gene picker (GS.4)."),
    "_page_controls.js.j2": (824, 1098, "Controls: live ORA, views, filters and live calls."),
}
assert lines[12] == "{% raw %}\n" and lines[-1] == "{% endraw %}\n"
assert lines[20] == "  }\n" and lines[62] == "\n" and lines[319:321] == ["\n", "\n"]
assert lines[1097] == "  }\n" and lines[1098] == "})();\n" and len(lines) == ORIGINAL_LINES

for name, (first, last, what) in PARTS.items():
    body = "".join(lines[first - 1 : last])
    (T / name).write_text(
        f"{{# {what} Included by _page.js.j2, inside its one function. #}}\n"
        f"{{% raw %}}{body}{{% endraw %}}\n"
    )

head = "".join(lines[:21])
middle = "".join(lines[63:321])  # line 63, blank, is the newline the reopening raw tag keeps
out = (
    head
    + "{% endraw %}\n"
    + '{% include "_page_theme.js.j2" %}\n'
    + "{% raw %}\n"
    + middle
    + "{% endraw %}\n"
    + "".join(f'{{% include "{name}" %}}\n' for name in list(PARTS)[1:])
    + "{% raw %}})();\n{% endraw %}\n"
)
(T / "_page.js.j2").write_text(out)
for name in [*PARTS, "_page.js.j2"]:
    print(name, len((T / name).read_text().splitlines()))
