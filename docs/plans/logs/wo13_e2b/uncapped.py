"""Run one ctx command with the page cap lifted, to measure what put a page over it.

python uncapped.py <ctx arguments ...>; the page is written whatever its size. Never for a
page that ships: the cap is the shipped rule.
"""

import sys

from ctx.report import builder

builder.MAX_BYTES = 10**9

from ctx.cli import app  # noqa: E402

app(sys.argv[1:], standalone_mode=False)
