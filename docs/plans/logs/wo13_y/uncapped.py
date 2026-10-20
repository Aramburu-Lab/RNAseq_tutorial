"""wo13-y: one ctx command with the page cap lifted, so a page over 15 MB is written and measured.

Usage: python uncapped.py <ctx arguments>. Prints which ctx it imported. Measurement only: the
cap stays in force in every committed path.
"""

import sys

import ctx
from ctx.cli import app
from ctx.report import builder

print("ctx from", ctx.__file__, flush=True)
builder.MAX_BYTES = 10**9
app(sys.argv[1:], standalone_mode=False)
