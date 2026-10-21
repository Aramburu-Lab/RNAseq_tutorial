# ruff: noqa: D100, SIM115 - a proof record, run once and kept as run
import sys
import time

from ctx.report.builder import plotly_js
from ctx.report.figures import pack_json

t = time.time()
js = plotly_js()
p = pack_json(js)
print("pack s", round(time.time() - t, 2), len(js.encode()), len(p))
open(sys.argv[1], "w").write(p)
