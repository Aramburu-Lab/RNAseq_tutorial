"""``python -m ctx``: the CLI, as `ctx run` starts each module."""

from ctx.cli import app

app(prog_name="ctx")
