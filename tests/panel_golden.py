"""A Reproduce panel's spec fixed in every field, and the hashes of the text main wrote for it.

Y.1 moved the steps the tabs share out of each tab; what a reader copies and downloads must not
change. ``GOLDEN`` holds the SHA-256 of ``panel.body`` and ``panel.script`` for each tab and
each place ctx can come from, written by main's code at bbb0252 on this spec:
``python tests/panel_golden.py`` with main's src first on PYTHONPATH prints them.
"""

import hashlib
from dataclasses import dataclass

from ctx.figures import panel
from ctx.figures.spec import FigureSpec, Input

COMMIT = "0123456789abcdef0123456789abcdef01234567"
WHERE = {"git": f"git {COMMIT}", "dirty": f"git {COMMIT} (uncommitted changes)",
         "wheel": "installed wheel"}  # fmt: skip


@dataclass(frozen=True)
class Params:
    """Three parameters whose comments the content file holds."""

    dataset: str
    contrast: str
    alpha: float


def spec(where: str) -> FigureSpec:
    """The fixed spec, drawn from ``WHERE[where]``; an input under each base."""
    return FigureSpec(
        figure_id="desynth_m5_volcano", dataset="desynth", module="de",
        params=Params("desynth", "a_case", 0.05),
        inputs=(Input("results", "2026-10-20_de_desynth_a_case_results_v1.parquet",
                      "outputs/05_de/2026-10-20_de_desynth_a_case_results_v1.parquet", "root",
                      "ab" * 32),
                Input("dataset_file", "desynth.toml", "conf/datasets/desynth.toml", "repo",
                      "cd" * 32),
                Input("qc_registry", "qc_metrics.toml", "~/shared files/qc_metrics.toml", "path",
                      "ef" * 32)),
        plotting="ctx.de.figures:volcano", plotting_file="src/ctx/de/figures.py",
        stem="2026-10-20_de_desynth_volcano_v1", data_sha256="12" * 32, root="~/ctx_demo",
        commit=WHERE[where], version="0.9.0", upstream="ctx --root ~/ctx_demo de --dataset desynth",
    )  # fmt: skip


def hashes() -> dict[str, str]:
    """SHA-256 of every tab's pasted block and script, by ``<where>/<tab>/<body|script>``."""
    found = {}
    for where in WHERE:
        for tab in panel.TABS:
            for kind, write in (("body", panel.body), ("script", panel.script)):
                found[f"{where}/{tab}/{kind}"] = hashlib.sha256(
                    write(spec(where), tab).encode()
                ).hexdigest()
    return found


GOLDEN = {
    "git/docker/body": "d715a025d7b7e45c32c55be6179de27ebad982e9d031e90742e2ab7b38a388dd",
    "git/docker/script": "d35fb3d4d042ea1ed28ebbf23ec702a5f93378a5c1ef39df22f75e439ccd31c4",
    "git/apptainer/body": "e8ed5ecaed548b71a6111e303506a7e79624fae7f83d07f7cbb74ee4038dbb7c",
    "git/apptainer/script": "e6eb955e3ad7bf51535652a7c3469e8cde1b9497edfd7dcb398fb2c91a8db187",
    "git/pixi/body": "a5108d9f9576240cfd43ce52802f728903dc31f6317b459af5051842092ef42e",
    "git/pixi/script": "981ae6fb02dcb1144eacfda7463daff8cc88c0a23d5ca8cf6cd753b44d0eb209",
    "dirty/docker/body": "d715a025d7b7e45c32c55be6179de27ebad982e9d031e90742e2ab7b38a388dd",
    "dirty/docker/script": "d35fb3d4d042ea1ed28ebbf23ec702a5f93378a5c1ef39df22f75e439ccd31c4",
    "dirty/apptainer/body": "e8ed5ecaed548b71a6111e303506a7e79624fae7f83d07f7cbb74ee4038dbb7c",
    "dirty/apptainer/script": "e6eb955e3ad7bf51535652a7c3469e8cde1b9497edfd7dcb398fb2c91a8db187",
    "dirty/pixi/body": "efc371318b3a0bf87f8957856179fd24f8a094bc8b1eca6419f9515aa65eb360",
    "dirty/pixi/script": "f33a01be9a63fd4e53703eb80ab222a206d07847eacff41d72fc957cdb9bfe07",
    "wheel/docker/body": "e9fcf32a4ff78ef7fc84efeeaa42de978c77047c37ebd390202ae0f5cf256181",
    "wheel/docker/script": "50c22a73615d14d30c5f6a4a3bd5046df6bd867e830b83c55d5dc0a858798866",
    "wheel/apptainer/body": "64f4143793e859a3f1de8728485cd745fa76af27cb831c87d1b0d80c5d5d2c3d",
    "wheel/apptainer/script": "125e3f5327b23bdfe797238b05f7368f17eef3102421e2597f68853fe8e1564c",
    "wheel/pixi/body": "21f6fc25fb1cda5380c224c7b707d656f1c5e8dc1c366d57d41c36afa3ab8e72",
    "wheel/pixi/script": "9f7aed34729f78cbe97cbfe9f5e3dd1db1f8a23d36800e9396a39cf0235756ad",
}

if __name__ == "__main__":
    panel._lock = lambda: None  # the lock of whichever checkout runs this has no say
    for key, value in hashes().items():
        print(f'    "{key}": "{value}",')
