"""R1: the guide's TCGA part, each line held to fetch_gdc.sbatch's Usage block or the CLI.

The laptop blocks run as pasted, behind test_cluster_scripts' stubs, against test_fetch_gdc's
made-up API, and gather's reader then takes the folders. No project ID is typed here. Each test
names the mutation it catches.
"""

import shlex
from pathlib import Path

from gdc_fixture import PROJECT, RELEASE
from test_cluster_scripts import _run, _stub, _usage
from test_fetch_gdc import CDR_NAME, api, site  # noqa: F401 - pytest fixtures, used by name
from test_guide import CLONE_PLACEHOLDER, GIT, SOURCE, SOURCE_PLACEHOLDER, codes, parse

from ctx.dataset import load
from ctx.gather import gdc
from ctx.report import guide, guide_gdc, text

ROOT = Path(__file__).resolve().parents[1]
FOUR = ROOT / "conf" / "datasets" / "Non-cancerous-blood-four-runs.toml"
# A pixi that installs nothing: the clone's fetch environment is the site fixture's stubs.
PIXI = '#!/bin/sh\necho "$@" >> "$STUB_LOG/pixi.calls"\n'


def _part(spec):
    """The TCGA part of a dataset without [gdc]: the guide's last group, after Part 3 and the
    mass-spec part (Part 4), so their numbers never meet."""
    part = guide.blocks(spec, spec.path)[-1]
    assert part.title == text.fill("guide", "gdc.group", n=5)
    return part


def _gdc_spec(tmp_path):
    """The GDC template copied with the made-up project and a slug of its own."""
    raw = load(guide_gdc.TEMPLATE)
    path = tmp_path / "conf" / "datasets" / "fake.toml"
    path.parent.mkdir(parents=True)
    body = guide_gdc.TEMPLATE.read_text().replace(raw.gdc.project, PROJECT)
    body = body.replace(raw.gdc.release, RELEASE)
    path.write_text(body.replace(f'slug = "{raw.slug}"', 'slug = "fake"'))
    return load(path)


def test_usage_lines_are_the_scripts():
    """Each launcher line of the TCGA part is a line of fetch_gdc.sbatch's Usage block, and each
    block holds exactly its lines, placeholders kept for a dataset without [gdc]. Catches a page
    line edited alone, and a Usage line edited alone."""
    blocks = codes([_part(load(FOUR))])
    usage = set(_usage(guide_gdc.SCRIPT))
    for key, (script, lines) in guide_gdc.USAGE.items():
        assert script == guide_gdc.SCRIPT and set(lines) <= usage, set(lines) - usage
        assert blocks[text.word("guide", f"gdc.{key}")].text == "\n".join(lines) + "\n"
    every = "".join(c.text for c in blocks.values())
    assert "<project>" in every and "<slug>" in every


def test_a_gdc_dataset_gets_its_part_filled(tmp_path):
    """A dataset file with [gdc] gets the TCGA part alone, every placeholder filled with its
    project and slug, and its own file named; one without gets every part, TCGA last. Catches
    the fetchngs parts shown for a TCGA project, and a placeholder left in a filled block."""
    spec = _gdc_spec(tmp_path)
    (only,) = guide.blocks(spec, spec.path)
    every = "".join(c.text for c in codes([only]).values())
    assert "<project>" not in every and "<slug>" not in every
    assert f' {PROJECT} "$RAW/fake"' in every
    groups = guide.blocks(load(FOUR), FOUR)
    # find, 1, 1b, 1c (counts, R2), 2, 3, 4 (mass spec, R3), 5 (TCGA)
    assert groups[-1] == _part(load(FOUR)) and len(groups) == 8


def test_ctx_lines_parse_and_name_the_template(tmp_path):
    """Every ctx line of the TCGA part is one the CLI takes, the run reads the three folders the
    script writes, and the dataset file starts as the GDC template without overwriting one.
    Catches an option renamed in the CLI alone, a folder the run does not name, and a copy
    that overwrites a student's edits."""
    spec = _gdc_spec(tmp_path)
    blocks = codes([guide_gdc.group(spec, 1)])
    every = [line for c in blocks.values() for line in c.text.splitlines()]
    lines = [line for line in every if " ctx --root " in line]
    parsed = [parse(line) for line in lines]
    assert [c.info_name for c in parsed] == ["setup", "run", "run"]
    run = parsed[-1].params
    raw = Path("~/ctx_fake/data/raw/fake").expanduser()
    assert [Path(p).expanduser() for p in run["gdc_dir"]] == [raw]
    assert Path(run["mirna_dir"]).expanduser() == Path(f"{raw}_mirna")
    assert Path(run["isoform_dir"]).expanduser() == Path(f"{raw}_isoform")
    assert parsed[1].params["check"] and not run["check"]
    start = blocks[text.word("guide", "gdc.start")].text
    template = guide_gdc.TEMPLATE.relative_to(ROOT).as_posix()
    assert f"cp -n {template} conf/datasets/fake.toml" in start
    assert (ROOT / "scripts" / "bash" / guide_gdc.SCRIPT).is_file()


def test_size_is_read_before_any_download():
    """Each route's metadata run and its download sit in separate blocks, the metadata first.
    Catches a block that pastes the size run and the download together."""
    blocks = codes([_part(load(FOUR))])
    for size, download in (("laptop", "laptop_download"), ("cluster", "cluster_download")):
        first = blocks[text.word("guide", f"gdc.{size}")].text
        second = blocks[text.word("guide", f"gdc.{download}")].text
        assert "DATA=0 " in first and "fetch_gdc.sbatch" in first
        assert "DATA=0" not in second and "read it, then download" in second
    titles = list(blocks)
    assert titles.index(text.word("guide", "gdc.laptop")) < titles.index(
        text.word("guide", "gdc.laptop_download"))  # fmt: skip


def test_laptop_blocks_fetch_what_gather_reads(site, api, tmp_path):  # noqa: F811 - fixtures
    """The laptop blocks, pasted in page order into one shell from a fresh folder, clone, size
    and download the three folders where the run block names them, and gather's reader takes
    them. Catches a block that only works after an unstated step, and a folder named apart
    from where the script writes it."""
    root, env = site
    _, url = api
    _stub(root / "stubs" / "git", GIT)
    _stub(root / "stubs" / "pixi", PIXI)
    spec = _gdc_spec(tmp_path)
    blocks = codes([guide_gdc.group(spec, 1)])
    script = "".join(blocks[text.word("guide", f"gdc.{k}")].text
                     for k in ("laptop", "laptop_download"))  # fmt: skip
    script = script.replace(CLONE_PLACEHOLDER, shlex.quote(env["CTX_REPO"]))
    (tmp_path / "paste.sh").write_text(script.replace(SOURCE_PLACEHOLDER, SOURCE))
    home = tmp_path / "home"
    home.mkdir()
    fresh = tmp_path / "fresh"
    fresh.mkdir()
    done = _run([tmp_path / "paste.sh"], fresh, env, HOME=str(home), GDC_API=url,
                GDC_WEB=f"{url}/web", CDR="0")  # fmt: skip
    assert done.returncode == 0, done.stdout[-3000:] + done.stderr[-2000:]
    assert "PASTE bytes of GDC files, star mirna isoform:" in done.stdout
    assert done.stdout.rstrip().endswith("every check passed")
    raw = home / "ctx_fake" / "data" / "raw" / "fake"
    star = gdc.read_folder(raw, spec.gdc)
    gdc.verify(star)
    assert Path(f"{raw}_mirna").is_dir() and Path(f"{raw}_isoform").is_dir()
    assert not (raw / CDR_NAME).exists()
