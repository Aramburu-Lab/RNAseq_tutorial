"""U.5 and U.9: the Code result kind, and the custom-data guide (`ctx guide`).

T1: each `ctx` line parses with the CLI's own click parser. T2: each launcher line is its
script's Usage line, and the laptop blocks run in page order from a fresh folder behind the
stubs of test_cluster_scripts. T3: the copy's includes are COPIED_FILES. T4: no number is typed
into the guide's text. No accession is typed here. Each test names the mutation it catches.
"""

import hashlib
import html
import re
import shlex
import tomllib
from pathlib import Path

import pytest
import typer
from test_cluster_scripts import BASH, FETCHNGS, _ids, _recorded, _run, _stub, _usage
from test_cluster_scripts import site as _cluster_site
from test_sites import DOCKER, UNAME

from ctx.cli import app
from ctx.dataset import load
from ctx.gather.fetchngs import FETCHNGS_VERSION
from ctx.gather.rnaseq import COPIED_FILES, RNASEQ_VERSION
from ctx.report import builder, course, guide, prose, text
from ctx.report.builder import (
    Code,
    Group,
    ModulePage,
    ParameterCard,
    Reproduction,
    ResultTable,
)


@pytest.fixture
def site(tmp_path):
    """test_cluster_scripts' site: the clone, stubs and job environment."""
    return _cluster_site.__wrapped__(tmp_path)


ROOT = Path(__file__).resolve().parents[1]
FOUR = ROOT / "conf" / "datasets" / "four_runs.toml"
CLONE_PLACEHOLDER = "<your clone of the course repository>"
# curl -fL --remote-name-all URL...: each file named by its URL, MD5SUMS listing the others.
CURL = """#!/bin/sh
echo "$@" >> "$STUB_LOG/curl.calls"
for a in "$@"; do
    case "$a" in http*) n=${a##*/}; [ "$n" = MD5SUMS ] || echo "$n" > "$n" ;; esac
done
for f in *.gz; do printf '%s  %s\\n' "$(md5sum < "$f" | cut -d' ' -f1)" "$f"; done > MD5SUMS
"""


def codes(groups) -> dict[str, Code]:
    """Every code block of the page, by its title."""
    return {i.title: i for g in groups for i in g.items if isinstance(i, Code)}


@pytest.fixture(scope="module")
def page_blocks():
    spec = load(FOUR)
    return spec, guide.blocks(spec, FOUR)


def _ctx_lines(blocks) -> list[str]:
    return [line for c in codes(blocks).values() for line in c.text.splitlines()
            if " ctx --root " in line]  # fmt: skip


def parse(line: str) -> object:
    """A `ctx` line through the CLI's own parser, as typing it would: unknown or missing
    options raise."""
    args = shlex.split(line[line.index("ctx --root") :])[1:]
    group = typer.main.get_command(app)
    at = next(i for i, a in enumerate(args) if a in group.commands)
    head = group.make_context("ctx", args[: at + 1])
    return group.commands[args[at]].make_context(args[at], args[at + 1 :], parent=head)


def test_ctx_lines_parse(page_blocks):
    """T1: every ctx line on the page, and the README's, is one the CLI takes. Catches an
    option renamed in the CLI and not on the page, and a required option left off."""
    lines = _ctx_lines(page_blocks[1]) + [
        line for line in (ROOT / "README.md").read_text().splitlines() if " ctx --root " in line
    ]
    assert len(lines) >= 8, lines
    parsed = [parse(line) for line in lines]
    assert {c.info_name for c in parsed} >= {"setup", "run", "gather", "guide"}
    # typer bundles its own click, so its errors are matched by their words.
    run = next(line for line in lines if " run --dataset four_runs " in line)
    with pytest.raises(Exception, match="No such option"):
        parse(run.replace("--fetchngs-dir", "--fetchngs-folder"))
    with pytest.raises(Exception, match="Missing parameter: dataset"):
        parse(run.replace("--dataset four_runs ", ""))


def test_usage_lines_are_the_scripts(page_blocks):
    """T2: each launcher line of the page is a line of its script's Usage block, filled.
    Catches a page line edited alone, and a Usage line edited alone."""
    spec, blocks = page_blocks
    shown = codes(blocks)
    for key, (script, lines) in guide.USAGE.items():
        assert set(lines) <= set(_usage(script)), (key, set(lines) - set(_usage(script)))
        (block,) = [c for c in shown.values() if c.title == text.word("guide", _titles()[key])]
        assert block.text == "".join(guide.filled(line, spec) + "\n" for line in lines)
    joined = " ".join(spec.accessions)
    assert f"printf '%s\\n' {joined} > ids.txt" in shown[text.word("guide", "part1.laptop")].text


def _titles() -> dict[str, str]:
    return {"fetch_laptop": "part1.laptop", "fetch_laptop_download": "part1.laptop_download",
            "setup": "part1.setup", "fetch_cluster": "part1.cluster",
            "fetch_cluster_size": "part1.cluster_size",
            "fetch_cluster_download": "part1.cluster_download", "count_laptop": "part2.laptop",
            "count_cluster": "part2.cluster"}  # fmt: skip


LAUNCH = re.compile(
    r"(?:^|\s)(?:CTX_SITE=laptop |SKIP_DOWNLOAD=1 |DRY_RUN=1 )*"
    r"(?:bash|sbatch -A \S+ \S+) \S*fetchngs\.sbatch",
    re.M,
)


def test_no_block_downloads_before_the_size_is_read(page_blocks):
    """Each route's metadata run and its download are in separate Copy blocks, the size read
    between them, so pasting one block never downloads unseen. Catches the laptop's two launches
    in one block (U.10 pasted it and downloaded at once), and the Kebnekaise grep and download
    in one."""
    for c in codes(page_blocks[1]).values():
        launches = [m.group(0) for m in LAUNCH.finditer(c.text)]
        downloads = [x for x in launches if "SKIP_DOWNLOAD" not in x and "DRY_RUN" not in x]
        assert not (downloads and len(launches) > len(downloads)), c.title
        assert not (downloads and "bytes of FASTQ at ENA" in c.text), c.title
    titles = list(codes(page_blocks[1]))
    order = [text.word("guide", k) for k in ("part1.laptop", "part1.laptop_download",
                                             "part1.cluster", "part1.cluster_size",
                                             "part1.cluster_download")]  # fmt: skip
    assert [t for t in titles if t in order] == order


def test_laptop_gtf_is_where_setup_finds_it(page_blocks):
    """Part 2's laptop block saves the GTF at the dataset's annotation_gtf dest under its ctx
    root, so ctx setup finds it in place, and the page lists every file setup downloads with its
    size. Catches the GTF fetched twice (U.10: 124,650,284 bytes again) and a download the page
    does not name."""
    spec, blocks = page_blocks
    lines = codes(blocks)[text.word("guide", "part2.laptop")].text.splitlines()
    (refs,) = [line.split("=", 1)[1].strip('"') for line in lines if line.startswith("REFS=")]
    count = lines[-1]
    gtf = re.search(r'"\$REFS/(\S+\.gtf\.gz)"', count).group(1)
    root = f"$HOME/ctx_{spec.slug}"
    assert f"{refs}/{gtf}" == f"{root}/{spec.resource('annotation_gtf').dest}"
    assert f" --root ~/ctx_{spec.slug} " in codes(blocks)[text.word("guide", "part3.laptop")].text
    (table,) = [i for g in blocks for i in g.items if isinstance(i, ResultTable)]
    assert [(r[0], r[1]) for r in table.rows] == [(r.name, f"{r.bytes:,}") for r in spec.resources]
    found = {r[0]: r[3] for r in table.rows}
    assert found["annotation_gtf"] == text.word("guide", "part3.found")
    assert found["msigdb_hallmark"] == text.word("guide", "part3.fetched")


def test_laptop_needs_from_one_place(page_blocks):
    """Part 2 states Docker's memory floor and the run folder's size from guide.NEEDS, the
    measured values. Catches a value typed into the page text instead."""
    note = next(i for g in page_blocks[1] for i in g.items if isinstance(i, str)
                and str(guide.NEEDS["index_gib"]) in i)  # fmt: skip
    assert note == text.fill("guide", "part2.needs", **guide.NEEDS)
    assert text.placeholders(text.raw("guide", "part2.needs")) == set(guide.NEEDS)


def test_reference_is_the_datasets(page_blocks):
    """The GTF Part 2 downloads is the dataset's annotation_gtf, Step-1's GENCODE release.
    Catches the Usage moved to another release alone."""
    spec, _ = page_blocks
    (base,) = [line[2:] for line in guide.USAGE["count_laptop"][1] if line.startswith("B=")]
    (curl,) = [line for line in guide.USAGE["count_laptop"][1] if line.startswith("curl ")]
    (gtf,) = re.findall(r'"\$B/(gencode[^"]+)"', curl)
    assert spec.resource("annotation_gtf").urls == (f"{base}/{gtf}",)


def test_laptop_blocks_run_in_page_order(site, tmp_path):
    """T2: the laptop blocks of Parts 1 and 2, pasted in order into one shell from a fresh
    folder, download, check the references and count, behind stubs. Catches a block that only
    works after an unstated step, and the pin, profile or parser lost on the way."""
    root, env = site
    _stub(root / "stubs" / "uname", UNAME)
    _stub(root / "stubs" / "docker", DOCKER)
    _stub(root / "stubs" / "curl", CURL)
    spec_path = tmp_path / "fixture.toml"
    ids = ", ".join(f'"{i}"' for i in _ids("download"))
    spec_path.write_text(FOUR.read_text().replace(
        tomllib.loads(FOUR.read_text())["dataset"]["accessions"].__repr__().replace("'", '"'),
        f"[{ids}]"))  # fmt: skip
    spec = load(spec_path)
    assert spec.accessions == tuple(_ids("download"))
    shown = codes(guide.blocks(spec, spec_path))
    part1 = shown[text.word("guide", "part1.laptop")].text
    part1 += shown[text.word("guide", "part1.laptop_download")].text
    part2 = shown[text.word("guide", "part2.laptop")].text
    clone = shlex.quote(env["CTX_REPO"])
    # Between the parts, the harness keeps what Part 1's last launch recorded.
    keep = 'cp "$STUB_LOG/nextflow.argv" "$STUB_LOG/part1.argv"\n'
    keep += 'cp "$STUB_LOG/nextflow.env" "$STUB_LOG/part1.env"\n'
    script = part1.replace(CLONE_PLACEHOLDER, clone) + keep + part2
    (tmp_path / "paste.sh").write_text("set -e\n" + script)
    home = tmp_path / "home"
    home.mkdir()
    fresh = tmp_path / "fresh"
    fresh.mkdir()
    done = _run([tmp_path / "paste.sh"], fresh, env, HOME=str(home), STUB_ARCH="arm64",
                STUB_FIXTURE=str(FETCHNGS / "download"))  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    log = root / "log"
    first = (log / "part1.argv").read_text().splitlines()
    assert first[:4] == ["run", "nf-core/fetchngs", "-r", FETCHNGS_VERSION]
    assert first[first.index("-profile") + 1] == "docker,arm"
    assert "NXF_SYNTAX_PARSER=v1" in (log / "part1.env").read_text().splitlines()
    argv, nf_env, _ = _recorded(root)
    assert argv[:4] == ["run", "nf-core/rnaseq", "-r", RNASEQ_VERSION]
    assert argv[argv.index("-profile") + 1] == "docker,arm64"
    assert nf_env["NXF_SYNTAX_PARSER"] == "v2"
    assert done.stdout.count("MD5 equal to") == 2, done.stdout
    run = home / "ctx_runs" / spec.slug
    assert (run / "nfcore_fetchngs_meta" / "results.complete").is_file()
    assert (run / "nfcore_rnaseq" / "results" / "salmon").is_dir()


def test_copy_includes_are_copied_files(page_blocks):
    """T3: the copy from Kebnekaise includes exactly COPIED_FILES, in their order, after the
    folder rule and before the exclusion. Catches a pattern dropped or added on the page."""
    block = codes(page_blocks[1])[text.word("guide", "part3.copy")].text
    rsync = block.splitlines()[1]
    assert re.findall(r"--include='([^']+)'", rsync) == ["*/", *COPIED_FILES]
    assert rsync.index("--exclude='*'") > rsync.rindex("--include=")


def test_no_number_typed_in_the_guide():
    """T4: every number the guide's words carry comes from a placeholder or sits in code.
    Catches a count or version typed into the page text."""
    raw = tomllib.loads((ROOT / "src" / "ctx" / "report" / "content" / "guide.toml").read_text())

    def strings(node):
        if isinstance(node, str):
            yield node
        elif isinstance(node, dict):
            for value in node.values():
                yield from strings(value)
        elif isinstance(node, list):
            for value in node:
                yield from strings(value)

    typed = [s for s in strings(raw) if prose.quantities(s)]
    assert typed == []


def _code_page(*blocks: Code) -> str:
    return builder.render(ModulePage(
        title="Your own data: test", **text.skeleton("guide"),
        parameters=(ParameterCard("--x", "x", "1", "x"),),
        results=(Group("Blocks", blocks),),
        reproduce=Reproduction("ctx guide", 1.0, 1, 1, "laptop", "none", {}),
    ))  # fmt: skip


def test_code_kept_as_written():
    """U.5: a code block prints its text byte for byte, escaped, with that text's SHA-256.
    Catches the text reflowed or stripped, and another block's hash shown."""
    first = Code("one", "a  b\n\tc <d> & 'e'\n\n  f\n", "f.sh")
    second = Code("two", "x\n")
    page = _code_page(first, second)
    shown = re.findall(r'<code id="ctx-code-\d+">(.*?)</code>', page, re.S)
    found = [html.unescape(t) for t in shown]
    assert found == [first.text, second.text]
    hashes = re.findall(r'<p class="sha">[^<]*<code>([0-9a-f]{64})</code>', page)
    assert hashes == [hashlib.sha256(c.text.encode()).hexdigest() for c in (first, second)]
    assert "<code>f.sh</code>" in page


def test_each_copy_button_copies_its_own_block():
    """U.5: each Copy button names its own block's id, unique on the page beside the Reproduce
    command's. Catches two buttons copying one block."""
    page = _code_page(Code("one", "a\n"), Code("two", "b\n"), Code("three", "c\n"))
    ids = re.findall(r'<code id="(ctx-code-\d+)">', page)
    buttons = re.findall(r'data-copy="(ctx-code-\d+)"', page)
    assert ids == buttons and len(set(ids)) == 3


def test_guide_written_and_linked_from_start(tmp_path):
    """`ctx guide` writes the page under outputs/09_report, and START.html links it. Catches
    the guide left off the course home page, and a page named outside the naming rule."""
    from conftest import ctx_run  # noqa: PLC0415

    done = ctx_run("--root", str(tmp_path), "guide", cwd=tmp_path)
    assert done.returncode == 0, done.stderr[-3000:]
    (page,) = (tmp_path / "outputs" / "09_report").glob("*_report_guide_v1.html")
    assert done.stdout.strip() == str(page)
    home = course.page(tmp_path, [(tmp_path, True)])
    assert home.guide == page.relative_to(tmp_path).as_posix()
    start = builder.render(home)
    assert f'<a href="{home.guide}">{text.word("course", "guide")}</a>' in start


def test_scripts_named_on_the_page_exist(page_blocks):
    """Every repository path the page names is a file of the repository. Catches a script or
    the template renamed and the page left behind."""
    every = "\n".join(c.text for c in codes(page_blocks[1]).values())
    path = r"(?:\$CTX_REPO/|^|\s)((?:scripts|conf)/[\w/.]+\.(?:sbatch|sh|py|toml))"
    named = set(re.findall(path, every, re.M))
    assert named and all((ROOT / p).is_file() for p in named), named
    assert all((BASH / f"site_{s}.sh").is_file() for s in ("kebnekaise", "laptop"))
