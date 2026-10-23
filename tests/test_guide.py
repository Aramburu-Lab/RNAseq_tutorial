"""U.5 and U.9: the Code result kind, and the custom-data guide (`ctx guide`).

T1: each `ctx` line parses with the CLI's own click parser. T2: each launcher line is its
script's Usage line, and the laptop blocks run in page order from a fresh folder behind the
stubs of test_cluster_scripts. T3: the copy's includes are CLUSTER_FILES. T4: no number is typed
into the guide's text. No accession is typed here. Each test names the mutation it catches.
"""

import csv
import gzip
import hashlib
import html
import re
import shlex
import shutil
import tomllib
from pathlib import Path

import pytest
import typer
from test_cluster_scripts import BASH, FETCHNGS, _ids, _recorded, _run, _stub, _usage
from test_cluster_scripts import site as _cluster_site
from test_sites import DOCKER, UNAME

from ctx.cli import app
from ctx.dataset import load
from ctx.gather import own_fastq
from ctx.gather.fetchngs import FETCHNGS_VERSION
from ctx.gather.rnaseq import CLUSTER_FILES, RNASEQ_VERSION
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
FOUR = ROOT / "conf" / "datasets" / "Non-cancerous-blood-four-runs.toml"
CLONE_PLACEHOLDER = "<your clone of the course repository>"
SOURCE_PLACEHOLDER = "<course repository>"
# The clone's source as a student types it; the git stub records the call and clones nothing,
# since the site fixture's clone is already in place.
SOURCE = "https://example.org/course.git"
GIT = """#!/bin/sh
echo "$@" >> "$STUB_LOG/git.calls"
"""
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
            if " ctx --root " in line or " ctx reads " in line]  # fmt: skip


def parse(line: str) -> object:
    """A `ctx` line through the CLI's own parser, as typing it would: unknown or missing
    options raise."""
    args = shlex.split(line[line.index(" ctx ") + 1 :])[1:]
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
    assert {c.info_name for c in parsed} >= {"setup", "run", "gather", "guide", "reads"}
    # typer bundles its own click, so its errors are matched by their words.
    run = next(line for line in lines if " run --dataset Non-cancerous-blood-four-runs " in line)
    with pytest.raises(Exception, match="No such option"):
        parse(run.replace("--fetchngs-dir", "--fetchngs-folder"))
    with pytest.raises(Exception, match="Missing parameter: dataset"):
        parse(run.replace("--dataset Non-cancerous-blood-four-runs ", ""))


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
    title = text.word("guide", "part3.downloads")
    tables = [i for g in blocks for i in g.items if isinstance(i, ResultTable)]
    (table,) = [t for t in tables if t.title == title]
    assert [(r[0], r[1]) for r in table.rows] == [(r.name, f"{r.bytes:,}") for r in spec.resources]
    found = {r[0]: r[3] for r in table.rows}
    assert found["annotation_gtf"] == text.word("guide", "part3.found")
    assert found["msigdb_hallmark"] == text.word("guide", "part3.fetched")


def test_accession_comes_before_the_download(page_blocks):
    """M1a on the page: the accession levels and the class exercise lead, before any launcher,
    one row per level from study to run. Catches the group moved after Part 1, a level dropped,
    and a level's prefixes left out of an archive's column."""
    first, *rest = page_blocks[1]
    assert first.title == text.word("guide", "find.group") and rest
    assert not any(isinstance(i, Code) for i in first.items)
    (levels,) = [i for i in first.items if isinstance(i, ResultTable)]
    assert [r[0] for r in levels.rows] == ["Study", "Sample", "Experiment", "Run"]
    assert [r[1:4] for r in levels.rows][-1] == ("SRR", "ERR", "DRR")
    assert all(all(cell for cell in row) for row in levels.rows)
    assert text.word("guide", "find.exercise") in first.items


def test_laptop_needs_from_one_place(page_blocks):
    """U.18: Part 2 states, before its first block, Docker's memory floor by both readings, the
    disk a laptop needs and the references' size, from guide.NEEDS, the laptop proof's report.
    The numbers are typed here, apart from the code. Catches a value typed into the page text,
    the Kebnekaise 16.8 GiB back, KiB read as kB, and the references left out of the disk."""
    (part2,) = [g for g in page_blocks[1] if g.title == text.fill("guide", "part2.group", n=2)]
    first = next(n for n, i in enumerate(part2.items) if isinstance(i, Code))
    needs, refs = part2.items[:first]
    assert needs == text.fill("guide", "part2.needs", **guide.needs())
    assert refs == text.fill("guide", "part2.refs", release=guide.release(), **guide.needs())
    for key in ("part2.needs", "part2.refs"):
        assert text.placeholders(text.raw("guide", key)) <= set(guide.needs()) | {"release"}
    for shown in ("17.020 GiB by docker stats", "17.5 GB by Nextflow's trace", "60.7 GB",
                  "59,298,436 KiB", "38.4 GB", "61.7 GB free", "2026-10-09", "20 GiB", "13 CPUs",
                  "Apple-silicon"):  # fmt: skip
        assert shown in needs, shown
    assert "970,288,386 B" in refs and "release 50" in refs and "16.8" not in needs


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
    _stub(root / "stubs" / "git", GIT)
    # The marker ends the metadata run's output, whose line before last must be the size.
    cut = "--- metadata run ends ---"
    comment = guide.USAGE["fetch_laptop_download"][1][0]
    assert "line before last is the size" in comment
    part1 = part1.replace(comment, f'echo "{cut}"\n{comment}')
    script = part1.replace(CLONE_PLACEHOLDER, clone).replace(SOURCE_PLACEHOLDER, SOURCE)
    script += keep + part2
    (tmp_path / "paste.sh").write_text("set -e\n" + script)
    home = tmp_path / "home"
    home.mkdir()
    fresh = tmp_path / "fresh"
    fresh.mkdir()
    done = _run([tmp_path / "paste.sh"], fresh, env, HOME=str(home), STUB_ARCH="arm64",
                STUB_FIXTURE=str(FETCHNGS / "download"))  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    log = root / "log"
    assert (log / "git.calls").read_text() == f"clone {SOURCE} {env['CTX_REPO']}\n"
    printed = done.stdout.split(cut)[0].splitlines()
    assert "bytes of FASTQ at ENA" in printed[-2] and printed[-1].startswith("complete: ")
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
    """T3: the copy from Kebnekaise includes exactly CLUSTER_FILES, in their order, after the
    folder rule and before the exclusion. Catches a pattern dropped or added on the page."""
    block = codes(page_blocks[1])[text.word("guide", "part3.copy")].text
    rsync = block.splitlines()[1]
    assert re.findall(r"--include='([^']+)'", rsync) == ["*/", *CLUSTER_FILES]
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
    path = r"(?:\$CTX_REPO/|^|\s)((?:scripts|conf|src)/[\w/.]+\.(?:sbatch|sh|py|toml|csv))"
    named = set(re.findall(path, every, re.M))
    assert named and all((ROOT / p).is_file() for p in named), named
    assert all((BASH / f"site_{s}.sh").is_file() for s in ("kebnekaise", "laptop"))


def _own_blocks(spec) -> tuple[str, str]:
    shown = codes(guide.blocks(spec, spec.path))
    return (shown[text.word("guide", "part1b.laptop")].text,
            shown[text.word("guide", "part1b.cluster")].text)  # fmt: skip


def test_own_fastq_lines_name_one_command(page_blocks):
    """Part 1b's two routes run one function: `ctx reads` on the laptop, the same file on
    Kebnekaise, each line parsed by its own parser to the same options; the table and its
    sample column are the dataset's own table and join. Catches an option renamed in one
    entry point, and a column other than the one curation joins on."""
    spec, _ = page_blocks
    laptop, cluster = _own_blocks(spec)
    (typed,) = [line for line in laptop.splitlines() if " ctx reads " in line]
    (scripted,) = [line for line in cluster.splitlines() if guide.OWN_SCRIPT in line]
    clicked = parse(typed).params
    argv = shlex.split(scripted.replace("$CTX_REPO", "REPO").replace("$RUN", "RUN"))
    assert argv[:2] == ["python3", f"REPO/{guide.OWN_SCRIPT}"]
    scripted_args = vars(own_fastq.parser().parse_args(argv[2:]))
    assert clicked["sample_column"] == scripted_args["sample_column"] == spec.sdrf.join[0]
    assert str(clicked["table"]) == spec.table().relative_to(ROOT).as_posix()
    assert str(scripted_args["table"]) == f"REPO/{clicked['table']}"
    assert str(clicked["fastq_dir"]) == str(scripted_args["fastq_dir"]) == guide.FASTQ_DIR


# `pixi run -e laptop ctx ...` runs this checkout's ctx; any other pixi call does nothing.
PIXI_RUN = """#!/bin/sh
echo "$@" >> "$STUB_LOG/pixi.calls"
if [ "$1" = run ]; then shift 4; PYTHONPATH="{src}" exec "{ctx}" "$@"; fi
"""
REFS_OWN = {"genome.fa.gz": b"genome", "annotation.gtf.gz": b"annotation"}


def _own_files(folder: Path, spec) -> Path:
    """One gzipped single-end file per run of the dataset, named by the run, from chr22_se."""
    data = gzip.compress((ROOT / "tests" / "fixtures" / "chr22_se" / "se_c.fq").read_bytes(),
                         mtime=0)  # fmt: skip
    folder.mkdir(parents=True)
    for run in spec.accessions:
        (folder / f"{run}.fastq.gz").write_bytes(data)
    return folder


def _cut(run: Path) -> list[dict[str, str]]:
    with (run / "nfcore_rnaseq" / "samplesheet.csv").open(newline="") as fh:
        return list(csv.DictReader(fh))


def test_own_fastq_then_part2_laptop(site, tmp_path):
    """Part 1b's laptop block, pasted in a fresh folder, then Part 2's laptop block unchanged:
    nf-core/rnaseq takes the samplesheet ctx reads wrote, each run's file where it is. Catches
    Part 1b leaving RUN or CTX_REPO unset for Part 2, the outdir where rnaseq.sbatch does not
    look, and a samplesheet naming a link instead of the file."""
    root, env = site
    _stub(root / "stubs" / "uname", UNAME)
    _stub(root / "stubs" / "docker", DOCKER)
    _stub(root / "stubs" / "curl", CURL)
    _stub(root / "stubs" / "pixi", PIXI_RUN.format(src=ROOT / "src", ctx=shutil.which("ctx")))
    spec = load(FOUR)
    clone = Path(env["CTX_REPO"])
    (clone / "conf" / "datasets").mkdir(parents=True)
    shutil.copy(spec.table(), clone / "conf" / "datasets")
    files = _own_files(tmp_path / "my_fastq", spec)
    laptop, _ = _own_blocks(spec)
    part2 = codes(guide.blocks(spec, FOUR))[text.word("guide", "part2.laptop")].text
    _stub(root / "stubs" / "git", GIT)
    script = laptop.replace(CLONE_PLACEHOLDER, shlex.quote(env["CTX_REPO"]))
    script = script.replace(SOURCE_PLACEHOLDER, SOURCE)
    script = script.replace(guide.FASTQ_DIR, str(files))
    (tmp_path / "paste.sh").write_text("set -e\n" + script + part2)
    home = tmp_path / "home"
    home.mkdir()
    fresh = tmp_path / "fresh"
    fresh.mkdir()
    done = _run([tmp_path / "paste.sh"], fresh, env, HOME=str(home), STUB_ARCH="arm64")
    assert done.returncode == 0, done.stderr[-3000:]
    argv, _, _ = _recorded(root)
    assert argv[:4] == ["run", "nf-core/rnaseq", "-r", RNASEQ_VERSION]
    run = home / "ctx_runs" / spec.slug
    assert [r["sample"] for r in _cut(run)] == list(spec.accessions)
    assert {r["fastq_1"] for r in _cut(run)} == {str(p.resolve()) for p in files.iterdir()}
    assert {r["strandedness"] for r in _cut(run)} == {"auto"}
    assert (run / own_fastq.METADATA / "results.complete").is_file()


def test_own_fastq_then_part2_kebnekaise(site, tmp_path):
    """Part 1b's Kebnekaise block runs the file on the fetch environment's Python, then
    rnaseq.sbatch, as Part 2's job runs it, takes that samplesheet. Catches the file run on
    another Python, and an outdir the cluster job does not find."""
    root, env = site
    spec = load(FOUR)
    clone = Path(env["CTX_REPO"])
    for relative in (guide.OWN_SCRIPT, spec.table().relative_to(ROOT).as_posix()):
        (clone / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(ROOT / relative, clone / relative)
    store = tmp_path / "store"
    store.mkdir()
    files = _own_files(tmp_path / "my_fastq", spec)
    _, cluster = _own_blocks(spec)
    fetch = clone / ".pixi" / "envs" / "fetch" / "bin"
    (tmp_path / "paste.sh").write_text("set -e\n" + cluster.replace(guide.FASTQ_DIR, str(files)))
    done = _run([tmp_path / "paste.sh"], tmp_path, env, CTX_STORE=str(store),
                PATH=f"{fetch}:{env['PATH']}")  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    assert (root / "log" / "python3.calls").read_text().split() == ["called"]
    run = store / spec.slug
    for name, data in REFS_OWN.items():
        (run / "refs" / name).parent.mkdir(parents=True, exist_ok=True)
        (run / "refs" / name).write_bytes(data)
    job = _run([BASH / "rnaseq.sbatch", run, *(run / "refs" / n for n in REFS_OWN), "true"],
               run, env)  # fmt: skip
    assert job.returncode == 0, job.stderr[-3000:]
    assert [r["sample"] for r in _cut(run)] == list(spec.accessions)
    assert {r["fastq_1"] for r in _cut(run)} == {str(p.resolve()) for p in files.iterdir()}


def test_every_install_is_frozen(page_blocks):
    """U.18: every pixi install on the page is --frozen, the page says why before Part 1's
    first block, and the laptop route clones the repository as Kebnekaise's does. Catches a
    non-frozen install, which rewrites pixi.lock and dirties the kit's commit, and the clone
    left off the laptop route or off Part 1b."""
    every = "\n".join(c.text for c in codes(page_blocks[1]).values())
    installs = re.findall(r'pixi"? install[^&\n]*', every)
    assert len(installs) >= 5 and all("--frozen" in i for i in installs), installs
    (part1,) = [g for g in page_blocks[1] if g.title == text.fill("guide", "part1.group", n=1)]
    assert part1.items[0] == text.word("guide", "part1.frozen")
    clone = f'git clone {SOURCE_PLACEHOLDER} "$CTX_REPO"'
    for key in ("part1.laptop", "part1b.laptop", "part1.setup"):
        lines = codes(page_blocks[1])[text.word("guide", key)].text.splitlines()
        install = next(n for n, line in enumerate(lines) if " install " in line)
        assert clone in lines[:install], key


def test_part3_says_why_the_metadata_outdir(page_blocks):
    """U.18: Part 3 says why --fetchngs-dir names the metadata run's outdir, before the first
    block that names it. Catches the reason dropped or moved below the commands."""
    (part3,) = [g for g in page_blocks[1] if g.title == text.fill("guide", "part3.group", n=3)]
    first = next(n for n, i in enumerate(part3.items) if isinstance(i, Code))
    assert text.word("guide", "part3.meta") in part3.items[:first]
    assert f'--fetchngs-dir "$RUN/{guide.META}"' in part3.items[first].text
