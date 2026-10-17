"""scripts/python/kit_datasets.py on a made-up repository history (K.12).

A rebuild command must reproduce the run: the dataset file it names holds the bytes the run
read. Each test names the mutation it catches.
"""

import hashlib
import json
import os
import subprocess
from pathlib import Path

import pytest
from test_kit import NEW, SCRIPT, gdc_root, load

from ctx.report import text

ROOT = Path(__file__).resolve().parents[1]
FILE = "conf/datasets/x.toml"
OTHER = "conf/datasets/y.toml"


@pytest.fixture(scope="module")
def kit_script():
    return load(SCRIPT, "build_kit_k12")


def sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


class History:
    """A git repository whose commits get one second apart, so log's newest-first is fixed."""

    def __init__(self, repo: Path) -> None:
        self.repo, self.second = repo, 0
        repo.mkdir(parents=True)
        self.git("init", "-q", "-b", "main")

    def git(self, *args: str) -> str:
        self.second += 1
        when = f"2099-01-01T00:00:{self.second:02d}Z"
        env = {**os.environ, "GIT_AUTHOR_DATE": when, "GIT_COMMITTER_DATE": when}
        return subprocess.run(["git", "-C", str(self.repo), "-c", "user.name=t", "-c",
                               "user.email=t@example.invalid", *args], env=env, check=True,
                              capture_output=True, text=True).stdout.strip()  # fmt: skip

    def commit(self, **files: str) -> str:
        for path, body in files.items():
            target = self.repo / path
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_text(body)
            self.git("add", path)
        self.git("commit", "-q", "-m", "c", "--allow-empty")
        return self.git("rev-parse", "HEAD")


@pytest.fixture
def history(tmp_path):
    """v1, v2, v1 again, v3; a side branch at v4 and main at v5, merged into v6, which exists
    in the merge alone; then v7 at HEAD. y.toml holds y1 from the first commit on."""
    h = History(tmp_path / "repo")
    commits = {"v1": h.commit(**{FILE: "v1\n", OTHER: "y1\n"}), "v2": h.commit(**{FILE: "v2\n"})}
    commits["v1 again"] = h.commit(**{FILE: "v1\n"})
    commits["v3"] = h.commit(**{FILE: "v3\n"})
    h.git("checkout", "-q", "-b", "side")
    commits["v4"] = h.commit(**{FILE: "v4\n"})
    h.git("checkout", "-q", "main")
    commits["v5"] = h.commit(**{FILE: "v5\n"})
    h.git("merge", "-q", "-s", "ours", "--no-commit", "side")
    commits["v6"] = h.commit(**{FILE: "v6\n"})
    commits["v7"] = h.commit(**{FILE: "v7\n"})
    return h.repo, commits


def test_each_version_names_the_newest_commit_holding_it(history, kit_script):
    """K.12: a version is named by the newest commit holding it, a merge's own version
    included, and by none where the kit commit's copy holds it; git show at that commit gives
    the very bytes. Mutations caught: the oldest commit kept (v1 at its first commit); merges
    skipped (v6 missing); the kit commit's copy not compared (v7 and y1 given a commit)."""
    repo, commits = history
    found = kit_script.codebase_datasets(repo)
    want = {f"v{n}": commits[f"v{n}"] for n in (2, 3, 4, 5, 6)} | {"v1": commits["v1 again"]}
    for version, commit in want.items():
        assert found[sha(f"{version}\n")] == (FILE, commit), version
        shown = subprocess.run(["git", "-C", str(repo), "show", f"{commit}:{FILE}"],
                               capture_output=True, check=True).stdout  # fmt: skip
        assert hashlib.sha256(shown).hexdigest() == sha(f"{version}\n")
    assert found[sha("v7\n")] == (FILE, None) and found[sha("y1\n")] == (OTHER, None)
    assert len(found) == 8, found


def test_a_rebuild_reads_the_bytes_the_run_read(tmp_path, history, kit_script, monkeypatch):
    """K.12: a gather that read an earlier version is rebuilt from that version: LEFT_OUT.json
    and README.md put `git show <commit>:conf/datasets/<file> > <file>` before the gather
    command, which names <file>, and LEFT_OUT.json names the commit; the written file has the
    recorded SHA-256. Mutations caught: the command still naming conf/datasets/<file>; the
    git show line dropped or put after the command; the commit left out of LEFT_OUT.json."""
    repo, commits = history
    clone = tmp_path / "builder" / "conf" / "datasets" / "x.toml"
    real, folder = gdc_root(tmp_path, (str(clone), sha("v2\n")))
    codebase = kit_script.codebase_datasets(repo)
    monkeypatch.setattr(kit_script, "codebase_datasets", lambda: codebase)
    out = tmp_path / "kit"
    shipped = kit_script.build(out, [("step", real, True)])
    kit_script.readme(out, [("step", True)], shipped, "abc1234")
    left = json.loads((out / "LEFT_OUT.json").read_text())["files"]
    restore = f"git show {commits['v2']}:{FILE} > x.toml"
    gather = f"ctx --root {real} gather --dataset tcgax --gdc-dir {folder} --dataset-file x.toml"
    assert [e["rebuild"][1:] for e in left] == [[restore, gather]] * 2, left
    assert {e["dataset_file_commit"] for e in left} == {commits["v2"]}
    readme = (out / "README.md").read_text()
    assert readme.index(restore) < readme.index(gather), readme
    assert f"{NEW}_gather_tcgax_counts_v1.csv" in readme
    subprocess.run(restore, shell=True, cwd=repo, check=True)  # noqa: S602 - the line a reader runs
    assert hashlib.sha256((repo / "x.toml").read_bytes()).hexdigest() == sha("v2\n")


@pytest.mark.parametrize("version", ["v2", "v7"])
def test_readme_explains_a_git_show_line_exactly_where_one_is(tmp_path, history, kit_script,
                                                              monkeypatch, version):  # fmt: skip
    """K.14: README.md says in one sentence what a `git show` line writes, under each left-out
    file whose commands hold one (v2, an earlier version), and nowhere when none does (v7, the
    kit commit's copy). Mutations caught: the sentence never written; written for every entry."""
    repo, commits = history
    clone = tmp_path / "builder" / "conf" / "datasets" / "x.toml"
    real, _ = gdc_root(tmp_path, (str(clone), sha(f"{version}\n")))
    codebase = kit_script.codebase_datasets(repo)
    monkeypatch.setattr(kit_script, "codebase_datasets", lambda: codebase)
    out = tmp_path / "kit"
    shipped = kit_script.build(out, [("step", real, True)])
    readme = (kit_script.readme(out, [("step", True)], shipped, "abc1234")).read_text()
    left = json.loads((out / "LEFT_OUT.json").read_text())["files"]
    shows = sum(c.startswith("git show ") for e in left for c in e["rebuild"])
    said = text.fill("course", "readme.left_out_earlier", name="x.toml", path=FILE,
                     commit=commits["v2"]).strip()  # fmt: skip
    assert shows == (2 if version == "v2" else 0), left
    assert readme.count(said) == shows, readme
    assert readme.count("an earlier version;") == shows
