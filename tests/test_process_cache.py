"""ctx.process.cache: a keyed folder is reused only while intact, and never overwritten (D-037)."""

import json

import pytest

from ctx.process import cache

INPUTS = {"tool": "STAR", "version": "2.7.11b", "genome_sha256": "aa", "sa_index_bases": 7}


class Builder:
    """A build that writes two files and counts its calls."""

    def __init__(self) -> None:
        self.calls = 0

    def __call__(self, dest):
        self.calls += 1
        (dest / "Genome").write_text("ACGT" * 100)
        (dest / "sub").mkdir()
        (dest / "sub" / "SA").write_text("x" * 50)


def test_a_second_call_reuses_the_folder(tmp_path):
    """The same inputs build nothing and touch nothing.

    Mutation caught: a manifest that does not record the key (every call rebuilds)."""
    build = Builder()
    folder, built = cache.keyed(tmp_path, "star", INPUTS, build)
    assert built and folder.name == f"star_{cache.key(INPUTS)[:12]}"
    stamps = {p: p.stat().st_mtime_ns for p in folder.rglob("*")}
    again, built = cache.keyed(tmp_path, "star", INPUTS, build)
    assert (again, built, build.calls) == (folder, False, 1)
    assert {p: p.stat().st_mtime_ns for p in folder.rglob("*")} == stamps
    manifest = json.loads((folder / cache.MANIFEST).read_text())
    assert manifest["inputs"] == INPUTS
    assert manifest["files"] == {"Genome": 400, "sub/SA": 50}


@pytest.mark.parametrize("name", sorted(INPUTS))
def test_each_input_is_in_the_key(tmp_path, name):
    """Mutation caught: a key that leaves an input out (two contents share one folder)."""
    build = Builder()
    first, _ = cache.keyed(tmp_path, "star", INPUTS, build)
    changed = {**INPUTS, name: f"{INPUTS[name]}1"}
    second, built = cache.keyed(tmp_path, "star", changed, build)
    assert built and second != first
    assert first.is_dir() and (first / cache.MANIFEST).is_file()


def test_a_truncated_file_forces_a_rebuild(tmp_path):
    """The damaged folder is set aside, not overwritten or deleted.

    Mutation caught: reuse without comparing sizes."""
    build = Builder()
    folder, _ = cache.keyed(tmp_path, "star", INPUTS, build)
    (folder / "sub" / "SA").write_text("x" * 10)
    again, built = cache.keyed(tmp_path, "star", INPUTS, build)
    assert (again, built, build.calls) == (folder, True, 2)
    assert (folder / "sub" / "SA").stat().st_size == 50
    aside = folder.with_name(f"{folder.name}.stale1")
    assert (aside / "sub" / "SA").stat().st_size == 10


def test_a_missing_file_forces_a_rebuild(tmp_path):
    """A listed file that is gone is damage like a short one: set aside, rebuilt.

    Mutation caught: the size read without checking the file is there (FileNotFoundError)."""
    build = Builder()
    folder, _ = cache.keyed(tmp_path, "star", INPUTS, build)
    (folder / "sub" / "SA").unlink()
    again, built = cache.keyed(tmp_path, "star", INPUTS, build)
    assert (again, built, build.calls) == (folder, True, 2)
    assert (folder / "sub" / "SA").stat().st_size == 50
    aside = folder.with_name(f"{folder.name}.stale1")
    assert (aside / "Genome").is_file() and not (aside / "sub" / "SA").exists()


def test_a_leftover_part_is_never_reused(tmp_path):
    """A crashed build's `.part`, even one holding a valid manifest, is rebuilt from nothing.

    Mutation caught: finishing a leftover `.part` by renaming it."""
    build = Builder()
    folder, _ = cache.keyed(tmp_path, "star", INPUTS, build)
    part = folder.with_name(f"{folder.name}.part")
    folder.rename(part)
    (part / "leftover").write_text("from the crashed build")
    again, built = cache.keyed(tmp_path, "star", INPUTS, build)
    assert (again, built, build.calls) == (folder, True, 2)
    assert not (folder / "leftover").exists()
    assert (folder.with_name(f"{part.name}.stale1") / "leftover").is_file()


def test_a_path_is_not_an_input(tmp_path):
    """A path names a place, not a content; inputs give a file by its SHA-256.

    Mutation caught: json.dumps with default=str (a moved file would keep its key)."""
    with pytest.raises(TypeError):
        cache.key({"genome": tmp_path})
