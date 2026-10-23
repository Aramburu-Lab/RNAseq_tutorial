"""R2.3: a dataset file that is not valid TOML is a DatasetError naming file, line and column."""

import pytest
from conftest import ctx_run

from ctx.dataset import DatasetError, load

BROKEN = '[dataset]\nslug = "x"\nlabel = \n'


def test_invalid_toml_names_line_and_column(tmp_path):
    """load raises DatasetError, not TOMLDecodeError, with the file and the place. Catches the
    decode error let through, and the file's name dropped from the message."""
    path = tmp_path / "broken.toml"
    path.write_text(BROKEN)
    with pytest.raises(DatasetError) as caught:
        load(path)
    message = str(caught.value)
    assert str(path) in message
    assert "not valid TOML" in message
    assert "line 3, column 9" in message


@pytest.mark.parametrize("command", [("run", "--check"), ("report",), ("de",)])
def test_commands_print_no_traceback(tmp_path, command):
    """ctx run, report and de on the broken file exit 1 with the problem, never a traceback.
    Catches report, which caught no DatasetError, printing its stack."""
    path = tmp_path / "broken.toml"
    path.write_text(BROKEN)
    done = ctx_run("--root", str(tmp_path / "root"), command[0], "--dataset", "x",
                   "--dataset-file", str(path), *command[1:], cwd=tmp_path)  # fmt: skip
    assert done.returncode == 1
    printed = done.stdout + done.stderr
    assert "line 3, column 9" in printed
    assert "Traceback" not in printed
    assert "TOMLDecodeError" not in printed


@pytest.mark.parametrize(
    "slug",
    [
        "Breast tissue (BRCA)",
        "Breast-tissue-(BRCA)",
        "BRCA-immune-cells:-macrophages",
        "-Breast",
        "Breast--BRCA",
    ],
)
def test_a_slug_outside_the_naming_rule_is_refused(tmp_path, slug):
    """A slug opens every output name, so one the naming lint refuses stops at load (D-114).
    Catches: a label typed as the slug passing until the first page is written."""
    path = tmp_path / "planted.toml"
    path.write_text(f'[dataset]\nslug = "{slug}"\n')
    with pytest.raises(DatasetError, match="dataset.slug"):
        load(path)


@pytest.mark.parametrize(
    "slug", ["Breast-tissue-BRCA", "Non-cancerous-blood-four-runs", "tcga_brca"]
)
def test_a_name_from_a_label_loads(tmp_path, slug):
    path = tmp_path / f"{slug}.toml"
    path.write_text(f'[dataset]\nslug = "{slug}"\n')
    assert load(path).slug == slug
