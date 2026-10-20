"""M9 gene cards: GTEx V11 files, the HPA and Ensembl cache, the card list, the page.

Each test names the mutation it catches; every one was planted and seen to fail (TODO.md M9).
GTEx values are made up (``genecard_root``); NPTX2's HPA and Ensembl responses are real.
"""

import gzip
import hashlib
import json
import os
import re
import subprocess
import sys
import urllib.error
from datetime import date
from pathlib import Path

import numpy as np
import polars as pl
import pytest
from genecard_root import ACTB, FIXTURES, NPTX2, OURS, SLUG, THEIRS, Stub, make, responses

from ctx import runtime
from ctx.genecard import GenecardError, annotation, cards, gtex, hpa
from ctx.genecard import run as module
from ctx.genecard.cache import CACHE, MANIFEST, SUPERSEDED, Cache
from ctx.genecard.summary import spread, tau
from ctx.report import genes
from ctx.report.builder import unpacked
from ctx.setup import MANIFEST as SETUP_MANIFEST
from ctx.setup import fetch

CTX = Path(sys.executable).parent / "ctx"
GENE = "ENSG00000106236"


def _xml() -> hpa.Xml:
    return hpa.parse_xml(gzip.decompress((FIXTURES / "hpa_v25_ENSG00000106236.xml.gz")
                                         .read_bytes()), "v25")  # fmt: skip


def _files(root: Path) -> gtex.Files:
    folder = root / "data" / "raw" / "gtex_v11"
    return gtex.Files(*(next(folder.glob(p)) for p in ("*gene_tpm.gct.gz", "*.parquet",
                                                        "*SampleAttributesDS.txt",
                                                        "*SubjectPhenotypesDS.txt")))  # fmt: skip


@pytest.fixture
def root(tmp_path):
    spec = make(tmp_path / "root")
    return tmp_path / "root", spec


# tau and quantiles: the one definition


def test_tau_logs_before_the_maximum():
    # Mutation: tau on linear nTPM gives 0.9575 for NPTX2's 51 consensus tissues.
    values = [v for _, v in _xml().tissues]
    assert len(values) == 51  # noqa: PLR2004 - HPA v25's consensus tissues
    assert tau(values) == pytest.approx(0.6710, abs=5e-5)


def test_tau_leaves_out_the_pooled_total_and_noise():
    # Mutation: the pooled PBMC total counted as a cell type; tau on values all below 1 nTPM.
    from ctx.genecard.figures import immune_tau  # noqa: PLC0415 - imports plotting contracts

    cells = (("A", "x", 100.0), ("B", "x", 0.0), (hpa.POOLED, hpa.POOLED, 50.0))
    xml = hpa.Xml("25", (), cells, (), (), "")
    assert immune_tau(xml) == pytest.approx(1.0)
    assert immune_tau(_xml()) is None  # NPTX2: all 18 cell types under 1 nTPM


def test_spread_refuses_an_empty_group():
    # Mutation: an empty group summarised as zero.
    with pytest.raises(ValueError, match="no values"):
        spread([])
    assert spread([1.0, 2.0, 3.0, 10.0]) == spread(np.array([10.0, 3.0, 2.0, 1.0]))
    assert spread([1.0, 2.0, 3.0, 10.0]).q3 == pytest.approx(4.75)  # linear rule, not nearest


# GTEx V11 files


def test_tpm_joins_by_number_and_checks_dimensions(root):
    # Mutations: the versioned ID as the key (our .6 never meets GTEx's .4); the _PAR_Y copy
    # read as the gene; a cut file accepted.
    path = _files(root[0]).tpm
    samples, rows = gtex.read_tpm(path, [NPTX2])
    assert rows[NPTX2].gene_id == THEIRS[NPTX2] != OURS[NPTX2]
    assert 999.0 not in rows[NPTX2].tpm
    body = gzip.decompress(path.read_bytes()).decode().rsplit("\n", 2)[0] + "\n"
    path.write_bytes(gzip.compress(body.encode()))
    with pytest.raises(GenecardError, match="declares"):
        gtex.read_tpm(path, [NPTX2])


def test_tissue_spread_by_sex_absent_not_zero(root):
    # Mutation: a sex with no sample in a tissue given a row of zeros.
    files = _files(root[0])
    samples, rows = gtex.read_tpm(files.tpm, [NPTX2])
    table = gtex.tissue_spread(rows[NPTX2].tpm, samples, gtex.donors(files))
    uterus = table.filter(pl.col("tissue") == "Uterus")
    assert uterus["sex"].to_list() == ["all", "female"]
    pituitary = table.filter(pl.col("tissue") == "Pituitary", pl.col("sex") == "all").row(
        0, named=True
    )
    assert (pituitary["n"], pituitary["median"]) == (4, 240.0)


def test_exons_read_in_batches_match_one_read(root):
    # Mutation: a batch's columns written one place off, or its last partial batch dropped.
    files = _files(root[0])
    samples, _ = gtex.read_tpm(files.tpm, [NPTX2])
    whole = gtex.read_exons(files.exons, [NPTX2, ACTB], samples, batch=len(samples))
    for batch in (1, 2, 4):
        parts = gtex.read_exons(files.exons, [NPTX2, ACTB], samples, batch=batch)
        assert np.array_equal(parts[NPTX2].reads, whole[NPTX2].reads)
    assert whole[NPTX2].exons == ("E1", "E2", "E3", "E4", "E5")
    assert whole[ACTB].reads.shape == (1, len(samples))


def test_exon_share_per_sample_before_the_median(root):
    # Mutation: the median of reads first, then the share; and a sample without reads kept.
    files = _files(root[0])
    samples, _ = gtex.read_tpm(files.tpm, [NPTX2])
    reads = gtex.ExonReads("g", ("E1", "E2"), np.array([[1.0, 9.0, 0.0], [9.0, 1.0, 0.0]]))
    three = samples[:3]
    who = pl.DataFrame({"sample": three, "tissue": ["T"] * 3, "sex": ["male"] * 3})
    shares, dropped = gtex.exon_shares(reads, three, who)
    assert dropped == 1
    assert shares["share"].to_list() == [0.5, 0.5]  # median of (0.1, 0.9); reads first gives 0.5
    reads = gtex.ExonReads("g", ("E1", "E2"), np.array([[1.0, 90.0, 5.0], [9.0, 10.0, 5.0]]))
    shares, _ = gtex.exon_shares(reads, three, who)
    assert shares["share"].to_list() == [0.5, 0.5]  # per sample (0.1, 0.9, 0.5): median 0.5
    medians_first = np.median(reads.reads, axis=1) / np.median(reads.reads, axis=1).sum()
    assert not np.allclose(medians_first, [0.5, 0.5])


# HPA, Ensembl and the GTF


def test_hpa_stamp_types_and_numbers():
    # Mutations: XML 25 filed under v26; a string concentration left a string; integer
    # against string; the expected NPTX2 numbers (HPA v25.1, plan section 3).
    body = gzip.decompress((FIXTURES / "hpa_v25_ENSG00000106236.xml.gz").read_bytes())
    with pytest.raises(GenecardError, match="version '25'"):
        hpa.parse_xml(body, "v26")
    xml = _xml()
    assert max(xml.single_cell, key=lambda c: c[1]) == ("Corticotrophs", 335.9)
    assert dict(xml.tissues)["Pituitary gland"] == 200.0  # noqa: PLR2004 - HPA v25
    facts = hpa.parse_json(gzip.decompress(
        (FIXTURES / "hpa_v25_ENSG00000106236.json.gz").read_bytes()))  # fmt: skip
    assert facts.blood_im == 2_490_000.0  # pg/L: 2.49 ug/L, TODO.md's 2.5
    assert facts.blood_ms == 220_000.0
    assert facts.categories["RNA blood cell specificity"] == "Not detected in immune cells"
    with pytest.raises(GenecardError, match="neither"):
        hpa.parse_json(json.dumps({hpa.BLOOD_IM: ["2490000"]}).encode())


def test_gtf_and_rest_name_the_same_mane_exons(root):
    # Mutation: introns one base off their exons, or exons read in file order on minus strand.
    gtf = root[0] / "resources" / "gencode" / "made_up.gtf.gz"
    mane = annotation.gtf_mane(gtf, [NPTX2])[NPTX2]
    rest = annotation.rest_mane(gzip.decompress(
        (FIXTURES / "ensembl_116_ENSG00000106236_lookup.json.gz").read_bytes()))  # fmt: skip
    assert rest is not None and (rest.transcript, rest.refseq) == ("ENST00000265634.4",
                                                                   "NM_002523.3")  # fmt: skip
    assert mane.exons == rest.exons and len(mane.exons) == 5  # noqa: PLR2004 - MANE, 5 exons
    assert mane.introns()[0] == (98617888, 98619642)


def test_rest_exons_in_transcription_order_on_minus_strand():
    # Mutation: a minus-strand gene's exons left in genomic order, so they never equal the GTF's.
    body = {
        "strand": -1,
        "seq_region_name": "7",
        "Transcript": [
            {
                "MANE": [
                    {"type": "MANE_Select", "id": "ENST1", "version": 2, "refseq_match": "NM_1.1"}
                ],
                "Exon": [
                    {"start": 10, "end": 20},
                    {"start": 50, "end": 60},
                    {"start": 30, "end": 40},
                ],
            }
        ],
    }
    mane = annotation.rest_mane(json.dumps(body).encode())
    assert mane is not None and mane.exons == ((50, 60), (30, 40), (10, 20))
    assert mane.introns() == ((21, 29), (41, 49))


def test_qpcr_pair_skips_exons_without_reads():
    # Mutation: a NaN partner passing, as max(0.05, nan) is 0.05 (found by this test); exons
    # ordered as text, E10 before E2.
    pair = annotation.qpcr_pair({"E1": float("nan"), "E2": 0.4, "E3": 0.2, "E10": 0.05,
                                 "E11": float("nan")})  # fmt: skip
    assert pair == annotation.Pair("E3", "E10", 0.2)
    assert annotation.qpcr_pair({"E1": float("nan"), "E2": float("nan")}) is None


# The cache


def test_cache_drift_offline_refresh(tmp_path):
    # Mutations: a changed body silently overwriting; a request under --offline; --refresh
    # deleting the old copy; a tampered cache file read; ctx setup's manifest touched.
    url = hpa.url("v25", GENE, "json")
    stub = Stub({url: b'{"a": 1}'})
    first = Cache(tmp_path, opener=stub).get("hpa", "v25", GENE, "json", url)
    assert first.body == b'{"a": 1}' and len(stub.calls) == 1
    again = Cache(tmp_path, opener=stub).get("hpa", "v25", GENE, "json", url)
    assert again.body == first.body and len(stub.calls) == 1
    offline = Cache(tmp_path, offline=True, opener=stub)
    missing = offline.get("hpa", "v25", GENE, "xml", hpa.url("v25", GENE, "xml"))
    assert missing.body is None and "--offline" in missing.error and len(stub.calls) == 1
    cached = tmp_path / CACHE / "hpa" / "v25" / GENE / "json.gz"
    cached.unlink()
    stub.bodies[url] = b'{"a": 2}'
    with pytest.raises(GenecardError, match="another body"):
        Cache(tmp_path, opener=stub).get("hpa", "v25", GENE, "json", url)
    fresh = Cache(tmp_path, refresh=True, opener=stub, today=date(2026, 10, 8))
    assert fresh.get("hpa", "v25", GENE, "json", url).body == b'{"a": 2}'
    cached.write_bytes(gzip.compress(b'{"a": 3}'))
    with pytest.raises(GenecardError, match="SHA-256"):
        Cache(tmp_path, opener=stub).get("hpa", "v25", GENE, "json", url)
    assert (tmp_path / SUPERSEDED / "2026-10-08" / "hpa" / "v25" / GENE).is_dir()
    recorded = json.loads((tmp_path / MANIFEST).read_text())
    assert recorded[f"hpa/v25/{GENE}/json.gz"]["sha256"] == hashlib.sha256(b'{"a": 2}').hexdigest()
    assert not (tmp_path / SETUP_MANIFEST).exists()


def test_failed_call_is_not_cached(tmp_path):
    # Mutation: an error body cached, so the next online run never asks again; a 404 is an
    # answer and is cached (test_genecard_absent).
    def broken(url: str) -> None:
        raise urllib.error.HTTPError(url, 500, "Server Error", None, None)  # type: ignore[arg-type]

    answer = Cache(tmp_path, opener=broken).get("hpa", "v25", GENE, "json", "https://x/y")
    assert answer.body is None and "500" in answer.error
    assert not (tmp_path / MANIFEST).exists()


# The card list


def test_symbols_ambiguous_or_unknown_stop_the_run(root):
    # Mutation: the first ID of an ambiguous symbol taken; an unknown name skipped silently.
    _, table = cards.gene_table(root[0], SLUG)
    with pytest.raises(GenecardError, match="TWIN names ENSG00000900001, ENSG00000900002"):
        cards.named("TWIN", table, "--genes")
    with pytest.raises(GenecardError, match="NDUFA4L2 names no gene"):
        cards.named("NPTX2 NDUFA4L2", table, "--genes")
    assert cards.named("nptx2\nENSG00000075624", table, "--genes") == [NPTX2, ACTB]


def test_default_list_order(root):
    # Mutation: top calls ordered by |log2 FC| first, or a gene listed twice.
    from ctx.dataset import load  # noqa: PLC0415

    spec = load(root[1])
    _, table = cards.gene_table(root[0], SLUG)
    listed = cards.card_list(spec, table, None, None, root[0])
    assert listed == [(NPTX2, "registered"), (900014, "up in a_case"), (900011, "up in a_case"),
                      (900012, "up in a_case"), (900013, "down in a_case")]  # fmt: skip


# The run, end to end


def _register(root: Path, spec: Path) -> None:
    from ctx.dataset import load  # noqa: PLC0415

    hits: list[str] = []
    fetch(root, load(spec).resources, opener=lambda url, **_: hits.append(url))
    assert hits == [], "files in place are registered, never fetched"


def _cards():
    import typer  # noqa: PLC0415

    from ctx.cli import app  # noqa: PLC0415
    from ctx.report.cards import parameter_cards  # noqa: PLC0415

    group = typer.main.get_command(app)
    return parameter_cards(group, group.commands["genecard"])  # type: ignore[attr-defined]


def _settings(root: Path) -> runtime.Settings:
    return runtime.Settings("laptop", root, 1, None, 0.0)


def test_run_degrades_and_stamps_the_footer(root):
    # Mutations: the whole card dropped when HPA fails; the footer without GTEx's dated
    # acknowledgement, the currency line, or HPA's versioned URL; a placeholder left unfilled.
    folder, spec = root
    _register(folder, spec)
    setup_manifest = (folder / SETUP_MANIFEST).read_bytes()
    options = module.Options(spec, genes=("NPTX2", "ACTB"))
    written = module.run(
        _settings(folder), _cards(), "ctx genecard", SLUG, options, Stub(responses())
    )
    page = unpacked(written[0].read_text())
    assert "Human Protein Atlas unavailable for ACTB" in page
    assert "Ensembl REST unavailable for ACTB" in page
    assert "Exon use of ACTB" in page and "NPTX2 across immune cell types" in page
    assert re.search(r"obtained from the GTEx Portal on \d{4}-\d{2}-\d{2}", page)
    assert "may not reflect the most current data" in page
    assert "https://v25.proteinatlas.org/ENSG00000106236" in page
    shown = re.sub(r"<(script|style)\b.*?</\1>", "", page, flags=re.S)
    assert not re.search(r"\{[a-z_]+(:[^}]*)?\}", re.sub(r"<[^>]+>", " ", shown))
    assert "same MANE Select exons" in page
    assert (folder / SETUP_MANIFEST).read_bytes() == setup_manifest
    assert len(page.encode()) < 15_000_000  # noqa: PLR2004 - the HTML report cap
    written_dirs = {p.parent.relative_to(folder).as_posix() for p in written}
    assert written_dirs == {"outputs/11_genecard"}
    tissues = pl.read_parquet(next(p for p in written if p.name.endswith("_tissues_v1.parquet")))
    assert set(tissues["symbol"]) == {"NPTX2", "ACTB"}


def test_cli_offline_reads_only_the_cache(root):
    # Mutation: --offline asking the network; the CLI's cards out of step with its options.
    folder, spec = root
    _register(folder, spec)
    module.run(
        _settings(folder),
        _cards(),
        "ctx genecard",
        SLUG,
        module.Options(spec, genes=("NPTX2",)),
        Stub(responses()),
    )
    env = {k: v for k, v in os.environ.items() if not k.startswith(("SLURM_", "CTX_"))}
    env["https_proxy"] = env["http_proxy"] = "http://127.0.0.1:9"  # any request would fail
    done = subprocess.run(
        [str(CTX), "--root", str(folder), "genecard", "--dataset", SLUG, "--dataset-file",
         str(spec), "--genes", "NPTX2,ACTB", "--offline"],
        cwd=folder, env=env, capture_output=True, text=True, check=False,
    )  # fmt: skip
    assert done.returncode == 0, done.stderr[-3000:]
    (page,) = (folder / "outputs" / "11_genecard").glob("*_overview_v1.html")
    html = unpacked(page.read_text())
    assert "NPTX2 in the Human Protein Atlas v25" in html
    assert "not cached, and --offline makes no request" in html  # ACTB


def test_unregistered_gtex_degrades(root):
    # Mutation: a missing GTEx registration crashing the run instead of naming the command.
    folder, spec = root
    options = module.Options(spec, genes=("NPTX2",), offline=True)
    written = module.run(_settings(folder), _cards(), "ctx genecard", SLUG, options)
    assert "GTEx unavailable for NPTX2" in unpacked(written[0].read_text())
    assert "ctx setup --dataset-file" in unpacked(written[0].read_text())


def test_number_joins_versions():
    # Mutation: identity by versioned ID.
    assert genes.number(OURS[NPTX2]) == genes.number(THEIRS[NPTX2]) == NPTX2
