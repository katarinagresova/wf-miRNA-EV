"""MultiQC custom content of the QC report (report_tables.py), with summary tables written by hand."""

import json

import pytest

from helpers import run_script, write_tsv

SAMPLES = ["S1", "S2"]
# Off-target check of each library: length, reads, result. S2 has no 17-nt reads.
OFF_TARGET = {
    "S1": [(17, 30, "miRNA_hairpin_other"), (17, 10, "rRNA"), (22, 80, "miRNA"), (22, 20, "tied")],
    "S2": [(22, 90, "miRNA"), (22, 10, "rRNA"), (23, 5, "unmapped")],
}


@pytest.fixture
def out(tmp_path):
    write_tsv(
        tmp_path / "library_summary.tsv",
        [(s, 1000, 50.0, 20.0, 3.0, 1.0, 20.0, 15.0, 12, 7, 75.0, 5.0, 10.0, 10.0) for s in SAMPLES],
        header=["sample", "kept_reads", "pct_miRNA", "pct_rRNA", "pct_mycoplasma", "pct_phiX", "mirbase_pct_tied",
                "mirbase_pct_tied_rRNA", "mirbase_detected_ge1", "mirbase_detected_ge10",
                "mirbase_off_target_pct_miRNA", "mirbase_off_target_pct_tied",
                "mirbase_off_target_pct_rRNA", "mirbase_off_target_pct_miRNA_hairpin_other"],
    )
    write_tsv(
        tmp_path / "composition.tsv",
        [(s, "mirbase", cls, "16-28", 100) for s in SAMPLES for cls in ["rRNA", "miRNA"]],
        header=["sample", "db", "class", "len_bin", "reads"],
    )
    write_tsv(tmp_path / "contaminants.tsv", [], header=["sample", "class", "accession", "species", "reads"])
    write_tsv(tmp_path / "insert_lengths.tsv", [(0, 1.0, 2.0), (1, 0.0, 0.0), (22, 60.0, 50.0), (23, 0.0, 0.0)],
              header=["length", *SAMPLES])
    off = []
    for s, rows in OFF_TARGET.items():
        d = tmp_path / "mirna" / s
        d.mkdir(parents=True)
        off.append(write_tsv(d / "off_target_mirbase.tsv.gz", rows, header=["length", "reads", "result"]))
    out = tmp_path / "out"
    out.mkdir()
    run_script(
        "report_tables.py",
        "--dbs", "mirbase",
        "--summary", tmp_path / "library_summary.tsv",
        "--composition", tmp_path / "composition.tsv",
        "--contaminants", tmp_path / "contaminants.tsv",
        "--insert-lengths", tmp_path / "insert_lengths.tsv",
        "--off-target", *off,
        "--out-dir", out,
    )
    return out


def read(out, name):
    with open(out / f"{name}_mqc.json") as fh:
        return json.load(fh)


def test_general_stats(out):
    data = read(out, "general_stats")["data"]["S1"]
    assert data["pct_contaminant"] == pytest.approx(4.0)
    # Another class only: rRNA and miRNA_hairpin_other, not miRNA or tied
    assert data["pct_off_target"] == pytest.approx(20.0)
    # The largest detected_min
    assert data["detected"] == 7


def test_mirna_counts_headers(out):
    headers = read(out, "mirna")["headers"]
    assert headers["mirbase_pct_tied"]["title"] == "% tied (mirbase)"
    assert headers["mirbase_pct_tied_rRNA"]["title"] == "% tied rRNA (mirbase)"


def test_classes_in_priority_order(out):
    section = read(out, "read_classes")
    assert section["categories"] == ["miRNA", "rRNA"]


def test_insert_lengths_end_at_last_length_with_reads(out):
    data = read(out, "insert_lengths")["data"]
    assert list(data["S1"]) == ["0", "1", "22"]


def test_no_contaminants(out):
    assert read(out, "contaminants")["plot_type"] == "html"


def test_off_target(out):
    results = read(out, "off_target")["data"][0]
    assert results["S1"] == {"miRNA_hairpin_other": 30, "rRNA": 10, "miRNA": 80, "tied": 20}
    lengths = read(out, "off_target_length")["data"][0]
    # Off-target: another class, not miRNA, tied or unmapped
    assert lengths["S1"] == {"17": 100.0, "22": 0.0, "23": None}
    # The same lengths for all libraries, a gap without reads
    assert lengths["S2"] == {"17": None, "22": 10.0, "23": 0.0}
