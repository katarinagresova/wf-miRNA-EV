"""Summary of step 6 (mirna_summary.py), with counter interface tables written by hand."""

import pandas as pd
import pytest

from helpers import run_script, write_tsv

SAMPLES = ["S1", "S2"]
# Interface table of the counter (mirna, reads, members): hsa-miR-199a-3p and hsa-miR-199b-3p have the same sequence
TABLE = [
    ("hsa-let-7a-5p", 50, "hsa-let-7a-5p"),
    ("hsa-miR-199a-3p/hsa-miR-199b-3p", 30, "hsa-miR-199a-3p;hsa-miR-199b-3p"),
    ("hsa-miR-21-5p", 5, "hsa-miR-21-5p"),
]
# Off-target check: database entry, reads, result, tied classes
OFF_TARGET = [
    ("hsa-miR-199b-3p", 10, "miRNA", ""),
    ("hsa-miR-199b-3p", 2, "rRNA", ""),
    ("hsa-let-7a-5p", 45, "miRNA", ""),
    ("hsa-let-7a-5p", 2, "tied", "rRNA,tRNA"),
    ("hsa-let-7a-5p", 1, "tied", "tRNA"),
]


@pytest.fixture
def out(tmp_path):
    write_tsv(tmp_path / "samples.tsv", [(s, f"{s}.fastq.gz", "run1") for s in SAMPLES], header=["sample", "fastq", "batch"])
    write_tsv(
        tmp_path / "composition.tsv",
        [(s, "mirbase", cls, len_bin, 100) for s in SAMPLES for cls, len_bin in [("miRNA", "16-28"), ("rRNA", "29-50")]],
        header=["sample", "db", "class", "len_bin", "reads"],
    )
    counts, off = [], []
    for s in SAMPLES:
        d = tmp_path / "mirna" / s
        d.mkdir(parents=True)
        counts.append(write_tsv(d / "bowtie_mirbase.tsv", TABLE, header=["mirna", "reads", "members"]))
        off.append(write_tsv(d / "off_target_mirbase.tsv.gz", OFF_TARGET,
                             header=["mirna", "reads", "result", "tied_classes"]))
    out = tmp_path / "out"
    run_script(
        "mirna_summary.py",
        "--samples", tmp_path / "samples.tsv",
        "--dbs", "mirbase",
        "--composition", tmp_path / "composition.tsv",
        "--counts", *counts,
        "--off-target", *off,
        "--detected-min", 1, 10,
        "--out-dir", out,
        "--out-summary", out / "library_summary.tsv",
        "--out-off-target", out / "off_target.tsv",
    )
    return out


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)


def test_counts_and_detected_mirnas(out):
    counts = read(out / "counts_mirbase.tsv", index_col="mirna")
    assert list(counts.index) == [m for m, *_ in TABLE]
    assert counts.loc["hsa-miR-199a-3p/hsa-miR-199b-3p", "S2"] == 30
    summary = read(out / "library_summary.tsv", index_col="sample")
    assert summary.loc["S1", "mirbase_mirna_reads"] == 85
    assert summary.loc["S1", "mirbase_detected_ge1"] == 3
    assert summary.loc["S1", "mirbase_detected_ge10"] == 2


def test_off_target_links_merged_entry(out):
    off = read(out / "off_target.tsv").set_index("mirna")
    assert off.loc["hsa-miR-199b-3p", "counter_mirna"] == "hsa-miR-199a-3p/hsa-miR-199b-3p"
    assert off.loc["hsa-miR-199b-3p", "reads"] == 2 * 12  # two libraries
    assert off.loc["hsa-miR-199b-3p", "off_target_frac"] == pytest.approx(2 / 12, abs=1e-4)
    summary = read(out / "library_summary.tsv", index_col="sample")
    assert summary.loc["S1", "mirbase_off_target_pct_rRNA"] == pytest.approx(100 * 2 / 60, abs=1e-4)


def test_tied_reads(out):
    summary = read(out / "library_summary.tsv", index_col="sample")
    cols = list(summary.columns)
    i = cols.index("mirbase_pct_of_kept")
    assert cols[i + 1 : i + 3] == ["mirbase_pct_tied", "mirbase_pct_tied_rRNA"]
    # Percent of the miRNA reads (85); a read is tied with rRNA if rRNA is one of its tied classes
    assert summary.loc["S1", "mirbase_pct_tied"] == pytest.approx(100 * 3 / 85, abs=1e-4)
    assert summary.loc["S1", "mirbase_pct_tied_rRNA"] == pytest.approx(100 * 2 / 85, abs=1e-4)
