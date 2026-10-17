"""Summary of step 6 (mirna_summary.py), with counter interface tables written by hand."""

import pandas as pd
import pytest

from helpers import run_script, write_tsv

SAMPLES = ["S1", "S2"]
# Interface tables (mirna, reads, members): miRge3 merges hsa-miR-103a-3p and hsa-miR-107 into one entry
TABLES = {
    "mirge3": [
        ("hsa-miR-103a-3p/107", 30, "hsa-miR-103a-3p;hsa-miR-107"),
        ("hsa-let-7a-5p", 50, "hsa-let-7a-5p"),
    ],
    "bowtie": [
        ("hsa-miR-103a-3p", 20, "hsa-miR-103a-3p"),
        ("hsa-miR-107", 12, "hsa-miR-107"),
        ("hsa-let-7a-5p", 48, "hsa-let-7a-5p"),
        ("hsa-miR-21-5p", 5, "hsa-miR-21-5p"),
    ],
}
# Off-target check: database entry, reads, result
OFF_TARGET = [("hsa-miR-107", 10, "miRNA"), ("hsa-miR-107", 2, "rRNA"), ("hsa-let-7a-5p", 48, "miRNA")]
OFF_TARGET_COUNTER = [("miRNA", 78), ("rRNA", 2)]


def summarize(tmp_path, counter, cross_check):
    tmp_path.mkdir(exist_ok=True)
    write_tsv(tmp_path / "samples.tsv", [(s, f"{s}.fastq.gz", "run1") for s in SAMPLES], header=["sample", "fastq", "batch"])
    write_tsv(
        tmp_path / "composition.tsv",
        [(s, "mirbase", cls, len_bin, 100) for s in SAMPLES for cls, len_bin in [("miRNA", "16-28"), ("rRNA", "29-50")]],
        header=["sample", "db", "class", "len_bin", "reads"],
    )
    counts, off, off_counter = [], [], []
    for s in SAMPLES:
        d = tmp_path / "mirna" / s
        d.mkdir(parents=True, exist_ok=True)
        counts += [write_tsv(d / f"{c}_mirbase.tsv", rows, header=["mirna", "reads", "members"]) for c, rows in TABLES.items()]
        off.append(write_tsv(d / "off_target_mirbase.tsv.gz", OFF_TARGET, header=["mirna", "reads", "result"]))
        off_counter.append(write_tsv(d / "off_target_mirbase.counter.tsv", OFF_TARGET_COUNTER, header=["result", "reads"]))
    out = tmp_path / "out"
    run_script(
        "mirna_summary.py",
        "--samples", tmp_path / "samples.tsv",
        "--dbs", "mirbase",
        "--counter", counter,
        "--cross-check", cross_check,
        "--composition", tmp_path / "composition.tsv",
        "--counts", *counts,
        "--off-target", *off,
        "--off-target-counter", *off_counter,
        "--detected-min", 1, 10,
        "--out-dir", out,
        "--out-summary", out / "library_summary.tsv",
        "--out-agreement", out / "method_agreement.tsv",
        "--out-off-target", out / "off_target.tsv",
    )
    return out


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)


@pytest.fixture
def out(tmp_path):
    return summarize(tmp_path, counter="mirge3", cross_check="bowtie")


def test_counters_joined_through_members(out):
    cross = read(out / "cross_check" / "bowtie_mirbase.tsv", index_col="mirna")
    # The bowtie counts are summed into the groups of the miRge3 entries
    assert cross.loc["hsa-miR-103a-3p/107", "S1"] == 32
    assert cross.loc["hsa-let-7a-5p", "S1"] == 48
    assert cross.loc["hsa-miR-21-5p", "S1"] == 5  # no main-counter entry: its own name
    agreement = read(out / "method_agreement.tsv").set_index("sample")
    assert agreement.loc["S1", "n_ge10_main"] == 2
    assert agreement.loc["S1", "n_ge10_cross_check"] == 2
    assert agreement.loc["S1", "n_ge10_both"] == 2


def test_off_target_links_merged_entry(out):
    off = read(out / "off_target.tsv").set_index("mirna")
    assert off.loc["hsa-miR-107", "counter_mirna"] == "hsa-miR-103a-3p/107"
    assert off.loc["hsa-miR-107", "reads"] == 2 * 12  # two libraries
    assert off.loc["hsa-miR-107", "off_target_frac"] == pytest.approx(2 / 12, abs=1e-4)
    summary = read(out / "library_summary.tsv", index_col="sample")
    assert summary.loc["S1", "mirbase_off_target_counter_pct_other_class"] == pytest.approx(100 * 2 / 80, abs=1e-4)


def test_change_of_main_counter(tmp_path, out):
    other = summarize(tmp_path / "other", counter="bowtie", cross_check="mirge3")
    a = read(out / "library_summary.tsv", index_col="sample")
    b = read(other / "library_summary.tsv", index_col="sample")
    # The main counter has the same column names and paths; only the cross-check columns change their prefix
    assert [c.replace("bowtie_", "mirge3_") for c in a.columns] == list(b.columns)
    assert a.loc["S1", "mirbase_mirna_reads"] == 80
    assert b.loc["S1", "mirbase_mirna_reads"] == 85
    assert b.loc["S1", "mirge3_mirbase_mirna_reads"] == 80
    assert set(read(other / "counts_mirbase.tsv", index_col="mirna").index) == {m for m, *_ in TABLES["bowtie"]}
