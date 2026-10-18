"""Off-target check (off_target.py): step 5 and the miRNA-only mapping, with hit tables written by hand."""

import pandas as pd

from helpers import COUNT, STEP5, classify, options, run_script, sequence, write_tsv

MIR_A = "mir-a|miRNA|mir-a"
RRNA = "KY962518.1|rRNA|45S_rDNA"


def off_target(step5, hits, counter_reads):
    """Run off_target.py on the outputs of classify(). hits: rows of the miRNA-only mapping; counter_reads: rows
    (sequence, mirna, reads) of the per-read table of the miRNA counter."""
    d = step5.dir
    run_script(
        "off_target.py",
        "--db", "mirbase",
        "--hits", write_tsv(d / "off_target.hits.tsv.gz", hits),
        "--arms", d / "mirbase_arms.tsv",
        "--reads", d / "collapsed.fa.gz",
        "--classes", d / "read_classes.tsv.gz",
        "--mirna-reads", d / "mirna_mirbase.reads.tsv.gz",
        "--counter-reads", write_tsv(d / "counter.reads.tsv", counter_reads, header=["sequence", "mirna", "reads"]),
        *options({k: v for k, v in STEP5.items() if k in ["min-len", "max-len", "start-tol", "end-tol"]}),
        "--out-reads", d / "off_target_mirbase.tsv.gz",
    )
    return pd.read_csv(d / "off_target_mirbase.tsv.gz", sep="\t").set_index("sequence")


def test_off_target_results(tmp_path):
    # All reads are on the arm miR-a-5p (5..27 of mir-a)
    reads = {1: 22, 2: 22, 3: 22, 4: 22, 5: 17}
    step5 = classify(
        tmp_path,
        reads=reads,
        hits={
            "mirbase": [(1, MIR_A, 5, 22, 1), (2, MIR_A, 5, 22, 0), (4, MIR_A, 5, 22, 0), (5, MIR_A, 5, 17, 1)],
            "rrna": [(1, RRNA, 300, 22, 0), (2, RRNA, 300, 22, 0)],
        },
    )
    seq = {rank: sequence(rank, length) for rank, length in reads.items()}
    reads = off_target(
        step5,
        hits=[
            (1, MIR_A, 5, 22, 1),  # 1 mismatch to the miRNA, 0 to the rRNA
            (2, MIR_A, 5, 22, 0),  # equal hits on the miRNA and on the rRNA
            (3, MIR_A, 5, 22, 2),  # 2 mismatches to the miRNA, no hit in step 5
            (4, MIR_A, 5, 22, 0),  # a miRNA read
            (5, MIR_A, 5, 17, 1),  # a 17-nt read with 1 mismatch: step 5 does not count it
        ],
        counter_reads=[(seq[2], "miR-a-5p", COUNT), (seq[4], "miR-a-5p", COUNT)],
    )
    assert reads.loc[seq[1], "result"] == "rRNA"
    assert reads.loc[seq[2], "result"] == "tied"
    assert reads.loc[seq[3], "result"] == "unmapped"
    assert reads.loc[seq[4], "result"] == "miRNA"
    assert reads.loc[seq[5], "result"] == "miRNA_hairpin_other"
    assert reads.loc[seq[2], "tied_classes"] == "rRNA"
    assert reads["tied_classes"].isna().sum() == 4
    assert (reads["mirna"] == "miR-a-5p").all()
    assert reads.loc[[seq[2], seq[4]], "counter_mirna"].eq("miR-a-5p").all()
    assert reads["counter_mirna"].isna().sum() == 3


def test_off_target_read_on_two_mature_mirnas(tmp_path):
    step5 = classify(tmp_path, reads={1: 22}, hits={"rrna": [(1, RRNA, 300, 22, 0)]})
    reads = off_target(
        step5,
        hits=[(1, MIR_A, 5, 22, 1), (1, "mir-b|miRNA|mir-b", 5, 22, 1)],
        counter_reads=[],
    )
    assert sorted(reads["mirna"]) == ["miR-a-5p", "miR-b-5p"]
    assert (reads["reads"] == COUNT / 2).all()  # divided equally, as in the bowtie counts
    assert (reads["result"] == "rRNA").all()
