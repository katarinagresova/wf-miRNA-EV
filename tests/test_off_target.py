"""Off-target check (off_target.py): step 5 and the miRNA-only mapping, with hit tables written by hand."""

import pandas as pd

from helpers import COUNT, STEP5, classify, options, run_script, sequence, write_tsv

MIR_A = "mir-a|miRNA|mir-a"
RRNA = "KY962518.1|rRNA|45S_rDNA"


def off_target(step5, hits, counter_reads):
    """Run off_target.py on the outputs of classify(). hits: rows of the miRNA-only mapping; counter_reads: rows
    (sequence, mirna, reads) of the main counter."""
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
        "--out-counter", d / "off_target_mirbase.counter.tsv",
    )
    reads = pd.read_csv(d / "off_target_mirbase.tsv.gz", sep="\t")
    counter = pd.read_csv(d / "off_target_mirbase.counter.tsv", sep="\t").set_index("result")["reads"]
    return reads.set_index("sequence"), counter


def test_off_target_results(tmp_path):
    # All reads are 22 nt on the arm miR-a-5p (5..27 of mir-a)
    step5 = classify(
        tmp_path,
        reads={1: 22, 2: 22, 3: 22, 4: 22},
        hits={
            "mirbase": [(1, MIR_A, 5, 22, 1), (2, MIR_A, 5, 22, 0), (4, MIR_A, 5, 22, 0)],
            "rrna": [(1, RRNA, 300, 22, 0), (2, RRNA, 300, 22, 0)],
        },
    )
    seq = {rank: sequence(rank, 22) for rank in [1, 2, 3, 4]}
    reads, counter = off_target(
        step5,
        hits=[
            (1, MIR_A, 5, 22, 1),  # 1 mismatch to the miRNA, 0 to the rRNA
            (2, MIR_A, 5, 22, 0),  # equal hits on the miRNA and on the rRNA
            (3, MIR_A, 5, 22, 2),  # 2 mismatches to the miRNA, no hit in step 5
            (4, MIR_A, 5, 22, 0),  # a miRNA read
        ],
        counter_reads=[(seq[1], "miR-a-5p", COUNT), (seq[4], "miR-a-5p", COUNT)],
    )
    assert reads.loc[seq[1], "result"] == "rRNA"
    assert reads.loc[seq[2], "result"] == "tied"
    assert reads.loc[seq[3], "result"] == "unmapped"
    assert reads.loc[seq[4], "result"] == "miRNA"
    assert (reads["mirna"] == "miR-a-5p").all()
    # The main counter counts read 1, which step 5 gives to the rRNA
    assert counter.to_dict() == {"rRNA": COUNT, "miRNA": COUNT}
    assert reads.loc[seq[1], "counter_mirna"] == "miR-a-5p"
    assert reads["counter_mirna"].isna().sum() == 2


def test_off_target_read_on_two_mature_mirnas(tmp_path):
    step5 = classify(tmp_path, reads={1: 22}, hits={"rrna": [(1, RRNA, 300, 22, 0)]})
    reads, _ = off_target(
        step5,
        hits=[(1, MIR_A, 5, 22, 1), (1, "mir-b|miRNA|mir-b", 5, 22, 1)],
        counter_reads=[],
    )
    assert sorted(reads["mirna"]) == ["miR-a-5p", "miR-b-5p"]
    assert (reads["reads"] == COUNT / 2).all()  # divided equally, as in the bowtie counts
    assert (reads["result"] == "rRNA").all()
