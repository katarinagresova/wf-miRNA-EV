"""rRNA fragment table (rrna_fragments.py): step 5 and the rRNA positions, with hit tables written by hand."""

import pandas as pd

from helpers import COUNT, STEP5, classify, run_script, write_tsv

RRNA = "KY962518.1|rRNA|45S_rDNA"
FEATURES = [(RRNA, 3657, 5526, "18S"), (RRNA, 6596, 6753, "5.8S"), (RRNA, 7920, 12971, "28S")]


def rrna_fragments(step5):
    """Run rrna_fragments.py on the outputs of classify()."""
    d = step5.dir
    run_script(
        "rrna_fragments.py",
        "--reads", d / "collapsed.fa.gz",
        "--classes", d / "read_classes.tsv.gz",
        "--hits", d / "rrna.hits.tsv.gz",
        "--contaminant", d / "contaminant.hits.tsv.gz",
        "--ends", d / "unmapped_ends.hits.tsv.gz",
        "--features", write_tsv(d / "rrna_features.bed", FEATURES),
        "--db", "mirbase",
        "--contaminant-min-len", STEP5["contaminant-min-len"],
        "--end-len", STEP5["end-len"],
        "--output", d / "fragments.tsv.gz",
    )
    return pd.read_csv(d / "fragments.tsv.gz", sep="\t")


def test_rrna_reads_with_a_position_from_an_end(tmp_path):
    step5 = classify(
        tmp_path,
        reads={1: 22, 2: 30, 3: 40, 4: 40},
        hits={"rrna": [(1, RRNA, 3700, 22, 0)]},
        ends=[
            # one rRNA end; the exact 18-nt end does not decide, because a 25-nt end has a hit
            ("2_5_25", "rrna", RRNA, 1, 6600),
            ("2_3_25", "genome", "21", 0, 1000),
            ("2_3_18", "rrna", RRNA, 0, 8000),
            # two joined fragments: both ends on rRNA
            ("3_5_25", "rrna", RRNA, 0, 3700),
            ("3_3_25", "rrna", RRNA, 0, 8000),
            # both ends on rRNA: the fewest mismatches
            ("4_5_25", "rrna", RRNA, 0, 3700),
            ("4_3_25", "rrna", RRNA, 1, 8000),
        ],
    )
    assert step5.classes.eq("rRNA").all()
    frag = rrna_fragments(step5)
    cols = ["rank", "length", "start", "end", "region", "mm", "n_pos", "from_end"]
    assert frag[cols].values.tolist() == [
        [1, 22, 3700, 3722, "18S", 0, 1, 0],
        [2, 30, 6600, 6625, "5.8S", 1, 1, 5],
        [3, 40, 3700, 3725, "18S", 0, 2, 5],
        [3, 40, 8000, 8025, "28S", 0, 2, 3],
        [4, 40, 3700, 3725, "18S", 0, 1, 5],
    ]
    assert (frag["count"] / frag["n_pos"]).sum() == 4 * COUNT
