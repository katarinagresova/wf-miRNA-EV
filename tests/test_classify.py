"""Classification rules of step 5 (quantify_mirna.py), with hit tables written by hand."""

import pytest

from helpers import COUNT, classify

RRNA = "KY962518.1|rRNA|45S_rDNA"
TRNA = "tRNA-Gly-GCC-1-1|tRNA|tRNA-Gly-GCC-1-1"
MYCO = "NC_014448.1|mycoplasma|M_hyorhinis"


def hairpin(name):
    return f"{name}|miRNA|{name}"


def test_tie_between_small_rna_classes(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 22, 2: 22},
        hits={
            "rrna": [(1, RRNA, 100, 22, 0), (2, RRNA, 100, 22, 1)],
            "trna": [(1, TRNA, 0, 22, 0), (2, TRNA, 0, 22, 0)],
        },
    )
    assert out.classes[1] == "rRNA"  # equal hits: the class priority decides
    assert out.classes[2] == "tRNA"  # the fewest mismatches come first


def test_contaminant_with_fewer_mismatches_than_human(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 20, 2: 17},
        hits={"rrna": [(1, RRNA, 100, 20, 1), (2, RRNA, 100, 17, 1)]},
        contaminant=[(1, MYCO, 5000, 20, 0, "-"), (2, MYCO, 5000, 17, 0, "-")],
    )
    assert out.classes[1] == "mycoplasma"
    assert out.classes[2] == "rRNA"  # shorter than contaminant_min_len, and a small RNA hit


@pytest.mark.parametrize(
    "genome_class, expected",
    [
        ("protein_coding", "mycoplasma"),
        ("lncRNA", "mycoplasma"),
        ("intron", "mycoplasma"),
        (".", "mycoplasma"),
        ("tRNA_locus", "tRNA_genomic"),  # not a weak human class
    ],
)
def test_contaminant_tie_with_genome_class(tmp_path, genome_class, expected):
    out = classify(
        tmp_path,
        reads={1: 20},
        genome=[(1, 0, genome_class)],
        contaminant=[(1, MYCO, 5000, 20, 0, "+")],
    )
    assert out.classes[1] == expected


def test_contaminant_tie_with_small_rna_stays_human(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 20},
        hits={"rrna": [(1, RRNA, 100, 20, 0)]},
        contaminant=[(1, MYCO, 5000, 20, 0, "+")],
    )
    assert out.classes[1] == "rRNA"


def test_short_read_with_contaminant_and_genome_hits(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 16, 2: 16, 3: 20},
        genome=[(1, 0, "intron"), (2, 1, "intron"), (3, 1, "intron")],
        contaminant=[(1, MYCO, 5000, 16, 0, "+"), (2, MYCO, 5000, 16, 1, "+"), (3, MYCO, 5000, 20, 1, "+")],
    )
    assert out.classes[1] == "mycoplasma"  # exact contaminant match, weak human evidence
    assert out.classes[2] == "intron"  # a short read needs an exact contaminant match
    assert out.classes[3] == "mycoplasma"  # a long read needs only a tie


def test_unmapped_reads_by_their_ends(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 40, 2: 40, 3: 40, 4: 40},
        ends=[
            # a contaminant end and no human end
            ("1_5_25", "contaminant", MYCO, 1),
            # a contaminant end and a human end
            ("2_5_25", "contaminant", MYCO, 0),
            ("2_3_25", "genome", "21", 0),
            # no 25-nt end hit: an exact 18-nt contaminant end
            ("3_5_18", "contaminant", MYCO, 0),
            # no 25-nt end hit: an 18-nt contaminant end with a mismatch does not count
            ("4_5_18", "contaminant", MYCO, 1),
        ],
    )
    assert out.classes[1] == "mycoplasma"
    assert out.classes[2] == "unmapped"
    assert out.classes[3] == "mycoplasma"
    assert out.classes[4] == "unmapped"


def test_hairpin_reads_and_arm_tolerances(tmp_path):
    # miR-a-5p is at 5..27 of mir-a; the tolerances are 3 nt at the 5' end and 5 nt at the 3' end
    reads = {
        1: (5, 22),  # the arm
        2: (8, 22),  # 5' end +3, 3' end +3
        3: (2, 22),  # 5' end -3, 3' end -3
        4: (5, 27),  # 3' end +5
        5: (9, 22),  # 5' end +4
        6: (5, 28),  # 3' end +6
        7: (30, 20),  # the loop
    }
    out = classify(
        tmp_path,
        reads={rank: length for rank, (_, length) in reads.items()},
        hits={"mirbase": [(rank, hairpin("mir-a"), offset, length, 0) for rank, (offset, length) in reads.items()]},
    )
    assert out.classes[[1, 2, 3, 4]].eq("miRNA").all()
    assert out.classes[[5, 6, 7]].eq("miRNA_hairpin_other").all()
    assert out.counts.loc["miR-a-5p", "reads"] == 4 * COUNT
    assert out.counts.loc["miR-a-5p", "canonical_reads"] == COUNT


def test_read_on_two_mature_mirnas(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 22, 2: 22},
        hits={
            "mirbase": [
                # two mature names: the count is divided equally
                (1, hairpin("mir-a"), 5, 22, 0),
                (1, hairpin("mir-b"), 5, 22, 0),
                # two paralogous hairpins with the same mature name: one mature name, not divided
                (2, hairpin("mir-c1"), 5, 22, 0),
                (2, hairpin("mir-c2"), 5, 22, 0),
            ]
        },
    )
    assert out.counts.loc["miR-a-5p", "reads"] == COUNT / 2
    assert out.counts.loc["miR-b-5p", "reads"] == COUNT / 2
    assert out.counts.loc["miR-c-5p", "reads"] == COUNT


def test_short_reads_need_an_exact_match(tmp_path):
    # miR-a-5p is at 5..27 of mir-a; reads of exact-max-len (18) nt or shorter count only with 0 mismatches
    out = classify(
        tmp_path,
        reads={1: 17, 2: 17, 3: 17, 4: 19},
        hits={
            "mirbase": [
                (1, hairpin("mir-a"), 5, 17, 0),
                (2, hairpin("mir-a"), 5, 17, 1),
                (3, hairpin("mir-a"), 5, 17, 1),
                (4, hairpin("mir-a"), 5, 19, 1),
            ],
            "rrna": [(3, RRNA, 100, 17, 1)],
        },
    )
    assert out.classes[1] == "miRNA"
    assert out.classes[2] == "miRNA_hairpin_other"
    assert out.classes[3] == "rRNA"  # the other hit with the same mismatches
    assert out.classes[4] == "miRNA"
    assert out.counts.loc["miR-a-5p", "reads"] == 2 * COUNT


def test_non_templated_3p_nucleotides(tmp_path):
    # Exact hairpin hits of the reads without 1-3 nt at the 3' end; they count only for reads without a hit in step 5
    out = classify(
        tmp_path,
        reads={1: 24, 2: 24, 3: 24, 4: 24},
        genome=[(3, 1, "intron")],
        nta=[
            (1, hairpin("mir-a"), 5, 22, 0),  # 2 nt removed: miR-a-5p
            (2, hairpin("mir-a"), 9, 21, 0),  # 3 nt removed: 5' end +4, outside the arm tolerances
            (3, hairpin("mir-a"), 5, 22, 0),  # a genome hit in step 5
            (4, hairpin("mir-a"), 45, 22, 0),  # 2 nt removed: miR-a-3p
            (4, hairpin("mir-b"), 5, 23, 0),  # 1 nt removed: miR-b-5p; the fewest removed nucleotides win
        ],
    )
    assert out.classes[[1, 4]].eq("miRNA").all()
    assert out.classes[2] == "unmapped"
    assert out.classes[3] == "intron"
    assert out.counts.loc["miR-a-5p", "nta_reads"] == COUNT
    assert out.counts.loc["miR-b-5p", "nta_reads"] == COUNT
    assert "miR-a-3p" not in out.counts.index
    assert out.counts["canonical_reads"].sum() == 0


def test_mature_names_with_the_same_sequence_are_one_entry(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 22},
        hits={"mirbase": [(1, hairpin("mir-d1"), 5, 22, 0), (1, hairpin("mir-d2"), 5, 22, 0)]},
    )
    assert out.counts.loc["miR-d1-5p/miR-d2-5p", "reads"] == COUNT  # not divided
    assert out.counts.loc["miR-d1-5p/miR-d2-5p", "members"] == "miR-d1-5p;miR-d2-5p"
