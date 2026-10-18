"""Classification rules of step 5 (quantify_mirna.py), with hit tables written by hand."""

import pandas as pd
import pytest

from helpers import COUNT, classify

RRNA = "KY962518.1|rRNA|45S_rDNA"
TRNA = "tRNA-Gly-GCC-1-1|tRNA|tRNA-Gly-GCC-1-1"
MYCO = "NC_014448.1|mycoplasma|M_hyorhinis"
ECOLI = "NC_000913.3|bacteria|E_coli"


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


def test_contaminant_reads_go_to_the_genomes_of_their_class(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 20},
        contaminant=[(1, MYCO, 5000, 20, 0, "+"), (1, ECOLI, 9000, 20, 0, "+")],
    )
    assert out.classes[1] == "mycoplasma"  # equal hits: the class priority decides
    cont = pd.read_csv(out.dir / "contaminants.tsv", sep="\t")
    assert cont[["class", "accession", "reads"]].values.tolist() == [["mycoplasma", "NC_014448.1", COUNT]]


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
            ("1_5_25", "contaminant", MYCO, 1, 5000),
            # a contaminant end and a human end
            ("2_5_25", "contaminant", MYCO, 0, 5000),
            ("2_3_25", "genome", "21", 0, 1000),
            # no 25-nt end hit: an exact 18-nt contaminant end
            ("3_5_18", "contaminant", MYCO, 0, 5000),
            # no 25-nt end hit: an 18-nt contaminant end with a mismatch does not count
            ("4_5_18", "contaminant", MYCO, 1, 5000),
        ],
    )
    assert out.classes[1] == "mycoplasma"
    assert out.classes[2] == "unmapped"
    assert out.classes[3] == "mycoplasma"
    assert out.classes[4] == "unmapped"


def test_unmapped_reads_by_their_human_ends(tmp_path):
    yrna = "ENST00000363120|Y_RNA|RNY1"
    out = classify(
        tmp_path,
        reads={1: 30, 2: 30, 3: 30, 4: 30, 5: 30, 6: 30, 7: 30, 8: 30, 9: 30, 10: 23},
        ends=[
            # a small RNA end with a mismatch
            ("1_3_25", "ncrna", yrna, 1, 0),
            # a 1-nt deletion: the 25-nt ends contain it and one has a genome hit; the exact 18-nt end gives the class
            ("2_5_25", "genome", "21", 1, 1000),
            ("2_3_18", "trna", TRNA, 0, 0),
            # the fewest mismatches, then the class priority
            ("3_5_25", "rrna", RRNA, 1, 100),
            ("3_3_25", "trna", TRNA, 0, 0),
            ("4_5_25", "rrna", RRNA, 0, 100),
            ("4_3_25", "trna", TRNA, 0, 0),
            # a 25-nt small RNA end wins over an exact 18-nt end
            ("5_5_25", "trna", TRNA, 1, 0),
            ("5_3_18", "rrna", RRNA, 0, 100),
            # genome and miRNA ends do not give a class
            ("6_5_25", "genome", "21", 0, 1000),
            ("7_5_25", "mirbase", hairpin("mir-a"), 0, 5),
            # a human end and a contaminant end
            ("8_5_25", "trna", TRNA, 0, 0),
            ("8_3_25", "contaminant", MYCO, 0, 5000),
            # an 18-nt end with a mismatch does not count
            ("9_5_18", "trna", TRNA, 1, 0),
            # a miRNA read by its non-templated 3' nucleotides keeps its class
            ("10_3_18", "trna", TRNA, 0, 0),
        ],
        nta=[(10, hairpin("mir-a"), 5, 22, 0)],
    )
    assert out.classes[1] == "Y_RNA"
    assert out.classes[2] == "tRNA"
    assert out.classes[3] == "tRNA"
    assert out.classes[4] == "rRNA"
    assert out.classes[5] == "tRNA"
    assert out.classes[[6, 7, 8, 9]].eq("unmapped").all()
    assert out.classes[10] == "miRNA"
    ends = pd.read_csv(out.dir / "unmapped_ends.tsv", sep="\t")
    assert ends.groupby("small_rna_class")["reads"].sum().to_dict() == {"Y_RNA": COUNT, "rRNA": COUNT, "tRNA": 5 * COUNT}


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
    # Exact hairpin hits of the reads without 1-3 nt at the 3' end
    out = classify(
        tmp_path,
        reads={1: 24, 2: 24, 4: 24},
        nta=[
            (1, hairpin("mir-a"), 5, 22, 0),  # 2 nt removed: miR-a-5p
            (2, hairpin("mir-a"), 9, 21, 0),  # 3 nt removed: 5' end +4, outside the arm tolerances
            (4, hairpin("mir-a"), 45, 22, 0),  # 2 nt removed: miR-a-3p
            (4, hairpin("mir-b"), 5, 23, 0),  # 1 nt removed: miR-b-5p; the fewest removed nucleotides win
        ],
    )
    assert out.classes[[1, 4]].eq("miRNA").all()
    assert out.classes[2] == "unmapped"
    assert out.counts.loc["miR-a-5p", "nta_reads"] == COUNT
    assert out.counts.loc["miR-b-5p", "nta_reads"] == COUNT
    assert "miR-a-3p" not in out.counts.index
    assert out.counts["canonical_reads"].sum() == 0


def test_non_templated_3p_nucleotides_against_other_hits(tmp_path):
    # An arm hit without 1 nt at the 3' end is a miRNA hit with 1 mismatch (miR-a-5p is at 5..27 of mir-a)
    out = classify(
        tmp_path,
        reads={1: 23, 2: 23, 3: 23, 4: 28, 5: 18},
        hits={
            "rrna": [(2, RRNA, 100, 23, 0), (3, RRNA, 100, 23, 1)],
            # 4: 3' end +6, outside the arm tolerances
            "mirbase": [(4, hairpin("mir-a"), 5, 28, 0), (5, hairpin("mir-a"), 5, 18, 1)],
        },
        # 1: the hairpin ends at the arm end; the added nucleotide is the mismatch of the genome hit
        genome=[(1, 1, "protein_coding")],
        nta=[(r, hairpin("mir-a"), 5, n - 1, 0) for r, n in [(1, 23), (2, 23), (3, 23), (4, 28), (5, 18)]],
    )
    reads = pd.read_csv(out.dir / "mirna_mirbase.reads.tsv.gz", sep="\t", index_col="rank")
    top = pd.read_csv(out.dir / "top_hits.tsv", sep="\t", index_col="rank")
    assert out.classes[1] == "miRNA"  # a small RNA class wins over a genome class
    assert reads.loc[1, "nta"] == 1
    assert top.loc[1, "mature"] == "miR-a-5p"
    assert out.classes[2] == "rRNA"  # an exact hit wins
    assert out.classes[3] == "miRNA"  # 1 mismatch each: the class priority decides
    assert reads.loc[3, "tied_classes"] == "rRNA"
    assert out.classes[4] == "miRNA_hairpin_other"  # an exact hairpin hit wins
    assert pd.isna(top.loc[4, "mature"])
    assert out.classes[5] == "miRNA"  # 18 nt with 1 mismatch: miRNA_hairpin_other without the NTA hit
    assert out.counts.loc["miR-a-5p", "nta_reads"] == 3 * COUNT
    assert out.counts.loc["miR-a-5p", "reads_tied_other_class"] == COUNT


def test_mature_names_with_the_same_sequence_are_one_entry(tmp_path):
    out = classify(
        tmp_path,
        reads={1: 22},
        hits={"mirbase": [(1, hairpin("mir-d1"), 5, 22, 0), (1, hairpin("mir-d2"), 5, 22, 0)]},
    )
    assert out.counts.loc["miR-d1-5p/miR-d2-5p", "reads"] == COUNT  # not divided
    assert out.counts.loc["miR-d1-5p/miR-d2-5p", "members"] == "miR-d1-5p;miR-d2-5p"
