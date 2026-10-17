"""Input of the miRNA counters for one library and miRNA database: its length-filtered reads without the contaminant
reads.

miRge3.0 accepts loose isomiRs (16-nt reads, 5' additions, mismatches). Given all reads, it counted Mesomycoplasma
hyorhinis fragments as human miRNAs in the SH-SY5Y libraries (e.g. GATTCCTGTCGGAGCC as a miR-181b-5p isomiR: 407 of
S18's 459 reported reads). Reads whose final class (read_classes.tsv.gz, items 4, 4b and 4c) is a contaminant are
therefore removed, so the counters see the human reads only, the denominator of the rRNA shares too.

Reads a FASTQ on stdin and writes the kept records to stdout unchanged.
"""

import argparse
import gzip
import sys

import pandas as pd

from quantify_mirna import CONTAMINANT_PRIORITY


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--classes", required=True, help="read_classes.tsv.gz: rank, count, length, class_<db>...")
    p.add_argument("--reads", required=True, help="collapsed.fa.gz (>s<rank>_x<count>)")
    p.add_argument("--db", required=True, help="miRNA database whose read classes are used")
    return p.parse_args()


def contaminant_sequences(classes_path, reads_path, db):
    cls = pd.read_csv(classes_path, sep="\t", usecols=["rank", f"class_{db}"])
    ranks = set(cls.loc[cls[f"class_{db}"].isin(CONTAMINANT_PRIORITY), "rank"])
    seqs = set()
    with gzip.open(reads_path, "rt") as fh:
        keep = False
        for line in fh:
            if line.startswith(">"):
                keep = int(line[2:].split("_x")[0]) in ranks
            elif keep:
                seqs.add(line.strip().encode())
    return seqs


def main():
    args = parse_args()
    drop = contaminant_sequences(args.classes, args.reads, args.db)
    n_in = n_out = 0
    out = sys.stdout.buffer
    for record in zip(*[iter(sys.stdin.buffer)] * 4):
        n_in += 1
        if record[1].rstrip(b"\n") in drop:
            continue
        n_out += 1
        out.writelines(record)
    print(
        f"class_{args.db}: {len(drop):,} contaminant sequences; {n_in:,} reads in, "
        f"{n_in - n_out:,} contaminant reads removed, {n_out:,} kept",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
