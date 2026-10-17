"""Interface tables of the miRNA counter bowtie, from the miRNA counts of step 5 (quantify_mirna.py).

Table: mirna (the entry: the mature names with the same sequence, joined by '/'), reads, members (the mature names,
joined by ';'), and the extra columns unique_seqs, canonical_reads and nta_reads.
Per-read table: sequence, mirna, reads; a read on several entries is split evenly between them.
"""

import argparse
import gzip

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True, help="step 5 counts: mirna, reads, members, unique_seqs, ...")
    p.add_argument("--mirna-reads", required=True, help="step 5 reads of the counts: rank, mirna, reads, ...")
    p.add_argument("--reads", required=True, help="collapsed.fa.gz (>s<rank>_x<count>)")
    p.add_argument("--out-table", required=True)
    p.add_argument("--out-reads", required=True)
    return p.parse_args()


def read_sequences(path, ranks):
    """rank -> sequence, for the given ranks."""
    seqs = {}
    with gzip.open(path, "rt") as fh:
        for header in fh:
            seq = next(fh)
            rank = int(header[2 : header.index("_x")])
            if rank in ranks:
                seqs[rank] = seq.rstrip()
    return seqs


def main():
    args = parse_args()
    table = pd.read_csv(args.counts, sep="\t")
    table[["mirna", "reads", "members", "unique_seqs", "canonical_reads", "nta_reads"]].to_csv(
        args.out_table, sep="\t", index=False
    )

    reads = pd.read_csv(args.mirna_reads, sep="\t", usecols=["rank", "mirna", "reads"])
    seqs = read_sequences(args.reads, set(reads["rank"]))
    reads = reads.assign(sequence=reads["rank"].map(seqs))
    assert reads["sequence"].notna().all()
    reads[["sequence", "mirna", "reads"]].to_csv(args.out_reads, sep="\t", index=False, compression="gzip")
    print(f"{len(table)} entries, {table['reads'].sum():,.0f} reads; {reads['sequence'].nunique():,} sequences")


if __name__ == "__main__":
    main()
