"""Interface tables of the miRNA counter bowtie, from the mature miRNA counts of step 5 (quantify_mirna.py).

Table: mirna, reads, members (the mature name itself), and the extra columns unique_seqs and canonical_reads.
Per-read table: sequence, mirna, reads; a read on several mature miRNAs is split evenly between them.
"""

import argparse
import gzip

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", required=True, help="step 5 counts: mature, reads, unique_seqs, canonical_reads, ...")
    p.add_argument("--mirna-reads", required=True, help="step 5 reads of the counts: rank, mature, reads, ...")
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
    counts = pd.read_csv(args.counts, sep="\t")
    table = counts.rename(columns={"mature": "mirna"}).assign(members=lambda t: t["mirna"])
    table[["mirna", "reads", "members", "unique_seqs", "canonical_reads"]].to_csv(args.out_table, sep="\t", index=False)

    reads = pd.read_csv(args.mirna_reads, sep="\t", usecols=["rank", "mature", "reads"])
    seqs = read_sequences(args.reads, set(reads["rank"]))
    reads = reads.assign(sequence=reads["rank"].map(seqs)).rename(columns={"mature": "mirna"})
    assert reads["sequence"].notna().all()
    reads[["sequence", "mirna", "reads"]].to_csv(args.out_reads, sep="\t", index=False, compression="gzip")
    print(f"{len(table)} miRNAs, {table['reads'].sum():,.0f} reads; {reads['sequence'].nunique():,} sequences")


if __name__ == "__main__":
    main()
