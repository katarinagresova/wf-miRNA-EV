"""Off-target check of one library and miRNA database: the reads that a miRNA-only mapping counts as miRNA, but that
step 5 gives to another class.

The miRNA-only mapping aligned the kept reads of min-len..max-len nt (contaminants included) to the miRNA references
only, with up to 3 mismatches (best stratum). A hit counts for a mature miRNA under the arm tolerances of the bowtie
counts (quantify_mirna.arm_assignments), also with mismatches in a short read. Each read with such a hit gets a
result from its class in step 5:

- miRNA: the class is miRNA, and no other class has an equal hit,
- tied: the class is miRNA, but another class has an equal hit (quantify_mirna.py, tied classes),
- unmapped: step 5 has no hit (the read has 2 or 3 mismatches to the miRNA),
- else the class of step 5: another reference has a better hit, or the read is short and has a mismatch
  (miRNA_hairpin_other; step 5 needs 0 mismatches for these reads).

A read on several mature miRNAs is split evenly between them, as in the bowtie counts. The check does not change a
count. The column counter_mirna gives the entry of the miRNA counter that counts the read (its per-read table).
"""

import argparse
import gzip

import numpy as np
import pandas as pd

from quantify_mirna import arm_assignments, read_hits


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", required=True)
    p.add_argument("--hits", required=True, help="hits of the miRNA-only mapping: rank, ref, offset, length, mismatches")
    p.add_argument("--arms", required=True, help="<db>_arms.tsv: precursor, mature, start, end, entry")
    p.add_argument("--reads", required=True, help="collapsed.fa.gz (>s<rank>_x<count>)")
    p.add_argument("--classes", required=True, help="read_classes.tsv.gz: rank, count, length, class_<db>...")
    p.add_argument("--mirna-reads", required=True, help="step 5 reads of the bowtie counts: rank, ..., tied_classes")
    p.add_argument("--counter-reads", required=True, help="per-read table of the miRNA counter: sequence, mirna, reads")
    p.add_argument("--min-len", type=int, required=True)
    p.add_argument("--max-len", type=int, required=True)
    p.add_argument("--start-tol", type=int, required=True)
    p.add_argument("--end-tol", type=int, required=True)
    p.add_argument("--out-reads", required=True)
    return p.parse_args()


def read_collapsed(path, ranks):
    """rank -> sequence for the given ranks."""
    by_rank = {}
    with gzip.open(path, "rt") as fh:
        for header in fh:
            seq = next(fh).rstrip()
            rank = int(header[2 : header.index("_x")])
            if rank in ranks:
                by_rank[rank] = seq
    return by_rank


def result(cls, tied):
    return np.where(cls == "miRNA", np.where(tied, "tied", "miRNA"), cls)


def main():
    args = parse_args()
    classes = pd.read_csv(args.classes, sep="\t", usecols=["rank", "count", "length", f"class_{args.db}"], index_col="rank")
    cls = classes[f"class_{args.db}"]
    mirna_reads = pd.read_csv(args.mirna_reads, sep="\t", usecols=["rank", "tied_classes"])
    tied = set(mirna_reads.loc[mirna_reads["tied_classes"].notna(), "rank"])
    counter = pd.read_csv(args.counter_reads, sep="\t")

    # Reads with a miRNA-only hit on a mature arm; best stratum only, so one mismatch count per read
    a = arm_assignments(read_hits(args.hits), pd.read_csv(args.arms, sep="\t"), args)
    a = a.drop_duplicates(["rank", "mature"])[["rank", "mature", "mm"]]
    by_rank = read_collapsed(args.reads, set(a["rank"]))

    count = classes.loc[a["rank"], "count"].to_numpy()
    a = a.assign(
        sequence=a["rank"].map(by_rank),
        length=classes.loc[a["rank"], "length"].to_numpy(),
        count=count,
        reads=count / a.groupby("rank")["mature"].transform("size").to_numpy(),
        **{"class": cls.loc[a["rank"]].to_numpy()},
    )
    is_tied = a["rank"].isin(tied).to_numpy()
    assert not (is_tied & (a["class"] != "miRNA").to_numpy()).any(), "tied reads must have the class miRNA"
    a["result"] = result(a["class"].to_numpy(), is_tied)
    a["counter_mirna"] = a["sequence"].map(counter.groupby("sequence")["mirna"].agg(lambda m: ";".join(sorted(set(m)))))
    out = a.rename(columns={"mature": "mirna"}).sort_values(["count", "sequence", "mirna"], ascending=[False, True, True])
    cols = ["sequence", "length", "count", "mm", "mirna", "reads", "class", "result", "counter_mirna"]
    out[cols].to_csv(args.out_reads, sep="\t", index=False, compression="gzip")

    reads = a.groupby("result")["reads"].sum()
    print(f"{args.db}: {a['rank'].nunique():,} sequences, {reads.sum():,.0f} reads with a miRNA-only hit; by result:")
    print(reads.sort_values(ascending=False).round(1).to_string())


if __name__ == "__main__":
    main()
