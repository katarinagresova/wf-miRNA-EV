"""Per-library table of the rRNA fragments (step 7).

The reads of class rRNA (quantify_mirna.py, primary database) with their positions on the single-copy rRNA reference
(bowtie, <= 1 mismatch, sense strand, all best hits). One row per read and position: a read with several equally good
positions has n_pos rows, so its weight per row is count / n_pos. `region` is the rRNA region at the middle of the
position. `contaminant_tie` marks reads of >= contaminant-min-len nt that match a contaminant genome with as few
mismatches as their rRNA position: the classification keeps them human, but they could also be contaminant (e.g.
mycoplasma) rRNA.

A read with the class rRNA from an end (quantify_mirna.py, reads without a full-length hit) has the positions of its
rRNA ends that gave the class: the deciding end length, then the fewest mismatches. One row per end and position
(from_end: 5 or 3; 0 for a full-length position), so a read with one position for each end has n_pos = 2. The row is
the end only (start, end, region, mm). The weight of the table is the rRNA reads of the class.
"""

import argparse
import gzip

import numpy as np
import pandas as pd

from quantify_mirna import END_CLASS_REFS, deciding_ends, read_end_hits


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads", required=True, help="collapsed.fa.gz (>s<rank>_x<count>)")
    p.add_argument("--classes", required=True, help="read_classes.tsv.gz: rank, count, length, class_<db>...")
    p.add_argument("--hits", required=True, help="rrna.hits.tsv.gz: rank, ref, offset, length, mismatches")
    p.add_argument("--contaminant", required=True, help="contaminant.hits.tsv.gz: as --hits, plus strand")
    p.add_argument("--ends", required=True,
                   help="unmapped_ends.hits.tsv.gz: end (<rank>_<5|3>_<end length>), index, reference, mismatches, "
                        "offset")
    p.add_argument("--end-len", type=int, required=True, help="end length aligned with <= mismatches")
    p.add_argument("--features", required=True, help="rRNA regions: ref, start, end, name")
    p.add_argument("--db", required=True, help="miRNA database whose read classes are used")
    p.add_argument("--contaminant-min-len", type=int, required=True)
    p.add_argument("--output", required=True)
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
    return pd.Series(seqs, name="sequence")


def region_at(ref, pos, features):
    """Name of the region containing each position (features: ref, start, end, name)."""
    out = np.full(len(pos), "", dtype=object)
    for f in features.itertuples(index=False):
        out[(ref == f.ref) & (pos >= f.start) & (pos < f.end)] = f.name
    return out


def main():
    args = parse_args()
    classes = pd.read_csv(args.classes, sep="\t", usecols=["rank", "count", "length", f"class_{args.db}"])
    reads = classes[classes[f"class_{args.db}"] == "rRNA"].drop(columns=f"class_{args.db}").set_index("rank")

    hits = pd.read_csv(
        args.hits, sep="\t", header=None, names=["rank", "ref", "offset", "length", "mm"],
        dtype={"rank": "int64", "ref": "category", "offset": "int32", "length": "int16", "mm": "int8"},
    ).drop(columns="length")
    hits = hits[hits["rank"].isin(reads.index)]
    # Reads without a full-length rRNA hit: the class rRNA from an end. Their rRNA end hits that gave the class
    # (quantify_mirna.py): the deciding end length, then the fewest mismatches
    from_end = reads.index.difference(hits["rank"].unique())
    ends = read_end_hits(args.ends)
    ends = deciding_ends(ends[ends["index"].isin(END_CLASS_REFS).to_numpy()], args.end_len)
    ends = ends[(ends["index"] == "rrna").to_numpy() & ends["rank"].isin(from_end).to_numpy()]
    ends = ends[ends["mm"].to_numpy() == ends.groupby("rank")["mm"].transform("min").to_numpy()]
    missing = from_end.difference(ends["rank"])
    assert not len(missing), f"{len(missing)} rRNA reads without an rRNA hit, e.g. rank {missing[0]}"

    # One row per position: the full-length hits, and the end hits (the position of the end)
    hits = hits.assign(end=hits["offset"] + reads.loc[hits["rank"], "length"].to_numpy(), from_end=0)
    ends = ends.assign(end=ends["offset"] + ends["end_len"], from_end=ends["side"])
    cols = ["rank", "ref", "offset", "end", "mm", "from_end"]
    frag = pd.concat([hits[cols], ends[cols]], ignore_index=True).rename(columns={"offset": "start"})
    # FASTA names are <id>|rRNA|<name>
    frag["ref"] = frag["ref"].astype(str).astype("category").cat.rename_categories(lambda r: r.split("|")[2])
    frag = frag.join(reads, on="rank")
    frag["n_pos"] = frag.groupby("rank")["rank"].transform("size").astype("int16")

    features = pd.read_csv(args.features, sep="\t", header=None, names=["ref", "start", "end", "name"])
    features["ref"] = features["ref"].str.split("|").str[2]
    frag["region"] = region_at(frag["ref"].astype(str).to_numpy(), ((frag["start"] + frag["end"]) // 2).to_numpy(), features)

    cont = pd.read_csv(args.contaminant, sep="\t", header=None, usecols=[0, 4], names=["rank", "mm"], dtype="int64")
    cont_mm = cont[cont["rank"].isin(reads.index)].groupby("rank")["mm"].min()
    frag["contaminant_tie"] = (
        (frag["length"] >= args.contaminant_min_len) & (frag["rank"].map(cont_mm).fillna(99) <= frag["mm"])
    ).astype("int8")

    frag = frag.join(read_sequences(args.reads, set(reads.index)), on="rank")
    cols = ["rank", "sequence", "count", "length", "ref", "start", "end", "region", "mm", "n_pos", "contaminant_tie",
            "from_end"]
    frag[cols].sort_values(["rank", "ref", "start"]).to_csv(args.output, sep="\t", index=False)

    w = frag["count"] / frag["n_pos"]
    end_reads = reads.loc[from_end, "count"].sum()
    print(f"rRNA reads with a full-length position: {reads['count'].sum() - end_reads:,} "
          f"({len(reads) - len(from_end):,} sequences)")
    print(f"rRNA reads with a position from an end: {end_reads:,} ({len(from_end):,} sequences)")
    print(f"rRNA reads: {reads['count'].sum():,} ({len(reads):,} sequences); rows: {len(frag):,}")
    assert round(w.sum()) == reads["count"].sum(), f"table weight {w.sum():,.1f}, rRNA reads {reads['count'].sum():,}"
    print("By region (% of rRNA reads):")
    print((100 * w.groupby(frag["region"]).sum() / w.sum()).sort_values(ascending=False).round(2).to_string())
    print(f"Equal contaminant match (reads >= {args.contaminant_min_len} nt): {100 * w[frag['contaminant_tie'] == 1].sum() / w.sum():.2f}%")


if __name__ == "__main__":
    main()
