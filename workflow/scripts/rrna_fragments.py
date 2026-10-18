"""Per-library table of the rRNA fragments (step 7).

The reads of class rRNA (quantify_mirna.py, primary database) with their positions on the single-copy rRNA reference
(bowtie, <= 1 mismatch, sense strand, all best hits). One row per read and position: a read with several equally good
positions has n_pos rows, so its weight per row is count / n_pos. `region` is the rRNA region at the read midpoint.
`contaminant_tie` marks reads of >= contaminant-min-len nt that match a contaminant genome with as few mismatches as
their rRNA position: the classification keeps them human, but they could also be contaminant (e.g. mycoplasma) rRNA.
The reads with the class rRNA from an end (quantify_mirna.py, reads without a full-length hit) have no position and are
not in the table.
"""

import argparse
import gzip

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads", required=True, help="collapsed.fa.gz (>s<rank>_x<count>)")
    p.add_argument("--classes", required=True, help="read_classes.tsv.gz: rank, count, length, class_<db>...")
    p.add_argument("--hits", required=True, help="rrna.hits.tsv.gz: rank, ref, offset, length, mismatches")
    p.add_argument("--contaminant", required=True, help="contaminant.hits.tsv.gz: as --hits, plus strand")
    p.add_argument("--ends", required=True,
                   help="unmapped_ends.hits.tsv.gz: end (<rank>_<5|3>_<end length>), index, reference, mismatches")
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
    # Reads without a full-length rRNA hit: the class rRNA from an end
    from_end = reads.index.difference(hits["rank"].unique())
    ends = pd.read_csv(args.ends, sep="\t", header=None, usecols=[0, 1], names=["end", "index"],
                       dtype={"end": str, "index": "category"})
    end_ranks = ends.loc[ends["index"] == "rrna", "end"].str.split("_").str[0].astype("int64")
    missing = from_end.difference(end_ranks)
    assert not len(missing), f"{len(missing)} rRNA reads without an rRNA hit, e.g. rank {missing[0]}"
    end_reads = reads.loc[from_end, "count"].sum()
    reads = reads.drop(from_end)
    # FASTA names are <id>|rRNA|<name>
    hits["ref"] = hits["ref"].cat.rename_categories(lambda r: r.split("|")[2])

    frag = hits.join(reads, on="rank")
    frag["n_pos"] = frag.groupby("rank")["rank"].transform("size").astype("int16")
    frag = frag.rename(columns={"offset": "start"})
    frag["end"] = frag["start"] + frag["length"]

    features = pd.read_csv(args.features, sep="\t", header=None, names=["ref", "start", "end", "name"])
    features["ref"] = features["ref"].str.split("|").str[2]
    frag["region"] = region_at(frag["ref"].astype(str).to_numpy(), (frag["start"] + frag["length"] // 2).to_numpy(), features)

    cont = pd.read_csv(args.contaminant, sep="\t", header=None, usecols=[0, 4], names=["rank", "mm"], dtype="int64")
    cont_mm = cont[cont["rank"].isin(reads.index)].groupby("rank")["mm"].min()
    frag["contaminant_tie"] = (
        (frag["length"] >= args.contaminant_min_len) & (frag["rank"].map(cont_mm).fillna(99) <= frag["mm"])
    ).astype("int8")

    frag = frag.join(read_sequences(args.reads, set(reads.index)), on="rank")
    cols = ["rank", "sequence", "count", "length", "ref", "start", "end", "region", "mm", "n_pos", "contaminant_tie"]
    frag[cols].sort_values(["rank", "ref", "start"]).to_csv(args.output, sep="\t", index=False)

    w = frag["count"] / frag["n_pos"]
    print(f"rRNA reads with a position: {reads['count'].sum():,} ({len(reads):,} sequences); rows: {len(frag):,}")
    print(f"rRNA reads with the class from an end (not in the table): {end_reads:,} ({len(from_end):,} sequences)")
    print("By region (% of rRNA reads):")
    print((100 * w.groupby(frag["region"]).sum() / w.sum()).sort_values(ascending=False).round(2).to_string())
    print(f"Equal contaminant match (reads >= {args.contaminant_min_len} nt): {100 * w[frag['contaminant_tie'] == 1].sum() / w.sum():.2f}%")


if __name__ == "__main__":
    main()
