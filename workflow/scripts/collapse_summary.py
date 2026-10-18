"""Summarize collapsed reads: library complexity per sample and the most abundant sequences shared across samples."""

import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--counts", nargs="+", required=True, help="counts.tsv.gz from rule collapse")
    p.add_argument("--samples", required=True, help="sample sheet")
    p.add_argument("--top-n", type=int, required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--top", required=True)
    return p.parse_args()


def complexity(counts):
    """counts: sequence counts sorted in descending order."""
    reads = int(counts["count"].sum())
    c = counts["count"]
    row = {
        "reads": reads,
        "unique_seqs": len(counts),
        "unique_pct": 100 * len(counts) / reads,
        "singletons_pct_of_unique": 100 * (c == 1).mean(),
        "singletons_pct_of_reads": 100 * c[c == 1].sum() / reads,
    }
    for n in (1, 10, 100, 1000):
        row[f"top{n}_pct_of_reads"] = 100 * c.iloc[:n].sum() / reads
    lengths = counts["sequence"].str.len()
    row["median_len_reads"] = weighted_median(lengths, c)
    row["median_len_unique"] = lengths.median()
    return row


def weighted_median(values, weights):
    df = pd.DataFrame({"v": values, "w": weights}).groupby("v")["w"].sum().sort_index()
    return float(df.index[df.cumsum().searchsorted(df.sum() / 2)])


def main():
    args = parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample").drop(columns="fastq")
    paths = {Path(p).parent.name: p for p in args.counts}

    rows, top_seqs, totals = {}, {}, {}
    for name, path in paths.items():
        counts = pd.read_csv(path, sep="\t", dtype={"sequence": str, "count": "int64"})
        rows[name] = complexity(counts)
        totals[name] = rows[name]["reads"]
        top_seqs[name] = counts.iloc[: args.top_n]

    summary = samples.join(pd.DataFrame.from_dict(rows, orient="index"), how="right")
    summary.index.name = "sample"
    summary.round(3).to_csv(args.summary, sep="\t")

    # Union of per-sample top sequences, quantified (CPM of kept reads) in every sample.
    # Second pass so sequences outside a sample's top N still get their count there.
    union = pd.Index(pd.concat(top_seqs.values())["sequence"].unique())
    cpm = {}
    for name, path in paths.items():
        found = pd.Series(0, index=union, dtype="int64")
        for chunk in pd.read_csv(path, sep="\t", dtype={"sequence": str, "count": "int64"}, chunksize=5_000_000):
            hit = chunk[chunk["sequence"].isin(union)]
            found[hit["sequence"]] = hit["count"].to_numpy()
        cpm[name] = found / totals[name] * 1e6
    top = pd.DataFrame(cpm)
    top.insert(0, "length", top.index.str.len())
    top_sets = [set(t["sequence"]) for t in top_seqs.values()]
    top.insert(1, "n_samples_in_top", [sum(s in t for t in top_sets) for s in top.index])
    top.insert(2, "mean_cpm", top[list(cpm)].mean(axis=1))
    top = top.sort_values("mean_cpm", ascending=False)
    top.index.name = "sequence"
    top.round(1).to_csv(args.top, sep="\t")


if __name__ == "__main__":
    main()
