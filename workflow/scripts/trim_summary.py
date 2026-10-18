"""Summarize adapter trimming per sample: read fates, insert length classes, 5' extra nt, insert length distribution."""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--json", nargs="+", required=True, help="cutadapt JSON reports")
    p.add_argument("--hist", nargs="+", required=True, help="hist.tsv from rule trim_hist")
    p.add_argument("--samples", required=True, help="sample sheet")
    p.add_argument("--min-len", type=int, required=True)
    p.add_argument("--max-len", type=int, required=True, help="mirna.max_len, the end of the miRNA insert length class")
    p.add_argument("--summary", required=True)
    p.add_argument("--lengths", required=True)
    p.add_argument("--plot", required=True)
    return p.parse_args()


def sample_name(path):
    # results/trim/{sample}/<file>
    return Path(path).parent.name


def read_hist(path):
    """Split hist.tsv into insert length counts per read set and the removed 5' nt counts."""
    hist = pd.read_csv(path, sep="\t", dtype={"value": str})
    lengths = hist[hist["metric"] == "length"].astype({"value": int})
    by_set = {r: g.set_index("value")["count"] for r, g in lengths.groupby("reads")}
    cut = hist[hist["metric"] == "cut_nt"].groupby("value")["count"].sum()
    empty = pd.Series(dtype=int)
    return by_set.get("adapter", empty), by_set.get("no_adapter", empty), cut


def summarize(report, ad, no, cut, min_len, max_len):
    rc, bp = report["read_counts"], report["basepair_counts"]
    raw = rc["input"]

    def n(series, lo, hi=None):
        idx = series.index
        return int(series[(idx >= lo) & ((idx <= hi) if hi is not None else True)].sum())

    counts = {
        "adapter_dimer": n(ad, 0, 0),
        f"insert_1_{min_len - 1}": n(ad, 1, min_len - 1),
        f"insert_{min_len}_{max_len}": n(ad, min_len, max_len),
        f"insert_{max_len + 1}_50": n(ad, max_len + 1, 50),
        "insert_51plus": n(ad, 51),
        # No adapter and < min_len left after NextSeq quality trimming: poly-G / dark-cycle or low-quality reads
        f"no_adapter_lt{min_len}": n(no, 0, min_len - 1),
        f"no_adapter_ge{min_len}": n(no, min_len),
    }
    assert sum(counts.values()) == raw, "read categories do not add up to the input reads"
    assert counts["adapter_dimer"] + n(ad, 1) == rc["read1_with_adapter"]

    row = {"raw_reads": raw, "kept_reads": n(ad, min_len)}
    row["kept_pct"] = 100 * row["kept_reads"] / raw
    row["with_adapter_pct"] = 100 * rc["read1_with_adapter"] / raw
    for k, v in counts.items():
        row[k] = v
        row[f"{k}_pct"] = 100 * v / raw
    row["quality_trimmed_bp_pct"] = 100 * bp["quality_trimmed"] / bp["input"]

    for nt in "ACGTN":
        row[f"cut_nt_{nt}_pct"] = 100 * cut.get(nt, 0) / raw
    return row, ad / raw * 100


def plot_lengths(lengths, samples, min_len, max_len, path):
    """One panel per batch (sample sheet column `batch`), or one panel for all libraries."""
    batch = samples["batch"].fillna("") if "batch" in samples else pd.Series("", index=samples.index)
    batches = sorted(batch.unique())
    fig, axes = plt.subplots(len(batches), 1, figsize=(9, 3 * len(batches)), sharex=True, squeeze=False)
    axes = axes[:, 0]
    for ax, b in zip(axes, batches):
        for name in batch.index[batch == b]:
            ax.plot(lengths.index, lengths[name], lw=1, label=name)
        for x in (min_len, max_len):
            ax.axvline(x, color="grey", lw=0.5, ls=":")
        ax.set_title(f"Batch {b}" if b else "All libraries")
        ax.set_ylabel("% of raw reads")
        ax.legend(fontsize=7, ncol=2)
    axes[-1].set_xlabel("Insert length after trimming (nt)")
    fig.tight_layout()
    fig.savefig(path)


def main():
    args = parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample").drop(columns="fastq")

    hists = {sample_name(p): read_hist(p) for p in args.hist}
    rows, lengths = {}, {}
    for path in args.json:
        name = sample_name(path)
        with open(path) as fh:
            report = json.load(fh)
        rows[name], lengths[name] = summarize(report, *hists[name], args.min_len, args.max_len)

    summary = samples.join(pd.DataFrame.from_dict(rows, orient="index"), how="right")
    summary.index.name = "sample"
    summary.round(3).to_csv(args.summary, sep="\t")

    lengths = pd.DataFrame(lengths).reindex(range(0, 100)).fillna(0)
    lengths.index.name = "length"
    lengths.round(4).to_csv(args.lengths, sep="\t")
    plot_lengths(lengths, samples.loc[lengths.columns], args.min_len, args.max_len, args.plot)


if __name__ == "__main__":
    main()
