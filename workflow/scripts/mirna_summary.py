"""Summarize step 6 over all libraries: the count matrices, the detected miRNAs, the off-target check per database
entry, and the library summary.

The step reads only the interface tables of the miRNA counter (results/mirna/<sample>/<counter>_<db>.tsv: mirna,
reads, members). The members link the database entries of the off-target check to the entries of the counter (for
example, hsa-miR-199b-3p to hsa-miR-199a-3p/hsa-miR-199b-3p).

In the library summary, the miRNA columns have no counter name, and the columns of the off-target check start with
the database and off_target. The off-target fraction is the share of the reads of the miRNA-only mapping that step 5
gives to another class (results other than miRNA, tied and unmapped). The tied percent is the share of the miRNA reads
with the result tied: step 5 gives them the class miRNA, but another class has an equal hit. The tied rRNA percent
counts only the tied reads with rRNA in their tied classes.
"""

import argparse
import re
from pathlib import Path

import pandas as pd

# Off-target results that are not another class of step 5
NOT_OFF_TARGET = ["miRNA", "tied", "unmapped"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True)
    p.add_argument("--dbs", nargs="+", required=True, help="database keys, the primary one first")
    p.add_argument("--composition", required=True, help="results/classify/composition.tsv")
    p.add_argument("--counts", nargs="+", required=True, help="results/mirna/<sample>/<counter>_<db>.tsv")
    p.add_argument("--off-target", nargs="+", required=True, help="results/mirna/<sample>/off_target_<db>.tsv.gz")
    p.add_argument("--detected-min", nargs="+", type=int, required=True)
    p.add_argument("--out-dir", required=True, help="counts_<db>.tsv")
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-off-target", required=True)
    return p.parse_args()


def parse_path(path, suffix):
    """results/mirna/<sample>/<prefix>_<db><suffix> -> sample, prefix, db."""
    p = Path(path)
    prefix, db = p.name.removesuffix(suffix).rsplit("_", 1)
    return p.parent.name, prefix, db


def matrix(cols, names):
    """Count matrix (entries x libraries), sorted by total reads."""
    mat = pd.DataFrame(cols).reindex(columns=names).fillna(0)
    mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index]
    mat.index.name = "mirna"
    return mat


def members(tables):
    """Entry -> its members, over all libraries."""
    t = pd.concat(tables.values(), ignore_index=True).drop_duplicates("mirna").set_index("mirna")["members"].fillna("")
    return {e: [m for m in ms.split(";") if m] for e, ms in t.items()}


def main():
    args = parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample").drop(columns="fastq")
    names = list(samples.index)
    out_dir = Path(args.out_dir)

    # --- kept reads and class shares (primary database, all lengths) from step 5
    comp = pd.read_csv(args.composition, sep="\t")
    comp = comp[comp["db"] == args.dbs[0]]
    summary = samples.copy()
    summary["kept_reads"] = comp.groupby("sample")["reads"].sum()
    # The miRNA length bin is the shortest one (mirna.min_len-mirna.max_len, e.g. 16-28 -> reads_16_28)
    mir_bin = min(comp["len_bin"].unique(), key=lambda b: int(re.split(r"[-+]", b)[0]))
    summary[f"reads_{mir_bin.replace('-', '_')}"] = comp[comp["len_bin"] == mir_bin].groupby("sample")["reads"].sum()
    cls_reads = comp.pivot_table(index="sample", columns="class", values="reads", aggfunc="sum")
    for c in cls_reads.columns:
        summary[f"pct_{c}"] = (100 * cls_reads[c] / summary["kept_reads"]).fillna(0)

    # --- counter tables: tables[db][sample]
    tables = {}
    for path in args.counts:
        sample, _, db = parse_path(path, ".tsv")
        tables.setdefault(db, {})[sample] = pd.read_csv(path, sep="\t", dtype={"members": str})

    entry_of = {}
    out_dir.mkdir(parents=True, exist_ok=True)
    for db in args.dbs:
        reads = {s: t.set_index("mirna")["reads"] for s, t in tables[db].items()}
        mat = matrix({s: r[r > 0] for s, r in reads.items()}, names)
        mat.round(3).to_csv(out_dir / f"counts_{db}.tsv", sep="\t")
        summary[f"{db}_mirna_reads"] = mat.sum()
        summary[f"{db}_pct_of_kept"] = 100 * mat.sum() / summary["kept_reads"]
        # From the off-target check below; 0 if the library has no miRNA reads
        summary[f"{db}_pct_tied"] = 0.0
        summary[f"{db}_pct_tied_rRNA"] = 0.0
        for n in args.detected_min:
            summary[f"{db}_detected_ge{n}"] = (mat >= n).sum()
        entry_of[db] = {m: e for e, ms in members(tables[db]).items() for m in ms}

    # --- off-target check: per library (summary columns) and per database entry (all libraries)
    off = []
    for path in args.off_target:
        sample, _, db = parse_path(path, ".tsv.gz")
        t = pd.read_csv(path, sep="\t", usecols=["mirna", "reads", "result", "tied_classes"],
                        dtype={"tied_classes": str}).assign(sample=sample, db=db)
        off.append(t)
        tied = t[t["result"] == "tied"]
        with_rrna = ["rRNA" in c.split(",") for c in tied["tied_classes"]]
        mirna_reads = summary.loc[sample, f"{db}_mirna_reads"]
        if mirna_reads > 0:
            summary.loc[sample, f"{db}_pct_tied"] = 100 * tied["reads"].sum() / mirna_reads
            summary.loc[sample, f"{db}_pct_tied_rRNA"] = 100 * tied.loc[with_rrna, "reads"].sum() / mirna_reads
        by_result = t.groupby("result")["reads"].sum()
        summary.loc[sample, f"{db}_off_target_reads"] = by_result.sum()
        for r, v in by_result.items():
            summary.loc[sample, f"{db}_off_target_pct_{r}"] = 100 * v / by_result.sum()
    pct_cols = [c for c in summary.columns if "_off_target_pct_" in c]
    summary[pct_cols] = summary[pct_cols].fillna(0)

    off = pd.concat(off, ignore_index=True)
    entry = off.pivot_table(index=["db", "mirna"], columns="result", values="reads", aggfunc="sum", fill_value=0)
    entry.columns = [f"reads_{r}" for r in entry.columns]
    for r in ["miRNA", "tied"]:
        if f"reads_{r}" not in entry:
            entry[f"reads_{r}"] = 0.0
    entry.insert(0, "reads", entry.sum(axis=1))
    entry.insert(1, "n_samples", off.groupby(["db", "mirna"])["sample"].nunique())
    is_other = [c for c in entry.columns if c.startswith("reads_") and c.removeprefix("reads_") not in NOT_OFF_TARGET]
    entry.insert(2, "tied_frac", entry["reads_tied"] / (entry["reads_miRNA"] + entry["reads_tied"]))
    entry.insert(3, "off_target_frac", entry[is_other].sum(axis=1) / entry["reads"])
    entry = entry.reset_index()
    entry.insert(2, "counter_mirna", [entry_of[db].get(m) for db, m in zip(entry["db"], entry["mirna"])])
    entry = entry.sort_values(["db", "reads"], ascending=[True, False])
    entry.round(4).to_csv(args.out_off_target, sep="\t", index=False)

    summary.index.name = "sample"
    summary.round(4).to_csv(args.out_summary, sep="\t")
    print(summary.T.to_string())


if __name__ == "__main__":
    main()
