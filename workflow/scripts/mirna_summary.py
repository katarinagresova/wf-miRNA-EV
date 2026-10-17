"""Summarize step 6 over all libraries: the count matrices, the detected miRNAs, the main counter against each
cross-check counter, the off-target check per database entry, and the library summary.

The step reads only the interface tables of the counters (results/mirna/<sample>/<counter>_<db>.tsv: mirna, reads,
members). Two counters can name their entries differently (for example, miRge3 merges hsa-miR-103a-3p and
hsa-miR-107 into one entry). Entries of two counters that share a member form one group. The cross-check counts are
summed into the groups; a group has the names of its main-counter entries (joined by '/' if more than one), or of
its cross-check entries if the main counter has none.

In the library summary, the miRNA columns of the main counter have no counter name, the columns of a cross-check
counter start with the counter name, and the columns of the off-target check start with the database and
off_target. The off-target fraction is the share of the reads of the miRNA-only mapping that step 5 gives to another
class (results other than miRNA, tied and unmapped).
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

# Off-target results that are not another class of step 5
NOT_OFF_TARGET = ["miRNA", "tied", "unmapped"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True)
    p.add_argument("--dbs", nargs="+", required=True, help="database keys, the primary one first")
    p.add_argument("--counter", required=True, help="main miRNA counter")
    p.add_argument("--cross-check", action="append", default=[], help="cross-check counter (repeat for more)")
    p.add_argument("--composition", required=True, help="results/classify/composition.tsv")
    p.add_argument("--counts", nargs="+", required=True, help="results/mirna/<sample>/<counter>_<db>.tsv")
    p.add_argument("--off-target", nargs="+", required=True, help="results/mirna/<sample>/off_target_<db>.tsv.gz")
    p.add_argument("--off-target-counter", nargs="+", required=True, help="results/mirna/<sample>/off_target_<db>.counter.tsv")
    p.add_argument("--detected-min", nargs="+", type=int, required=True)
    p.add_argument("--agreement-min", type=int, default=10, help="min reads in both counters for the correlation")
    p.add_argument("--out-dir", required=True, help="counts_<db>.tsv and cross_check/<counter>_<db>.tsv")
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-agreement", required=True)
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


def group_labels(main_members, cross_members):
    """Group name of each main and each cross-check entry; entries that share a member form one group."""
    parent = {}

    def find(x):
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        parent[find(a)] = find(b)

    main_of = {}
    for e, ms in main_members.items():
        find(("main", e))
        for m in ms:
            if m in main_of:
                union(("main", e), main_of[m])
            main_of[m] = ("main", e)
    for e, ms in cross_members.items():
        find(("cross", e))
        for m in ms:
            if m in main_of:
                union(("cross", e), main_of[m])
    groups = {}
    for node in list(parent):
        groups.setdefault(find(node), []).append(node)
    label = {}
    for nodes in groups.values():
        main = sorted(e for k, e in nodes if k == "main")
        name = "/".join(main) if main else "/".join(sorted(e for k, e in nodes if k == "cross"))
        for node in nodes:
            label[node] = name
    return {e: label[("main", e)] for e in main_members}, {e: label[("cross", e)] for e in cross_members}


def spearman(a, b):
    return pd.Series(a).rank().corr(pd.Series(b).rank())


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
    # The miRNA length bin is the shortest one (e.g. 16-28)
    mir_bin = min(comp["len_bin"].unique(), key=lambda b: int(re.split(r"[-+]", b)[0]))
    summary["reads_16_28"] = comp[comp["len_bin"] == mir_bin].groupby("sample")["reads"].sum()
    cls_reads = comp.pivot_table(index="sample", columns="class", values="reads", aggfunc="sum")
    for c in cls_reads.columns:
        summary[f"pct_{c}"] = (100 * cls_reads[c] / summary["kept_reads"]).fillna(0)

    # --- counter tables: tables[counter][db][sample]
    tables = {}
    for path in args.counts:
        sample, counter, db = parse_path(path, ".tsv")
        tables.setdefault(counter, {}).setdefault(db, {})[sample] = pd.read_csv(path, sep="\t", dtype={"members": str})

    def add_counts(prefix, db, mat):
        summary[f"{prefix}{db}_mirna_reads"] = mat.sum()
        summary[f"{prefix}{db}_pct_of_kept"] = 100 * mat.sum() / summary["kept_reads"]
        for n in args.detected_min:
            summary[f"{prefix}{db}_detected_ge{n}"] = (mat >= n).sum()

    agreement, main_of = [], {}
    (out_dir / "cross_check").mkdir(parents=True, exist_ok=True)
    for db in args.dbs:
        main_tables = tables[args.counter][db]
        main_reads = {s: t.set_index("mirna")["reads"] for s, t in main_tables.items()}
        mat = matrix({s: r[r > 0] for s, r in main_reads.items()}, names)
        mat.round(3).to_csv(out_dir / f"counts_{db}.tsv", sep="\t")
        add_counts("", db, mat)
        main_members = members(main_tables)
        main_of[db] = {m: e for e, ms in main_members.items() for m in ms}

        for counter in args.cross_check:
            cross_tables = tables[counter][db]
            cross_reads = {s: t.set_index("mirna")["reads"] for s, t in cross_tables.items()}
            add_counts(f"{counter}_", db, matrix({s: r[r > 0] for s, r in cross_reads.items()}, names))
            main_label, cross_label = group_labels(main_members, members(cross_tables))
            grouped = {s: r.groupby(r.index.map(cross_label)).sum() for s, r in cross_reads.items()}
            matrix(grouped, names).round(3).to_csv(out_dir / "cross_check" / f"{counter}_{db}.tsv", sep="\t")
            for s in names:
                m = main_reads[s][main_reads[s] > 0]
                m = m.groupby(m.index.map(main_label)).sum()
                both = pd.concat({"main": m, "cross_check": grouped[s]}, axis=1).fillna(0)
                det = both >= args.agreement_min
                shared = both[det.all(axis=1)]
                n = args.agreement_min
                agreement.append({
                    "sample": s,
                    "db": db,
                    "cross_check": counter,
                    f"n_ge{n}_main": int(det["main"].sum()),
                    f"n_ge{n}_cross_check": int(det["cross_check"].sum()),
                    f"n_ge{n}_both": len(shared),
                    "spearman": spearman(shared["main"], shared["cross_check"]) if len(shared) > 2 else np.nan,
                    "pearson_log10": np.log10(shared + 1).corr().iloc[0, 1] if len(shared) > 2 else np.nan,
                    "cross_check_over_main_reads": both["cross_check"].sum() / both["main"].sum(),
                })
    pd.DataFrame(agreement).round(4).to_csv(args.out_agreement, sep="\t", index=False)

    # --- off-target check: per library (summary columns) and per database entry (all libraries)
    off = []
    for path in args.off_target:
        sample, _, db = parse_path(path, ".tsv.gz")
        t = pd.read_csv(path, sep="\t", usecols=["mirna", "reads", "result"]).assign(sample=sample, db=db)
        off.append(t)
        by_result = t.groupby("result")["reads"].sum()
        summary.loc[sample, f"{db}_off_target_reads"] = by_result.sum()
        for r, v in by_result.items():
            summary.loc[sample, f"{db}_off_target_pct_{r}"] = 100 * v / by_result.sum()
    for path in args.off_target_counter:
        sample, _, db = parse_path(path, ".counter.tsv")
        c = pd.read_csv(path, sep="\t")
        if not len(c):
            continue
        summary.loc[sample, f"{db}_off_target_counter_reads"] = c["reads"].sum()
        other = c.loc[~c["result"].isin(NOT_OFF_TARGET), "reads"].sum()
        summary.loc[sample, f"{db}_off_target_counter_pct_other_class"] = 100 * other / c["reads"].sum()
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
    entry.insert(2, "counter_mirna", [main_of[db].get(m) for db, m in zip(entry["db"], entry["mirna"])])
    entry = entry.sort_values(["db", "reads"], ascending=[True, False])
    entry.round(4).to_csv(args.out_off_target, sep="\t", index=False)

    summary.index.name = "sample"
    summary.round(4).to_csv(args.out_summary, sep="\t")
    print(summary.T.to_string())


if __name__ == "__main__":
    main()
