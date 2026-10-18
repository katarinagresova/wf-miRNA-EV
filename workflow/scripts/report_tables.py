"""Write the MultiQC custom content of the QC report (step 8): the general statistics, the insert lengths, the read
classes, the contaminants, the miRNA counts and the off-target check of each library.

The read classes and the general statistics use the primary miRNA database (the first one). The miRNA counts and the
off-target check have one data set for each database.
"""

import argparse
import json
from pathlib import Path

import pandas as pd

from mirna_summary import NOT_OFF_TARGET
from quantify_mirna import ALL_CLASSES, CONTAMINANT_PRIORITY

# Off-target results, in plot order
RESULTS = ["miRNA", "tied"] + [c for c in ALL_CLASSES if c != "miRNA"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dbs", nargs="+", required=True, help="database keys, the primary one first")
    p.add_argument("--summary", required=True, help="results/library_summary.tsv")
    p.add_argument("--composition", required=True, help="results/classify/composition.tsv")
    p.add_argument("--contaminants", required=True, help="results/classify/contaminants.tsv")
    p.add_argument("--insert-lengths", required=True, help="results/trim/insert_lengths.tsv")
    p.add_argument("--off-target", nargs="+", required=True, help="results/mirna/<sample>/off_target_<db>.tsv.gz")
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def write(out_dir, name, section):
    with open(Path(out_dir) / f"{name}_mqc.json", "w") as fh:
        json.dump(section, fh, indent=1)


def records(df):
    """Data frame (libraries x columns) -> {library: {column: value}}, without missing values."""
    return {s: {k: v for k, v in row.items() if pd.notna(v)} for s, row in df.iterrows()}


def lines(df):
    """Data frame (lengths x libraries) -> {library: {length: value}}, a gap for a missing value. MultiQC 1.25 sorts
    the x values of a JSON line plot as text, so the line plots use categories: all libraries have the same lengths,
    in order."""
    return {s: {str(k): (None if pd.isna(v) else v) for k, v in col.items()} for s, col in df.items()}


def general_stats(summary, db, detected_max):
    """The main values of each library."""
    contaminant = summary.filter(items=[f"pct_{c}" for c in CONTAMINANT_PRIORITY]).sum(axis=1)
    off_target = summary.filter(regex=f"^{db}_off_target_pct_").rename(columns=lambda c: c.split("_pct_", 1)[1])
    other = off_target.drop(columns=NOT_OFF_TARGET, errors="ignore").sum(axis=1)
    detected = f"{db}_detected_ge{detected_max}"
    stats = pd.DataFrame({
        "kept_reads": summary["kept_reads"],
        "pct_miRNA": summary.get("pct_miRNA", 0.0),
        "pct_rRNA": summary.get("pct_rRNA", 0.0),
        "pct_contaminant": contaminant,
        "detected": summary[detected],
        "pct_off_target": other,
    })
    pct = {"suffix": "%", "min": 0, "max": 100, "format": "{:,.1f}"}
    headers = {
        "kept_reads": {"title": "Kept reads",
                       "description": "Reads after the adapter trim, the length filter and the 5' adapter filter",
                       "format": "{:,.0f}", "scale": "Blues"},
        "pct_miRNA": {"title": "% miRNA", "description": f"Kept reads with the class miRNA ({db})",
                      "scale": "Greens", **pct},
        "pct_rRNA": {"title": "% rRNA", "description": f"Kept reads with the class rRNA ({db})",
                     "scale": "Oranges", **pct},
        "pct_contaminant": {"title": "% contaminant", "description": "Kept reads with a contaminant class",
                            "scale": "Reds", **pct},
        "detected": {"title": f"miRNAs ≥{detected_max}", "description": f"miRNA entries with {detected_max} or more "
                     f"reads ({db})", "format": "{:,.0f}", "scale": "Purples"},
        "pct_off_target": {"title": "% off-target", "description": f"Reads of the miRNA-only mapping that step 5 "
                           f"gives to another class ({db})", "scale": "Reds", **pct},
    }
    return {
        "id": "wf_general_stats",
        "plot_type": "generalstats",
        "namespace": "wf-miRNA-EV",
        "headers": headers,
        "data": records(stats),
    }


def insert_lengths(path):
    """Insert length distribution of each library; length 0 is the adapter dimers."""
    lengths = pd.read_csv(path, sep="\t", index_col="length")
    lengths = lengths.loc[: lengths[lengths.sum(axis=1) > 0].index.max()]
    return {
        "id": "wf_insert_lengths",
        "section_name": "Insert lengths",
        "description": "Insert length after the adapter trim. Length 0 is the adapter dimers. The reads without an "
                       "adapter are not in the plot.",
        "plot_type": "linegraph",
        "pconfig": {"id": "wf_insert_lengths_plot", "title": "wf-miRNA-EV: insert lengths",
                    "xlab": "Insert length (nt)", "ylab": "% of raw reads", "ymin": 0, "categories": True},
        "data": lines(lengths),
    }


def read_classes(path, db):
    """Kept reads of each class (primary database)."""
    comp = pd.read_csv(path, sep="\t")
    comp = comp[comp["db"] == db].pivot_table(index="sample", columns="class", values="reads", aggfunc="sum")
    order = [c for c in ALL_CLASSES if c in comp.columns]
    return {
        "id": "wf_read_classes",
        "section_name": "Read classes",
        "description": f"Class of each kept read (step 5, {db}). The hit with the fewest mismatches gives the class, "
                       "and a fixed class priority breaks ties. Genome classes are reads without a small RNA hit.",
        "plot_type": "bargraph",
        "categories": order,
        "pconfig": {"id": "wf_read_classes_plot", "title": "wf-miRNA-EV: read classes", "ylab": "Reads",
                    "cpswitch_c_active": False},
        "data": records(comp[order]),
    }


def contaminants(path, libraries):
    """Kept reads of each contaminant genome."""
    section = {
        "id": "wf_contaminants",
        "section_name": "Contaminants",
        "description": "Kept reads of each contaminant genome (step 5).",
    }
    cont = pd.read_csv(path, sep="\t")
    if cont.empty:
        return {**section, "plot_type": "html", "data": "<p>No library has contaminant reads.</p>"}
    cont = cont.pivot_table(index="sample", columns="species", values="reads", aggfunc="sum")
    cont = cont.reindex(libraries).fillna(0)
    return {
        **section,
        "plot_type": "bargraph",
        "pconfig": {"id": "wf_contaminants_plot", "title": "wf-miRNA-EV: contaminants", "ylab": "Reads"},
        "data": records(cont),
    }


def mirna_counts(summary, dbs):
    """miRNA reads and detected miRNA entries of each database."""
    headers, cols = {}, []
    for db in dbs:
        for col in summary.columns:
            if not col.startswith(f"{db}_") or "_off_target_" in col:
                continue
            what = col.removeprefix(f"{db}_")
            if what == "mirna_reads":
                h = {"title": f"miRNA reads ({db})", "format": "{:,.0f}"}
            elif what == "pct_of_kept":
                h = {"title": f"% of kept ({db})", "suffix": "%", "min": 0, "max": 100, "format": "{:,.1f}"}
            elif what == "pct_tied":
                h = {"title": f"% tied ({db})", "description": "miRNA reads with an equal hit on another class",
                     "suffix": "%", "min": 0, "max": 100, "format": "{:,.1f}"}
            elif what == "pct_tied_rRNA":
                h = {"title": f"% tied rRNA ({db})", "description": "miRNA reads with an equal hit on rRNA",
                     "suffix": "%", "min": 0, "max": 100, "format": "{:,.1f}"}
            else:
                n = what.removeprefix("detected_ge")
                h = {"title": f"miRNAs ≥{n} ({db})", "description": f"miRNA entries with {n} or more reads",
                     "format": "{:,.0f}"}
            headers[col] = h
            cols.append(col)
    return {
        "id": "wf_mirna",
        "section_name": "miRNA counts",
        "description": "miRNA reads, tied miRNA reads and detected miRNA entries of each library (step 6).",
        "plot_type": "table",
        "headers": headers,
        "pconfig": {"id": "wf_mirna_table", "title": "wf-miRNA-EV: miRNA counts"},
        "data": records(summary[cols]),
    }


def off_target(paths, dbs):
    """Reads of the miRNA-only mapping: for each library and result, and the off-target percent by read length."""
    by_result, by_length = {db: {} for db in dbs}, {db: {} for db in dbs}
    for path in paths:
        p = Path(path)
        sample, db = p.parent.name, p.name.removesuffix(".tsv.gz").removeprefix("off_target_")
        t = pd.read_csv(path, sep="\t", usecols=["length", "reads", "result"])
        by_result[db][sample] = t.groupby("result")["reads"].sum().to_dict()
        other = t["reads"].where(~t["result"].isin(NOT_OFF_TARGET), 0)
        by_length[db][sample] = (100 * other.groupby(t["length"]).sum() / t.groupby("length")["reads"].sum()).round(4)
    present = {r for d in by_result.values() for s in d.values() for r in s}
    labels = [{"name": db} for db in dbs]
    results = {
        "id": "wf_off_target",
        "section_name": "miRNA off-target check",
        "description": "Reads with a hit in the miRNA-only mapping, by their result in step 5. <code>tied</code>: "
                       "miRNA, but another class has an equal hit. Another class: a better hit on another reference. "
                       "<code>unmapped</code>: no hit in step 5.",
        "plot_type": "bargraph",
        "categories": [r for r in RESULTS if r in present],
        "pconfig": {"id": "wf_off_target_plot", "title": "wf-miRNA-EV: miRNA off-target check", "ylab": "Reads",
                    "cpswitch_c_active": False, "data_labels": labels},
        "data": [by_result[db] for db in dbs],
    }
    lengths = {
        "id": "wf_off_target_length",
        "section_name": "miRNA off-target check by read length",
        "description": "Percent of the reads of the miRNA-only mapping that step 5 gives to another class, for each "
                       "read length. Short reads match a miRNA by chance more often.",
        "plot_type": "linegraph",
        "pconfig": {"id": "wf_off_target_length_plot", "title": "wf-miRNA-EV: off-target reads by read length",
                    "xlab": "Read length (nt)", "ylab": "% off-target", "ymin": 0, "categories": True,
                    "data_labels": labels},
        "data": [lines(pd.DataFrame(by_length[db]).sort_index()) for db in dbs],
    }
    return results, lengths


def main():
    args = parse_args()
    db = args.dbs[0]
    summary = pd.read_csv(args.summary, sep="\t", dtype={"sample": str}).set_index("sample")
    detected_max = max(int(c.rsplit("_ge", 1)[1]) for c in summary.columns if c.startswith(f"{db}_detected_ge"))

    write(args.out_dir, "general_stats", general_stats(summary, db, detected_max))
    write(args.out_dir, "insert_lengths", insert_lengths(args.insert_lengths))
    write(args.out_dir, "read_classes", read_classes(args.composition, db))
    write(args.out_dir, "contaminants", contaminants(args.contaminants, summary.index))
    write(args.out_dir, "mirna", mirna_counts(summary, args.dbs))
    results, lengths = off_target(args.off_target, args.dbs)
    write(args.out_dir, "off_target", results)
    write(args.out_dir, "off_target_length", lengths)


if __name__ == "__main__":
    main()
