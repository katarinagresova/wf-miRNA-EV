"""Summarize step 4 over all samples: read-type composition, miRNA count matrices, detected miRNAs, miRge3 vs bowtie.

miRge3 merges near-identical miRNAs into one name (e.g. hsa-miR-103a-3p/107), and its bundled MirGeneDB library
is an older release than ours, with renumbered families (miR-93 is Hsa-Mir-17-P1c_5p there, Hsa-Mir-17-P4d_5p in
MirGeneDB 3.0). To compare the two methods, each bowtie mature name is mapped to the miRge3 name with the identical
mature sequence (by name if there is none), then to miRge3's merged name, and the bowtie counts are summed per name. All samples are
summarized; the `exclude` column of the sample sheet is carried along so that downstream figures can drop them.

Also summarizes the contaminant reads per genome, and annotates the most abundant sequences (top_sequences.tsv) with
their read class and a short label from their best hits (rRNA region, tRNA fragment type, contaminant gene, ...).
"""

import argparse
import re

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--samples", required=True)
    p.add_argument("--dbs", nargs="+", required=True, help="database keys, e.g. mirbase mirgenedb")
    p.add_argument("--composition", nargs="+", required=True, help="results/mirna/<sample>/composition.tsv")
    p.add_argument("--bowtie-counts", nargs="+", required=True, help="results/mirna/<sample>/<db>_counts.tsv")
    p.add_argument("--mirge3-counts", nargs="+", required=True, help="results/mirge3/<db>/<sample>/miR.Counts.csv")
    p.add_argument("--mirge3-report", nargs="+", required=True, help="results/mirge3/<db>/<sample>/annotation.report.csv")
    p.add_argument("--merges", nargs="+", required=True, help="miRge3 merges file per --dbs, same order")
    p.add_argument("--mirge3-mature", nargs="+", required=True, help="miRge3 library mature FASTA per --dbs, same order")
    p.add_argument("--mature", nargs="+", required=True, help="FASTAs of the mature sequences used for bowtie (all dbs)")
    p.add_argument("--detected-min", nargs="+", type=int, required=True)
    p.add_argument("--contaminants", nargs="+", required=True, help="results/mirna/<sample>/contaminants.tsv")
    p.add_argument("--top-hits", nargs="+", required=True, help="results/mirna/<sample>/top_hits.tsv")
    p.add_argument("--top", required=True, help="results/qc/top_sequences.tsv")
    p.add_argument("--rrna-features", required=True, help="rRNA feature BED (reference, start, end, name)")
    p.add_argument("--contaminant-features", required=True, help="contaminant feature BED (contig, start, end, type:product, 0, strand)")
    p.add_argument("--ref-fastas", nargs="+", required=True, help="small RNA reference FASTAs, for the target lengths")
    p.add_argument("--agreement-min", type=int, default=10, help="min reads in both methods for the correlation")
    p.add_argument("--out-dir", required=True, help="count matrices: <out-dir>/<method>_<db>.tsv, method = mirge3, bowtie, bowtie_mirge3names")
    p.add_argument("--out-summary", required=True)
    p.add_argument("--out-composition", required=True)
    p.add_argument("--out-agreement", required=True)
    p.add_argument("--out-tied", required=True)
    p.add_argument("--out-contaminants", required=True)
    p.add_argument("--out-top-annotation", required=True)
    return p.parse_args()


def by_sample(paths, samples):
    """Map each sample to the one path that contains its name as a path component."""
    out = {}
    for s in samples:
        hit = [p for p in paths if f"/{s}/" in p]
        assert len(hit) == 1, (s, hit)
        out[s] = hit[0]
    return out


def by_sample_db(paths, samples, dbs):
    return {db: by_sample([p for p in paths if f"/{db}/" in p or p.endswith(f"/{db}_counts.tsv")], samples) for db in dbs}


def read_merges(path):
    """Member name -> miRge3 merged name. Members carry SNP variants (.SNPA, ...) that are dropped."""
    member_of = {}
    for line in open(path):
        fields = [f for f in line.rstrip("\n").split(",") if f]
        if not fields:
            continue
        merged = fields[0]
        for m in fields[1:]:
            member_of.setdefault(re.sub(r"\.SNP[A-Z]$", "", m), merged)
    return member_of


def read_fasta(path):
    """Name (first word) -> DNA sequence."""
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.upper().replace("U", "T")
    return seqs


def mirge3_names(our_seqs, mirge3_fasta, member_of):
    """Our mature name -> miRge3 (merged) name, via identical sequence, else via the name itself."""
    by_seq = {}
    for name, seq in read_fasta(mirge3_fasta).items():
        by_seq.setdefault(seq, re.sub(r"\.SNP[A-Z]$", "", name))
    out = {}
    for name, seq in our_seqs.items():
        m = by_seq.get(seq, name)
        out[name] = member_of.get(m, m)
    return out


def ref_lengths(paths):
    """Full reference name (<id>|<class>|<name>) -> length."""
    out, name = {}, None
    for path in paths:
        for line in open(path):
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            else:
                out[name] += len(line.strip())
    return out


def fragment_type(start, end, length, half_min=29):
    """Position of a fragment on a (tRNA-like) RNA: 5'/3' half or fragment, internal, full length."""
    at5, at3 = start <= 1, end >= length - 1
    if at5 and at3:
        return "full length"
    if at5:
        return "5′ half" if end - start >= half_min else "5′ fragment"
    if at3:
        return "3′ half" if end - start >= half_min else "3′ fragment"
    return "internal"


def short_species(name):
    """Mesomycoplasma_hyorhinis_HUB_1 -> M. hyorhinis; phiX174 stays."""
    parts = name.split("_")
    return f"{parts[0][0]}. {parts[1]}" if len(parts) > 1 else name


def trna_isotype(name):
    """tRNA-Gly-GCC-2-1 -> tRNA-Gly-GCC; MT-TL1 stays."""
    m = re.match(r"^(tRNA-[A-Za-z]+-[A-Z]{3})", name)
    return m.group(1) if m else name


def join_names(names, n=2):
    names = list(dict.fromkeys(names))
    return "/".join(names[:n]) + (f" (+{len(names) - n})" if len(names) > n else "")


# Short names of rRNA features and of contaminant RNA products
FEATURE_SHORT = {"RNA5S1": "5S", "5S_5'flank": "5′ flank", "5S_3'flank": "3′ flank", "5'ETS": "5′ETS", "3'ETS/IGS": "3′ETS/IGS"}
PRODUCT_SHORT = {
    "signal recognition particle sRNA small type": "SRP RNA",
    "signal recognition particle sRNA large type": "SRP RNA",
    "transfer-messenger RNA": "tmRNA",
    "RNase P RNA component class A": "RNase P RNA",
    "RNase P RNA component class B": "RNase P RNA",
}


class TopAnnotator:
    """Short label for a sequence from its best hits (all of its final read class)."""

    def __init__(self, rrna_bed, cont_bed, lengths):
        self.rrna = pd.read_csv(rrna_bed, sep="\t", header=None, names=["ref", "start", "end", "name"])
        self.cont = pd.read_csv(cont_bed, sep="\t", header=None, names=["ref", "start", "end", "name", "score", "strand"])
        self.cont = self.cont[self.cont["name"].str.split(":").str[0].isin(["tRNA", "rRNA", "tmRNA", "ncRNA", "CDS"])]
        self.cont = {r: g for r, g in self.cont.groupby("ref")}
        self.lengths = lengths

    def rrna_label(self, h):
        """Region and 1-based start within the region, e.g. 28S:1595; '→ 3′ flank' if the read runs into the next one."""
        labels = []
        for _, x in h.iterrows():
            s, e = int(x["offset"]), int(x["offset"] + x["length"])
            f = self.rrna[self.rrna["ref"] == x["ref"]]
            start = f[(f["start"] <= s) & (f["end"] > s)]
            end = f[(f["start"] < e) & (f["end"] >= e)]["name"]
            if not len(start):
                labels.append("?")
                continue
            a = start["name"].iloc[0]
            b = end.iloc[0] if len(end) else a
            label = f"{FEATURE_SHORT.get(a, a)}:{s - start['start'].iloc[0] + 1}"
            labels.append(label + (f" → {FEATURE_SHORT.get(b, b)}" if b != a else ""))
        return join_names(labels, 1)

    def contaminant_label(self, h):
        species = join_names([short_species(r.split("|")[2]) for r in h["ref"]], 1)
        feats = []
        for _, x in h.iterrows():
            f = self.cont.get(x["ref"])
            if f is None:
                continue
            s, e = int(x["offset"]), int(x["offset"] + x["length"])
            ov = f[(f["start"] < e) & (f["end"] > s)]
            # RNA genes before CDS
            ov = ov.assign(rna=~ov["name"].str.startswith("CDS:")).sort_values("rna", ascending=False)
            if not len(ov):
                feats.append("intergenic")
                continue
            g = ov.iloc[0]
            kind, product = g["name"].split(":", 1)
            if g["strand"] != x["strand"]:
                feats.append(f"antisense {kind}")
            elif kind == "tRNA":
                rs, re_ = (s - g["start"], e - g["start"]) if g["strand"] == "+" else (g["end"] - e, g["end"] - s)
                feats.append(f"{product} {fragment_type(rs, re_, g['end'] - g['start'])}")
            elif kind in ("rRNA", "tmRNA", "ncRNA"):
                feats.append(PRODUCT_SHORT.get(product, product.replace(" ribosomal RNA", " rRNA")))
            else:
                feats.append("CDS")
        return f"{species} {join_names(feats, 1)}".strip()

    def small_label(self, h, cls):
        # Prefer genes over pseudogenes (Ensembl names ending in P<n>)
        h = h.assign(name=h["ref"].str.split("|").str[2])
        h = h.assign(pseudo=h["name"].str.contains(r"P\d+$")).sort_values(["pseudo", "name"])
        x = h.iloc[0]
        length = self.lengths.get(x["ref"], 0)
        where = fragment_type(int(x["offset"]), int(x["offset"] + x["length"]), length) if length else ""
        if cls in ("tRNA", "mt_tRNA"):
            return f"{join_names(h['name'].map(trna_isotype))} {where}"
        return f"{x['name']} {where}".replace("half", "fragment")

    def __call__(self, cls, h):
        h = h.dropna(subset=["ref"])
        if cls == "miRNA":
            return join_names(h["mature"].dropna())
        if cls == "unmapped" or not len(h):
            return cls.replace("_", " ")
        if cls in ("mycoplasma", "bacteria", "phiX"):
            return self.contaminant_label(h)
        if cls == "rRNA":
            return self.rrna_label(h)
        if cls == "miRNA_hairpin_other":
            return join_names(h["ref"].str.split("|").str[2]) + " hairpin"
        return self.small_label(h, cls)


def spearman(a, b):
    return pd.Series(a).rank().corr(pd.Series(b).rank())


def main():
    args = parse_args()
    samples = pd.read_csv(args.samples, sep="\t", dtype=str).set_index("sample")
    names = list(samples.index)
    meta = samples[["sid", "run", "cell_line", "cleanup", "rd", "analysis", "exclude"]]

    # --- composition (long), % of kept reads
    comp = pd.concat(
        pd.read_csv(p, sep="\t").assign(sample=s) for s, p in by_sample(args.composition, names).items()
    )
    kept = comp[comp["db"] == args.dbs[0]].groupby("sample")["reads"].sum()
    comp["pct_of_kept"] = 100 * comp["reads"] / comp["sample"].map(kept)
    comp = comp.merge(meta[["sid", "exclude"]], left_on="sample", right_index=True)
    comp[["sample", "sid", "exclude", "db", "class", "len_bin", "reads", "unique_seqs", "pct_of_kept"]].round(4).to_csv(
        args.out_composition, sep="\t", index=False
    )

    summary = meta.copy()
    summary["kept_reads"] = kept
    # The miRNA length bin is the shortest one (e.g. 16-28)
    mir_bin = min(comp["len_bin"].unique(), key=lambda b: int(re.split(r"[-+]", b)[0]))
    mir_range = comp[(comp["db"] == args.dbs[0]) & (comp["len_bin"] == mir_bin)]
    summary["reads_16_28"] = mir_range.groupby("sample")["reads"].sum()
    # Class shares for the primary database, all lengths
    cls_pct = comp[comp["db"] == args.dbs[0]].pivot_table(index="sample", columns="class", values="pct_of_kept", aggfunc="sum")
    for c in cls_pct.columns:
        summary[f"pct_{c}"] = cls_pct[c].fillna(0)

    bowtie_paths = by_sample_db(args.bowtie_counts, names, args.dbs)
    mirge_paths = by_sample_db(args.mirge3_counts, names, args.dbs)
    report_paths = by_sample_db(args.mirge3_report, names, args.dbs)

    our_seqs = {}
    for path in args.mature:
        our_seqs.update(read_fasta(path))

    agreement, tied_rows = [], []
    for db, merges_path, mirge3_fasta in zip(args.dbs, args.merges, args.mirge3_mature):
        member_of = read_merges(merges_path)
        to_mirge3 = mirge3_names(our_seqs, mirge3_fasta, member_of)
        mats = {"mirge3": {}, "bowtie": {}, "bowtie_mirge3names": {}}
        for s in names:
            b = pd.read_csv(bowtie_paths[db][s], sep="\t", index_col=0)
            mats["bowtie"][s] = b["reads"]
            m = pd.read_csv(mirge_paths[db][s], index_col=0).iloc[:, 0]
            mats["mirge3"][s] = m[m > 0]
            rep = pd.read_csv(report_paths[db][s]).iloc[0]
            summary.loc[s, f"mirge3_{db}_all_mirna_reads"] = rep["All miRNA Reads"]
            summary.loc[s, f"bowtie_{db}_tied_other_class_reads"] = b["reads_tied_other_class"].sum()
            summary.loc[s, f"bowtie_{db}_canonical_reads"] = b["canonical_reads"].sum()
            tied_rows.append(b[["reads", "reads_tied_other_class", "tied_classes"]].fillna("").assign(db=db, sample=s))

            # bowtie counts summed into miRge3's merged names
            bm = b["reads"].groupby(b.index.map(lambda n: to_mirge3.get(n, member_of.get(n, n)))).sum()
            mats["bowtie_mirge3names"][s] = bm
            both = pd.concat({"mirge3": m, "bowtie": bm}, axis=1).fillna(0)
            shared = both[(both >= args.agreement_min).all(axis=1)]
            det = both >= args.agreement_min
            agreement.append(
                {
                    "sample": s,
                    "db": db,
                    f"n_ge{args.agreement_min}_mirge3": int(det["mirge3"].sum()),
                    f"n_ge{args.agreement_min}_bowtie": int(det["bowtie"].sum()),
                    f"n_ge{args.agreement_min}_both": len(shared),
                    "spearman": spearman(shared["mirge3"], shared["bowtie"]) if len(shared) > 2 else np.nan,
                    "pearson_log10": np.log10(shared + 1).corr().iloc[0, 1] if len(shared) > 2 else np.nan,
                    "bowtie_over_mirge3_reads": both["bowtie"].sum() / both["mirge3"].sum(),
                }
            )

        for method, cols in mats.items():
            mat = pd.DataFrame(cols).fillna(0)[names]
            mat = mat.loc[mat.sum(axis=1).sort_values(ascending=False).index]
            mat.index.name = "mirna"
            mat.round(3).to_csv(f"{args.out_dir}/{method}_{db}.tsv", sep="\t")
            if method == "bowtie_mirge3names":
                continue
            summary[f"{method}_{db}_mirna_reads"] = mat.sum()
            summary[f"{method}_{db}_pct_of_kept"] = 100 * mat.sum() / summary["kept_reads"]
            for n in args.detected_min:
                summary[f"{method}_{db}_detected_ge{n}"] = (mat >= n).sum()

    summary.index.name = "sample"
    summary.round(4).to_csv(args.out_summary, sep="\t")
    agr = pd.DataFrame(agreement).merge(meta[["sid", "exclude"]], left_on="sample", right_index=True)
    agr.round(4).to_csv(args.out_agreement, sep="\t", index=False)

    # miRNA entries whose reads match another class equally well (candidate artefacts, e.g. rRNA-derived miRBase entries)
    tied = pd.concat(tied_rows).reset_index().rename(columns={"index": "mature"})
    tied = tied.merge(meta[["exclude"]], left_on="sample", right_index=True)
    tied = tied[tied["exclude"] != "yes"]
    out = (
        tied.groupby(["db", "mature"])
        .agg(
            reads=("reads", "sum"),
            reads_tied=("reads_tied_other_class", "sum"),
            tied_classes=("tied_classes", lambda t: ",".join(sorted(set(",".join(t).split(",")) - {""}))),
            n_samples_detected=("sample", "nunique"),
            n_samples_tied=("reads_tied_other_class", lambda r: int((r > 0).sum())),
        )
        .reset_index()
    )
    out = out[out["reads_tied"] > 0]
    out["tied_frac"] = out["reads_tied"] / out["reads"]
    out.sort_values(["db", "reads_tied"], ascending=[True, False]).round(4).to_csv(args.out_tied, sep="\t", index=False)

    # --- contaminant reads per genome, % of kept reads
    cont = pd.concat(pd.read_csv(p, sep="\t").assign(sample=s) for s, p in by_sample(args.contaminants, names).items())
    cont = cont[cont["db"] == args.dbs[0]].drop(columns="db")
    cont = cont.groupby(["sample", "class", "accession", "species"], as_index=False)[["reads", "unique_seqs"]].sum()
    cont["pct_of_kept"] = 100 * cont["reads"] / cont["sample"].map(kept)
    cont = cont.merge(meta[["sid", "exclude"]], left_on="sample", right_index=True)
    cont[["sample", "sid", "exclude", "class", "accession", "species", "reads", "unique_seqs", "pct_of_kept"]].round(4).to_csv(
        args.out_contaminants, sep="\t", index=False
    )

    # --- annotation of the most abundant sequences (the class does not depend on the sample, so any sample will do)
    hits = pd.concat(pd.read_csv(p, sep="\t").assign(src=p) for p in args.top_hits)
    hits = hits[hits["db"] == args.dbs[0]]
    annotate = TopAnnotator(args.rrna_features, args.contaminant_features, ref_lengths(args.ref_fastas))
    top = pd.read_csv(args.top, sep="\t", index_col=0)
    rows = []
    for seq, h in hits.groupby("sequence"):
        h = h[h["src"] == h["src"].iloc[0]]
        cls = h["class"].iloc[0]
        rows.append({"sequence": seq, "class": cls, "mismatches": h["mm"].min(), "label": annotate(cls, h),
                     "n_hits": h["ref"].notna().sum()})
    ann = pd.DataFrame(rows).set_index("sequence")
    missing = top.index.difference(ann.index)
    assert not len(missing), f"{len(missing)} top sequences without hits, e.g. {missing[0]}"
    ann = top[["length"]].join(ann)
    ann.to_csv(args.out_top_annotation, sep="\t")
    print(summary.T.to_string())


if __name__ == "__main__":
    main()
