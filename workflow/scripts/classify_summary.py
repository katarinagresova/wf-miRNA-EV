"""Summarize step 5 over all libraries: the read classes (with the percent of the kept and of the human reads), the
contaminant reads per genome, and the most abundant sequences with their class.

Each sequence of top_sequences.tsv gets its read class and a short label from its best hits (rRNA region, tRNA fragment
type, contaminant gene, ...). Human reads are the kept reads minus the contaminant reads.
"""

import argparse
import re

import pandas as pd

from quantify_mirna import CONTAMINANT_PRIORITY


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dbs", nargs="+", required=True, help="miRNA database keys, the primary one first")
    p.add_argument("--composition", nargs="+", required=True, help="results/classify/<sample>/composition.tsv")
    p.add_argument("--contaminants", nargs="+", required=True, help="results/classify/<sample>/contaminants.tsv")
    p.add_argument("--top-hits", nargs="+", required=True, help="results/classify/<sample>/top_hits.tsv")
    p.add_argument("--top", required=True, help="results/collapse/top_sequences.tsv")
    p.add_argument("--rrna-features", required=True, help="rRNA feature BED (reference, start, end, name)")
    p.add_argument("--contaminant-features", required=True, help="contaminant feature BED (contig, start, end, type:product, 0, strand)")
    p.add_argument("--ref-fastas", nargs="+", required=True, help="small RNA reference FASTAs, for the target lengths")
    p.add_argument("--out-composition", required=True)
    p.add_argument("--out-contaminants", required=True)
    p.add_argument("--out-top", required=True)
    return p.parse_args()


def sample_name(path):
    # results/classify/{sample}/<file>
    return path.rstrip("/").split("/")[-2]


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


def main():
    args = parse_args()

    # --- composition (long), % of kept and of human reads
    comp = pd.concat((pd.read_csv(p, sep="\t").assign(sample=sample_name(p)) for p in args.composition), ignore_index=True)
    kept = comp[comp["db"] == args.dbs[0]].groupby("sample")["reads"].sum()
    comp["pct_of_kept"] = 100 * comp["reads"] / comp["sample"].map(kept)
    is_cont = comp["class"].isin(CONTAMINANT_PRIORITY)
    human = comp[~is_cont].groupby(["sample", "db"])["reads"].sum()
    human_reads = human.reindex(pd.MultiIndex.from_frame(comp[["sample", "db"]])).to_numpy()
    comp["pct_of_human"] = (100 * comp["reads"] / human_reads).where(~is_cont)
    comp[["sample", "db", "class", "len_bin", "reads", "unique_seqs", "pct_of_kept", "pct_of_human"]].round(4).to_csv(
        args.out_composition, sep="\t", index=False
    )

    # --- contaminant reads per genome, % of kept reads
    cont = pd.concat(pd.read_csv(p, sep="\t").assign(sample=sample_name(p)) for p in args.contaminants)
    cont = cont[cont["db"] == args.dbs[0]].drop(columns="db")
    cont = cont.groupby(["sample", "class", "accession", "species"], as_index=False)[["reads", "unique_seqs"]].sum()
    cont["pct_of_kept"] = 100 * cont["reads"] / cont["sample"].map(kept)
    cont[["sample", "class", "accession", "species", "reads", "unique_seqs", "pct_of_kept"]].round(4).to_csv(
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
    ann.to_csv(args.out_top, sep="\t")


if __name__ == "__main__":
    main()
