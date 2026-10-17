"""Assign every collapsed read to one RNA class and count mature miRNAs (bowtie1 cross-check), once per miRNA database.

Each read was aligned independently to every small RNA reference (best stratum per reference). A read is assigned
to the class of its fewest-mismatch hit over all references; ties are broken by CLASS_PRIORITY. Reads without any
small RNA hit were aligned to the genome and are classed by the same-strand annotation they overlap.

All reads were also aligned to contaminant genomes (both strands). A read of at least contaminant-min-len nt is
assigned to a contaminant class if it matches a contaminant with fewer mismatches than any human reference (small RNA
or genome); ties stay human. Reads with weak human evidence need less: at least contaminant-min-len nt, a contaminant
hit as good as the best human one (ties go to the contaminant); shorter, an exact contaminant match. Weak evidence is
no small RNA hit and genome hits only in genes, introns or intergenic regions (WEAK_HUMAN), or no human hit at all.
Short reads match the genome by chance (a 16-mer has ~70 expected hits with <= 1 mismatch, and the class priority picks
any exon among them), so a short contaminant read would otherwise land in one of these classes. A tie is ~300x more
likely for a contaminant read than for a human read (the genome sizes).

Reads longer than end-len-exact nt without any full-length hit (human or contaminant) were aligned by their ends. Such
a read is a contaminant if an end-len-nt end matches a contaminant (<= the cascade's mismatches) and no end-len-nt end
matches a human reference; or, when no end-len-nt end matches anything, if an end-len-exact-nt end matches a
contaminant exactly and none matches a human reference exactly. These are contaminant reads with more errors than the
full-length alignment allows (or an indel), and contaminant fragments joined to another fragment by the ligation.

A hairpin hit counts for a mature arm if the read is min-len..max-len nt and its 5'/3' ends lie within
start-tol/end-tol nt of the arm's ends. A read assigned to several mature names (paralogous hairpins, identical
arms) is split evenly between them.
"""

import argparse
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from pandas.api.types import union_categoricals

# Small RNA classes, in priority order for reads with equally good hits on several references
CLASS_PRIORITY = [
    "miRNA",
    "rRNA",
    "tRNA",
    "mt_tRNA",
    "Y_RNA",
    "snoRNA",
    "snRNA",
    "7SK",
    "7SL_SRP",
    "vault_RNA",
    "other_ncRNA",
    "miRNA_hairpin_other",
]
# Contaminant classes (names in the contaminant FASTA), in priority order for reads matching several genomes
CONTAMINANT_PRIORITY = ["mycoplasma", "bacteria", "phiX"]
# Genome classes (reads without a small RNA hit), in priority order
GENOME_PRIORITY = [
    "rRNA_genomic",
    "mt_other",
    "tRNA_genomic",
    "mRNA",
    "lncRNA",
    "other_gene",
    "intron",
    "intergenic",
]
GENOME_CLASS = {
    "rDNA_locus": "rRNA_genomic",
    "rRNA": "rRNA_genomic",
    "rRNA_pseudogene": "rRNA_genomic",
    "Mt_rRNA": "rRNA_genomic",
    "mt_genome": "mt_other",
    "tRNA_locus": "tRNA_genomic",
    "protein_coding": "mRNA",
    "lncRNA": "lncRNA",
    "intron": "intron",
    ".": "intergenic",
}
ALL_CLASSES = CLASS_PRIORITY + GENOME_PRIORITY + CONTAMINANT_PRIORITY + ["unmapped"]
# Classes with weak human evidence: ties with a contaminant go to the contaminant, and short reads can be contaminants.
# Genome classes that cover most of the genome, where a read lands by chance; not the rDNA copies, tRNA loci or the
# mt genome, which a random read rarely hits and a conserved contaminant rRNA/tRNA read may really match.
WEAK_HUMAN = ["mRNA", "lncRNA", "other_gene", "intron", "intergenic", "unmapped"]
# Human references of the end hits, in priority order (the human reference reported for a read with a human end)
END_HUMAN_PRIORITY = ["rrna", "trna", "ncrna", "mirbase", "mirgenedb", "genome"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads", required=True, help="collapsed.fa.gz (>s<rank>_x<count>)")
    p.add_argument("--hits", nargs="+", required=True, help="<ref>.hits.tsv.gz: rank, ref, offset, length, mismatches")
    p.add_argument("--contaminant", required=True, help="contaminant.hits.tsv.gz: as --hits, plus strand")
    p.add_argument("--genome", required=True, help="genome.classes.tsv.gz: rank, mismatches, classes")
    p.add_argument("--ends", required=True,
                   help="unmapped_ends.hits.tsv.gz: end (<rank>_<5|3>_<end length>), index, reference, mismatches")
    p.add_argument("--arms", nargs="+", required=True, help="<db>_arms.tsv: precursor, mature, start, end")
    p.add_argument("--min-len", type=int, required=True)
    p.add_argument("--max-len", type=int, required=True)
    p.add_argument("--start-tol", type=int, required=True)
    p.add_argument("--end-tol", type=int, required=True)
    p.add_argument("--contaminant-min-len", type=int, required=True)
    p.add_argument("--end-len", type=int, required=True, help="end length aligned with <= mismatches")
    p.add_argument("--end-len-exact", type=int, required=True, help="end length that must match exactly")
    p.add_argument("--top-n", type=int, required=True, help="most abundant sequences whose best hits go to --out-top")
    p.add_argument("--out-reads", required=True)
    p.add_argument("--out-composition", required=True)
    p.add_argument("--out-counts", nargs="+", required=True, help="<db>_counts.tsv, one per --arms")
    p.add_argument("--out-contaminants", required=True, help="reads per contaminant genome and length bin")
    p.add_argument("--out-contaminant-rules", required=True,
                   help="contaminant reads by the class they would have without contaminants, the rule and the read length")
    p.add_argument("--out-ends", required=True, help="reads tested by their ends, by outcome and read length")
    p.add_argument("--out-top", required=True, help="best hits of the --top-n most abundant sequences")
    return p.parse_args()


def ref_name(path, suffix):
    return Path(path).name.removesuffix(suffix)


def read_collapsed(path, top_n):
    """rank -> count, length; and the sequences of the top_n most abundant ranks."""
    ranks, counts, lengths, top = [], [], [], {}
    with gzip.open(path, "rt") as fh:
        for header in fh:
            seq = next(fh)
            name = header[2:].rstrip()
            rank, count = name.split("_x")
            ranks.append(int(rank))
            counts.append(int(count))
            lengths.append(len(seq) - 1)
            if int(rank) <= top_n:
                top[int(rank)] = seq.rstrip()
    return pd.DataFrame({"count": counts, "length": lengths}, index=pd.Index(ranks, name="rank")), top


def read_hits(path, strand=False):
    names = ["rank", "ref", "offset", "length", "mm"] + (["strand"] if strand else [])
    hits = pd.read_csv(
        path, sep="\t", header=None, names=names,
        dtype={"rank": "int64", "ref": "category", "offset": "int32", "length": "int16", "mm": "int8", "strand": "category"},
    )
    if not strand:
        hits["strand"] = pd.Categorical(["+"] * len(hits))
    # Reference names are <id>|<class>|<name>; parse once per distinct reference
    refs = hits["ref"].cat.categories
    codes = hits["ref"].cat.codes.to_numpy()
    hits["target"] = pd.Categorical(refs.str.split("|").str[0].to_numpy()[codes])
    hits["cls"] = refs.str.split("|").str[1].to_numpy()[codes]
    return hits


def concat_hits(frames):
    """Concatenate hit tables, keeping the categorical columns categorical."""
    frames = [f.copy() for f in frames]
    for c in ["ref", "target", "strand"]:
        cats = union_categoricals([f[c] for f in frames if c in f]).categories
        for f in frames:
            if c in f:
                f[c] = f[c].cat.set_categories(cats)
    return pd.concat(frames, ignore_index=True)


def arm_assignments(hits, arms, args):
    """Hairpin hits that match a mature arm: rank, mm, mature."""
    m = hits.merge(arms, left_on="target", right_on="precursor")
    ok = (
        m["length"].between(args.min_len, args.max_len)
        & ((m["offset"] - m["start"]).abs() <= args.start_tol)
        & ((m["offset"] + m["length"] - m["end"]).abs() <= args.end_tol)
    )
    m = m[ok]
    m = m.assign(canonical=(m["offset"] == m["start"]) & (m["offset"] + m["length"] == m["end"]) & (m["mm"] == 0))
    return m[["rank", "target", "offset", "mm", "mature", "canonical"]]


def best_class(cand, priority):
    """Per read: fewest mismatches, then class priority. cand: rank, mm, cls."""
    order = {c: i for i, c in enumerate(priority)}
    cand = cand.assign(prio=cand["cls"].map(order).astype("int16"))
    cand = cand.sort_values(["rank", "mm", "prio"])
    best = cand.drop_duplicates("rank").set_index("rank")
    return best, cand


def genome_classes(path):
    g = pd.read_csv(path, sep="\t", header=None, names=["rank", "mm", "classes"], dtype={"rank": "int64", "mm": "int8"})
    g = g.assign(cls=g["classes"].str.split(",")).explode("cls")
    g["cls"] = g["cls"].map(lambda c: GENOME_CLASS.get(c, "other_gene"))
    # A hit overlapping both an exon and the gene body is exonic: the priority order takes care of that
    best, _ = best_class(g[["rank", "mm", "cls"]], GENOME_PRIORITY)
    return best


def end_evidence(path, end_len, exact_len):
    """Per read with an end hit: the deciding end length, the outcome, the best contaminant class and human reference;
    and the contaminant candidates (best stratum over the deciding ends).

    The end_len-nt ends decide when any of them has a hit; otherwise the exact hits of the exact_len-nt ends."""
    e = pd.read_csv(path, sep="\t", header=None, names=["end", "index", "ref", "mm"],
                    dtype={"end": str, "index": "category", "ref": "category", "mm": "int8"})
    parts = e["end"].str.split("_", expand=True)
    e = e.assign(rank=parts[0].astype("int64"), end_len=parts[2].astype("int16")).drop(columns="end")
    e = e[(e["end_len"] == end_len) | (e["mm"] == 0)]
    e = e[e["end_len"].to_numpy() == e.groupby("rank")["end_len"].transform("max").to_numpy()]

    is_cont = (e["index"] == "contaminant").to_numpy()
    human = e[~is_cont]
    prio = human["index"].map({r: i for i, r in enumerate(END_HUMAN_PRIORITY)}).astype("int16")
    human_best = human.assign(prio=prio).sort_values(["rank", "mm", "prio"]).drop_duplicates("rank").set_index("rank")
    cont = e[is_cont]
    refs = cont["ref"].cat.categories
    cont = cont.assign(cls=refs.str.split("|").str[1].to_numpy()[cont["ref"].cat.codes.to_numpy()])
    cont_best, cont_cand = best_class(cont[["rank", "ref", "mm", "cls"]], CONTAMINANT_PRIORITY)
    cont_cand = cont_cand[cont_cand["mm"].to_numpy() == cont_best.loc[cont_cand["rank"], "mm"].to_numpy()]

    ev = e.groupby("rank")["end_len"].first().to_frame()
    ev["contaminant"] = cont_best["cls"]
    ev["human_ref"] = human_best["index"].astype(str)
    c, h = ev["contaminant"].notna().to_numpy(), ev["human_ref"].notna().to_numpy()
    ev["outcome"] = np.select([c & ~h, h & ~c], ["contaminant_end", "human_end"], default="both_ends")
    return ev, cont_cand


def length_bin(length, min_len, max_len):
    bins = [min_len - 1, max_len, 50, np.inf]
    labels = [f"{min_len}-{max_len}", f"{max_len + 1}-50", "51+"]
    return pd.cut(length, bins=bins, labels=labels)


def main():
    args = parse_args()
    reads, top_seqs = read_collapsed(args.reads, args.top_n)
    print(f"{len(reads):,} unique sequences, {reads['count'].sum():,} reads")

    hits = {ref_name(p, ".hits.tsv.gz"): read_hits(p) for p in args.hits}
    arms = {ref_name(p, "_arms.tsv"): pd.read_csv(p, sep="\t") for p in args.arms}
    cols = ["rank", "ref", "target", "offset", "length", "strand", "mm", "cls"]
    other = concat_hits([h[cols] for ref, h in hits.items() if ref not in arms])
    genome = genome_classes(args.genome)

    # Contaminant hits of reads long enough, and exact hits of shorter reads; best stratum per read, then class priority
    cont = read_hits(args.contaminant, strand=True)
    no_hit = reads["length"].to_numpy() > args.end_len_exact
    for ranks in [h["rank"] for h in hits.values()] + [genome.index, cont["rank"]]:
        no_hit &= ~reads.index.isin(ranks)
    cont_long = reads.loc[cont["rank"], "length"].to_numpy() >= args.contaminant_min_len
    cont = cont[cont_long | (cont["mm"].to_numpy() == 0)][cols]
    cont_best, cont_cand = best_class(cont, CONTAMINANT_PRIORITY)
    cont_cand = cont_cand[cont_cand["mm"].to_numpy() == cont_best.loc[cont_cand["rank"], "mm"].to_numpy()]
    cont_best_long = reads.loc[cont_best.index, "length"].to_numpy() >= args.contaminant_min_len

    # Reads tested by their ends: longer than the exact end length and without a hit on any reference
    ends, end_cand = end_evidence(args.ends, args.end_len, args.end_len_exact)
    tested = reads.index[no_hit]
    assert ends.index.isin(tested).all()

    composition, read_classes, contaminants, rules, top_hits = [], {}, [], [], []
    for (db, db_arms), out_counts in zip(arms.items(), args.out_counts):
        mir = hits[db]
        assigned = arm_assignments(mir, db_arms, args)
        # Hairpin hits that match an arm are miRNA; other hairpin hits (loop, long or offset reads) are not
        key = ["rank", "target", "offset"]
        on_arm = pd.MultiIndex.from_frame(mir[key]).isin(pd.MultiIndex.from_frame(assigned[key]))
        mir_cand = mir[cols].assign(cls=np.where(on_arm, "miRNA", "miRNA_hairpin_other"))
        best, cand = best_class(concat_hits([mir_cand, other]), CLASS_PRIORITY)

        cls = pd.Series("unmapped", index=reads.index, dtype=object)
        cls[genome.index] = genome["cls"]
        cls[best.index] = best["cls"]
        # Contaminant with fewer mismatches than the best human hit; with weak human evidence, also ties (long reads)
        # and exact matches (short reads, which have only exact contaminant hits here)
        human_cls = cls.copy()
        human_mm = pd.Series(np.inf, index=reads.index)
        human_mm[genome.index] = genome["mm"]
        human_mm[best.index] = best["mm"]
        c_mm, h_mm = cont_best["mm"].to_numpy(), human_mm[cont_best.index].to_numpy()
        weak = human_cls[cont_best.index].isin(WEAK_HUMAN).to_numpy()
        rule = np.select(
            [cont_best_long & (c_mm < h_mm), cont_best_long & weak & (c_mm == h_mm), ~cont_best_long & weak],
            ["fewer_mm", "tie_weak", "exact_short"], default="",
        )
        wins = cont_best.index[rule != ""]
        cls[wins] = cont_best.loc[wins, "cls"]
        # Reads without any full-length hit: a contaminant end and no human end
        end_wins = ends.index[ends["outcome"].eq("contaminant_end").to_numpy() & cls[ends.index].eq("unmapped").to_numpy()]
        cls[end_wins] = ends.loc[end_wins, "contaminant"]
        read_classes[db] = cls
        moved = wins.append(end_wins)
        end_rule = ("end_" + ends.loc[end_wins, "end_len"].astype(str)).to_numpy()
        r = pd.DataFrame({"human_class": human_cls[moved].to_numpy(), "contaminant_class": cls[moved].to_numpy(),
                          "rule": np.concatenate([rule[rule != ""], end_rule]), "length": reads.loc[moved, "length"].to_numpy(),
                          "count": reads.loc[moved, "count"].to_numpy()})
        rules.append(r.groupby(["human_class", "contaminant_class", "rule", "length"])
                     .agg(reads=("count", "sum"), unique_seqs=("count", "size")).reset_index().assign(db=db))
        print(f"{db}: {reads.loc[moved, 'count'].sum():,} contaminant reads ({len(moved):,} sequences); by rule:")
        print(r.groupby("rule")["count"].sum().to_string())
        best = best[cls[best.index].to_numpy() == best["cls"].to_numpy()]

        # Reads per contaminant genome, split evenly between equally good genomes
        cc = pd.concat([cont_cand.loc[cont_cand["rank"].isin(wins), ["rank", "ref", "cls"]],
                        end_cand.loc[end_cand["rank"].isin(end_wins), ["rank", "ref", "cls"]].astype({"ref": str})])
        cc = cc.astype({"ref": str}).drop_duplicates(["rank", "ref"])
        cc = cc.assign(
            reads=reads.loc[cc["rank"], "count"].to_numpy() / cc.groupby("rank")["ref"].transform("size").to_numpy(),
            len_bin=length_bin(reads.loc[cc["rank"], "length"].to_numpy(), args.min_len, args.max_len),
        )
        contaminants.append(
            cc.groupby(["cls", "ref", "len_bin"], observed=True)
            .agg(reads=("reads", "sum"), unique_seqs=("rank", "nunique"))
            .reset_index()
            .assign(db=db)
        )

        # Best hits of the most abundant sequences: the hits of their final class, at the best mismatch count
        top_ranks = [r for r in top_seqs if r in reads.index]
        t = concat_hits([cand[cand["rank"].isin(top_ranks)], cont_cand[cont_cand["rank"].isin(top_ranks)]])
        t = t[t["cls"].to_numpy() == cls[t["rank"]].to_numpy()].drop(columns="prio")
        t = t[t["mm"] == t.groupby("rank")["mm"].transform("min")]
        t = t.merge(assigned[["rank", "target", "offset", "mature"]].drop_duplicates(), how="left", on=["rank", "target", "offset"]) if len(t) else t
        top = pd.DataFrame({"rank": top_ranks})
        top = top.assign(
            sequence=top["rank"].map(top_seqs), count=reads.loc[top_ranks, "count"].to_numpy(),
            cls=cls[top_ranks].to_numpy(), db=db,
        )
        top_hits.append(top.merge(t.drop(columns="cls"), how="left", on="rank"))

        # Mature counts from reads whose best class is miRNA, using only their best-stratum arm hits
        best_mm = best["mm"]
        a = assigned[assigned["rank"].map(best["cls"]).eq("miRNA")]
        a = a[a["mm"].to_numpy() == best_mm.loc[a["rank"]].to_numpy()]
        a = a.groupby(["rank", "mature"], as_index=False)["canonical"].any()
        a["reads"] = reads.loc[a["rank"], "count"].to_numpy() / a.groupby("rank")["mature"].transform("size")
        # Other classes hit equally well by the same read (e.g. miRBase entries derived from rRNA or tRNA)
        c = concat_hits([cand, cont_cand])
        c = c[c["rank"].isin(a["rank"]) & ~c["cls"].isin(["miRNA", "miRNA_hairpin_other"])]
        c = c[c["mm"].to_numpy() == best_mm.loc[c["rank"]].to_numpy()]
        tied = c.drop_duplicates(["rank", "cls"]).groupby("rank")["cls"].agg(lambda s: ",".join(sorted(s)))
        a["tied"] = a["rank"].map(tied)
        a["canonical_reads"] = a["reads"] * a["canonical"]
        a["tied_reads"] = a["reads"] * a["tied"].notna()
        counts = a.groupby("mature").agg(
            reads=("reads", "sum"),
            unique_seqs=("rank", "size"),
            canonical_reads=("canonical_reads", "sum"),
            reads_tied_other_class=("tied_reads", "sum"),
            tied_classes=("tied", lambda t: ",".join(sorted(set(",".join(t.dropna()).split(",")) - {""}))),
        ).sort_values("reads", ascending=False)
        counts.round(3).to_csv(out_counts, sep="\t")
        print(f"{db}: {counts['reads'].sum():,.0f} miRNA reads, {len(counts)} mature miRNAs with >= 1 read")

        comp = reads.assign(cls=cls, len_bin=length_bin(reads["length"], args.min_len, args.max_len))
        comp = comp.groupby(["cls", "len_bin"], observed=True).agg(reads=("count", "sum"), unique_seqs=("count", "size"))
        composition.append(comp.reset_index().assign(db=db))

    comp = pd.concat(composition)
    comp["cls"] = pd.Categorical(comp["cls"], categories=ALL_CLASSES, ordered=True)
    comp = comp.sort_values(["db", "cls", "len_bin"])[["db", "cls", "len_bin", "reads", "unique_seqs"]]
    comp.rename(columns={"cls": "class"}).to_csv(args.out_composition, sep="\t", index=False)
    assert all(comp.groupby("db")["reads"].sum() == reads["count"].sum())

    out = reads.assign(**{f"class_{db}": c for db, c in read_classes.items()})
    out.to_csv(args.out_reads, sep="\t", compression="gzip")

    cont_out = pd.concat(contaminants)
    cont_out = cont_out.assign(species=cont_out["ref"].str.split("|").str[2], accession=cont_out["ref"].str.split("|").str[0])
    cont_out.rename(columns={"cls": "class"})[["db", "class", "accession", "species", "len_bin", "reads", "unique_seqs"]].round(3).to_csv(
        args.out_contaminants, sep="\t", index=False
    )
    ev = ends.reindex(tested)
    e = pd.DataFrame({"outcome": ev["outcome"].fillna("no_end_hit").to_numpy(), "end_len": ev["end_len"].astype("Int16").to_numpy(),
                      "human_ref": ev["human_ref"].to_numpy(), "length": reads.loc[tested, "length"].to_numpy(),
                      "count": reads.loc[tested, "count"].to_numpy()})
    e = e.groupby(["outcome", "end_len", "human_ref", "length"], dropna=False).agg(reads=("count", "sum"), unique_seqs=("count", "size"))
    e.reset_index().to_csv(args.out_ends, sep="\t", index=False)
    print(f"tested by their ends: {e['reads'].sum():,} reads ({e['unique_seqs'].sum():,} sequences); by outcome:")
    print(e.groupby("outcome")["reads"].sum().to_string())

    rules_out = pd.concat(rules)
    rules_out[["db", "human_class", "contaminant_class", "rule", "length", "reads", "unique_seqs"]].to_csv(
        args.out_contaminant_rules, sep="\t", index=False
    )
    top_out = pd.concat(top_hits).rename(columns={"cls": "class"})
    top_out[["db", "rank", "sequence", "count", "class", "mm", "ref", "offset", "length", "strand", "mature"]].to_csv(
        args.out_top, sep="\t", index=False
    )


if __name__ == "__main__":
    main()
