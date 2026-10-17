"""Interface tables of the miRNA counter mirge3, from the miRge3.0 output of one library and miRNA database.

Table: mirna (the miRge3 entry of miR.Counts.csv, all entries of the library), reads, members (the names of the
workflow's miRNA database release that the entry includes, separated by ';').

miRge3 merges near-identical miRNAs into one name (e.g. hsa-miR-103a-3p/107), and its bundled MirGeneDB library
is an older release than ours, with renumbered families (miR-93 is Hsa-Mir-17-P1c_5p there, Hsa-Mir-17-P4d_5p in
MirGeneDB 3.0). Each mature name of the workflow's release is mapped to the miRge3 name with the identical mature
sequence (by name if there is none), then to miRge3's merged name. The members of an entry are the names mapped to it.

Per-read table: sequence, mirna, reads. The sequences of mapped.csv with an exact or isomiR miRNA, after the count
filter of miRge3 (mirge/libs/summary.py, mirge_can): an entry counts only if it has 2 or more exact reads and its
exact reads are more than 10% of its isomiR reads. The isomiR GFF cannot give this table: it does not have the
exact/isomiR calls of the count, and it leaves out the entries without an annotation. The script checks that the
per-read table gives the counts of miR.Counts.csv.
"""

import argparse
import re

import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True, help="the count column of mapped.csv")
    p.add_argument("--counts", required=True, help="miRge3 miR.Counts.csv")
    p.add_argument("--mapped", required=True, help="miRge3 mapped.csv")
    p.add_argument("--merges", required=True, help="miRge3 library merges file of the database")
    p.add_argument("--mirge3-mature", required=True, help="miRge3 library mature FASTA of the database")
    p.add_argument("--mature", nargs="+", required=True, help="mature FASTAs of the workflow's database release")
    p.add_argument("--arms", required=True, help="<db>_arms.tsv: the mature names of the workflow's release")
    p.add_argument("--exact-ratio", type=float, default=0.1, help="miRge3 -ex (crThreshold)")
    p.add_argument("--out-table", required=True)
    p.add_argument("--out-reads", required=True)
    return p.parse_args()


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


def counted_reads(mapped, sample, merges, ratio_min):
    """The reads that miRge3 counts, with the entry that counts them (as in mirge_can and summarize of miRge3)."""
    m = pd.read_csv(mapped, usecols=["Sequence", "exact miRNA", "isomiR miRNA", sample], dtype=str)
    m[sample] = m[sample].astype(float)
    exact = m[m["exact miRNA"].notna()].assign(name=lambda d: d["exact miRNA"], exact=True)
    iso = m[m["isomiR miRNA"].notna()].assign(name=lambda d: d["isomiR miRNA"], exact=False)
    x = exact.groupby("name")[sample].sum()
    y = iso.groupby("name")[sample].sum().reindex(x.index, fill_value=0)
    ratio = (x / y).where(y > 0, x)
    keep = x.index[(x >= 2) & (ratio > ratio_min)]
    reads = pd.concat([exact, iso])
    reads = reads[reads["name"].isin(keep)]
    # miRge3 takes the merged name of the last merges line that lists the name, with the SNP suffix
    merged_of = {}
    for line in open(merges):
        fields = line.strip().split(",")
        for item in fields[1:]:
            merged_of[item] = fields[0]
    return pd.DataFrame({
        "sequence": reads["Sequence"].to_numpy(),
        "mirna": reads["name"].map(lambda n: merged_of.get(n, n)).to_numpy(),
        "reads": reads[sample].to_numpy(),
    })


def main():
    args = parse_args()
    counts = pd.read_csv(args.counts, index_col=0).iloc[:, 0]

    reads = counted_reads(args.mapped, args.sample, args.merges, args.exact_ratio)
    check = pd.concat({"reads": reads.groupby("mirna")["reads"].sum(), "counts": counts}, axis=1).fillna(0)
    bad = check[(check["reads"] - check["counts"]).abs() > 1e-6]
    if len(bad):
        raise SystemExit(f"per-read table does not give miR.Counts.csv for {len(bad)} entries:\n{bad.head(10)}")
    reads.sort_values(["reads", "sequence"], ascending=[False, True]).to_csv(
        args.out_reads, sep="\t", index=False, compression="gzip"
    )

    names = set(pd.read_csv(args.arms, sep="\t")["mature"])
    our_seqs = {}
    for path in args.mature:
        our_seqs.update({n: s for n, s in read_fasta(path).items() if n in names})
    to_mirge3 = mirge3_names(our_seqs, args.mirge3_mature, read_merges(args.merges))
    members = pd.Series(list(to_mirge3), index=list(to_mirge3.values())).groupby(level=0).agg(lambda n: ";".join(sorted(n)))
    table = pd.DataFrame({"mirna": counts.index, "reads": counts.to_numpy(), "members": counts.index.map(members)})
    table.sort_values("reads", ascending=False, kind="stable").to_csv(args.out_table, sep="\t", index=False)

    unlinked = sorted(n for n, m in to_mirge3.items() if m not in counts.index)
    print(f"{len(table)} miRge3 entries, {table['reads'].sum():,.0f} reads, {len(reads):,} sequences; "
          f"{table['members'].notna().sum()} entries with members; {len(unlinked)} names of the release map to no "
          f"miRge3 entry (e.g. {unlinked[:5]})")


if __name__ == "__main__":
    main()
