"""Build the small RNA references for step 4 (and the single-copy rRNA reference shared with Part 2), and the
contaminant genome reference.

All output FASTAs are DNA (U -> T), one line per sequence. Reference names are `<id>|<class>|<name>`,
so the class of every alignment can be read from the reference name alone.
"""

import argparse
import gzip
import re
import tarfile
from collections import Counter


def read_fasta(path):
    """Yield (header, sequence) pairs; handles gzip and multi-line records."""
    opener = gzip.open if str(path).endswith(".gz") else open
    header, seq = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header, seq = line[1:], []
            elif line:
                seq.append(line)
    if header is not None:
        yield header, "".join(seq)


def dna(seq):
    return seq.upper().replace("U", "T")


def revcomp(seq):
    return seq[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def genome_region(path, chrom, start, end):
    """Sequence of chrom[start:end] (0-based) from an uncompressed FASTA, reading only up to that chromosome."""
    seq, found = [], False
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if found:
                    break
                found = line[1:].split()[0] == chrom
            elif found:
                seq.append(line.rstrip())
    if not found:
        raise SystemExit(f"{chrom} not in {path}")
    return "".join(seq)[start:end].upper()


def write_fasta(records, path):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def arms_table(precursors, matures):
    """Locate every mature sequence in every precursor that contains it (0-based, half-open).

    A mature sequence can sit on several precursors (e.g. hsa-let-7a-5p on let-7a-1/-2/-3), and in
    MirGeneDB the same sequence can be listed under several names, so all pairs are kept.
    """
    rows = []
    for pname, pseq in precursors.items():
        for mname, mseq in matures.items():
            start = pseq.find(mseq)
            while start != -1:
                rows.append((pname, mname, start, start + len(mseq)))
                start = pseq.find(mseq, start + 1)
    return rows


def write_arms(rows, matures, path):
    """Write the arms table: precursor, mature, start, end, entry. The entry is the counted miRNA: the mature names
    with the same sequence, joined by '/' (e.g. hsa-miR-199a-3p/hsa-miR-199b-3p)."""
    names = {}
    for name, seq in matures.items():
        names.setdefault(seq, []).append(name)
    entry = {n: "/".join(sorted(group)) for group in names.values() for n in group}
    with open(path, "w") as fh:
        fh.write("precursor\tmature\tstart\tend\tentry\n")
        for row in sorted(rows):
            fh.write("\t".join(map(str, [*row, entry[row[1]]])) + "\n")
    merged = [g for g in names.values() if len(g) > 1]
    print(f"{len(merged)} entries of 2 or more mature names with the same sequence ({sum(map(len, merged))} names)")


def mirbase(args):
    """miRBase: the hairpins of one species as the reference, its mature sequences located on them."""
    prefix = f"{args.species}-"
    hairpins = {h.split()[0]: dna(s) for h, s in read_fasta(args.hairpin) if h.startswith(prefix)}
    matures = {h.split()[0]: dna(s) for h, s in read_fasta(args.mature) if h.startswith(prefix)}
    write_fasta(((f"{n}|miRNA|{n}", s) for n, s in hairpins.items()), args.fasta)
    rows = arms_table(hairpins, matures)
    write_arms(rows, matures, args.arms)
    located = {r[1] for r in rows}
    print(f"miRBase: {len(hairpins)} {args.species} hairpins, {len(matures)} {args.species} mature, {len(located)} located on a hairpin")
    missing = sorted(set(matures) - located)
    if missing:
        raise SystemExit(f"{len(missing)} mature sequences not found on any hairpin: {missing[:10]}")


def mirgenedb(args):
    """MirGeneDB: precursors with 30 nt flanks (`_pri`) as the reference, mature and star arms located on them."""
    pri = {h.split()[0].removesuffix("_pri"): dna(s) for h, s in read_fasta(args.pri)}
    matures = {}
    for path in (args.mature, args.star):
        for h, s in read_fasta(path):
            matures[h.split()[0]] = dna(s)
    write_fasta(((f"{n}|miRNA|{n}", s) for n, s in pri.items()), args.fasta)
    rows = arms_table(pri, matures)
    # Keep only arms located on their own precursor (Hsa-Let-7-P1b_5p on Hsa-Let-7-P1b); a sequence shared
    # between paralogues is still found through its own entry on each paralogue.
    own = [r for r in rows if r[1].rsplit("_", 1)[0] == r[0]]
    write_arms(own, matures, args.arms)
    missing = sorted(set(matures) - {r[1] for r in own})
    print(f"MirGeneDB: {len(pri)} precursors, {len(matures)} mature/star arms, {len(missing)} not on own precursor")
    if missing:
        raise SystemExit(f"arms not found on their own precursor: {missing[:10]}")


def ensembl_attrs(header):
    """Parse `key:value` fields of an Ensembl FASTA header (description may contain spaces)."""
    fields = header.split()
    attrs = {"id": fields[0].split(".")[0], "location": fields[2]}
    for f in fields[3:]:
        key, sep, val = f.partition(":")
        if key == "description":
            break
        if sep:
            attrs[key] = val
    return attrs


def ncrna_class(a):
    """Small ncRNA class of an Ensembl ncRNA transcript, or None if it belongs to another reference or is not small."""
    biotype, symbol = a["transcript_biotype"], a.get("gene_symbol", "")
    if biotype in {"snRNA", "misc_RNA", "snoRNA", "scaRNA", "vault_RNA", "ribozyme", "sRNA", "scRNA"}:
        if re.match(r"(RNY\d|RNY$|Y_RNA$)", symbol):
            return "Y_RNA"
        if re.match(r"(RN7SL|Metazoa_SRP$)", symbol):
            return "7SL_SRP"
        if re.match(r"(RN7SK|7SK$)", symbol):
            return "7SK"
        if biotype == "vault_RNA" or re.match(r"(VTRNA|Vault$)", symbol):
            return "vault_RNA"
        if biotype in {"snoRNA", "scaRNA"}:
            return "snoRNA"
        if biotype == "snRNA":
            return "snRNA"
        return "other_ncRNA"
    # miRNA, rRNA and Mt_tRNA have their own references; lncRNA and the rest are annotated from the genome
    return None


def ncrna(args):
    """Ensembl small ncRNAs other than miRNA, rRNA and tRNA (Y RNA, 7SL, 7SK, vault, snoRNA, snRNA, other)."""
    records, classes = [], Counter()
    for h, s in read_fasta(args.ncrna):
        a = ensembl_attrs(h)
        cls = ncrna_class(a)
        if cls:
            records.append((f"{a['id']}|{cls}|{a.get('gene_symbol', a['gene'].split('.')[0])}", dna(s)))
            classes[cls] += 1
    write_fasta(records, args.fasta)
    print("ncRNA transcripts per class:", dict(classes))


def trna(args):
    """Mature tRNAs: GtRNAdb nuclear (high confidence) plus Ensembl mt-tRNAs, all with the 3' CCA added."""
    records = []
    for h, s in read_fasta(args.gtrnadb):
        # Genus_species_tRNA-Gly-GCC-2-1 -> tRNA-Gly-GCC-2-1
        name = re.sub(r"^[A-Z][a-z]+_[a-z]+_", "", h.split()[0])
        records.append((f"{name}|tRNA|{name}", dna(s) + "CCA"))
    n_nuc = len(records)
    for h, s in read_fasta(args.ncrna):
        a = ensembl_attrs(h)
        if a["transcript_biotype"] == "Mt_tRNA":
            records.append((f"{a['id']}|mt_tRNA|{a['gene_symbol']}", dna(s) + "CCA"))
    write_fasta(records, args.fasta)
    print(f"tRNA: {n_nuc} nuclear (GtRNAdb), {len(records) - n_nuc} mitochondrial (Ensembl)")

    # Genomic tRNA loci (pre-tRNA, leader/trailer, pseudogenes) for annotating genome hits, Ensembl chromosome names
    with tarfile.open(args.gtrnadb_tar) as tar, open(args.bed, "w") as out:
        member = next(m for m in tar.getmembers() if m.name.endswith("-tRNAs.bed"))
        n = 0
        for line in tar.extractfile(member).read().decode().splitlines():
            f = line.split("\t")
            chrom = f[0].removeprefix("chr")
            chrom = "MT" if chrom == "M" else chrom
            out.write("\t".join([chrom, f[1], f[2], f"{f[3]}|tRNA_locus", "0", f[5]]) + "\n")
            n += 1
    print(f"tRNA loci: {n}")


def rrna(args):
    """Single-copy rRNA reference: rDNA repeat (45S unit + IGS), 5S, and mt-rRNA."""
    records = []
    (h, s), *rest = list(read_fasta(args.rdna))
    if rest:
        raise SystemExit(f"expected one rDNA record in {args.rdna}")
    rdna_id = h.split()[0]
    records.append((f"{rdna_id}|rRNA|45S_rDNA", dna(s)))
    extra = {}
    loc_5s = None
    for h, s in read_fasta(args.ncrna):
        a = ensembl_attrs(h)
        if a["id"] == args.rna5s or a.get("gene_symbol") in args.mt_rrna:
            name = a.get("gene_symbol")
            extra[name] = (f"{a['id']}|rRNA|{name}", dna(s))
            if a["id"] == args.rna5s:
                loc_5s = a["location"]
    for name in ["RNA5S1", *args.mt_rrna]:
        if name not in extra:
            raise SystemExit(f"{name} not found in {args.ncrna}")

    # 5S with genomic flanks: EV reads run past the annotated ends (e.g. the 3' oligo-U of the Pol III
    # transcript), and would not align end-to-end to the bare 119 nt gene
    ref5s, core = extra["RNA5S1"]
    _, _, chrom, start, end, strand = loc_5s.split(":")
    start, end, f = int(start) - 1, int(end), args.rna5s_flank
    locus = genome_region(args.genome, chrom, start - f, end + f)
    locus = revcomp(locus) if strand == "-1" else locus
    if locus[f : f + len(core)] != core:
        raise SystemExit("RNA5S1 transcript does not match the genome at its annotated location")
    extra["RNA5S1"] = (ref5s, locus)
    records.extend(extra[name] for name in ["RNA5S1", *args.mt_rrna])
    write_fasta(records, args.fasta)
    print("rRNA reference:", ", ".join(f"{r.split('|')[2]} ({len(s)} nt)" for r, s in records))


def rrna_features(args):
    """Feature BED of the rRNA reference (0-based half-open, reference coordinates): the rDNA regions given, mature 5S
    and its genomic flanks, and the mt-rRNAs as whole references."""
    records = list(read_fasta(args.fasta))
    (rdna, rdna_seq), (ref5s, seq5s), *mt = records
    f, n5s = args.rna5s_flank, args.rna5s_length
    if f + n5s > len(seq5s):
        raise SystemExit(f"5S of {n5s} nt does not fit in {ref5s} ({len(seq5s)} nt, {f} nt 5' flank)")
    with open(args.bed, "w") as out:
        for feat, start, end in args.features:
            out.write(f"{rdna}\t{start}\t{min(int(end), len(rdna_seq))}\t{feat}\n")
        out.write(f"{ref5s}\t0\t{f}\t5S_5'flank\n{ref5s}\t{f}\t{f + n5s}\tRNA5S1\n")
        out.write(f"{ref5s}\t{f + n5s}\t{len(seq5s)}\t5S_3'flank\n")
        for ref, seq in mt:
            out.write(f"{ref}\t0\t{len(seq)}\t{ref.split('|')[2]}\n")
    print(f"Mature 5S: {seq5s[f:f + n5s]} ({n5s} nt)")


def read_feature_table(path):
    """NCBI feature table -> (type, product, strand, [(start, end), ...]) with 0-based half-open intervals."""
    feats, cur = [], None
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if line.startswith(">"):
            continue
        if f[0]:
            a, b = (int(x.strip("<>")) for x in f[:2])
            iv = (min(a, b) - 1, max(a, b))
            if len(f) > 2 and f[2]:
                cur = {"type": f[2], "product": "", "strand": "+" if b >= a else "-", "iv": [iv]}
                feats.append(cur)
            elif cur is not None:
                cur["iv"].append(iv)
        elif cur is not None and len(f) >= 5 and f[3] in ("product", "gene") and not cur["product"]:
            cur["product"] = f[4]
    return feats


def contaminant(args):
    """Contaminant genomes, named <accession>|<class>|<species>, plus their RNA and CDS features as BED."""
    species = {acc: (cls, re.sub(r"\W+", "_", sp).strip("_")) for acc, cls, sp in args.species}
    records = []
    with open(args.bed, "w") as out:
        for fasta, ft in zip(sorted(args.genomes), sorted(args.feature_tables)):
            (h, s), *rest = list(read_fasta(fasta))
            acc = h.split()[0]
            if rest or acc not in species or not ft.endswith(f"/{acc}.ft"):
                raise SystemExit(f"unexpected contaminant input: {fasta}, {ft}")
            cls, sp = species[acc]
            ref = f"{acc}|{cls}|{sp}"
            records.append((ref, dna(s)))
            n = Counter()
            for feat in read_feature_table(ft):
                # gene rows duplicate their CDS/tRNA/rRNA rows
                if feat["type"] == "gene":
                    continue
                n[feat["type"]] += 1
                for start, end in feat["iv"]:
                    out.write(f"{ref}\t{start}\t{end}\t{feat['type']}:{feat['product']}\t0\t{feat['strand']}\n")
            print(f"{ref}: {len(s):,} nt, " + ", ".join(f"{v} {k}" for k, v in sorted(n.items())))
    write_fasta(records, args.fasta)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("mirbase")
    s.add_argument("--species", required=True, help="miRBase name prefix, for example hsa")
    s.add_argument("--mature", required=True)
    s.add_argument("--hairpin", required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--arms", required=True)
    s.set_defaults(func=mirbase)

    s = sub.add_parser("mirgenedb")
    s.add_argument("--mature", required=True)
    s.add_argument("--star", required=True)
    s.add_argument("--pri", required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--arms", required=True)
    s.set_defaults(func=mirgenedb)

    s = sub.add_parser("ncrna")
    s.add_argument("--ncrna", required=True)
    s.add_argument("--fasta", required=True)
    s.set_defaults(func=ncrna)

    s = sub.add_parser("trna")
    s.add_argument("--gtrnadb", required=True)
    s.add_argument("--gtrnadb-tar", required=True)
    s.add_argument("--ncrna", required=True)
    s.add_argument("--fasta", required=True)
    s.add_argument("--bed", required=True)
    s.set_defaults(func=trna)

    s = sub.add_parser("rrna")
    s.add_argument("--rdna", required=True)
    s.add_argument("--ncrna", required=True)
    s.add_argument("--rna5s", required=True)
    s.add_argument("--rna5s-flank", type=int, required=True, help="genomic flank added on both sides of 5S (nt)")
    s.add_argument("--genome", required=True, help="uncompressed genome FASTA, for the 5S flanks")
    s.add_argument("--mt-rrna", nargs="+", required=True)
    s.add_argument("--fasta", required=True)
    s.set_defaults(func=rrna)

    s = sub.add_parser("rrna-features")
    s.add_argument("--fasta", required=True, help="rRNA reference: rDNA repeat, 5S with flanks, mt-rRNAs")
    s.add_argument("--rna5s-flank", type=int, required=True)
    s.add_argument("--rna5s-length", type=int, required=True, help="mature 5S (nt)")
    s.add_argument("--features", nargs=3, action="append", default=[], metavar=("NAME", "START", "END"))
    s.add_argument("--bed", required=True)
    s.set_defaults(func=rrna_features)

    s = sub.add_parser("contaminant")
    s.add_argument("--genomes", nargs="+", required=True)
    s.add_argument("--feature-tables", nargs="+", required=True, help="NCBI feature tables, <accession>.ft")
    s.add_argument("--species", nargs=3, action="append", required=True, metavar=("ACCESSION", "CLASS", "SPECIES"))
    s.add_argument("--fasta", required=True)
    s.add_argument("--bed", required=True)
    s.set_defaults(func=contaminant)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
