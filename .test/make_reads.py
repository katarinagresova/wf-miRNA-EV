"""Make the synthetic test libraries (.test/data/) from the test references (.test/refs/).

Each read is a sequence of known origin (the insert), with one extra 5' base and the RealSeq adapter, as in a
RealSeq-Biofluids library. truth.tsv gives the insert, its origin, the read count in each library and the class that
step 5 must give it (expected_<db>), or "removed" for inserts shorter than trim.min_len.

The inserts are kept only if their class is certain: an insert of one class must not match a reference of another
class exactly; genome inserts and short contaminant inserts must not match a small RNA reference with <= 1 mismatch
(step 5 aligns to the genome only the reads without a small RNA hit). The miRNA inserts are isomiRs with known
changes of miRNAs that have the same sequence in miRBase and MirGeneDB; an 18-nt isomiR with a substitution must not
count (miRNA_hairpin_other). Some rRNA inserts match a miRBase or MirGeneDB
mature miRNA with 0-3 mismatches: the off-target check must give them the result rRNA or tied (off_target_<db>).

Run in the environment of workflow/envs/py.yaml:

    python .test/make_reads.py --refs .test/refs --out .test/data
"""

import argparse
import gzip
import random
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "workflow" / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from make_refs import REGIONS  # noqa: E402
from prepare_refs import dna, ensembl_attrs, ncrna_class, read_fasta, read_feature_table, revcomp  # noqa: E402

ADAPTER = "TGGAATTCTCGGGTGCCAAGG"
TAIL = "AACTCCAGTCACATCTCGTATGCCGTCTTCTGCTTA"  # after the 3' adapter
READ_LEN = 75
# Libraries: sample, batch, and the read count factor of each origin
LIBRARIES = {
    "ev_1": ("run1", {"miRNA": 1, "rRNA": 5, "tRNA": 4, "ncRNA": 2, "contaminant": 3, "genome": 1, "other": 1}),
    "cells_1": ("run2", {"miRNA": 5, "rRNA": 1, "tRNA": 1, "ncRNA": 1, "contaminant": 0, "genome": 2, "other": 1}),
}
# Common EV miRNAs with the same sequence in miRBase and MirGeneDB
MIRNAS = [
    "hsa-let-7a-5p", "hsa-miR-16-5p", "hsa-miR-92a-3p", "hsa-miR-181a-5p", "hsa-miR-191-5p", "hsa-miR-26a-5p",
    "hsa-miR-103a-3p", "hsa-miR-107", "hsa-miR-126-3p", "hsa-miR-423-5p", "hsa-miR-451a", "hsa-miR-486-5p",
]
# KY962518.1 regions (0-based), as in config/config.yaml
RDNA_REGIONS = {"18S": (3657, 5526), "5.8S": (6596, 6753), "28S": (7920, 12971)}
NCRNA_SOURCES = {"Y_RNA": "RNY1", "snoRNA": "SNORD3A", "snRNA": "RNU2-1", "7SK": "RN7SK", "7SL_SRP": "RN7SL1", "vault_RNA": "VTRNA1-1"}
SMALL_CLASSES = ["miRNA", "rRNA", "tRNA", "mt_tRNA", "Y_RNA", "snoRNA", "snRNA", "7SK", "7SL_SRP", "vault_RNA", "other_ncRNA"]


class Refs:
    """The test references as the workflow builds them, as one string per class for substring searches."""

    def __init__(self, d):
        self.mirbase_hairpin = {h.split()[0]: dna(s) for h, s in read_fasta(d / "mirbase/hairpin.fa")}
        self.mirbase_mature = {h.split()[0]: dna(s) for h, s in read_fasta(d / "mirbase/mature.fa")}
        self.mirgenedb_pri = [dna(s) for _, s in read_fasta(d / "mirgenedb/pri.fa")]
        self.mirgenedb_mature = {
            h.split()[0]: dna(s) for kind in ["mature", "star"] for h, s in read_fasta(d / f"mirgenedb/{kind}.fa")
        }
        self.rdna = dna(next(read_fasta(d / "ncbi/KY962518.1.fa"))[1])
        self.ncrna = {}  # gene symbol -> (class, sequence) of the first transcript
        small = {c: [] for c in SMALL_CLASSES}
        small["miRNA"] = [*self.mirbase_hairpin.values(), *self.mirgenedb_pri]
        small["rRNA"].append(self.rdna)
        for h, s in read_fasta(d / "ensembl/ncrna.fa.gz"):
            a = ensembl_attrs(h)
            symbol, biotype = a.get("gene_symbol", ""), a["transcript_biotype"]
            cls = ncrna_class(a) or {"Mt_tRNA": "mt_tRNA", "Mt_rRNA": "rRNA", "rRNA": "rRNA"}.get(biotype)
            seq = dna(s) + ("CCA" if cls == "mt_tRNA" else "")
            small[cls].append(seq)
            self.ncrna.setdefault(symbol, (cls, seq))
        self.trna = {h.split()[0]: dna(s) + "CCA" for h, s in read_fasta(d / "gtrnadb/mature-tRNAs.fa")}
        small["tRNA"] += self.trna.values()
        self.small = {c: "|".join(seqs) for c, seqs in small.items()}
        self.rrna = [self.rdna, self.ncrna["RNA5S1"][1], self.ncrna["MT-RNR1"][1], self.ncrna["MT-RNR2"][1]]

        genome = {h.split()[0]: s.upper() for h, s in read_fasta(d / "ensembl/genome.fa.gz")}
        self.genome = genome
        regions = [genome[c][s:e] for c, s, e in REGIONS]
        self.genome_both = "|".join(regions + [revcomp(r) for r in regions])
        self.contaminants = {}
        for acc in ["NC_014448.1", "NC_001422.1"]:
            self.contaminants[acc] = dna(next(read_fasta(d / f"ncbi/{acc}.fa"))[1])
        self.contaminant_both = "|".join([*self.contaminants.values(), *map(revcomp, self.contaminants.values())])
        self.myco_features = read_feature_table(d / "ncbi/NC_014448.1.ft")
        self.gtf = d / "ensembl/gtf.gz"

    def exact_classes(self, seq):
        """Small RNA classes, 'genome' and 'contaminant' with an exact match (small RNAs: sense strand only)."""
        hits = {c for c, s in self.small.items() if seq in s}
        hits |= {"genome"} if seq in self.genome_both else set()
        hits |= {"contaminant"} if seq in self.contaminant_both else set()
        return hits

    def near_small(self, seq):
        """Small RNA classes with a match with <= 1 mismatch."""
        variants = [seq] + [seq[:i] + b + seq[i + 1 :] for i in range(len(seq)) for b in "ACGT" if b != seq[i]]
        return {c for c, s in self.small.items() if any(v in s for v in variants)}


def gtf_features(path):
    """GTF rows: chrom, start (0-based), end, strand, feature, gene name, gene biotype."""
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            attrs = dict(x.strip().split(" ", 1) for x in f[8].split(";") if x.strip())
            name = attrs.get("gene_name", attrs["gene_id"]).strip('"')
            rows.append((f[0], int(f[3]) - 1, int(f[4]), f[6], f[2], name, attrs["gene_biotype"].strip('"')))
    return rows


class Inserts:
    def __init__(self, refs, rng):
        self.refs, self.rng, self.rows, self.seen = refs, rng, [], set()

    def add(self, seq, origin, detail, expected, weight, off_target=(), expected_mirgenedb=None):
        if seq in self.seen:
            return False
        self.seen.add(seq)
        self.rows.append({
            "sequence": seq, "origin": origin, "detail": detail,
            "expected_mirbase": expected, "expected_mirgenedb": expected_mirgenedb or expected,
            "off_target": ",".join(off_target), "weight": weight,
        })
        return True

    def fragments(self, seq, n, lengths, origin, detail, expected, weight, clean):
        """n fragments of seq with lengths from `lengths` that pass clean(fragment)."""
        added, tries = 0, 0
        while added < n and tries < 1000:
            tries += 1
            length = self.rng.choice(lengths)
            if len(seq) < length:
                continue
            start = self.rng.randrange(len(seq) - length + 1)
            frag = seq[start : start + length]
            if "N" not in frag and clean(frag) and self.add(frag, origin, f"{detail}:{start}", expected, weight):
                added += 1
        if added < n:
            raise SystemExit(f"only {added} of {n} clean fragments of {detail}")


def mirnas(ins):
    """isomiRs of MIRNAS: the mature sequence, templated 5'/3' changes, a substitution, 1-3 non-templated 3' A or T,
    and the 18-nt 5' part with and without a substitution (step 5 counts reads of 16-18 nt only with 0 mismatches)."""
    r = ins.refs
    for name in MIRNAS:
        mature = r.mirbase_mature[name]
        if mature not in r.mirgenedb_mature.values():
            raise SystemExit(f"{name}: not the same sequence in MirGeneDB")
        hairpin = next(h for h in r.mirbase_hairpin.values() if mature in h)
        s = hairpin.index(mature)
        e = s + len(mature)
        sub = mature[:10] + ("A" if mature[10] != "A" else "C") + mature[11:]
        nta = "".join("T" if hairpin[e + i : e + i + 1] == "A" else "A" for i in range(3))
        weight = ins.rng.choice([5, 10, 20, 40])
        variants = {
            "canonical": (mature, weight * 10),
            "5p_ext1": (hairpin[s - 1 : e], weight),
            "5p_trim1": (mature[1:], weight),
            "3p_ext1": (hairpin[s : e + 1], weight * 2),
            "3p_ext2": (hairpin[s : e + 2], weight),
            "3p_trim1": (mature[:-1], weight * 2),
            "snv_10": (sub, weight),
            f"nta_{nta[0]}": (mature + nta[0], weight),
            "5p_18nt": (mature[:18], weight),
        }
        for kind, (seq, w) in variants.items():
            ins.add(seq, "miRNA", f"{name}:{kind}", "miRNA", w)
        # Without a hit in step 5 (2 or 3 mismatches to the hairpin): counted by the rule for non-templated
        # nucleotides. With 1 mismatch in 18 nt: not counted
        for n in [2, 3]:
            seq = mature + nta[:n]
            if not r.exact_classes(seq) and not r.near_small(seq):
                ins.add(seq, "miRNA", f"{name}:nta_{nta[:n]}", "miRNA", weight)
        seq = sub[:18]
        if not r.exact_classes(seq) and r.near_small(seq) == {"miRNA"}:
            ins.add(seq, "miRNA", f"{name}:5p_18nt_snv_10", "miRNA_hairpin_other", weight)


def hamming_windows(ref, query):
    """Mismatches of query against every window of ref."""
    a = np.frombuffer(ref.encode(), dtype=np.uint8)
    q = np.frombuffer(query.encode(), dtype=np.uint8)
    if len(a) < len(q):
        return np.array([], dtype=int)
    windows = np.lib.stride_tricks.sliding_window_view(a, len(q))
    return (windows != q).sum(axis=1)


def off_target_rrna(ins, per_mm=2):
    """rRNA windows that match a mature miRNA with 0-3 mismatches (per_mm per mismatch count)."""
    r = ins.refs
    found = {}
    for name, mature in sorted(r.mirbase_mature.items()):
        for ref in r.rrna:
            d = hamming_windows(ref, mature)
            if len(d) and d.min() <= 3:
                pos = int(d.argmin())
                found.setdefault(int(d.min()), []).append((name, ref[pos : pos + len(mature)]))
    for mm in range(4):
        added = 0
        for name, seq in found.get(mm, []):
            mgdb = min(int(hamming_windows(seq, m)[0]) if len(m) == len(seq) else 99 for m in r.mirgenedb_mature.values())
            # An exact hairpin match away from an equal-length mature: the class depends on the arm rules, so skip.
            # A better hit elsewhere on the hairpins: the miRNA-only mapping keeps only that hit, so skip.
            in_mirbase = any(seq in h for h in r.mirbase_hairpin.values())
            in_mirgenedb = any(seq in p for p in r.mirgenedb_pri)
            best = min(int(d.min()) for h in r.mirbase_hairpin.values() if len(d := hamming_windows(h, seq)))
            if added == per_mm or (in_mirbase and mm > 0) or (in_mirgenedb and mgdb > 0) or best < mm:
                continue
            dbs = ["mirbase"] + (["mirgenedb"] if mgdb <= 3 else [])
            if ins.add(seq, "rRNA_off_target", f"{name}:{mm}mm", "miRNA" if mm == 0 else "rRNA", 4, off_target=dbs,
                       expected_mirgenedb="miRNA" if mgdb == 0 else "rRNA"):
                added += 1


def small_rnas(ins):
    r = ins.refs

    def clean(cls):
        return lambda f: r.exact_classes(f) - {"genome"} == {cls}

    for region, (s, e) in RDNA_REGIONS.items():
        ins.fragments(r.rdna[s:e], 4, range(18, 41), "rRNA", region, "rRNA", 10, clean("rRNA"))
    for symbol in ["RNA5S1", "MT-RNR2"]:
        ins.fragments(r.ncrna[symbol][1], 2, range(20, 36), "rRNA", symbol, "rRNA", 10, clean("rRNA"))
    for name in ["tRNA-Gly-GCC-1-1", "tRNA-Glu-CTC-1-1", "tRNA-Val-CAC-1-1", "tRNA-Asp-GTC-2-1"]:
        trna = next(s for n, s in r.trna.items() if n.endswith(name))
        for length in [30, 33]:  # 5' halves
            if clean("tRNA")(trna[:length]):
                ins.add(trna[:length], "tRNA", f"{name}:5'half", "tRNA", 10)
        if clean("tRNA")(trna[-20:]):  # 3' fragment with CCA
            ins.add(trna[-20:], "tRNA", f"{name}:3'CCA", "tRNA", 10)
    ins.fragments(r.ncrna["MT-TF"][1], 2, range(20, 30), "mt_tRNA", "MT-TF", "mt_tRNA", 5, clean("mt_tRNA"))
    for cls, symbol in NCRNA_SOURCES.items():
        assert r.ncrna[symbol][0] == cls, (symbol, r.ncrna[symbol][0])
        ins.fragments(r.ncrna[symbol][1], 2, range(22, 34), "ncRNA", symbol, cls, 5, clean(cls))


def contaminants(ins):
    r = ins.refs
    myco = r.contaminants["NC_014448.1"]

    def clean(f):
        return r.exact_classes(f) == {"contaminant"}

    # tRNA fragments, the most common mycoplasma reads in the miRNA-PD libraries
    trnas = [f for f in r.myco_features if f["type"] == "tRNA"][:6]
    for f in trnas:
        (s, e), = f["iv"][:1]
        seq = myco[s:e] if f["strand"] == "+" else revcomp(myco[s:e])
        ins.fragments(seq, 1, range(25, 36), "contaminant", f"M_hyorhinis:{f['product']}", "mycoplasma", 10, clean)
    for acc, cls in [("NC_014448.1", "mycoplasma"), ("NC_001422.1", "phiX")]:
        g = r.contaminants[acc]
        for strand_seq in [g, revcomp(g)]:
            ins.fragments(strand_seq[:200_000], 2, range(25, 41), "contaminant", acc, cls, 5, clean)
    # Short exact matches with no human hit (rule exact_short)
    ins.fragments(myco[:200_000], 2, [16, 17], "contaminant", "NC_014448.1:short", "mycoplasma", 5,
                  lambda f: clean(f) and not r.near_small(f) and f not in r.genome_both)
    # A contaminant 25-nt end joined to 12 random nt: no full-length hit; class from the ends (item 4c)
    for i in range(2):
        start = ins.rng.randrange(100_000)
        seq = myco[start : start + 25] + "".join(ins.rng.choices("ACGT", k=12))
        if not r.exact_classes(seq) and not r.near_small(seq):
            ins.add(seq, "contaminant", f"NC_014448.1:end:{start}", "mycoplasma", 5)


def genome_reads(ins):
    r = ins.refs
    feats = gtf_features(r.gtf)

    def covered(chrom, s, e, strand, kinds=("exon", "gene")):
        return [f for f in feats if f[0] == chrom and f[1] < e and f[2] > s and f[3] == strand and f[4] in kinds]

    def clean(f):
        return r.exact_classes(f) == {"genome"} and not r.near_small(f)

    def sample(cls, chrom, intervals, strand, n, ok):
        added, tries = 0, 0
        while added < n and tries < 2000:
            tries += 1
            s0, e0 = intervals[ins.rng.randrange(len(intervals))]
            length = ins.rng.choice(range(25, 31))
            if e0 - s0 < length:
                continue
            s = ins.rng.randrange(s0, e0 - length + 1)
            seq = r.genome[chrom][s : s + length]
            seq = seq if strand == "+" else revcomp(seq)
            if ok(s, s + length) and clean(seq) and ins.add(seq, "genome", f"{chrom}:{s}:{strand}", cls, 3):
                added += 1
        if added < n:
            raise SystemExit(f"only {added} of {n} {cls} reads")

    def exons(gene):
        return [(f[1], f[2]) for f in feats if f[5] == gene and f[4] == "exon"]

    sod1 = next(f for f in feats if f[5] == "SOD1" and f[4] == "gene")
    sample("mRNA", "21", exons("SOD1"), "+", 3, lambda s, e: True)
    sample("intron", "21", [(sod1[1], sod1[2])], "+", 3, lambda s, e: not covered("21", s, e, "+", ("exon",)))
    sample("lncRNA", "21", exons("SOD1-DT"), "-", 3,
           lambda s, e: all(f[6] == "lncRNA" for f in covered("21", s, e, "-", ("exon",))))
    sample("intergenic", "21", [(31_600_000, 31_700_000)], "+", 3,
           lambda s, e: not covered("21", s, e, "+") and not covered("21", s, e, "-"))
    sample("mt_other", "MT", [(5903, 7445)], "+", 2, lambda s, e: True)  # MT-CO1
    # Pre-tRNA: the 3' end of tRNA-Gly-GCC-1-5 (minus strand, 21:17454788-17454859) and 10 nt of its trailer
    seq = revcomp(r.genome["21"][17_454_778:17_454_803])
    assert clean(seq)
    ins.add(seq, "genome", "tRNA-Gly-GCC-1-5:trailer", "tRNA_genomic", 3)


def others(ins):
    for _ in range(5):
        seq = "".join(ins.rng.choices("ACGT", k=ins.rng.randrange(25, 41)))
        if not ins.refs.exact_classes(seq) and not ins.refs.near_small(seq):
            ins.add(seq, "random", "random", "unmapped", 1)
    ins.add("", "adapter_dimer", "empty insert", "removed", 20)
    for length in [8, 12, 15]:
        ins.add("".join(ins.rng.choices("ACGT", k=length)), "too_short", f"{length} nt", "removed", 2)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--refs", type=Path, default=Path(".test/refs"))
    p.add_argument("--out", type=Path, default=Path(".test/data"))
    p.add_argument("--seed", type=int, default=1)
    args = p.parse_args()
    rng = random.Random(args.seed)
    ins = Inserts(Refs(args.refs), rng)
    mirnas(ins)
    off_target_rrna(ins)
    small_rnas(ins)
    contaminants(ins)
    genome_reads(ins)
    others(ins)

    args.out.mkdir(parents=True, exist_ok=True)
    group = {"miRNA": "miRNA", "rRNA": "rRNA", "rRNA_off_target": "rRNA", "tRNA": "tRNA", "mt_tRNA": "tRNA",
             "ncRNA": "ncRNA", "contaminant": "contaminant", "genome": "genome"}
    for row in ins.rows:
        for sample, (_, factor) in LIBRARIES.items():
            row[sample] = row["weight"] * factor[group.get(row["origin"], "other")]
    for sample in LIBRARIES:
        reads = [row["sequence"] for row in ins.rows for _ in range(row[sample])]
        rng.shuffle(reads)
        with gzip.open(args.out / f"{sample}.fastq.gz", "wt", compresslevel=9) as fh:
            for i, insert in enumerate(reads):
                seq = (rng.choice("ACGT") + insert + ADAPTER + TAIL)[:READ_LEN]
                fh.write(f"@{sample}:{i + 1}\n{seq}\n+\n{'I' * len(seq)}\n")
        print(f"{sample}: {len(reads):,} reads, {sum(row[sample] > 0 for row in ins.rows)} inserts")
    cols = ["sequence", "origin", "detail", "expected_mirbase", "expected_mirgenedb", "off_target", *LIBRARIES]
    with open(args.out / "truth.tsv", "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for row in ins.rows:
            fh.write("\t".join(str(row[c]) for c in cols) + "\n")
    with open(args.out.parent / "config/samples.tsv", "w") as fh:
        fh.write("sample\tfastq\tbatch\n")
        for sample, (batch, _) in LIBRARIES.items():
            fh.write(f"{sample}\tdata/{sample}.fastq.gz\t{batch}\n")
    print(f"{len(ins.rows)} inserts; by origin:")
    for origin in dict.fromkeys(row["origin"] for row in ins.rows):
        print(f"  {origin}: {sum(row['origin'] == origin for row in ins.rows)}")


if __name__ == "__main__":
    main()
