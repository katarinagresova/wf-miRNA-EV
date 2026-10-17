"""Cut the small test references in .test/refs/ from the full references of a workflow run.

Run once, in the environment of workflow/envs/mirge3.yaml (bowtie 1.3.1; the bowtie-inspect wrapper does not work
with the Python 3.15 of align.yaml). The output goes into the repository:

    python .test/make_refs.py --resources resources --out .test/refs

- miRBase: the human entries. MirGeneDB and GtRNAdb: unchanged (human only, small).
- Ensembl ncRNA: the small ncRNA biotypes, the mt-tRNAs and mt-rRNAs, and RNA5S1.
- Genome: the REGIONS at their real positions; all other bases are N, and each chromosome ends after its last region.
  So the GTF and the ncRNA headers keep their coordinates. GTF: the lines that overlap a region.
- NCBI records: the rDNA and two contaminant genomes (M. hyorhinis, phiX).
- miRge3.0 library: the files that miRge3.0 reads without -trf, -ai and -nmir. The miRNA, hairpin, tRNA, snoRNA and
  rRNA indices are complete; the ncRNA and mRNA indices have only the transcripts in the REGIONS. The indices are
  built with a small lookup table (--ftabchars 6), so that the files are small.
"""

import argparse
import gzip
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

# Genome regions (chromosome, 0-based start, end) of the test genome
REGIONS = [
    ("1", 228_609_000, 228_612_000),  # RNA5S1 (5S flanks of the rRNA reference)
    ("21", 17_452_000, 17_457_000),  # tRNA-Gly-GCC-1-5 locus
    ("21", 31_600_000, 31_700_000),  # SOD1 (mRNA, intron), SOD1-DT (lncRNA), intergenic
    ("MT", 0, 16_569),  # mt genome
]
NCRNA_BIOTYPES = {"snRNA", "misc_RNA", "snoRNA", "scaRNA", "vault_RNA", "ribozyme", "sRNA", "scRNA", "Mt_tRNA", "Mt_rRNA"}
RNA5S1 = "ENST00000362482"
NCBI = ["KY962518.1", "NC_014448.1", "NC_001422.1"]
MIRGE3_FILES = [
    "annotation.Libs/human_merges_miRBase.csv",
    "annotation.Libs/human_merges_MirGeneDB.csv",
    "annotation.Libs/human_miRBase.gff3",
    "annotation.Libs/human_MirGeneDB.gff3",
    "fasta.Libs/human_mature_miRBase.fa",
    "fasta.Libs/human_mature_MirGeneDB.fa",
]
MIRGE3_INDICES = [
    "mirna_miRBase", "mirna_MirGeneDB", "hairpin_miRBase", "hairpin_MirGeneDB",
    "mature_trna", "pre_trna", "snorna", "rrna", "ncrna_others", "mrna",
]


def read_fasta(path):
    """Yield (header, sequence); handles gzip and multi-line records."""
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


def write_fasta(records, path, width=None):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            step = width or len(seq) or 1
            fh.writelines(seq[i : i + step] + "\n" for i in range(0, len(seq), step))


def overlaps(chrom, start, end):
    return any(c == chrom and start < e and end > s for c, s, e in REGIONS)


def genome(src, dst):
    """The REGIONS of the genome FASTA (sorted by position within a chromosome); other bases N."""
    chroms = list(dict.fromkeys(c for c, *_ in REGIONS))
    records = {}
    for header, seq in read_fasta(src):
        chrom = header.split()[0]
        if chrom in chroms:
            parts, pos = [], 0
            for c, s, e in REGIONS:
                if c == chrom:
                    parts += ["N" * (s - pos), seq[s:e]]
                    pos = e
            records[chrom] = "".join(parts)
    write_fasta(((c, records[c]) for c in chroms), dst, width=60)


def gtf(src, dst):
    with gzip.open(src, "rt") as fh, gzip.open(dst, "wt") as out:
        for line in fh:
            f = line.split("\t")
            if line.startswith("#") or overlaps(f[0], int(f[3]) - 1, int(f[4])):
                out.write(line)


def ncrna(src, dst):
    def keep(header):
        biotype = re.search(r"transcript_biotype:(\S+)", header).group(1)
        return biotype in NCRNA_BIOTYPES or header.startswith(RNA5S1)

    write_fasta(((h, s) for h, s in read_fasta(src) if keep(h)), dst, width=60)


def mirge3_index_names(header):
    """Chromosome and intervals of a miRge3.0 ncRNA/mRNA name: ENST... chr1 segs:1-9,10-981 cds:+:65565-65573,..."""
    f = header.split()
    chrom = f[1].removeprefix("chr")
    chrom = "MT" if chrom == "M" else chrom
    intervals = [tuple(map(int, iv.split("-"))) for iv in f[-1].split(":")[-1].split(",")]
    return chrom, intervals


def mirge3_lib(src, dst):
    for name in MIRGE3_FILES:
        (dst / name).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(src / name, dst / name)
    (dst / "index.Libs").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        for name in MIRGE3_INDICES:
            index = f"human_{name}"
            fasta = subprocess.run(
                ["bowtie-inspect", "-a", "1000000", src / "index.Libs" / index],
                check=True, capture_output=True, text=True,
            ).stdout.splitlines()
            records = list(zip([h[1:] for h in fasta[0::2]], fasta[1::2]))
            if name in ["ncrna_others", "mrna"]:
                records = [
                    (h, s) for h, s in records
                    if any(overlaps(c, a - 1, b) for c, ivs in [mirge3_index_names(h)] for a, b in ivs)
                ]
            write_fasta(records, f"{tmp}/{index}.fa")
            subprocess.run(
                ["bowtie-build", "--ftabchars", "6", f"{tmp}/{index}.fa", dst / "index.Libs" / index],
                check=True, capture_output=True,
            )
            print(f"miRge3.0 {index}: {len(records)} sequences")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--resources", type=Path, default=Path("resources"), help="resources/ of a workflow run")
    p.add_argument("--out", type=Path, default=Path(".test/refs"))
    args = p.parse_args()
    res, out = args.resources, args.out
    for d in ["mirbase", "mirgenedb", "gtrnadb", "ensembl", "ncbi", "mirge3/human"]:
        (out / d).mkdir(parents=True, exist_ok=True)

    for kind in ["mature", "hairpin"]:
        records = [(h, s) for h, s in read_fasta(res / f"download/mirbase/{kind}.fa") if h.startswith("hsa-")]
        write_fasta(records, out / f"mirbase/{kind}.fa")
    for kind in ["mature", "star", "pri"]:
        shutil.copy(res / f"download/mirgenedb/{kind}.fa", out / f"mirgenedb/{kind}.fa")
    shutil.copy(res / "download/gtrnadb/hg38-mature-tRNAs.fa", out / "gtrnadb/mature-tRNAs.fa")
    shutil.copy(res / "download/gtrnadb/hg38-tRNAs.tar.gz", out / "gtrnadb/tRNAs.tar.gz")

    genome(res / "refs/genome.fa", out / "ensembl/genome.fa.gz")
    gtf(res / "download/ensembl/gtf.gz", out / "ensembl/gtf.gz")
    ncrna(res / "download/ensembl/ncrna.gz", out / "ensembl/ncrna.fa.gz")

    shutil.copy(res / "download/rdna.fa", out / f"ncbi/{NCBI[0]}.fa")
    for acc in NCBI[1:]:
        for ext in ["fa", "ft"]:
            shutil.copy(res / f"download/contaminants/{acc}.{ext}", out / f"ncbi/{acc}.{ext}")

    mirge3_lib(res / "mirge3/human", out / "mirge3/human")


if __name__ == "__main__":
    main()
