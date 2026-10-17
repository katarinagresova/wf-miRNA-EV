# Configuration

The workflow needs two files: a sample sheet and a config file. `config/config.yaml` is an example config for human
libraries made with the RealSeq-Biofluids kit. The workflow checks both files against the schemas in
`workflow/schemas/` before it starts.

## Sample sheet

A TSV file with one row for each library. Set its path in `samples`. `config/samples.tsv` is an example.

| Column | Required | Content |
|---|---|---|
| `sample` | yes | Library name. It is part of the output paths: letters, digits, `.`, `_` and `-` only. |
| `fastq` | yes | Path to the single-end FASTQ file (gzip). |
| `batch` | no | Sequencing run. The insert length plot (`results/trim/insert_lengths.pdf`) has one panel for each batch. |
| other columns | no | Copied into the summary tables, for example cell line, material or condition. |

The workflow processes every row. To leave out a library, remove its row.

## Library kit and trimming

`kit` selects a preset of trim settings. The values in `trim` override the preset.

| Kit | Adapter (3′) | `cut_5p` | Note |
|---|---|---|---|
| `realseq_biofluids` | `TGGAATTCTCGGGTGCCAAGG` | 1 | One extra nucleotide at the 5′ end. No UMIs. |
| `custom` | set `trim.adapter` | set `trim.cut_5p` | |

| Key | Default | Content |
|---|---|---|
| `trim.adapter` | from the kit | 3′ adapter sequence. |
| `trim.cut_5p` | from the kit | Nucleotides to remove at the 5′ end. |
| `trim.nextseq_trim` | 0 | Quality cutoff for `cutadapt --nextseq-trim` (two-colour chemistry, removes poly-G). 0 turns it off. |
| `trim.min_len` | 16 | Minimum insert length. Shorter reads are removed. |

Kits with UMIs (for example, QIAseq miRNA) are not supported yet.

## References

The workflow downloads the references into `resources/`. Each source in `refs` is a URL or a local path. Use local
paths on a cluster without internet access. A source with `url` and `md5` is checked against the md5 sum.

| Key | Content |
|---|---|
| `species.mirbase`, `species.mirgenedb` | Species prefix of the miRNA names (human: `hsa`, `Hsa`). |
| `refs.mirbase` | miRBase mature and hairpin FASTA (default 22.1). |
| `refs.mirgenedb` | MirGeneDB mature, star and precursor FASTA (default 3.0). MirGeneDB has no fixed download for each release, so the workflow writes the md5 of the download to the log. For a fixed version, give local copies. |
| `refs.gtrnadb` | GtRNAdb mature tRNAs and the tRNA gene loci. |
| `refs.ensembl` | Ensembl genome, gene annotation (GTF) and ncRNA FASTA (default release 113, GRCh38). |
| `refs.rrna` | The single-copy rRNA reference: the rDNA repeat (NCBI accession), 5S with genomic flanks, the mitochondrial rRNAs, and the rRNA regions (`features`). |
| `refs.rdna_loci` | rDNA copies in the genome that the annotation does not mark as rRNA. Genome hits in these loci get the class `rRNA_genomic`. |
| `refs.contaminants` | Contaminant genomes: NCBI accession, read class and species. The default list has 6 *Mycoplasma* species, *E. coli* and phiX. |
| `refs.ncbi_dir` | Optional directory with local copies of the NCBI records: `<accession>.fa`, and `<accession>.ft` (feature table) for the contaminants. |
| `genome_index` | Optional prefix of a bowtie 1 index of the genome that already exists. The index build needs 32 GB of memory and several hours. |

The default references are for human. For another species, change `species` and the references. The code has no
species names.

## Read classification

| Key | Default | Content |
|---|---|---|
| `classify.mismatches` | 1 | Mismatches for the alignment to each reference (bowtie 1 `-v`). |
| `classify.max_hits` | 20 | Hits for each read on the multi-copy references (ncRNA, genome). |
| `classify.arm_start_tol` | 3 | A hairpin read counts for a mature miRNA if its 5′ end is this many nucleotides or less from the 5′ end of the mature sequence. |
| `classify.arm_end_tol` | 5 | The same for the 3′ end. |
| `classify.contaminant_min_len` | 18 | A read of this length or longer is a contaminant read if its contaminant hit has fewer mismatches than its human hits. Shorter reads need an exact contaminant match. |
| `classify.end_len` | 25 | Length of the read ends that align to the references, for reads without a full-length hit. |
| `classify.end_len_exact` | 18 | The same, for ends that must match exactly. |

[Read classification](https://katarinagresova.github.io/wf-miRNA-EV/classification/) explains these rules.

## miRNA counts

| Key | Default | Content |
|---|---|---|
| `mirna.databases` | `[mirbase, mirgenedb]` | miRNA databases. The first one is the primary database: the read classes and `library_summary.tsv` use it. |
| `mirna.min_len`, `mirna.max_len` | 16, 28 | Read lengths that can count for a miRNA. |
| `mirna.exact_max_len` | 18 | Reads of this length or shorter count for a miRNA only with 0 mismatches. |
| `mirna.nta_max` | 3 | Maximum of non-templated 3′ nucleotides (for example, adenylation). |
| `mirna.detected_min_reads` | `[1, 5, 10]` | Read thresholds for the detected miRNAs. |
| `mirna.off_target_mismatches` | 3 | Mismatches of the miRNA-only mapping of the off-target check (maximum 3). |

[miRNA counts](https://katarinagresova.github.io/wf-miRNA-EV/mirna_counters/) explains these rules.

## Other keys

| Key | Default | Content |
|---|---|---|
| `collapse.top_n` | 50 | The most abundant sequences of each library that go into `results/collapse/top_sequences.tsv`. |
