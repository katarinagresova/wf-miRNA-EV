# Outputs

All outputs are in `results/`. Each step has one directory. `results/<step>/{sample}/` has the files of one library,
and `results/<step>/` has the tables for all libraries. The tables are tab-separated. In the tables for all libraries,
the extra columns of the sample sheet come after `sample`.

This page lists the outputs that other analyses can use. Their paths and columns change only with a new minor
version. Other files in `results/` and `resources/` (for example, the alignment hits) are intermediate files.

`{db}` is a miRNA database of `mirna.databases` (`mirbase`, `mirgenedb`). The primary database is the first one.

## Library summary

`results/library_summary.tsv` is the main table. It has one row for each library.

| Column | Content |
|---|---|
| `kept_reads` | Reads with an insert of `trim.min_len` nt or longer, without the 5′ adapter artefacts. All percentages below use this number. |
| `reads_16_28` | Kept reads in the miRNA length range. The name changes with `mirna.min_len` and `mirna.max_len` (`reads_<min_len>_<max_len>`). |
| `pct_<class>` | Percent of the kept reads in each [read class](classification.md#read-classes) (primary database). |
| `{db}_mirna_reads` | Reads of all miRNA entries. |
| `{db}_pct_of_kept` | `{db}_mirna_reads` as a percent of the kept reads. |
| `{db}_pct_tied` | Tied miRNA reads as a percent of `{db}_mirna_reads`. A tied read has the class miRNA, but another class has an equally good hit (`tied` in the [off-target check](classification.md#mirna-off-target-check)). |
| `{db}_pct_tied_rRNA` | Tied miRNA reads with an equally good rRNA hit, as a percent of `{db}_mirna_reads`. |
| `{db}_detected_ge<n>` | miRNA entries with `<n>` or more reads, for each `<n>` of `mirna.detected_min_reads`. |
| `{db}_off_target_reads` | Reads of the miRNA-only mapping of the [off-target check](classification.md#mirna-off-target-check). |
| `{db}_off_target_pct_<result>` | Percent of these reads with each result of the off-target check. |

## QC report

`results/report/multiqc_report.html` is a MultiQC report with these sections:

| Section | Content |
|---|---|
| General statistics | Kept reads, percent miRNA, rRNA and contaminant, detected miRNAs, percent off-target (primary database), and the FastQC and cutadapt values. |
| Insert lengths | Insert length distribution after the adapter trim, as a percent of the raw reads. Length 0 is the adapter dimers. |
| Read classes | Kept reads of each read class (primary database). |
| Contaminants | Kept reads of each contaminant genome. |
| miRNA counts | miRNA reads, tied miRNA reads and detected miRNAs of each database. |
| miRNA off-target check | Reads of the miRNA-only mapping by result, for each library and database. |
| miRNA off-target check by read length | Percent off-target reads for each read length. Short reads match a miRNA by chance more often. |
| FastQC (raw reads), cutadapt, FastQC (trimmed reads) | The standard MultiQC sections. |

## Step 1: raw reads

| File | Content |
|---|---|
| `results/qc/raw_reads.tsv` | `reads`, `min_len`, `max_len` and `Q30_pct` of the raw FASTQ file of each library. |
| `results/qc/fastqc_raw/{sample}_fastqc.html` | FastQC report of the raw reads. |
| `results/qc/fastqc_trimmed/{sample}_fastqc.html` | FastQC report of the trimmed reads. |

## Step 2: adapter trimming

| File | Content |
|---|---|
| `results/trim/summary.tsv` | Read fates of each library (see below). |
| `results/trim/insert_lengths.tsv` | Insert length distribution: one row for each length (0–99 nt), one column for each library, as a percent of the raw reads. Only reads with an adapter. |
| `results/trim/insert_lengths.pdf` | Plot of `insert_lengths.tsv`, one panel for each batch. |
| `results/trim/{sample}/trimmed.fastq.gz` | Kept reads, after the adapter trim, the length filter and the 5′ adapter filter. |
| `results/trim/{sample}/no_adapter.fastq.gz` | Reads without an adapter. |

Columns of `summary.tsv`. Each count has a second column with the suffix `_pct`, the percent of the raw reads. The
names with 15 and 16 change with `trim.min_len`. The names with 28 and 29 change with `mirna.max_len`.

| Column | Content |
|---|---|
| `raw_reads` | Reads in the FASTQ file. |
| `kept_reads`, `kept_pct` | Reads with an adapter and an insert of `trim.min_len` nt or longer, without `adapter_5p_artefact`. |
| `with_adapter_pct` | Percent of the raw reads with an adapter. |
| `adapter_dimer` | Reads with an empty insert. |
| `insert_1_15`, `insert_16_28`, `insert_29_50`, `insert_51plus` | Reads with an adapter, by insert length. |
| `no_adapter_lt16`, `no_adapter_ge16` | Reads without an adapter, shorter or longer than `trim.min_len` after the quality trim. Short reads without an adapter are often poly-G reads. |
| `adapter_5p_artefact` | Reads with an insert of `trim.min_len` nt or longer that has 15 nt or more of the 5′ adapter (`trim.adapter_5p`). The workflow removes them. They are also in the `insert_*` columns. |
| `quality_trimmed_bp_pct` | Percent of the bases that the quality trim removed. |
| `cut_nt_<A,C,G,T,N>_pct` | The nucleotide that `trim.cut_5p` removed at the 5′ end, as a percent of the raw reads. |

## Step 3: collapse

| File | Content |
|---|---|
| `results/collapse/{sample}/counts.tsv.gz` | Each sequence (`sequence`) and its read count (`count`), most abundant first. |
| `results/collapse/{sample}/collapsed.fa.gz` | The same sequences as FASTA. The name `s<rank>_x<count>` gives the rank in `counts.tsv.gz` and the count. |
| `results/collapse/summary.tsv` | Library complexity (see below). |
| `results/collapse/top_sequences.tsv` | The `collapse.top_n` most abundant sequences of each library: `sequence`, `length`, `n_samples_in_top` (libraries in which the sequence is in the top), `mean_cpm`, and the CPM in each library. |

Columns of `summary.tsv`: `reads`, `unique_seqs`, `unique_pct` (unique sequences as a percent of the reads),
`singletons_pct_of_unique` and `singletons_pct_of_reads` (sequences with 1 read), `top<n>_pct_of_reads` (reads of the
1, 10, 100 and 1000 most abundant sequences), `median_len_reads` and `median_len_unique`.

## Step 4: references

| File | Content |
|---|---|
| `resources/refs/rrna.fa` | The single-copy rRNA reference: the 45S rDNA repeat, 5S with genomic flanks, and the mitochondrial rRNAs. |
| `resources/refs/rrna_features.bed` | The rRNA regions on this reference (for example, 18S, ITS1, 28S). |

## Step 5: read classification

| File | Content |
|---|---|
| `results/classify/{sample}/read_classes.tsv.gz` | The class of each collapsed sequence: `rank` (the rank in `collapsed.fa.gz`), `count`, `length`, and `class_{db}` for each database. |
| `results/classify/composition.tsv` | Reads for each library, database, class and length bin (see below). |
| `results/classify/contaminants.tsv` | Reads for each library and contaminant genome (primary database): `class`, `accession`, `species`, `reads`, `unique_seqs`, `pct_of_kept`. |
| `results/classify/top_sequences_classes.tsv` | The sequences of `top_sequences.tsv` with their `length`, their `class` (primary database), the `mismatches` and the number (`n_hits`) of their best hits, and a short `label` (for example, the miRNA name, the rRNA region and position, or the tRNA fragment type). |
| `results/classify/{sample}/composition.tsv` | Reads of one library for each database, class and length bin: `db`, `class`, `len_bin`, `reads`, `unique_seqs`. |
| `results/classify/{sample}/contaminants.tsv` | Reads of one library for each database, contaminant genome and length bin: `db`, `class`, `accession`, `species`, `len_bin`, `reads`, `unique_seqs`. |

Columns of `composition.tsv`: `db`, `class`, `len_bin` (16–28, 29–50 or 51+ nt), `reads`, `unique_seqs`,
`pct_of_kept`, and `pct_of_human`. `pct_of_human` divides by the human reads: the kept reads without the contaminant
reads. Use `pct_of_human` to compare libraries, because the contaminant share can be very different between
libraries.

## Step 6: miRNA counts and off-target check

| File | Content |
|---|---|
| `results/mirna/counts_{db}.tsv` | Count matrix: one row for each miRNA entry (`mirna`), one column for each library. The most abundant entries are first. |
| `results/mirna/{sample}/bowtie_{db}.tsv` | The counts of one library (see below). |
| `results/mirna/{sample}/bowtie_{db}.reads.tsv.gz` | The counted reads: `sequence`, `mirna` (entry) and `reads`. |
| `results/mirna/{sample}/off_target_{db}.tsv.gz` | Each read of the miRNA-only mapping (see below). |
| `results/mirna/off_target.tsv` | The off-target check for each database entry, summed over all libraries (see below). |

An entry is a mature miRNA. Mature names with the same sequence are one entry, with the names joined by `/` (for
example, `hsa-miR-199a-3p/hsa-miR-199b-3p`). [miRNA counts](mirna_counters.md) gives the counting rules.

Columns of `bowtie_{db}.tsv`:

| Column | Content |
|---|---|
| `mirna` | The entry. |
| `reads` | The count. A read on two entries counts 0.5 for each. |
| `members` | The mature names of the entry, joined by `;`. |
| `unique_seqs` | Different sequences that count for the entry. |
| `canonical_reads` | Reads with exactly the mature sequence. |
| `nta_reads` | Reads with non-templated 3′ nucleotides. |

Columns of `off_target_{db}.tsv.gz`: `sequence`, `length`, `count`, `mm` (mismatches of the miRNA-only hit), `mirna`
(the mature name of the hit), `reads` (the count, divided between the hits of the read), `class` (the class in step
5), `result` (the [result](classification.md#mirna-off-target-check) of the off-target check), `tied_classes` (the
other classes with an equally good hit, joined by `,`, empty if the result is not `tied`) and `counter_mirna` (the
entry that counts the read, empty if no entry counts it).

Columns of `off_target.tsv`:

| Column | Content |
|---|---|
| `db`, `mirna` | Database and mature name. |
| `counter_mirna` | The entry that contains this mature name. |
| `reads`, `n_samples` | Reads of the miRNA-only mapping, and the libraries with such reads. |
| `tied_frac` | `tied` reads, divided by the `miRNA` and `tied` reads. |
| `off_target_frac` | Reads with another class in step 5, divided by `reads`. |
| `reads_<result>` | Reads with each result. |

## Step 7: rRNA fragments

`results/rrna/{sample}/fragments.tsv.gz` lists the reads with the class rRNA (primary database) and their positions on
the single-copy rRNA reference.

| Column | Content |
|---|---|
| `rank`, `sequence`, `count`, `length` | The read. |
| `ref`, `start`, `end` | Position on the rRNA reference (0-based start, end not included). |
| `region` | rRNA region at the middle of the read. |
| `mm` | Mismatches. |
| `n_pos` | Positions with the same number of mismatches. A read has one row for each position. Use `count / n_pos` as the weight of a row. |
| `contaminant_tie` | 1 if the read is `classify.contaminant_min_len` nt or longer, and matches a contaminant genome with the same number of mismatches or fewer. The read can be a contaminant rRNA (for example, *Mycoplasma*). |

## Version

`results/wf_version.txt` gives the workflow version, its git commit, and the miRNA counter with its version. The
workflow writes this file again in each run that changes an output, so the commit is that of the last such run.
