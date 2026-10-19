# Changelog

The notable changes of each release. The output paths and columns in [docs/outputs.md](docs/outputs.md) change only
with a new minor version, and this file lists each such change.

## 0.1.0 (2026-10-19)

First release.

- Steps 1–8 for single-end small RNA libraries (gzipped FASTQ): raw read QC (FastQC, SeqKit), adapter trimming and
  length filter (cutadapt), collapse of identical reads, references, read classification (bowtie 1), miRNA counts
  with an off-target check, rRNA fragments, and a QC report (MultiQC).
- One class for each read: miRNA, rRNA, tRNA, other small RNA, genome or contaminant (for example, *Mycoplasma*).
  Reads without a full-length hit can get a class from the hits of their ends.
- miRNA counts for miRBase and MirGeneDB. The off-target check shows the reads that a miRNA-only mapping would count,
  but that match another class better or equally well.
- The kit preset `realseq_biofluids`. For other kits without UMIs, set the trim values in the config.
- References from URLs or local paths. A source with an md5 sum in the config is checked against it.
- Use as a Snakemake module, with the settings under the config key `wf_mirna_ev`.
- A SLURM profile (`profiles/slurm`).
- Unit tests and a synthetic test data set with the known class and miRNA of each read. GitHub Actions runs them at
  each push.
- A documentation site: <https://katarinagresova.github.io/wf-miRNA-EV/>.

Known limitations:

- Human only. The default references are for human, and two parts of the code work only for human (see
  [config/README.md](config/README.md)).
- Kits with UMIs (for example, QIAseq miRNA) are not supported.
