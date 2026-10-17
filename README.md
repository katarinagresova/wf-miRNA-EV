# wf-miRNA-EV

A Snakemake workflow for small RNA sequencing data from extracellular vesicles (EVs). It also accepts whole-cell
libraries, so that you can compare EVs with their source cells.

In EV libraries, most reads are not miRNA: rRNA fragments of miRNA length and contaminant RNA (for example,
*Mycoplasma*) can be most of the library. Thus the workflow gives each read one class (miRNA, rRNA, tRNA, other small
RNA, genome classes, contaminants), not only the miRNA reads. It marks the miRNA entries that also match another class
equally well.

**Status: in development.** Steps 1–7 are ported from the exploratory repository miRNA-PD and give the same results
on its data. The config in `config/` is still the miRNA-PD one, and there is no example sample sheet yet.

## What it does

The input is a set of single-end FASTQ files (gzip) and a sample sheet. The steps:

1. QC of the raw reads (FastQC, seqkit, MultiQC).
2. Adapter trimming and length filter (cutadapt). All reads of 16 nt or longer are kept.
3. Collapse of identical reads.
4. Download and build of the references and indices (miRBase, MirGeneDB, rRNA, tRNA, Ensembl ncRNA and genome,
   contaminant genomes).
5. Classification of each read (bowtie 1), contaminants included.
6. miRNA counts.
7. rRNA fragment tables.
8. QC report (MultiQC).

## Setup

Create the environment that runs Snakemake. Strict channel priority makes the conda solve much faster.

```bash
CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yaml
conda activate wf-miRNA-EV
```

The rules use their own conda environments (`workflow/envs/`). Snakemake creates them at the first run.

## Run

Local (laptop or workstation):

```bash
snakemake -n                          # dry run
snakemake --use-conda --cores 8
```

SLURM cluster: first, set the partition (and the account, if necessary) in `profiles/slurm/config.yaml`. Then
create the rule environments on the login node and submit the Snakemake process as a SLURM job. The process submits
one SLURM job for each rule.

```bash
CONDA_CHANNEL_PRIORITY=strict snakemake --profile profiles/slurm --conda-create-envs-only
mkdir -p logs
sbatch --job-name=wf-miRNA-EV --time=2-00:00:00 --mem=4G --output=logs/snakemake_%j.log \
    --wrap "snakemake --profile profiles/slurm"
```

## Repository layout

```
├── environment.yaml     # Snakemake and the SLURM executor plugin
├── config/              # example config and sample sheet
├── workflow/
│   ├── Snakefile
│   ├── rules/           # one file per step
│   ├── scripts/
│   ├── envs/            # conda environment of each rule
│   └── schemas/         # config and sample sheet schemas
├── profiles/slurm/      # example SLURM profile
├── docs/
├── tests/               # unit tests (pytest)
├── .test/               # synthetic test data
└── .github/workflows/   # CI
```

The outputs go to `results/`, the references to `resources/` and the logs to `logs/`.

## License

Apache License 2.0. See [LICENSE](LICENSE).
