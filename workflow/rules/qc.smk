# Step 1: raw read QC and read counts.


rule fastqc_raw:
    input:
        lambda wc: SAMPLES.loc[wc.sample, "fastq"],
    output:
        html="results/qc/fastqc_raw/{sample}_fastqc.html",
        zip="results/qc/fastqc_raw/{sample}_fastqc.zip",
    log:
        "logs/fastqc_raw/{sample}.log",
    threads: 2
    resources:
        mem_mb=4000,
        runtime=240,
    conda:
        "../envs/qc.yaml"
    shadow:
        "minimal"
    shell:
        # Symlink so the report is named after the sample, not the FASTQ file
        """
        ln -s $(realpath {input}) {wildcards.sample}.fastq.gz
        fastqc --threads {threads} --memory 2048 --outdir results/qc/fastqc_raw {wildcards.sample}.fastq.gz &> {log}
        """


rule fastqc_trimmed:
    input:
        "results/trim/{sample}/trimmed.fastq.gz",
    output:
        html="results/qc/fastqc_trimmed/{sample}_fastqc.html",
        zip="results/qc/fastqc_trimmed/{sample}_fastqc.zip",
    log:
        "logs/fastqc_trimmed/{sample}.log",
    threads: 2
    resources:
        mem_mb=4000,
        runtime=240,
    conda:
        "../envs/qc.yaml"
    shadow:
        "minimal"
    shell:
        # Symlink so the report is named after the sample, not the FASTQ file
        """
        ln -s $(realpath {input}) {wildcards.sample}.fastq.gz
        fastqc --threads {threads} --memory 2048 --outdir results/qc/fastqc_trimmed {wildcards.sample}.fastq.gz &> {log}
        """


rule seqkit_stats_raw:
    input:
        lambda wc: SAMPLES.loc[wc.sample, "fastq"],
    output:
        "results/qc/seqkit_raw/{sample}.tsv",
    log:
        "logs/seqkit_raw/{sample}.log",
    threads: 2
    resources:
        runtime=120,
    conda:
        "../envs/qc.yaml"
    shell:
        "seqkit stats --all --tabular --threads {threads} {input} > {output} 2> {log}"


rule read_counts:
    input:
        script=workflow.source_path("../scripts/read_counts.py"),
        stats=expand("results/qc/seqkit_raw/{sample}.tsv", sample=SAMPLES.index),
        samples=config["samples"],
    output:
        "results/qc/raw_reads.tsv",
    log:
        "logs/read_counts.log",
    conda:
        "../envs/py.yaml"
    shell:
        # shell rather than script: the script preamble puts Snakemake's site-packages first on sys.path
        """
        python {input.script} \
            --stats {input.stats} \
            --samples {input.samples} \
            --output {output} &> {log}
        """
