# Step 8: QC report. MultiQC with FastQC (raw and trimmed reads), cutadapt and the sections of report_tables.py.

REPORT_TABLES = [
    "general_stats",
    "insert_lengths",
    "read_classes",
    "contaminants",
    "mirna",
    "off_target",
    "off_target_length",
]


rule report_tables:
    input:
        script=workflow.source_path("../scripts/report_tables.py"),
        summary="results/library_summary.tsv",
        composition="results/classify/composition.tsv",
        contaminants="results/classify/contaminants.tsv",
        insert_lengths="results/trim/insert_lengths.tsv",
        off_target=expand("results/mirna/{sample}/off_target_{db}.tsv.gz", sample=SAMPLES.index, db=MIRNA_DBS),
    output:
        expand("results/report/custom/{name}_mqc.json", name=REPORT_TABLES),
    params:
        dbs=MIRNA_DBS,
    log:
        "logs/report_tables.log",
    resources:
        mem_mb=16000,
        runtime=30,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
            --dbs {params.dbs} \
            --summary {input.summary} \
            --composition {input.composition} \
            --contaminants {input.contaminants} \
            --insert-lengths {input.insert_lengths} \
            --off-target {input.off_target} \
            --out-dir results/report/custom &> {log}
        """


rule multiqc:
    input:
        config=workflow.source_path("../report/multiqc_config.yaml"),
        fastqc_raw=expand("results/qc/fastqc_raw/{sample}_fastqc.zip", sample=SAMPLES.index),
        fastqc_trimmed=expand("results/qc/fastqc_trimmed/{sample}_fastqc.zip", sample=SAMPLES.index),
        cutadapt=expand("results/trim/{sample}/{sample}.cutadapt.log", sample=SAMPLES.index),
        custom=rules.report_tables.output,
    output:
        "results/report/multiqc_report.html",
    log:
        "logs/multiqc.log",
    resources:
        mem_mb=8000,
        runtime=60,
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc --force --no-version-check --config {input.config} --filename multiqc_report.html \
            --outdir results/report {input.fastqc_raw} {input.fastqc_trimmed} {input.cutadapt} {input.custom} &> {log}
        """
