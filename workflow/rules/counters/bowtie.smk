# miRNA counter bowtie: the mature miRNA counts of step 5 (quantify_mirna.py). This file only makes the interface
# tables from the step 5 tables.

PER_READ_COUNTERS.add("bowtie")


rule bowtie_table:
    input:
        script=workflow.source_path("../../scripts/bowtie_table.py"),
        counts="results/classify/{sample}/mirna_{db}.tsv",
        mirna_reads="results/classify/{sample}/mirna_{db}.reads.tsv.gz",
        reads="results/collapse/{sample}/collapsed.fa.gz",
    output:
        table="results/mirna/{sample}/bowtie_{db}.tsv",
        reads="results/mirna/{sample}/bowtie_{db}.reads.tsv.gz",
    log:
        "logs/bowtie_table/{db}/{sample}.log",
    resources:
        mem_mb=8000,
        runtime=30,
    conda:
        "../../envs/py.yaml"
    shell:
        """
        python {input.script} --counts {input.counts} --mirna-reads {input.mirna_reads} --reads {input.reads} \
            --out-table {output.table} --out-reads {output.reads} &> {log}
        """


rule bowtie_version:
    output:
        "results/mirna/bowtie_version.txt",
    log:
        "logs/bowtie_version.log",
    conda:
        "../../envs/align.yaml"
    shell:
        "(bowtie --version | awk 'NR == 1 {{print \"bowtie\", $3}}' > {output}) 2> {log}"
