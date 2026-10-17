# Step 3: collapse identical trimmed reads into sequence counts.


rule collapse:
    input:
        "results/trim/{sample}/trimmed.fastq.gz",
    output:
        # sequence<TAB>count, sorted by count (descending)
        "results/collapse/{sample}/counts.tsv.gz",
    log:
        "logs/collapse/{sample}.log",
    threads: 8
    resources:
        mem_mb=24000,
        runtime=480,
    conda:
        "../envs/trim.yaml"
    shell:
        """
        tmp=$(dirname {output})/sort_tmp
        mkdir -p $tmp
        (
        printf 'sequence\\tcount\\n'
        pigz -dc -p 2 {input} \
            | awk 'NR % 4 == 2' \
            | LC_ALL=C sort -S 16G --parallel {threads} -T $tmp \
            | uniq -c \
            | LC_ALL=C sort -k1,1nr -k2,2 -S 4G --parallel {threads} -T $tmp \
            | awk -v OFS='\\t' '{{print $2, $1}}'
        ) 2> {log} | pigz -p {threads} > {output}
        rm -rf $tmp
        """


rule collapse_fasta:
    input:
        "results/collapse/{sample}/counts.tsv.gz",
    output:
        # >s<rank>_x<count>, input for the aligners
        "results/collapse/{sample}/collapsed.fa.gz",
    log:
        "logs/collapse_fasta/{sample}.log",
    threads: 2
    resources:
        runtime=120,
    conda:
        "../envs/trim.yaml"
    shell:
        """
        pigz -dc {input} \
            | awk 'NR > 1 {{printf ">s%d_x%d\\n%s\\n", NR - 1, $2, $1}}' \
            | pigz -p {threads} > {output} 2> {log}
        """


rule collapse_summary:
    input:
        script="workflow/scripts/collapse_summary.py",
        counts=expand("results/collapse/{sample}/counts.tsv.gz", sample=SAMPLES.index),
        samples=config["samples"],
    output:
        summary="results/qc/collapse_summary.tsv",
        top="results/qc/top_sequences.tsv",
    params:
        top_n=config["collapse"]["top_n"],
    log:
        "logs/collapse_summary.log",
    resources:
        mem_mb=16000,
        runtime=120,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
            --counts {input.counts} \
            --samples {input.samples} \
            --top-n {params.top_n} \
            --summary {output.summary} \
            --top {output.top} &> {log}
        """
