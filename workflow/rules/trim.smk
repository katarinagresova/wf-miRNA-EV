# Step 2: adapter trimming, length filter and trimming statistics.

# Library kit presets; the values of config["trim"] override them
KITS = {
    # RealSeq-Biofluids: one extra nt at the 5' end, no UMIs
    "realseq_biofluids": {"adapter": "TGGAATTCTCGGGTGCCAAGG", "cut_5p": 1},
    "custom": {},
}
TRIM = {**KITS[config["kit"]], **config["trim"]}


rule cutadapt:
    input:
        lambda wc: SAMPLES.loc[wc.sample, "fastq"],
    output:
        # All reads with an adapter, including empty inserts (adapter dimers), so they can be counted
        adapter=temp("results/trim/{sample}/adapter.fastq.gz"),
        no_adapter="results/trim/{sample}/no_adapter.fastq.gz",
        json="results/trim/{sample}/cutadapt.json",
        report="results/trim/{sample}/{sample}.cutadapt.log",
    params:
        adapter=TRIM["adapter"],
        cut_5p=TRIM["cut_5p"],
        nextseq_trim=f"--nextseq-trim {TRIM['nextseq_trim']}" if TRIM["nextseq_trim"] else "",
    log:
        "logs/cutadapt/{sample}.log",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=480,
    conda:
        "../envs/trim.yaml"
    shell:
        # --rename keeps the removed 5' nt in the header (cut=<nt>) for the first-nt statistics
        """
        cutadapt --cores {threads} \
            --cut {params.cut_5p} \
            {params.nextseq_trim} \
            --adapter {params.adapter} \
            --rename '{{id}} cut={{cut_prefix}}' \
            --untrimmed-output {output.no_adapter} \
            --output {output.adapter} \
            --json {output.json} \
            {input} > {output.report} 2> {log}
        """


rule length_filter:
    input:
        "results/trim/{sample}/adapter.fastq.gz",
    output:
        "results/trim/{sample}/trimmed.fastq.gz",
    params:
        min_len=TRIM["min_len"],
    log:
        "logs/length_filter/{sample}.log",
    threads: 4
    resources:
        runtime=240,
    conda:
        "../envs/trim.yaml"
    shell:
        "seqkit seq --min-len {params.min_len} --threads {threads} {input} -o {output} 2> {log}"


rule trim_hist:
    input:
        adapter="results/trim/{sample}/adapter.fastq.gz",
        no_adapter="results/trim/{sample}/no_adapter.fastq.gz",
    output:
        "results/trim/{sample}/hist.tsv",
    log:
        "logs/trim_hist/{sample}.log",
    threads: 2
    resources:
        runtime=240,
    conda:
        "../envs/trim.yaml"
    shell:
        # Long format: reads (adapter|no_adapter), metric (length|cut_nt), value, count
        """
        (
        printf 'reads\\tmetric\\tvalue\\tcount\\n'
        for f in adapter no_adapter; do
            [ $f = adapter ] && fq={input.adapter} || fq={input.no_adapter}
            pigz -dc -p {threads} $fq | awk -v f=$f -v OFS='\\t' '
                NR % 4 == 1 {{ n = split($2, a, "="); cut[a[n]]++ }}
                NR % 4 == 2 {{ len[length($0)]++ }}
                END {{
                    for (l in len) print f, "length", l, len[l]
                    for (c in cut) print f, "cut_nt", c, cut[c]
                }}'
        done
        ) > {output} 2> {log}
        """


rule trim_summary:
    input:
        script=workflow.source_path("../scripts/trim_summary.py"),
        json=expand("results/trim/{sample}/cutadapt.json", sample=SAMPLES.index),
        hist=expand("results/trim/{sample}/hist.tsv", sample=SAMPLES.index),
        samples=config["samples"],
    output:
        summary="results/trim/summary.tsv",
        lengths="results/trim/insert_lengths.tsv",
        plot="results/trim/insert_lengths.pdf",
    params:
        min_len=TRIM["min_len"],
    log:
        "logs/trim_summary.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
            --json {input.json} \
            --hist {input.hist} \
            --samples {input.samples} \
            --min-len {params.min_len} \
            --summary {output.summary} \
            --lengths {output.lengths} \
            --plot {output.plot} &> {log}
        """
