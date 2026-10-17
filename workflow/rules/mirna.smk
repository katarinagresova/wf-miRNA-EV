# Step 6: miRNA counts, summary and off-target check. Each miRNA counter has its rule file in counters/; it writes
# results/mirna/{sample}/{counter}_{db}.tsv (mirna, reads, members) and, optionally, the per-read table
# results/mirna/{sample}/{counter}_{db}.reads.tsv.gz (sequence, mirna, reads).

if MIRNA["counter"] in MIRNA["cross_check"]:
    raise WorkflowError(f"mirna.cross_check contains the main counter {MIRNA['counter']}")
COUNTERS = [MIRNA["counter"], *MIRNA["cross_check"]]
# Counters that write a per-read table; each counter file adds its name
PER_READ_COUNTERS = set()

for counter in COUNTERS:
    include: f"counters/{counter}.smk"

MIRNA_SUMMARY = {
    k: v
    for k, v in {
        "summary": "results/library_summary.tsv",
        "counts": expand("results/mirna/counts_{db}.tsv", db=MIRNA_DBS),
        "cross_check": expand("results/mirna/cross_check/{counter}_{db}.tsv", counter=MIRNA["cross_check"], db=MIRNA_DBS),
        "agreement": "results/mirna/method_agreement.tsv",
        "off_target": "results/mirna/off_target.tsv",
    }.items()
    if v
}


rule mirna_reads:
    # Input of the counters: the reads in the miRNA length range, without the reads classed as a contaminant (human
    # reads only; miRge3 had counted mycoplasma fragments as miRNA isomiRs)
    input:
        script=workflow.source_path("../scripts/mirna_reads.py"),
        reads="results/trim/{sample}/trimmed.fastq.gz",
        collapsed="results/collapse/{sample}/collapsed.fa.gz",
        classes="results/classify/{sample}/read_classes.tsv.gz",
    output:
        temp("results/mirna/{sample}/reads_{db}.fastq.gz"),
    params:
        min_len=MIRNA["min_len"],
        max_len=MIRNA["max_len"],
    log:
        "logs/mirna_reads/{db}/{sample}.log",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240,
    conda:
        "../envs/py.yaml"
    shell:
        """
        seqkit seq -m {params.min_len} -M {params.max_len} -j {threads} {input.reads} 2> {log} \
            | python {input.script} --classes {input.classes} --reads {input.collapsed} --db {wildcards.db} 2>> {log} \
            | pigz -p {threads} > {output}
        """


rule off_target_align:
    # The miRNA-only mapping of the off-target check: all kept reads in the miRNA length range (contaminants
    # included), on the miRNA references only, with up to off_target_mismatches mismatches
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        index="resources/bowtie/{db}/.done",
    output:
        # read rank, reference, 0-based offset, read length, mismatches; best stratum only, sense strand only
        temp("results/mirna/{sample}/off_target_{db}.hits.tsv.gz"),
    params:
        prefix="resources/bowtie/{db}/{db}",
        min_len=MIRNA["min_len"],
        max_len=MIRNA["max_len"],
        mm=MIRNA["off_target_mismatches"],
    log:
        "logs/off_target_align/{db}/{sample}.log",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=240,
    conda:
        "../envs/align.yaml"
    shell:
        """
        seqkit seq -m {params.min_len} -M {params.max_len} {input.reads} 2> {log} \
            | bowtie -f -v {params.mm} -a --best --strata --norc -p {threads} -x {params.prefix} - 2>> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print substr($1, 2, index($1, "_") - 2), $3, $4, length($5), mm}}' \
            | pigz -p 2 > {output}
        """


rule off_target:
    input:
        script=workflow.source_path("../scripts/off_target.py"),
        hits="results/mirna/{sample}/off_target_{db}.hits.tsv.gz",
        arms="resources/refs/{db}_arms.tsv",
        reads="results/collapse/{sample}/collapsed.fa.gz",
        classes="results/classify/{sample}/read_classes.tsv.gz",
        mirna_reads="results/classify/{sample}/mirna_{db}.reads.tsv.gz",
        counter=lambda wc: (
            f"results/mirna/{wc.sample}/{MIRNA['counter']}_{wc.db}.reads.tsv.gz"
            if MIRNA["counter"] in PER_READ_COUNTERS
            else []
        ),
    output:
        # Each read with a miRNA-only hit: its mismatches, its result and the entry of the main counter
        reads="results/mirna/{sample}/off_target_{db}.tsv.gz",
        # The reads of the main counter by their class in step 5 (empty if the counter has no per-read table)
        counter="results/mirna/{sample}/off_target_{db}.counter.tsv",
    params:
        counter=lambda wc, input: f"--counter-reads {input.counter}" if input.counter else "",
        min_len=MIRNA["min_len"],
        max_len=MIRNA["max_len"],
        start_tol=CLASSIFY["arm_start_tol"],
        end_tol=CLASSIFY["arm_end_tol"],
    log:
        "logs/off_target/{db}/{sample}.log",
    resources:
        mem_mb=32000,
        runtime=60,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} --db {wildcards.db} --hits {input.hits} --arms {input.arms} --reads {input.reads} \
            --classes {input.classes} --mirna-reads {input.mirna_reads} {params.counter} \
            --min-len {params.min_len} --max-len {params.max_len} \
            --start-tol {params.start_tol} --end-tol {params.end_tol} \
            --out-reads {output.reads} --out-counter {output.counter} &> {log}
        """


rule mirna_summary:
    input:
        script=workflow.source_path("../scripts/mirna_summary.py"),
        samples=config["samples"],
        composition="results/classify/composition.tsv",
        counts=expand("results/mirna/{sample}/{counter}_{db}.tsv", sample=SAMPLES.index, counter=COUNTERS, db=MIRNA_DBS),
        off_target=expand("results/mirna/{sample}/off_target_{db}.tsv.gz", sample=SAMPLES.index, db=MIRNA_DBS),
        off_target_counter=expand("results/mirna/{sample}/off_target_{db}.counter.tsv", sample=SAMPLES.index, db=MIRNA_DBS),
    output:
        **MIRNA_SUMMARY,
    params:
        dbs=MIRNA_DBS,
        counter=MIRNA["counter"],
        cross_check=" ".join(f"--cross-check {c}" for c in MIRNA["cross_check"]),
        detected_min=MIRNA["detected_min_reads"],
    log:
        "logs/mirna_summary.log",
    resources:
        mem_mb=16000,
        runtime=30,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
            --samples {input.samples} \
            --dbs {params.dbs} \
            --counter {params.counter} {params.cross_check} \
            --composition {input.composition} \
            --counts {input.counts} \
            --off-target {input.off_target} \
            --off-target-counter {input.off_target_counter} \
            --detected-min {params.detected_min} \
            --out-dir $(dirname {output.agreement}) \
            --out-summary {output.summary} \
            --out-agreement {output.agreement} \
            --out-off-target {output.off_target} &> {log}
        """
