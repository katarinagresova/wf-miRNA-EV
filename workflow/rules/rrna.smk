# Step 7: rRNA fragment tables.


rule rrna_fragments:
    input:
        script=workflow.source_path("../scripts/rrna_fragments.py"),
        reads="results/collapse/{sample}/collapsed.fa.gz",
        classes="results/classify/{sample}/read_classes.tsv.gz",
        hits="results/classify/{sample}/hits/rrna.hits.tsv.gz",
        contaminant="results/classify/{sample}/hits/contaminant.hits.tsv.gz",
        features="resources/refs/rrna_features.bed",
    output:
        "results/rrna/{sample}/fragments.tsv.gz",
    params:
        # The read classes of the primary miRNA database
        db=MIRNA_DBS[0],
        min_len=CLASSIFY["contaminant_min_len"],
    log:
        "logs/rrna_fragments/{sample}.log",
    resources:
        mem_mb=16000,
        runtime=60,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} --reads {input.reads} --classes {input.classes} --hits {input.hits} \
            --contaminant {input.contaminant} --features {input.features} --db {params.db} \
            --contaminant-min-len {params.min_len} --output {output} &> {log}
        """
