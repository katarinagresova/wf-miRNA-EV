# miRNA counter mirge3: miRge3.0 on the reads of rule mirna_reads. All miRge3.0-specific code is here: the library,
# the run, and the interface tables (mirge3_table.py: the merges, the translation of the library names to the miRNA
# database releases of the workflow, the per-read table).

PER_READ_COUNTERS.add("mirge3")

MIRGE3_ORGANISM = config["species"]["mirge3"]
MIRGE3_LIB_DIR = ref_source(REFS["miRge3"], f"resources/mirge3/{MIRGE3_ORGANISM}")
# miRge3 database names, and their spelling in the library file names
MIRGE3_DB = {"mirbase": "miRBase", "mirgenedb": "miRGeneDB"}
MIRGE3_LIB = {"mirbase": "miRBase", "mirgenedb": "MirGeneDB"}
# Mature sequences of the workflow's database releases, for the name translation
MATURE_FASTAS = {
    "mirbase": [ref_source(REFS["mirbase"]["mature"], "resources/download/mirbase/mature.fa")],
    "mirgenedb": [
        ref_source(REFS["mirgenedb"]["mature"], "resources/download/mirgenedb/mature.fa"),
        ref_source(REFS["mirgenedb"]["star"], "resources/download/mirgenedb/star.fa"),
    ],
}


localrules:
    download_mirge3_lib,


rule download_mirge3_lib:
    output:
        directory(f"resources/mirge3/{MIRGE3_ORGANISM}"),
    params:
        # Only for a URL; a local library in the config is not downloaded
        url=lambda wc: REFS["miRge3"]["url"],
        md5=lambda wc: REFS["miRge3"]["md5"],
        organism=MIRGE3_ORGANISM,
    log:
        "logs/download_mirge3_lib.log",
    conda:
        "../../envs/download.yaml"
    shell:
        """
        cd resources/mirge3
        [ -s {params.organism}.tar.gz ] && echo "{params.md5}  {params.organism}.tar.gz" | md5sum -c - &> /dev/null \
            || curl -sSfL --retry 3 -o {params.organism}.tar.gz "{params.url}" 2> ../../{log}
        echo "{params.md5}  {params.organism}.tar.gz" | md5sum -c - >> ../../{log} 2>&1
        tar xzf {params.organism}.tar.gz 2>> ../../{log}
        rm {params.organism}.tar.gz
        """


rule mirge3:
    # On the 16-28 nt human reads (already trimmed, so no adapter given)
    input:
        reads="results/mirna/{sample}/reads_{db}.fastq.gz",
        lib=MIRGE3_LIB_DIR,
    output:
        counts="results/mirna/{sample}/mirge3_{db}/miR.Counts.csv",
        rpm="results/mirna/{sample}/mirge3_{db}/miR.RPM.csv",
        annotation="results/mirna/{sample}/mirge3_{db}/annotation.report.csv",
        gff="results/mirna/{sample}/mirge3_{db}/sample_miRge3.gff",
        # Each annotated sequence with its exact or isomiR miRNA and its count
        mapped=temp("results/mirna/{sample}/mirge3_{db}/mapped.csv"),
    params:
        db=lambda wc: MIRGE3_DB[wc.db],
        organism=MIRGE3_ORGANISM,
        outdir="results/mirna/{sample}/mirge3_{db}",
    log:
        "logs/mirge3/{db}/{sample}.log",
    threads: 8
    resources:
        # Peak memory in miRNA-PD: 6.2 GB for a 1.4 GB input
        mem_mb=lambda wc, input: 4000 + 6 * int(input.size_mb),
        runtime=720,
    conda:
        "../../envs/mirge3.yaml"
    shadow:
        "minimal"
    shell:
        # miRge3 names its output columns after the input file, so the input is linked as <sample>.fastq.gz. It
        # writes into a time-stamped miRge.<date> subdirectory; its results are moved up one level. -q 0: no 3'
        # quality trim (default 10), so that miRge3 counts the reads as given; its per-read table must contain only
        # sequences of the collapsed reads
        """
        ln -s $(realpath {input.reads}) {wildcards.sample}.fastq.gz
        miRge3.0 -s {wildcards.sample}.fastq.gz -lib $(dirname $(realpath {input.lib})) -on {params.organism} \
            -db {params.db} -q 0 -gff -cpu {threads} -o mirge3_out &> {log}
        mkdir -p {params.outdir}
        mv mirge3_out/miRge.*/* {params.outdir}/
        """


rule mirge3_table:
    input:
        script=workflow.source_path("../../scripts/mirge3_table.py"),
        counts="results/mirna/{sample}/mirge3_{db}/miR.Counts.csv",
        mapped="results/mirna/{sample}/mirge3_{db}/mapped.csv",
        lib=MIRGE3_LIB_DIR,
        mature=lambda wc: MATURE_FASTAS[wc.db],
        arms="resources/refs/{db}_arms.tsv",
    output:
        table="results/mirna/{sample}/mirge3_{db}.tsv",
        reads="results/mirna/{sample}/mirge3_{db}.reads.tsv.gz",
    params:
        merges=lambda wc: f"annotation.Libs/{MIRGE3_ORGANISM}_merges_{MIRGE3_LIB[wc.db]}.csv",
        mirge3_mature=lambda wc: f"fasta.Libs/{MIRGE3_ORGANISM}_mature_{MIRGE3_LIB[wc.db]}.fa",
    log:
        "logs/mirge3_table/{db}/{sample}.log",
    resources:
        mem_mb=8000,
        runtime=30,
    conda:
        "../../envs/py.yaml"
    shell:
        """
        python {input.script} --sample {wildcards.sample} --counts {input.counts} --mapped {input.mapped} \
            --merges {input.lib}/{params.merges} --mirge3-mature {input.lib}/{params.mirge3_mature} \
            --mature {input.mature} --arms {input.arms} \
            --out-table {output.table} --out-reads {output.reads} &> {log}
        """


rule mirge3_version:
    output:
        "results/mirna/mirge3_version.txt",
    log:
        "logs/mirge3_version.log",
    conda:
        "../../envs/mirge3.yaml"
    shell:
        """python -c "import importlib.metadata as m; print('miRge3.0', m.version('mirge3'))" > {output} 2> {log}"""
