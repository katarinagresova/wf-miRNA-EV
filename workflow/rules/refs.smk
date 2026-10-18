# Step 4: download and build the references (miRNA, rRNA, tRNA, other small ncRNA, contaminant genomes, genome) and their bowtie indices.

import re
import shlex

REFS = config["refs"]
SPECIES = config["species"]

# Bowtie1 references; each FASTA has names `<id>|<class>|<name>`, except the genome
BOWTIE_REFS = ["mirbase", "mirgenedb", "rrna", "trna", "ncrna", "contaminant", "genome"]
CONTAMINANTS = {acc: (cls, species) for acc, cls, species in REFS["contaminants"]}
# Bowtie index of the genome: the one in the config, or built by the workflow
GENOME_INDEX = config["genome_index"] or "resources/bowtie/genome/genome"
GENOME_INDEX_FILE = f"{GENOME_INDEX}.1.ebwt" if config["genome_index"] else "resources/bowtie/genome/.done"


rule download_mirbase:
    localrule: True
    output:
        "resources/download/mirbase/{kind}.fa",
    wildcard_constraints:
        kind="mature|hairpin",
    params:
        url=lambda wc: REFS["mirbase"][wc.kind]["url"],
        md5=lambda wc: REFS["mirbase"][wc.kind]["md5"],
    log:
        "logs/download_mirbase/{kind}.log",
    conda:
        "../envs/download.yaml"
    shell:
        """
        curl -sSfL --retry 3 -o {output} "{params.url}" 2> {log}
        echo "{params.md5}  {output}" | md5sum -c - >> {log} 2>&1
        """


rule download_mirgenedb:
    localrule: True
    output:
        "resources/download/mirgenedb/{kind}.fa",
    wildcard_constraints:
        kind="mature|star|pri",
    params:
        url=lambda wc: REFS["mirgenedb"][wc.kind],
        species=SPECIES["mirgenedb"],
    log:
        "logs/download_mirgenedb/{kind}.log",
    conda:
        "../envs/download.yaml"
    shell:
        """
        curl -sSfL --retry 3 -o {output} "{params.url}" 2> {log}
        grep -q '^>{params.species}-' {output}
        md5sum {output} >> {log}
        """


rule download_gtrnadb:
    localrule: True
    output:
        mature="resources/download/gtrnadb/mature-tRNAs.fa",
        tar="resources/download/gtrnadb/tRNAs.tar.gz",
    params:
        mature=REFS["gtrnadb"]["mature"],
        tar=REFS["gtrnadb"]["tar"],
    log:
        "logs/download_gtrnadb.log",
    conda:
        "../envs/download.yaml"
    shell:
        """
        curl -sSfL --retry 3 -o {output.mature} "{params.mature}" 2> {log}
        curl -sSfL --retry 3 -o {output.tar} "{params.tar}" 2>> {log}
        md5sum {output} >> {log}
        """


rule download_ensembl:
    localrule: True
    output:
        "resources/download/ensembl/{kind}.gz",
    wildcard_constraints:
        kind="genome|gtf|ncrna",
    params:
        url=lambda wc: REFS["ensembl"][wc.kind],
    log:
        "logs/download_ensembl/{kind}.log",
    conda:
        "../envs/download.yaml"
    shell:
        # Ensembl publishes CHECKSUMS (BSD `sum`) next to every file
        """
        curl -sSfL --retry 3 -o {output} "{params.url}" 2> {log}
        expected=$(curl -sSfL "$(dirname {params.url})/CHECKSUMS" | awk -v f=$(basename {params.url}) '$3 == f {{print $1, $2}}')
        observed=$(sum {output} | awk '{{print $1, $2}}')
        echo "CHECKSUMS: $expected, observed: $observed" >> {log}
        [ "$expected" = "$observed" ]
        """


rule download_rdna:
    localrule: True
    output:
        "resources/download/rdna.fa",
    params:
        id=REFS["rrna"]["rdna"],
    log:
        "logs/download_rdna.log",
    conda:
        "../envs/download.yaml"
    shell:
        """
        curl -sSf "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={params.id}&rettype=fasta&retmode=text" \
            -o {output} 2> {log}
        grep -q '^>{params.id}' {output}
        """


rule download_contaminant:
    localrule: True
    output:
        fasta="resources/download/contaminants/{acc}.fa",
        # NCBI feature table (genes, tRNA, rRNA, CDS), to say which genes the contaminant reads come from
        ft="resources/download/contaminants/{acc}.ft",
    wildcard_constraints:
        acc="|".join(re.escape(a) for a in CONTAMINANTS),
    log:
        "logs/download_contaminant/{acc}.log",
    conda:
        "../envs/download.yaml"
    shell:
        """
        url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={wildcards.acc}&retmode=text"
        curl -sSf --retry 3 "$url&rettype=fasta" -o {output.fasta} 2> {log}
        curl -sSf --retry 3 "$url&rettype=ft" -o {output.ft} 2>> {log}
        grep -q '^>{wildcards.acc}' {output.fasta}
        md5sum {output} >> {log}
        """


rule contaminant_reference:
    input:
        script=workflow.source_path("../scripts/prepare_refs.py"),
        fasta=[ncbi_source(acc, "fa", f"resources/download/contaminants/{acc}.fa") for acc in CONTAMINANTS],
        ft=[ncbi_source(acc, "ft", f"resources/download/contaminants/{acc}.ft") for acc in CONTAMINANTS],
    output:
        fasta="resources/refs/contaminant.fa",
        # contig, start, end, <type>:<product>, 0, strand
        bed="resources/refs/contaminant_features.bed",
    params:
        species=" ".join(f"--species {acc} {cls} {shlex.quote(sp)}" for acc, (cls, sp) in CONTAMINANTS.items()),
    log:
        "logs/refs/contaminant.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} contaminant \
            --genomes {input.fasta} --feature-tables {input.ft} {params.species} \
            --fasta {output.fasta} --bed {output.bed} &> {log}
        """


rule mirbase_reference:
    input:
        script=workflow.source_path("../scripts/prepare_refs.py"),
        mature=ref_source(REFS["mirbase"]["mature"], "resources/download/mirbase/mature.fa"),
        hairpin=ref_source(REFS["mirbase"]["hairpin"], "resources/download/mirbase/hairpin.fa"),
    output:
        fasta="resources/refs/mirbase.fa",
        arms="resources/refs/mirbase_arms.tsv",
    params:
        species=SPECIES["mirbase"],
    log:
        "logs/refs/mirbase.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} mirbase --species {params.species} \
            --mature {input.mature} --hairpin {input.hairpin} \
            --fasta {output.fasta} --arms {output.arms} &> {log}
        """


rule mirgenedb_reference:
    input:
        script=workflow.source_path("../scripts/prepare_refs.py"),
        mature=ref_source(REFS["mirgenedb"]["mature"], "resources/download/mirgenedb/mature.fa"),
        star=ref_source(REFS["mirgenedb"]["star"], "resources/download/mirgenedb/star.fa"),
        pri=ref_source(REFS["mirgenedb"]["pri"], "resources/download/mirgenedb/pri.fa"),
    output:
        fasta="resources/refs/mirgenedb.fa",
        arms="resources/refs/mirgenedb_arms.tsv",
    log:
        "logs/refs/mirgenedb.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} mirgenedb \
            --mature {input.mature} --star {input.star} --pri {input.pri} \
            --fasta {output.fasta} --arms {output.arms} &> {log}
        """


rule rrna_reference:
    input:
        script=workflow.source_path("../scripts/prepare_refs.py"),
        rdna=ncbi_source(REFS["rrna"]["rdna"], "fa", "resources/download/rdna.fa"),
        ncrna=ref_source(REFS["ensembl"]["ncrna"], "resources/download/ensembl/ncrna.gz"),
        genome="resources/refs/genome.fa",
    output:
        "resources/refs/rrna.fa",
    params:
        rna5s=REFS["rrna"]["rna5s"],
        rna5s_flank=REFS["rrna"]["rna5s_flank"],
        mt_rrna=REFS["rrna"]["mt_rrna"],
    log:
        "logs/refs/rrna.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} rrna \
            --rdna {input.rdna} --ncrna {input.ncrna} \
            --rna5s {params.rna5s} --rna5s-flank {params.rna5s_flank} --genome {input.genome} --mt-rrna {params.mt_rrna} \
            --fasta {output} &> {log}
        """


# Separate from the reference, so that a change of the region boundaries does not trigger realignment
rule rrna_features:
    input:
        script=workflow.source_path("../scripts/prepare_refs.py"),
        fasta="resources/refs/rrna.fa",
    output:
        "resources/refs/rrna_features.bed",
    params:
        rna5s_flank=REFS["rrna"]["rna5s_flank"],
        rna5s_length=REFS["rrna"]["rna5s_length"],
        features=" ".join(
            f"--features {shlex.quote(name)} {start} {end}"
            for name, start, end in REFS["rrna"]["features"][REFS["rrna"]["rdna"]]
        ),
    log:
        "logs/refs/rrna_features.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} rrna-features --fasta {input.fasta} --rna5s-flank {params.rna5s_flank} \
            --rna5s-length {params.rna5s_length} {params.features} --bed {output} &> {log}
        """


rule trna_reference:
    input:
        script=workflow.source_path("../scripts/prepare_refs.py"),
        gtrnadb=ref_source(REFS["gtrnadb"]["mature"], "resources/download/gtrnadb/mature-tRNAs.fa"),
        tar=ref_source(REFS["gtrnadb"]["tar"], "resources/download/gtrnadb/tRNAs.tar.gz"),
        ncrna=ref_source(REFS["ensembl"]["ncrna"], "resources/download/ensembl/ncrna.gz"),
    output:
        fasta="resources/refs/trna.fa",
        bed="resources/refs/trna_loci.bed",
    log:
        "logs/refs/trna.log",
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} trna \
            --gtrnadb {input.gtrnadb} --gtrnadb-tar {input.tar} --ncrna {input.ncrna} \
            --fasta {output.fasta} --bed {output.bed} &> {log}
        """


rule ncrna_reference:
    input:
        script=workflow.source_path("../scripts/prepare_refs.py"),
        ncrna=ref_source(REFS["ensembl"]["ncrna"], "resources/download/ensembl/ncrna.gz"),
    output:
        "resources/refs/ncrna.fa",
    log:
        "logs/refs/ncrna.log",
    conda:
        "../envs/py.yaml"
    shell:
        "python {input.script} ncrna --ncrna {input.ncrna} --fasta {output} &> {log}"


rule genome_reference:
    input:
        ref_source(REFS["ensembl"]["genome"], "resources/download/ensembl/genome.gz"),
    output:
        "resources/refs/genome.fa",
    log:
        "logs/refs/genome.log",
    resources:
        runtime=60,
    conda:
        "../envs/align.yaml"
    shell:
        # Keep only the sequence name in the headers
        "(zcat {input} | cut -d ' ' -f 1 > {output}) 2> {log}"


rule gene_annotation:
    input:
        gtf=ref_source(REFS["ensembl"]["gtf"], "resources/download/ensembl/gtf.gz"),
        trna="resources/refs/trna_loci.bed",
    output:
        # chrom, start, end, <gene_name>|<class>, 0, strand; class = gene biotype of exons, or intron
        "resources/refs/genome_annotation.bed",
    log:
        "logs/refs/gene_annotation.log",
    conda:
        "../envs/align.yaml"
    shell:
        """
        (
        zcat {input.gtf} | awk -F '\\t' -v OFS='\\t' '
            $3 == "exon" || $3 == "gene" {{
                match($9, /gene_biotype "[^"]+"/); bt = substr($9, RSTART + 14, RLENGTH - 15)
                name = $9; if (match($9, /gene_name "[^"]+"/)) name = substr($9, RSTART + 11, RLENGTH - 12)
                else {{ match($9, /gene_id "[^"]+"/); name = substr($9, RSTART + 9, RLENGTH - 10) }}
                print $1, $4 - 1, $5, name "|" ($3 == "exon" ? bt : "intron"), 0, $7
            }}'
        cat {input.trna}
        ) | sort -k1,1 -k2,2n -S 4G > {output} 2> {log}
        """


rule genome_annotation_merged:
    input:
        "resources/refs/genome_annotation.bed",
    output:
        # Same-strand intervals merged within each class (name = class), so exons and gene bodies (intron) stay
        # separate; plus rDNA loci and the mt genome on both strands
        "resources/refs/genome_annotation_merged.bed",
    params:
        extra="\\n".join(
            f"{chrom}\\t{start}\\t{end}\\t{cls}\\t0\\t{strand}"
            for chrom, start, end, cls in [*[(*r, "rDNA_locus") for r in REFS["rdna_loci"]], ("MT", 0, 16569, "mt_genome")]
            for strand in "+-"
        ),
    log:
        "logs/refs/genome_annotation_merged.log",
    conda:
        "../envs/align.yaml"
    shell:
        # Merge per class by making the class part of the chromosome name (<class>@<chrom>) during the merge
        """
        (
        (awk -F '\\t' -v OFS='\\t' '{{split($4, a, "|"); $4 = a[2]; print}}' {input}; printf '{params.extra}\\n') \
            | awk -F '\\t' -v OFS='\\t' '{{print $4 "@" $1, $2, $3, $4, $5, $6}}' \
            | sort -k1,1 -k2,2n -S 4G \
            | bedtools merge -s -c 6 -o distinct \
            | awk -F '\\t' -v OFS='\\t' '{{split($1, k, "@"); print k[2], $2, $3, k[1], 0, $4}}' \
            | sort -k1,1 -k2,2n -S 4G
        ) > {output} 2> {log}
        """


rule bowtie_index:
    input:
        "resources/refs/{ref}.fa",
    output:
        touch("resources/bowtie/{ref}/.done"),
    wildcard_constraints:
        ref="|".join(BOWTIE_REFS),
    params:
        prefix="resources/bowtie/{ref}/{ref}",
    log:
        "logs/bowtie_index/{ref}.log",
    threads: lambda wc: 16 if wc.ref == "genome" else 2
    resources:
        mem_mb=lambda wc: 32000 if wc.ref == "genome" else 4000,
        runtime=lambda wc: 480 if wc.ref == "genome" else 60,
    conda:
        "../envs/align.yaml"
    shell:
        "bowtie-build --threads {threads} {input} {params.prefix} &> {log}"

