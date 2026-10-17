# Steps 5-7: read classification (bowtie cascade on the collapsed reads, incl. contaminants), miRNA
# quantification (miRge3.0, bowtie cross-check) and rRNA fragment tables.

MIRNA = config["mirna"]
MIRNA_DBS = ["mirbase", "mirgenedb"]
SMALL_REFS = MIRNA_DBS + ["rrna", "trna", "ncrna"]
# Multi-copy references: class-level assignment only, so the number of hits per read is capped
MULTICOPY_REFS = {"ncrna", "genome"}

MIRNA_SUMMARY = {
    "summary": "results/qc/mirna_summary.tsv",
    "composition": "results/qc/mirna_composition.tsv",
    "agreement": "results/qc/mirna_method_agreement.tsv",
    "tied": "results/qc/mirna_tied_other_class.tsv",
    "contaminants": "results/qc/contaminants.tsv",
    "top_annotation": "results/qc/top_sequences_annotated.tsv",
    "matrices": expand("results/mirna/matrices/{method}_{db}.tsv", method=["mirge3", "bowtie", "bowtie_mirge3names"], db=MIRNA_DBS),
}


def bowtie_hits_opts(ref):
    if ref in MULTICOPY_REFS:
        return f"-k {MIRNA['max_hits']} --best --strata"
    return "-a --best --strata"


rule bowtie_small:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        index="resources/bowtie/{ref}/.done",
    output:
        # read rank (s<rank>_x<count>), reference, 0-based offset, read length, mismatches; best stratum only, sense strand only
        "results/align/{sample}/{ref}.hits.tsv.gz",
    wildcard_constraints:
        ref="|".join(SMALL_REFS),
    params:
        prefix="resources/bowtie/{ref}/{ref}",
        mm=MIRNA["mismatches"],
        hits=lambda wc: bowtie_hits_opts(wc.ref),
    log:
        "logs/bowtie_small/{sample}_{ref}.log",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=240,
    conda:
        "../envs/align.yaml"
    shell:
        """
        bowtie -f -v {params.mm} {params.hits} --norc -p {threads} -x {params.prefix} {input.reads} 2> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print substr($1, 2, index($1, "_") - 2), $3, $4, length($5), mm}}' \
            | pigz -p 2 > {output}
        """


rule bowtie_contaminant:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        index="resources/bowtie/contaminant/.done",
    output:
        # As bowtie_small, plus the strand: genes of the contaminant genomes lie on both strands
        "results/align/{sample}/contaminant.hits.tsv.gz",
    params:
        prefix="resources/bowtie/contaminant/contaminant",
        mm=MIRNA["mismatches"],
        hits=f"-k {MIRNA['max_hits']} --best --strata",
    log:
        "logs/bowtie_contaminant/{sample}.log",
    threads: 8
    resources:
        mem_mb=8000,
        runtime=240,
    conda:
        "../envs/align.yaml"
    shell:
        """
        bowtie -f -v {params.mm} {params.hits} -p {threads} -x {params.prefix} {input.reads} 2> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print substr($1, 2, index($1, "_") - 2), $3, $4, length($5), mm, $2}}' \
            | pigz -p 2 > {output}
        """


rule genome_unassigned:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        hits=expand("results/align/{{sample}}/{ref}.hits.tsv.gz", ref=SMALL_REFS),
    output:
        temp("results/align/{sample}/unassigned_small.fa.gz"),
    log:
        "logs/genome_unassigned/{sample}.log",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    conda:
        "../envs/align.yaml"
    shell:
        # Reads without a hit on any small RNA reference go on to the genome
        """
        (
        pigz -dc {input.hits} | cut -f 1 \
            | awk 'NR == FNR {{hit[$1] = 1; next}}
                   /^>/ {{keep = !(substr($1, 3, index($1, "_") - 3) in hit)}}
                   keep' - <(pigz -dc {input.reads}) \
            | pigz -p {threads} > {output}
        ) 2> {log}
        """


rule bowtie_genome:
    input:
        reads="results/align/{sample}/unassigned_small.fa.gz",
        index="resources/bowtie/genome/.done",
        annotation="resources/refs/genome_annotation_merged.bed",
    output:
        # read rank, mismatches, comma-separated classes of the same-strand annotation it overlaps ('.' = none)
        "results/align/{sample}/genome.classes.tsv.gz",
    params:
        prefix="resources/bowtie/genome/genome",
        mm=MIRNA["mismatches"],
        hits=bowtie_hits_opts("genome"),
    log:
        "logs/bowtie_genome/{sample}.log",
    threads: 16
    resources:
        mem_mb=16000,
        runtime=480,
    conda:
        "../envs/align.yaml"
    shell:
        # Both strands; annotation overlaps are same-strand only (antisense hits fall to intron/intergenic)
        """
        bowtie -f -v {params.mm} {params.hits} -p {threads} -x {params.prefix} {input.reads} 2> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print $3, $4, $4 + length($5), substr($1, 2, index($1, "_") - 2), mm, $2}}' \
            | bedtools intersect -a stdin -b {input.annotation} -s -wao 2>> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{print $4, $5, $10}}' \
            | LC_ALL=C sort -u -S 8G --parallel 4 \
            | pigz -p 4 > {output}
        """


rule unmapped_ends:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        hits=expand("results/align/{{sample}}/{ref}.hits.tsv.gz", ref=SMALL_REFS),
        contaminant="results/align/{sample}/contaminant.hits.tsv.gz",
        genome="results/align/{sample}/genome.classes.tsv.gz",
        index=expand("resources/bowtie/{ref}/.done", ref=SMALL_REFS + ["genome", "contaminant"]),
    output:
        # end (<rank>_<5|3>_<end length>), index, reference, mismatches: best hits (best stratum, up to max_hits for
        # the contaminants, one for the human references)
        "results/align/{sample}/unmapped_ends.hits.tsv.gz",
    params:
        small=SMALL_REFS,
        mm=MIRNA["mismatches"],
        max_hits=MIRNA["max_hits"],
        end_len=MIRNA["end_len"],
        exact_len=MIRNA["end_len_exact"],
    log:
        "logs/unmapped_ends/{sample}.log",
    threads: 8
    resources:
        mem_mb=24000,
        runtime=240,
    conda:
        "../envs/align.yaml"
    shell:
        # Item 4c: the ends of reads longer than exact_len nt without a hit on any reference, in a temporary FASTA
        # next to the output (node-local /tmp is too small for the run 24 libraries). The exact_len-nt ends are
        # aligned with <= mm mismatches too; quantify_mirna.py keeps their exact hits only
        """
        ends=$(mktemp -p $(dirname {output}) unmapped_ends.XXXXXX.fa)
        trap 'rm -f $ends' EXIT
        (
        (pigz -dc {input.hits} {input.contaminant} | cut -f 1; pigz -dc {input.genome} | cut -f 1) \
            | awk -v s={params.exact_len} -v l={params.end_len} \
                'NR == FNR {{hit[$1] = 1; next}}
                 /^>/ {{rank = substr($1, 3, index($1, "_") - 3); next}}
                 !(rank in hit) && length($0) > s {{
                     n = length($0)
                     print ">" rank "_5_" s; print substr($0, 1, s); print ">" rank "_3_" s; print substr($0, n - s + 1)
                     if (n > l) {{print ">" rank "_5_" l; print substr($0, 1, l); print ">" rank "_3_" l; print substr($0, n - l + 1)}}
                 }}' - <(pigz -dc {input.reads}) > $ends
        echo "$(grep -c '^>' $ends) ends"
        for ref in {params.small} genome contaminant; do
            case $ref in
                genome) opts="-k 1 --best" ;;
                contaminant) opts="-k {params.max_hits} --best --strata" ;;
                *) opts="-k 1 --best --norc" ;;
            esac
            bowtie -f -v {params.mm} $opts -p {threads} -x resources/bowtie/$ref/$ref $ends \
                | awk -F '\\t' -v OFS='\\t' -v ref=$ref '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print $1, ref, $3, mm}}'
        done | pigz -p 2 > {output}
        ) &> {log}
        """


rule quantify_mirna:
    input:
        script="workflow/scripts/quantify_mirna.py",
        reads="results/collapse/{sample}/collapsed.fa.gz",
        hits=expand("results/align/{{sample}}/{ref}.hits.tsv.gz", ref=SMALL_REFS),
        contaminant="results/align/{sample}/contaminant.hits.tsv.gz",
        genome="results/align/{sample}/genome.classes.tsv.gz",
        ends="results/align/{sample}/unmapped_ends.hits.tsv.gz",
        arms=expand("resources/refs/{db}_arms.tsv", db=MIRNA_DBS),
    output:
        reads="results/mirna/{sample}/read_classes.tsv.gz",
        composition="results/mirna/{sample}/composition.tsv",
        counts=expand("results/mirna/{{sample}}/{db}_counts.tsv", db=MIRNA_DBS),
        # Reads per contaminant genome, and the best hits of the most abundant sequences (for the report)
        contaminants="results/mirna/{sample}/contaminants.tsv",
        top="results/mirna/{sample}/top_hits.tsv",
        # Contaminant reads by the class they would have had and the rule that assigned them (item 4b check)
        contaminant_rules="results/mirna/{sample}/contaminant_rules.tsv",
        # Unmapped reads tested by their ends (item 4c), by outcome and read length
        ends="results/mirna/{sample}/unmapped_ends.tsv",
    params:
        min_len=MIRNA["min_len"],
        max_len=MIRNA["max_len"],
        start_tol=MIRNA["arm_start_tol"],
        end_tol=MIRNA["arm_end_tol"],
        contaminant_min_len=MIRNA["contaminant_min_len"],
        end_len=MIRNA["end_len"],
        exact_len=MIRNA["end_len_exact"],
        top_n=config["collapse"]["top_n"],
    log:
        "logs/quantify_mirna/{sample}.log",
    resources:
        mem_mb=64000,
        runtime=240,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
            --reads {input.reads} \
            --hits {input.hits} \
            --contaminant {input.contaminant} \
            --genome {input.genome} \
            --ends {input.ends} \
            --arms {input.arms} \
            --min-len {params.min_len} --max-len {params.max_len} \
            --start-tol {params.start_tol} --end-tol {params.end_tol} \
            --contaminant-min-len {params.contaminant_min_len} \
            --end-len {params.end_len} --end-len-exact {params.exact_len} \
            --top-n {params.top_n} \
            --out-reads {output.reads} \
            --out-composition {output.composition} \
            --out-counts {output.counts} \
            --out-contaminants {output.contaminants} \
            --out-contaminant-rules {output.contaminant_rules} \
            --out-ends {output.ends} \
            --out-top {output.top} &> {log}
        """


# miRge3.0: the primary miRNA quantification, on the 16-28 nt reads (already trimmed, so no adapter given)
MIRGE3_DB = {"mirbase": "miRBase", "mirgenedb": "miRGeneDB"}
# File names in the miRge3 library use a different spelling
MIRGE3_LIB = {"mirbase": "miRBase", "mirgenedb": "MirGeneDB"}


rule mirge3_reads:
    # miRge3 input: the reads in the miRNA length range, without the reads classed as a contaminant (human reads
    # only; miRge3 had counted mycoplasma fragments as miRNA isomiRs)
    input:
        script="workflow/scripts/mirge3_reads.py",
        reads="results/trim/{sample}/trimmed.fastq.gz",
        collapsed="results/collapse/{sample}/collapsed.fa.gz",
        classes="results/mirna/{sample}/read_classes.tsv.gz",
    output:
        # Named after the sample: miRge3 names its output columns after the file
        temp("results/mirge3/reads/{db}/{sample}.fastq.gz"),
    wildcard_constraints:
        db="|".join(MIRNA_DBS),
    params:
        min_len=MIRNA["min_len"],
        max_len=MIRNA["max_len"],
    log:
        "logs/mirge3_reads/{db}/{sample}.log",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=240,
    conda:
        "../envs/mirge3.yaml"
    shell:
        """
        seqkit seq -m {params.min_len} -M {params.max_len} -j {threads} {input.reads} 2> {log} \
            | python {input.script} --classes {input.classes} --reads {input.collapsed} --db {wildcards.db} 2>> {log} \
            | pigz -p {threads} > {output}
        """


rule mirge3:
    input:
        reads="results/mirge3/reads/{db}/{sample}.fastq.gz",
        lib="resources/mirge3/human",
    output:
        counts="results/mirge3/{db}/{sample}/miR.Counts.csv",
        rpm="results/mirge3/{db}/{sample}/miR.RPM.csv",
        annotation="results/mirge3/{db}/{sample}/annotation.report.csv",
        gff="results/mirge3/{db}/{sample}/sample_miRge3.gff",
    wildcard_constraints:
        db="|".join(MIRNA_DBS),
    params:
        db=lambda wc: MIRGE3_DB[wc.db],
        outdir="results/mirge3/{db}/{sample}",
    log:
        "logs/mirge3/{db}/{sample}.log",
    threads: 8
    resources:
        mem_mb=lambda wc: 96000 if SAMPLES.loc[wc.sample, "run"] == "24" else 32000,
        runtime=720,
    conda:
        "../envs/mirge3.yaml"
    shadow:
        "minimal"
    shell:
        # miRge3 writes into a time-stamped miRge.<date> subdirectory; its results are moved up one level
        """
        miRge3.0 -s {input.reads} -lib $(dirname $(realpath {input.lib})) -on human -db {params.db} \
            -gff -cpu {threads} -o mirge3_out &> {log}
        mkdir -p {params.outdir}
        mv mirge3_out/miRge.*/* {params.outdir}/
        """


rule mirna_summary:
    input:
        samples=config["samples"],
        composition=expand("results/mirna/{sample}/composition.tsv", sample=SAMPLES.index),
        bowtie=expand("results/mirna/{sample}/{db}_counts.tsv", sample=SAMPLES.index, db=MIRNA_DBS),
        mirge3=expand("results/mirge3/{db}/{sample}/miR.Counts.csv", sample=SAMPLES.index, db=MIRNA_DBS),
        mirge3_report=expand("results/mirge3/{db}/{sample}/annotation.report.csv", sample=SAMPLES.index, db=MIRNA_DBS),
        lib="resources/mirge3/human",
        mature=["resources/download/mirbase/mature.fa", "resources/download/mirgenedb/mature.fa", "resources/download/mirgenedb/star.fa"],
        contaminants=expand("results/mirna/{sample}/contaminants.tsv", sample=SAMPLES.index),
        top_hits=expand("results/mirna/{sample}/top_hits.tsv", sample=SAMPLES.index),
        top="results/qc/top_sequences.tsv",
        rrna_features="resources/refs/rrna_features.bed",
        contaminant_features="resources/refs/contaminant_features.bed",
        ref_fastas=expand("resources/refs/{ref}.fa", ref=["rrna", "trna", "ncrna"]),
    output:
        **MIRNA_SUMMARY,
    params:
        dbs=MIRNA_DBS,
        merges=[f"annotation.Libs/human_merges_{MIRGE3_LIB[db]}.csv" for db in MIRNA_DBS],
        mirge3_mature=[f"fasta.Libs/human_mature_{MIRGE3_LIB[db]}.fa" for db in MIRNA_DBS],
        detected_min=MIRNA["detected_min_reads"],
    log:
        "logs/mirna_summary.log",
    resources:
        mem_mb=8000,
        runtime=30,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python workflow/scripts/mirna_summary.py \
            --samples {input.samples} \
            --dbs {params.dbs} \
            --composition {input.composition} \
            --bowtie-counts {input.bowtie} \
            --mirge3-counts {input.mirge3} \
            --mirge3-report {input.mirge3_report} \
            --merges $(printf '{input.lib}/%s ' {params.merges}) \
            --mirge3-mature $(printf '{input.lib}/%s ' {params.mirge3_mature}) \
            --mature {input.mature} \
            --detected-min {params.detected_min} \
            --contaminants {input.contaminants} \
            --top-hits {input.top_hits} \
            --top {input.top} \
            --rrna-features {input.rrna_features} \
            --contaminant-features {input.contaminant_features} \
            --ref-fastas {input.ref_fastas} \
            --out-dir $(dirname {output.matrices[0]}) \
            --out-summary {output.summary} \
            --out-composition {output.composition} \
            --out-agreement {output.agreement} \
            --out-tied {output.tied} \
            --out-contaminants {output.contaminants} \
            --out-top-annotation {output.top_annotation} &> {log}
        """


rule rrna_fragments:
    input:
        script="workflow/scripts/rrna_fragments.py",
        reads="results/collapse/{sample}/collapsed.fa.gz",
        classes="results/mirna/{sample}/read_classes.tsv.gz",
        hits="results/align/{sample}/rrna.hits.tsv.gz",
        contaminant="results/align/{sample}/contaminant.hits.tsv.gz",
        features="resources/refs/rrna_features.bed",
    output:
        "results/rrna/{sample}/fragments.tsv.gz",
    params:
        min_len=MIRNA["contaminant_min_len"],
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
            --contaminant {input.contaminant} --features {input.features} --contaminant-min-len {params.min_len} \
            --output {output} &> {log}
        """
