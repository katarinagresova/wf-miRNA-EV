# Step 5: read classification (bowtie cascade on the collapsed reads, incl. contaminants) and the miRNA counts.

CLASSIFY = config["classify"]
MIRNA = config["mirna"]
MIRNA_DBS = MIRNA["databases"]
SMALL_REFS = MIRNA_DBS + ["rrna", "trna", "ncrna"]
# bowtie reports all hits of the best stratum (-a). With -k, the reported subset of the hits depends on the read
# name, which contains the rank: the class of a sequence would then depend on the other reads of the library

wildcard_constraints:
    db="|".join(MIRNA_DBS),

CLASSIFY_SUMMARY = {
    "composition": "results/classify/composition.tsv",
    "contaminants": "results/classify/contaminants.tsv",
    "top": "results/classify/top_sequences_classes.tsv",
}


rule bowtie_small:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        index="resources/bowtie/{ref}/.done",
    output:
        # read rank (s<rank>_x<count>), reference, 0-based offset, read length, mismatches; best stratum only, sense strand only
        "results/classify/{sample}/hits/{ref}.hits.tsv.gz",
    wildcard_constraints:
        ref="|".join(SMALL_REFS),
    params:
        prefix="resources/bowtie/{ref}/{ref}",
        mm=CLASSIFY["mismatches"],
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
        bowtie -f -v {params.mm} -a --best --strata --norc -p {threads} -x {params.prefix} {input.reads} 2> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print substr($1, 2, index($1, "_") - 2), $3, $4, length($5), mm}}' \
            | pigz -p 2 > {output}
        """


rule bowtie_nta:
    # Reads with non-templated 3' nucleotides (quantify_mirna.py): the reads without n nt at the 3' end (n = 1 to
    # nta_max), exact hits on the miRNA hairpins. Only the reads that are max_len nt or shorter, and min_len nt or
    # longer after the removal
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        index="resources/bowtie/{db}/.done",
    output:
        # As bowtie_small; the read length is the length after the removal
        "results/classify/{sample}/hits/{db}_nta.hits.tsv.gz",
    params:
        prefix="resources/bowtie/{db}/{db}",
        min_len=MIRNA["min_len"],
        max_len=MIRNA["max_len"],
        nta_max=MIRNA["nta_max"],
    log:
        "logs/bowtie_nta/{sample}_{db}.log",
    threads: 4
    resources:
        mem_mb=4000,
        runtime=120,
    conda:
        "../envs/align.yaml"
    shell:
        """
        for n in $(seq 1 {params.nta_max}); do
            seqkit seq -m $(({params.min_len} + n)) -M {params.max_len} {input.reads} \
                | bowtie -f -v 0 -3 $n -a --norc -p {threads} -x {params.prefix} -
        done 2> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{print substr($1, 2, index($1, "_") - 2), $3, $4, length($5), 0}}' \
            | pigz -p 2 > {output}
        """


rule bowtie_contaminant:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        index="resources/bowtie/contaminant/.done",
    output:
        # As bowtie_small, plus the strand: genes of the contaminant genomes lie on both strands
        "results/classify/{sample}/hits/contaminant.hits.tsv.gz",
    params:
        prefix="resources/bowtie/contaminant/contaminant",
        mm=CLASSIFY["mismatches"],
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
        bowtie -f -v {params.mm} -a --best --strata -p {threads} -x {params.prefix} {input.reads} 2> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print substr($1, 2, index($1, "_") - 2), $3, $4, length($5), mm, $2}}' \
            | pigz -p 2 > {output}
        """


rule genome_unassigned:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        hits=expand("results/classify/{{sample}}/hits/{ref}.hits.tsv.gz", ref=["rrna", "trna", "ncrna"]),
    output:
        temp("results/classify/{sample}/hits/unassigned_small.fa.gz"),
    log:
        "logs/genome_unassigned/{sample}.log",
    threads: 4
    resources:
        mem_mb=16000,
        runtime=120,
    conda:
        "../envs/align.yaml"
    shell:
        # Reads without a hit on rRNA, tRNA or ncRNA go on to the genome. Reads with a hit on only one miRNA database
        # need the genome class in the classification of the other database
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
        reads="results/classify/{sample}/hits/unassigned_small.fa.gz",
        index=GENOME_INDEX_FILE,
        annotation="resources/refs/genome_annotation_merged.bed",
    output:
        # read rank, mismatches, class of a same-strand annotation that a hit overlaps ('.' = none); one line for
        # each different class
        "results/classify/{sample}/hits/genome.classes.tsv.gz",
    params:
        prefix=GENOME_INDEX,
        mm=CLASSIFY["mismatches"],
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
        bowtie -f -v {params.mm} -a --best --strata -p {threads} -x {params.prefix} {input.reads} 2> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print $3, $4, $4 + length($5), substr($1, 2, index($1, "_") - 2), mm, $2}}' \
            | bedtools intersect -a stdin -b {input.annotation} -s -wao 2>> {log} \
            | awk -F '\\t' -v OFS='\\t' '{{print $4, $5, $10}}' \
            | LC_ALL=C sort -u -S 8G --parallel 4 \
            | pigz -p 4 > {output}
        """


rule unmapped_ends:
    input:
        reads="results/collapse/{sample}/collapsed.fa.gz",
        hits=expand("results/classify/{{sample}}/hits/{ref}.hits.tsv.gz", ref=SMALL_REFS),
        contaminant="results/classify/{sample}/hits/contaminant.hits.tsv.gz",
        genome="results/classify/{sample}/hits/genome.classes.tsv.gz",
        index=expand("resources/bowtie/{ref}/.done", ref=SMALL_REFS + ["contaminant"]),
        genome_index=GENOME_INDEX_FILE,
    output:
        # end (<rank>_<5|3>_<end length>), index, reference, mismatches, 0-based offset: best hits (best stratum; one
        # hit for the genome, of which only the presence and the mismatches count)
        "results/classify/{sample}/hits/unmapped_ends.hits.tsv.gz",
    params:
        small=SMALL_REFS,
        genome_index=GENOME_INDEX,
        mm=CLASSIFY["mismatches"],
        end_len=CLASSIFY["end_len"],
        exact_len=CLASSIFY["end_len_exact"],
    log:
        "logs/unmapped_ends/{sample}.log",
    threads: 8
    resources:
        mem_mb=24000,
        runtime=240,
    conda:
        "../envs/align.yaml"
    shell:
        # The ends of reads longer than exact_len nt without a hit on any reference, in a temporary FASTA next to
        # the output (node-local /tmp can be too small for large libraries). The exact_len-nt ends are aligned with
        # <= mm mismatches too; quantify_mirna.py keeps their exact hits only
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
                genome) opts="-k 1 --best"; prefix={params.genome_index} ;;
                contaminant) opts="-a --best --strata"; prefix=resources/bowtie/$ref/$ref ;;
                *) opts="-a --best --strata --norc"; prefix=resources/bowtie/$ref/$ref ;;
            esac
            bowtie -f -v {params.mm} $opts -p {threads} -x $prefix $ends \
                | awk -F '\\t' -v OFS='\\t' -v ref=$ref '{{mm = ($8 == "") ? 0 : gsub(/,/, ",", $8) + 1; print $1, ref, $3, mm, $4}}'
        done | pigz -p 2 > {output}
        ) &> {log}
        """


rule quantify_mirna:
    input:
        script=workflow.source_path("../scripts/quantify_mirna.py"),
        reads="results/collapse/{sample}/collapsed.fa.gz",
        hits=expand("results/classify/{{sample}}/hits/{ref}.hits.tsv.gz", ref=SMALL_REFS),
        contaminant="results/classify/{sample}/hits/contaminant.hits.tsv.gz",
        genome="results/classify/{sample}/hits/genome.classes.tsv.gz",
        ends="results/classify/{sample}/hits/unmapped_ends.hits.tsv.gz",
        arms=expand("resources/refs/{db}_arms.tsv", db=MIRNA_DBS),
        nta=expand("results/classify/{{sample}}/hits/{db}_nta.hits.tsv.gz", db=MIRNA_DBS),
    output:
        reads="results/classify/{sample}/read_classes.tsv.gz",
        composition="results/classify/{sample}/composition.tsv",
        # miRNA counts and their reads; the counter file counters/bowtie.smk makes the interface tables from them
        counts=expand("results/classify/{{sample}}/mirna_{db}.tsv", db=MIRNA_DBS),
        mirna_reads=expand("results/classify/{{sample}}/mirna_{db}.reads.tsv.gz", db=MIRNA_DBS),
        # Reads per contaminant genome, and the best hits of the most abundant sequences (for the report)
        contaminants="results/classify/{sample}/contaminants.tsv",
        top="results/classify/{sample}/top_hits.tsv",
        # Contaminant reads by the class they would have had and the rule that assigned them
        contaminant_rules="results/classify/{sample}/contaminant_rules.tsv",
        # Unmapped reads tested by their ends, by outcome and read length
        ends="results/classify/{sample}/unmapped_ends.tsv",
    params:
        min_len=MIRNA["min_len"],
        max_len=MIRNA["max_len"],
        exact_max_len=MIRNA["exact_max_len"],
        start_tol=CLASSIFY["arm_start_tol"],
        end_tol=CLASSIFY["arm_end_tol"],
        contaminant_min_len=CLASSIFY["contaminant_min_len"],
        end_len=CLASSIFY["end_len"],
        exact_len=CLASSIFY["end_len_exact"],
        top_n=config["collapse"]["top_n"],
    log:
        "logs/quantify_mirna/{sample}.log",
    resources:
        # Peak memory: approximately 1.1 GB per million unique sequences (unique_seqs in results/collapse/summary.tsv)
        mem_mb=16000,
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
            --nta-hits {input.nta} \
            --min-len {params.min_len} --max-len {params.max_len} --exact-max-len {params.exact_max_len} \
            --start-tol {params.start_tol} --end-tol {params.end_tol} \
            --contaminant-min-len {params.contaminant_min_len} \
            --end-len {params.end_len} --end-len-exact {params.exact_len} \
            --top-n {params.top_n} \
            --out-reads {output.reads} \
            --out-composition {output.composition} \
            --out-counts {output.counts} \
            --out-mirna-reads {output.mirna_reads} \
            --out-contaminants {output.contaminants} \
            --out-contaminant-rules {output.contaminant_rules} \
            --out-ends {output.ends} \
            --out-top {output.top} &> {log}
        """


rule classify_summary:
    input:
        script=workflow.source_path("../scripts/classify_summary.py"),
        composition=expand("results/classify/{sample}/composition.tsv", sample=SAMPLES.index),
        contaminants=expand("results/classify/{sample}/contaminants.tsv", sample=SAMPLES.index),
        top_hits=expand("results/classify/{sample}/top_hits.tsv", sample=SAMPLES.index),
        top="results/collapse/top_sequences.tsv",
        rrna_features="resources/refs/rrna_features.bed",
        contaminant_features="resources/refs/contaminant_features.bed",
        ref_fastas=expand("resources/refs/{ref}.fa", ref=["rrna", "trna", "ncrna"]),
        samples=config["samples"],
    output:
        **CLASSIFY_SUMMARY,
    params:
        dbs=MIRNA_DBS,
    log:
        "logs/classify_summary.log",
    resources:
        mem_mb=8000,
        runtime=30,
    conda:
        "../envs/py.yaml"
    shell:
        """
        python {input.script} \
            --dbs {params.dbs} \
            --composition {input.composition} \
            --contaminants {input.contaminants} \
            --top-hits {input.top_hits} \
            --top {input.top} \
            --rrna-features {input.rrna_features} \
            --contaminant-features {input.contaminant_features} \
            --ref-fastas {input.ref_fastas} \
            --samples {input.samples} \
            --out-composition {output.composition} \
            --out-contaminants {output.contaminants} \
            --out-top {output.top} &> {log}
        """
