# miRNA counts

The miRNA counts come from the [read classification](classification.md) of step 5. Only reads with the class miRNA
count. Thus the miRNA reads in the count tables and the miRNA class in `composition.tsv` are the same reads. All
filters are in the config.

## Counting rules

A read counts for a mature miRNA if all of these conditions are true:

- Its class in step 5 is `miRNA`: no other reference has a better hit.
- It is 16–28 nt long (`mirna.min_len`, `mirna.max_len`).
- It has 0 or 1 mismatch to the hairpin (`classify.mismatches`). A read of 18 nt or shorter must have 0 mismatches
  (`mirna.exact_max_len`).
- Its 5′ end is 3 nt or less from the 5′ end of the mature miRNA, and its 3′ end is 5 nt or less from the 3′ end
  (`classify.arm_start_tol`, `classify.arm_end_tol`). The 3′ end varies more in isomiRs.

A read that counts for two or more entries is divided equally between them.

A read without the first 4 nt or more of the mature miRNA does not count, because it has no seed. Its class is
`miRNA_hairpin_other`.

### Non-templated 3′ nucleotides

Many miRNAs have 1–3 extra nucleotides at the 3′ end that are not in the genome (for example, adenylation:
hsa-miR-92a-3p + AA). With these nucleotides, the read has too many mismatches and no hit in step 5.

A read that is still `unmapped` after step 5 counts for a mature miRNA if:

- without 1–3 nt at its 3′ end (`mirna.nta_max`), it matches the hairpin with 0 mismatches,
- the shorter read is 16 nt or longer and agrees with the arm tolerances.

The match with the fewest removed nucleotides wins. The read gets the class `miRNA`. The column `nta_reads` of
`bowtie_{db}.tsv` gives these reads.

### Entries

Mature names with the same sequence become one entry, because no read can separate them. The entry name is the names,
joined by `/` (for example, `hsa-miR-199a-3p/hsa-miR-199b-3p`). The column `members` gives the names. miRBase 22.1 has
16 such entries for human, and MirGeneDB 3.0 has 69.

## Why these rules

**Short reads with a mismatch.** A short read matches a miRNA hairpin by chance more often than a long read. We
aligned EV reads to a decoy reference: the miRBase hairpins, reversed. Most decoy hits were reads of 16–18 nt with 1
mismatch. Thus a read of 18 nt or shorter needs 0 mismatches. This rule removes more decoy hits, and fewer real miRNA
reads, than a minimum length of 18 nt.

**Only reads with the class miRNA.** In EV data, many fragments of rRNA, tRNA and introns match a miRNA with 1 or 2
mismatches. A counter that maps only to the miRNA references counts them. The [off-target
check](classification.md#mirna-off-target-check) shows these reads for each library.

**Why not miRge3.0.** An earlier development version used miRge3.0 (0.1.4). Its last alignment step aligns the reads without an
annotation to the miRNA index with 2 mismatches, after the removal of 1 nt at the 5′ end and 2 nt at the 3′ end. Thus
a 16-nt read aligns with 13 nt and 2 mismatches. By our estimate, approximately 40% of random 16-nt sequences get a
hit. Its count filter applies to each library separately: it sets the count of a miRNA to 0 if the exact reads are
10% or less of the isomiR reads. With many fragment reads, a miRNA then has all its reads in some libraries and 0 in
the others. The counting rules above keep the isomiR type that miRge3.0 found and the strict rules rejected
(non-templated 3′ nucleotides), and merge identical mature sequences as miRge3.0 does.

## Counter interface

The summary step (`mirna_summary.py`) reads only the tables of the counter. Thus another counter can replace the
bowtie counts, and only the code of the counter changes.

The counter writes these tables for each library and database:

| File | Columns |
|---|---|
| `results/mirna/{sample}/{counter}_{db}.tsv` | `mirna` (the entry), `reads` (the count after all filters), `members` (the names of the entry in the miRBase or MirGeneDB release of the config, joined by `;`). Other columns are allowed. The summary step does not use them. |
| `results/mirna/{sample}/{counter}_{db}.reads.tsv.gz` | The counted reads: `sequence`, `mirna`, `reads`. The off-target check uses this table. |

The bowtie counts add the columns `unique_seqs`, `canonical_reads` and `nta_reads`. Their counter file is
`workflow/rules/counters/bowtie.smk`. It only makes the tables above from the tables of step 5.

## Add a counter

1. Write a rule file `workflow/rules/counters/<counter>.smk`. Its rules write the two tables above. Put the scripts in
   `workflow/scripts/` and the conda environment in `workflow/envs/`.
2. In `workflow/rules/mirna.smk`, include the new rule file, and use `<counter>_{db}.tsv` and
   `<counter>_{db}.reads.tsv.gz` as the inputs of `off_target` and `mirna_summary`. In `workflow/Snakefile`, give
   the counter name and version to the rule `wf_version`.
3. Run the test data set in `.test/`. It has isomiRs with known classes (see the README, Tests).

A new counter can change the entry names in `counts_{db}.tsv`. Thus it needs a new minor version.
