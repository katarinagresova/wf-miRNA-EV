# Read classification

Step 5 gives each read one class. The miRNA counts, the read shares and the rRNA fragments all use these classes.
Thus a read counts only once.

## Why EV data need it

EV libraries are different from whole-cell libraries. A pipeline that maps the reads only to miRNA references gives
incorrect miRNA counts and read shares for EV data.

| EV property | Design decision |
|---|---|
| Most reads are not miRNA. rRNA is often the largest class. | Classify every read, not only the miRNA reads. |
| rRNA fragments have the length of miRNAs, and some match a miRNA with 1 or 2 mismatches. | Align each read to all small RNA references, and use the best hit over all references. Report the reads that a miRNA-only mapping counts as miRNA ([off-target check](#mirna-off-target-check)). |
| The genome has many rDNA copies, and most have no rRNA annotation. | Count rRNA on a single-copy rRNA reference (45S rDNA, 5S, mitochondrial rRNA). |
| The RNA input is low, so contaminant RNA (for example, *Mycoplasma*) can be a large share. | Align all reads to contaminant genomes. Report the read shares also as a percent of the human reads. |
| Short contaminant reads match the human genome by chance. | Use [looser contaminant rules](#contaminants) for reads with weak human evidence. |
| Fragments longer than miRNAs carry information (for example, tRNA halves and Y RNA fragments). | Keep all reads of 16 nt or longer. Apply the miRNA length range to the miRNA counts only. |

Each of these rules is also correct for whole-cell libraries.

## Read classes

| Group | Classes | Reference |
|---|---|---|
| miRNA | `miRNA`, `miRNA_hairpin_other` | miRBase or MirGeneDB hairpins. `miRNA` reads count for a mature miRNA. `miRNA_hairpin_other` reads are on a hairpin, but do not count (for example, the loop, or a mature miRNA without its 5′ end). |
| rRNA | `rRNA` | The single-copy rRNA reference. |
| tRNA | `tRNA`, `mt_tRNA` | GtRNAdb mature tRNAs, Ensembl mitochondrial tRNAs. |
| Other small RNA | `Y_RNA`, `snoRNA`, `snRNA`, `7SK`, `7SL_SRP`, `vault_RNA`, `other_ncRNA` | Ensembl ncRNA. |
| Genome | `rRNA_genomic`, `mt_other`, `tRNA_genomic`, `mRNA`, `lncRNA`, `other_gene`, `intron`, `intergenic` | The genome, with the Ensembl gene annotation. Only reads without a small RNA hit. |
| Contaminant | `mycoplasma`, `bacteria`, `phiX` | The contaminant genomes of `refs.contaminants`. |
| | `unmapped` | No hit. |

Each database of `mirna.databases` has its own classification, because the miRNA hairpins are different. The other
classes change only for reads with a miRNA hit.

## Rules

1. **Small RNA references.** bowtie 1 aligns each read to each small RNA reference (sense strand, up to
   `classify.mismatches` mismatches). The hit with the fewest mismatches over all references gives the class. If
   classes have equally good hits, this order decides: `miRNA`, `rRNA`, `tRNA`, `mt_tRNA`, `Y_RNA`, `snoRNA`,
   `snRNA`, `7SK`, `7SL_SRP`, `vault_RNA`, `other_ncRNA`, `miRNA_hairpin_other`. The
   [off-target check](#mirna-off-target-check) reports the miRNA reads with such a tie as `tied`.
2. **miRNA.** A hairpin hit is `miRNA` only if it follows the [counting rules](mirna_counters.md#counting-rules):
   16–28 nt, the read ends near the ends of the mature miRNA, and 0 mismatches for reads of 18 nt or shorter. Other
   hairpin hits are `miRNA_hairpin_other`. A read without such a hit can have a hit with
   [non-templated 3′ nucleotides](mirna_counters.md#non-templated-3-nucleotides): it matches a mature miRNA after the
   removal of 1–3 nt at the 3′ end. In rule 1, this hit is a `miRNA` hit with 1 mismatch.
3. **Genome.** Reads without an rRNA, tRNA or ncRNA hit align to the genome (both strands). A read without a small
   RNA hit gets its class from the annotation on the same strand. The rDNA copies of `refs.rdna_loci` and the rRNA
   genes give `rRNA_genomic`. Other hits in the mitochondrial genome give `mt_other`. An exon hit wins over an intron
   hit.
4. **Contaminants.** See [Contaminants](#contaminants).

## Contaminants

All reads align to the contaminant genomes (both strands). A random short sequence matches a large genome by chance:
a random 16-nt sequence matches the 10 Mb of the default contaminant genomes with 1 mismatch or less with a
probability of approximately 21%. A random 18-nt sequence has a probability of approximately 1.6%. Thus the rules
depend on the read length and on the human evidence.

- A read of `classify.contaminant_min_len` (18) nt or longer is a contaminant read if its contaminant hit has fewer
  mismatches than its best human hit. For a read with a small RNA hit, this is the best small RNA hit (rule 1); its
  genome hits do not count. For other reads, it is the best genome hit. If the hits are equally good, the read stays
  human.
- **Weak human evidence**: no small RNA hit, and genome hits only in `mRNA`, `lncRNA`, `other_gene`, `intron` or
  `intergenic`, or no human hit. A short read lands in these classes by chance, because they cover most of the
  genome. For these reads, an equally good contaminant hit is sufficient (18 nt or longer). A shorter read needs an
  exact contaminant match.
- **Reads without a full-length hit.** Some contaminant reads have more errors than the alignment allows, or two
  fragments joined in the ligation. The workflow aligns the 25-nt ends (`classify.end_len`) of these reads. The read
  is a contaminant read if an end matches a contaminant and no end matches a human reference. If no 25-nt end has a
  hit, the 18-nt ends (`classify.end_len_exact`) must match exactly.

The contaminant reads are not in the miRNA counts. `pct_of_human` in `results/classify/composition.tsv` divides by
the human reads (the kept reads without the contaminant reads).

## miRNA off-target check

A pipeline that maps the reads only to miRNA references counts reads of other RNAs as miRNA. The off-target check
shows these reads for each library. It does not change a count.

The check compares two mappings:

- **miRNA-only mapping.** bowtie 1 aligns all kept reads of 16–28 nt to the miRNA hairpins only, with up to
  `mirna.off_target_mismatches` (3) mismatches. The contaminant reads are included. The rule "0 mismatches for 18 nt
  or shorter" does not apply.
- **Classification.** The class of the read in step 5.

Each read with a miRNA-only hit on a mature miRNA gets one result:

| Result | Meaning |
|---|---|
| `miRNA` | Step 5 gives the class miRNA. No other class has an equally good hit. |
| `tied` | Step 5 gives the class miRNA, but another class has an equally good hit. |
| `miRNA_hairpin_other` | The read is 18 nt or shorter and has a mismatch. Step 5 does not count it. |
| another class (for example, `rRNA`, `tRNA`, `intron`, `mycoplasma`) | Another reference has a better hit. |
| `unmapped` | The read has 2 or 3 mismatches to the miRNA, and step 5 has no hit. |

For example, a read with 1 mismatch to a miRNA and 0 mismatches to the rRNA is probably an rRNA fragment. But it can
also be a miRNA with a SNP or an edited base. Thus the check reports and does not remove.

The results are in `library_summary.tsv` (percent of each result), in `results/mirna/off_target.tsv` (for each
miRNA, over all libraries) and in the QC report (by library and by read length). See [Outputs](outputs.md).
