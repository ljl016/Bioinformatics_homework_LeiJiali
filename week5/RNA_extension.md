# RNA-specific Extension: Small-RNA Analysis Plan

## Scope

The supplied matrix is treated as bulk RNA-seq gene-level count data. It does not include library-preparation details, read lengths, adapter information, or validated miRNA identifiers. Therefore, it cannot support a completed miRNA analysis. The following is a proposed workflow for a separate small-RNA library, not a claim about the current dataset.

## Proposed workflow

1. Confirm that the libraries used a small-RNA protocol and preserve the original FASTQ files.
2. Inspect read quality, adapter content, and length distributions. Remove the documented small-RNA adapter and retain the biologically appropriate read-length range.
3. Align reads to the relevant reference genome and a versioned mature/precursor miRNA reference. Record how multi-mapping reads are handled.
4. Generate a raw integer count matrix with stable miRNA identifiers. Remove contaminants and low-count features using a stated rule.
5. Match count columns to metadata rows exactly. Fit `~ batch + condition` in DESeq2 when the same balanced design applies.
6. Extract treated versus control, shrink effect sizes, and apply `padj < 0.05` together with an effect-size threshold.
7. Report mapped-read proportions, retained miRNAs, PCA, a differential-expression plot, and the complete result table.

## Main limitations

Small-RNA conclusions depend on the library protocol, adapter sequence, reference version, read-length choice, and multi-mapping policy. Gene-level bulk RNA-seq counts cannot be relabeled as miRNA measurements, and target-pathway claims would require separately validated annotations and enrichment analysis.
