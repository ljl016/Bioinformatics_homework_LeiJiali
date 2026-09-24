# Week 5 Homework 1 AI Verification Log

## AI use

- Tool: OpenAI Codex.
- Date: 2026-09-22.
- Purpose: interpret the assignment, review and repair the starter workflow, run the analysis locally, and draft the interpretation and verification record.

## Expanded prompt used for final analysis and quality control

> Act as a bioinformatics analysis assistant and complete Week 5 Homework 1 using `Week5_Homework_Count_Matrix.csv`, `Week5_Homework_Sample_Metadata.csv`, and the supplied starter code. First verify that counts are non-negative integers, sample names and metadata rows match exactly, condition and batch are factors, control is the reference, and the model matrix is full rank. Use DESeq2 with `design = ~ batch + condition`; retain genes with at least 10 counts in at least 3 samples; inspect `resultsNames(dds)`; extract treated versus control; and apply `apeglm` shrinkage. Define significance as `padj < 0.05` and `|shrunken log2FC| >= 1`. Export the complete result table, a PCA colored by condition and shaped by batch, a thresholded volcano plot, the fitted object, design checks, and `sessionInfo()`. Write a 100–150 word interpretation based only on computed results and a short small-RNA extension plan. Preserve the original inputs, do not use the instructor key as analysis evidence, do not invent gene functions or pathways, and document accepted advice, rejected advice, exact independent checks, and one unresolved uncertainty.

This expanded prompt was used for the final quality-control and packaging pass on 2026-09-22. It restates the assignment-specific decisions that were actually applied.

## Earlier short prompts preserved verbatim

1. `那你先给我homework1的作业方案`
2. `你看下这个跟homework1有没有关系，有关系就补充`
3. `ok，就按你这个方案完成homework1`

The second prompt accompanied the course rubric image. The third prompt authorized completion according to the immediately preceding supplemented plan. These short prompts are retained for an honest audit trail rather than presented as the detailed prompt above.

## Advice accepted

- Keep raw integer counts unchanged and use `design = ~ batch + condition`.
- Relevel `condition` so that `control` is the reference, inspect `resultsNames(dds)`, and verify `condition_treated_vs_control` before extraction.
- Retain genes with at least 10 counts in at least 3 samples.
- Use `apeglm` shrinkage and define significance with both `padj < 0.05` and `|shrunken log2FC| >= 1`.
- Replace the starter's `vst()` call with `varianceStabilizingTransformation()` because only 989 genes remain after filtering, fewer than the default 1,000-gene subset requested by `vst()`.
- Export the complete 989-row result table, while sorting by adjusted p value so the top 20 features are immediately visible.

## Advice or possibilities rejected

- Exporting only the top 20 rows was rejected because the data-specific instructions require the complete shrunken result table.
- The instructor-key CSV was identified during the initial file audit as a teacher reference and excluded from the analysis. Its truth values were not used to select genes, tune thresholds, interpret results, or validate the differential-expression calls.
- Pathway and gene-function stories were rejected because the dataset supplies generic gene IDs and no validated functional or biotype annotation.
- Samples were not deleted based only on the PCA plot.
- Raw p values were not used to call significance.

## Exact checks performed independently

- The count matrix contains 1,000 genes and 12 samples, no missing values, no duplicate gene or sample IDs, and only non-negative integers.
- Count-matrix columns exactly match metadata row names in the same order.
- Metadata contains six control and six treated samples; batches A, B, and C each contain two samples per condition.
- The model matrix for `~ batch + condition` has full column rank.
- Filtering retained 989 genes and removed 11.
- `resultsNames(dds)` returned `Intercept`, `batch_B_vs_A`, `batch_C_vs_A`, and `condition_treated_vs_control`.
- The output table has 989 unique genes and no missing adjusted p values.
- The significance totals were recomputed from the exported CSV: 60 total, 36 up in treated, and 24 down in treated.
- PCA variance labels were checked from the fitted object: PC1 = 24.1% and PC2 = 9.3%.
- Both PNG files were opened and visually checked for readable axes, legends, comparison direction, and thresholds.
- The final script completed successfully when launched from the original H-drive directory containing Chinese characters.
- The copied inputs matched the originals by SHA-256: counts `2FD8B4D68D2184D6548C7AECC28E1C3DD5391AF4FE91759642961A3CBB1BE91F`; metadata `34B467DF46ACABD121911853C8B1F84D6D9C4794C6A5DE8EB70EEBFBF62E722C`.

## Unresolved uncertainty

The generic rubric names `week5_homework.R`, whereas the data-specific Homework 1 instructions require `week5_deseq2_analysis.R`. The data-specific filename was retained. More importantly, the absence of validated gene symbols and biotypes prevents gene-specific functional interpretation.

## Methodological resource

Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. arXiv:2609.00065. https://doi.org/10.48550/arXiv.2609.00065
