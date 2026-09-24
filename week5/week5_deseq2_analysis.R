# Week 5 Homework 1: bulk RNA-seq differential expression with DESeq2

if (.Platform$OS.type == "windows") {
  invisible(suppressWarnings(Sys.setlocale("LC_CTYPE", ".UTF-8")))
}

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) == 1) {
  setwd(dirname(normalizePath(sub("^--file=", "", script_arg))))
}

required_packages <- c("DESeq2", "apeglm", "ggplot2", "ggrepel")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required packages before running: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
})

input_names <- c(
  "Week5_Homework_Count_Matrix.csv",
  "Week5_Homework_Sample_Metadata.csv"
)
input_dir <- if (all(file.exists(input_names))) {
  "."
} else if (all(file.exists(file.path("..", "Source_Materials", input_names)))) {
  file.path("..", "Source_Materials")
} else if (all(file.exists(file.path("..", input_names)))) {
  ".."
} else {
  stop("Input CSV files were not found beside the script or in Source_Materials.")
}
count_file <- file.path(input_dir, input_names[1])
metadata_file <- file.path(input_dir, input_names[2])
output_files <- c(
  "week5_deseq2_analysis.R", "week5_deseq2_results.csv", "week5_pca.png",
  "week5_de_plot.png", "week5_interpretation.md",
  "week5_AI_verification_log.md", "week5_deseq2_object.rds",
  "session_info.txt", "design_and_checks.tsv", "RNA_extension.md"
)

counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
coldata <- read.csv(metadata_file, row.names = 1, check.names = FALSE)
count_matrix <- as.matrix(counts)

stopifnot(
  !anyDuplicated(rownames(counts)),
  !anyDuplicated(colnames(counts)),
  !anyDuplicated(rownames(coldata)),
  ncol(counts) == nrow(coldata),
  identical(colnames(counts), rownames(coldata)),
  !anyNA(count_matrix),
  all(count_matrix >= 0),
  all(count_matrix == round(count_matrix)),
  setequal(coldata$condition, c("control", "treated")),
  setequal(coldata$batch, c("A", "B", "C"))
)

coldata$condition <- relevel(factor(coldata$condition), ref = "control")
coldata$batch <- factor(coldata$batch)
design_formula <- ~ batch + condition
design_matrix <- model.matrix(design_formula, coldata)
stopifnot(qr(design_matrix)$rank == ncol(design_matrix))

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = coldata,
  design = design_formula
)

genes_before <- nrow(dds)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
genes_after <- nrow(dds)
dds <- DESeq(dds)

coefficient_names <- resultsNames(dds)
target_coef <- "condition_treated_vs_control"
stopifnot(target_coef %in% coefficient_names)

res <- results(dds, name = target_coef, alpha = 0.05)
res_shrunk <- lfcShrink(dds, coef = target_coef, res = res, type = "apeglm")
res_df <- data.frame(gene_id = rownames(res_shrunk), as.data.frame(res_shrunk), check.names = FALSE)
res_df$significant <- !is.na(res_df$padj) &
  res_df$padj < 0.05 & abs(res_df$log2FoldChange) >= 1
res_df$direction <- ifelse(
  res_df$significant,
  ifelse(res_df$log2FoldChange > 0, "Up in treated", "Down in treated"),
  "Not significant"
)
res_df <- res_df[order(is.na(res_df$padj), res_df$padj, -abs(res_df$log2FoldChange)), ]
write.csv(res_df, "week5_deseq2_results.csv", row.names = FALSE)

set.seed(20260922)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
pca_df <- plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE)
percent_var <- round(100 * attr(pca_df, "percentVar"))

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = condition, shape = batch, label = name)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3, max.overlaps = Inf, seed = 20260922) +
  labs(
    title = "PCA of variance-stabilized counts",
    subtitle = "Color: condition; shape: batch",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw(base_size = 12)
ggsave("week5_pca.png", p_pca, width = 7, height = 5, dpi = 300)

plot_df <- res_df[!is.na(res_df$padj), ]
plot_df$neg_log10_padj <- -log10(pmax(plot_df$padj, 1e-300))
p_volcano <- ggplot(
  plot_df,
  aes(log2FoldChange, neg_log10_padj, color = direction, shape = direction)
) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_color_manual(values = c(
    "Up in treated" = "#C0392B",
    "Down in treated" = "#2F6DB3",
    "Not significant" = "grey70"
  )) +
  scale_shape_manual(values = c(
    "Up in treated" = 17,
    "Down in treated" = 15,
    "Not significant" = 16
  )) +
  labs(
    title = "Differential expression: treated versus control",
    subtitle = "Significant: padj < 0.05 and |shrunken log2FC| >= 1",
    x = "Shrunken log2 fold change (treated / control)",
    y = "-log10 adjusted p value",
    color = NULL,
    shape = NULL
  ) +
  theme_bw(base_size = 12)
ggsave("week5_de_plot.png", p_volcano, width = 7, height = 5, dpi = 300)

checks <- data.frame(
  item = c(
    "rna_class", "library_value_type", "count_dimensions", "sample_alignment",
    "nonnegative_integer_counts", "biological_unit", "condition_levels",
    "batch_levels", "reference_level", "design_formula", "model_full_rank",
    "filter_rule", "genes_before_filter", "genes_after_filter",
    "coefficient", "contrast_direction", "significance_rule", "random_seed",
    "input_preservation"
  ),
  status = c(
    "documented", "verified", "verified", "pass", "pass", "documented",
    "pass", "pass", "pass", "verified", "pass", "applied", "recorded",
    "recorded", "verified", "verified", "applied", "recorded", "pass"
  ),
  value = c(
    "bulk RNA-seq", "raw non-negative integer counts",
    paste(nrow(counts), "genes x", ncol(counts), "samples"),
    "count columns exactly match metadata rows", "TRUE", "independent sample",
    paste(levels(coldata$condition), collapse = ","),
    paste(levels(coldata$batch), collapse = ","), levels(coldata$condition)[1],
    "~ batch + condition", "TRUE", "count >= 10 in at least 3 samples",
    genes_before, genes_after, target_coef, "treated versus control",
    "padj < 0.05 and abs(shrunken log2FC) >= 1", "20260922",
    "original CSV inputs read only"
  ),
  stringsAsFactors = FALSE
)
write.table(checks, "design_and_checks.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(dds, "week5_deseq2_object.rds")
session_lines <- c(
  "Week 5 Homework 1 reproducibility record",
  paste("Working directory:", normalizePath(getwd(), winslash = "/")),
  paste("Input files preserved:", paste(c(count_file, metadata_file), collapse = ", ")),
  paste("Output files:", paste(output_files, collapse = ", ")),
  "Random seed: 20260922 (plot label placement only)",
  "",
  capture.output(sessionInfo())
)
writeLines(session_lines, "session_info.txt")

cat("Genes before/after filtering:", genes_before, "/", genes_after, "\n")
cat("Coefficients:", paste(coefficient_names, collapse = ", "), "\n")
cat("Significant genes:", sum(res_df$significant), "\n")
print(table(res_df$direction))
