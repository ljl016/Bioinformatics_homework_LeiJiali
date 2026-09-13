options(timeout = 300)
suppressPackageStartupMessages(library(DESeq2))

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
output_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1])))
} else {
  getwd()
}
data_dir <- file.path(dirname(output_dir), "week3_homework_data")
dir.create(data_dir, showWarnings = FALSE)

download_once <- function(url, name) {
  path <- file.path(data_dir, name)
  if (!file.exists(path)) download.file(url, path, mode = "wb", quiet = FALSE)
  path
}

field_value <- function(lines, prefix) {
  hit <- lines[startsWith(lines, prefix)]
  if (!length(hit)) return(NA_character_)
  trimws(sub("^[^=]+=[[:space:]]*", "", hit[1]))
}

field_values <- function(lines, prefix) {
  trimws(sub("^[^=]+=[[:space:]]*", "", lines[startsWith(lines, prefix)]))
}

characteristic_value <- function(lines, key) {
  prefix <- "!Sample_characteristics_ch1 = "
  values <- sub(paste0("^", prefix), "", lines[startsWith(lines, prefix)])
  hit <- values[startsWith(tolower(values), paste0(tolower(key), ":"))]
  if (!length(hit)) return(NA_character_)
  trimws(sub("^[^:]+:[[:space:]]*", "", hit[1]))
}

read_soft_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  starts <- grep("^\\^SAMPLE = ", lines)
  ends <- c(starts[-1] - 1L, length(lines))
  rows <- Map(function(a, b) {
    block <- lines[a:b]
    relation <- block[startsWith(block, "!Sample_relation = SRA:")]
    srx <- if (length(relation)) sub(".*term=", "", relation[1]) else NA_character_
    data.frame(
      gsm = sub("^\\^SAMPLE = ", "", block[1]),
      title = field_value(block, "!Sample_title = "),
      platform = field_value(block, "!Sample_platform_id = "),
      srx = srx,
      organism = field_value(block, "!Sample_organism_ch1 = "),
      molecule = field_value(block, "!Sample_molecule_ch1 = "),
      library_strategy = field_value(block, "!Sample_library_strategy = "),
      extract_protocol = field_value(block, "!Sample_extract_protocol_ch1 = "),
      cohort = characteristic_value(block, "cohort"),
      tissue = characteristic_value(block, "tissue"),
      age = characteristic_value(block, "age"),
      bmi = characteristic_value(block, "bmi"),
      gender = characteristic_value(block, "gender"),
      patient_id = characteristic_value(block, "patientid"),
      recorded_timepoint = characteristic_value(block, "timepoint"),
      stage = characteristic_value(block, "transplant stage"),
      iri = characteristic_value(block, "ischemia reperfusion injury (iri)"),
      stringsAsFactors = FALSE
    )
  }, starts, ends)
  do.call(rbind, rows)
}

read_soft_series <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  list(
    title = field_value(lines, "!Series_title = "),
    summary = paste(field_values(lines, "!Series_summary = "), collapse = " "),
    type = field_value(lines, "!Series_type = "),
    supplementary_file = field_value(lines, "!Series_supplementary_file = ")
  )
}

read_counts <- function(path, metadata_samples) {
  x <- read.delim(gzfile(path), check.names = FALSE)
  annotation_cols <- ncol(x) - metadata_samples
  expression <- x[, seq.int(annotation_cols + 1L, ncol(x)), drop = FALSE]
  is_count <- all(vapply(expression, function(v) {
    is.numeric(v) && all(is.finite(v)) && all(v >= 0) && all(v == floor(v))
  }, logical(1)))
  list(rows = nrow(x), annotation_cols = annotation_cols,
       samples = ncol(expression), is_count = is_count,
       sample_names = names(expression), expression = expression)
}

normalize_sample_id <- function(x) {
  x <- sub("^s_", "", x)
  x <- sub("_GeneCount$", "", x)
  x <- sub("^Sample_", "", x)
  x <- sub("\\.bam$", "", x)
  tolower(gsub("-", "", x))
}

geo_root <- "https://ftp.ncbi.nlm.nih.gov/geo/series"

# 1. OBTAIN: read the official files and calculate the answers from them.
g146_soft <- download_once(
  paste0(geo_root, "/GSE146nnn/GSE146853/soft/GSE146853_family.soft.gz"),
  "GSE146853_family.soft.gz"
)
g146_counts <- download_once(
  paste0(geo_root, "/GSE146nnn/GSE146853/suppl/GSE146853_GeneCount_raw.txt.gz"),
  "GSE146853_GeneCount_raw.txt.gz"
)
g874_soft <- download_once(
  paste0(geo_root, "/GSE87nnn/GSE87487/soft/GSE87487_family.soft.gz"),
  "GSE87487_family.soft.gz"
)
g874_counts <- download_once(
  paste0(geo_root, "/GSE87nnn/GSE87487/suppl/GSE87487_counts.20samples.txt.gz"),
  "GSE87487_counts.20samples.txt.gz"
)

g146 <- read_soft_samples(g146_soft)
g146$subject <- sub("-[12]$", "", g146$title)
g146$timepoint <- ifelse(grepl("-1$", g146$title), "T1",
                         ifelse(grepl("-2$", g146$title), "T2", NA_character_))
g146$cohort[g146$cohort == "Diarrheaiarrhea"] <- "Diarrhea"
g146_check <- read_counts(g146_counts, nrow(g146))

subject_frequency_146 <- table(g146$subject)
visit_frequency_counts_146 <- table(subject_frequency_146)
cohort_counts_146 <- table(g146$cohort, useNA = "ifany")
timepoint_counts_146 <- table(g146$timepoint, useNA = "ifany")
first_visit_146 <- subset(g146, timepoint == "T1" & !is.na(cohort) & cohort != "NA")
first_visit_146$analysis_group <- ifelse(first_visit_146$cohort == "Healthy", "Healthy", "IBS")

# Select an example by biological/design criteria, not by a known accession.
example_row <- subset(g146, cohort == "Healthy" & timepoint == "T1")
example_row <- example_row[order(example_row$subject, example_row$title), , drop = FALSE][1, ]
example_gsm <- example_row$gsm
example_srx <- example_row$srx  # obtained from !Sample_relation = SRA: in GEO SOFT
runinfo_path <- download_once(
  sprintf("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=%s",
          URLencode(example_srx)),
  paste0(example_srx, "_runinfo.csv")
)
runinfo <- read.csv(runinfo_path, check.names = FALSE)
example_srr <- runinfo$Run[1]
example_project <- runinfo$BioProject[1]
example_study <- runinfo$SRAStudy[1]

g874 <- read_soft_samples(g874_soft)
g874_series <- read_soft_series(g874_soft)
g874$subject <- sub("(?i)-?bx[12]$", "", g874$title, perl = TRUE)
g874_check <- read_counts(g874_counts, nrow(g874))
subject_frequency_874 <- table(g874$subject)
stage_counts_874 <- table(g874$stage)
iri_counts_874 <- table(g874$iri)
iri_subject_counts_874 <- tapply(g874$subject, g874$iri,
                                 function(x) length(unique(x)))
stages_by_subject_874 <- table(g874$subject, g874$stage)
g874_complete_pairs <- all(stages_by_subject_874[, c("pre-reperfusion", "post-reperfusion")] == 1L)
g874_samples_independent <- all(subject_frequency_874 == 1L)

g874_organisms <- unique(na.omit(g874$organism))
g874_tissues <- unique(na.omit(g874$tissue))
g874_library_strategies <- unique(na.omit(g874$library_strategy))
g874_molecules <- unique(na.omit(g874$molecule))
g874_is_liver_allograft_biopsy <- identical(g874_tissues, "Liver") &&
  grepl("allograft", g874_series$summary, ignore.case = TRUE) &&
  any(grepl("biops", g874$extract_protocol, ignore.case = TRUE))
g874_is_bulk_tissue <- g874_is_liver_allograft_biopsy &&
  identical(g874_molecules, "total RNA") &&
  identical(g874_library_strategies, "RNA-Seq")
g874_count_size_mib <- unname(file.info(g874_counts)$size / 1024^2)

# Discover the GSE87487 SRA accessions from SOFT, then prove raw runs exist for all samples.
example_row_874 <- subset(g874, stage == "pre-reperfusion")
example_row_874 <- example_row_874[order(example_row_874$gsm), , drop = FALSE][1, ]
example_gsm_874 <- example_row_874$gsm
example_srx_874 <- example_row_874$srx
runinfo_path_874 <- download_once(
  sprintf("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=%s",
          URLencode(example_srx_874)),
  paste0(example_srx_874, "_runinfo.csv")
)
runinfo_874 <- read.csv(runinfo_path_874, check.names = FALSE)
example_srr_874 <- sort(unique(runinfo_874$Run))[1]
study_874 <- unique(runinfo_874$SRAStudy)
project_874 <- unique(runinfo_874$BioProject)

study_runinfo_path_874 <- download_once(
  sprintf("https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=%s",
          URLencode(study_874)),
  paste0(study_874, "_runinfo.csv")
)
study_runinfo_874 <- read.csv(study_runinfo_path_874, check.names = FALSE)
study_srr_874 <- unique(study_runinfo_874$Run)
study_srx_874 <- unique(study_runinfo_874$Experiment)
study_gsm_874 <- unique(study_runinfo_874$SampleName)
g874_raw_reprocessing_possible <- setequal(g874$srx, study_srx_874) &&
  setequal(g874$gsm, study_gsm_874) &&
  all(!is.na(study_runinfo_874$download_path) & nzchar(study_runinfo_874$download_path))

comparison_answer_874 <- if (!g874_samples_independent && g874_complete_pairs) {
  c("Samples are not independent.",
    "Use paired post-vs-pre comparison blocking on subject.")
} else {
  "Pairing could not be established from the metadata."
}

paired_design_874 <- g874[, c("gsm", "title", "subject", "stage", "iri", "platform", "srx")]
paired_design_874$stage <- factor(paired_design_874$stage,
                                 levels = c("pre-reperfusion", "post-reperfusion"))
design_matrix_874 <- model.matrix(~ subject + stage, data = paired_design_874)

# Integrity, cleaning, and a focused QC figure. No differential expression is run.
g146_count_ids <- normalize_sample_id(g146_check$sample_names)
g874_count_ids <- normalize_sample_id(g874_check$sample_names)
g146_metadata_ids <- normalize_sample_id(g146$title)
g874_metadata_ids <- normalize_sample_id(g874$title)

g146_aligned <- setequal(g146_count_ids, g146_metadata_ids)
g874_aligned <- setequal(g874_count_ids, g874_metadata_ids)
g146_duplicate_gsm <- sum(duplicated(g146$gsm))
g874_duplicate_gsm <- sum(duplicated(g874$gsm))
g146_duplicate_design <- sum(duplicated(paste(g146$subject, g146$timepoint)))
g874_duplicate_design <- sum(duplicated(paste(g874$subject, g874$stage)))

g146_analysis_index <- match(normalize_sample_id(first_visit_146$title), g146_count_ids)
g146_analysis_counts <- as.matrix(g146_check$expression[, g146_analysis_index, drop = FALSE])
colnames(g146_analysis_counts) <- first_visit_146$title
min_samples_146 <- min(table(first_visit_146$analysis_group))
keep_146 <- rowSums(g146_analysis_counts >= 10) >= min_samples_146

g874_analysis_index <- match(g874_metadata_ids, g874_count_ids)
g874_analysis_counts <- as.matrix(g874_check$expression[, g874_analysis_index, drop = FALSE])
colnames(g874_analysis_counts) <- g874$title
min_samples_874 <- min(stage_counts_874)
keep_874 <- rowSums(g874_analysis_counts >= 10) >= min_samples_874

col_data_146 <- data.frame(
  analysis_group = factor(first_visit_146$analysis_group,
                          levels = c("Healthy", "IBS")),
  cohort = factor(first_visit_146$cohort),
  row.names = first_visit_146$title
)
dds_146 <- DESeqDataSetFromMatrix(
  countData = g146_analysis_counts[keep_146, , drop = FALSE],
  colData = col_data_146,
  design = ~ analysis_group
)
vst_146 <- vst(dds_146, blind = TRUE)
pca_146 <- prcomp(t(assay(vst_146)))
pca_percent_146 <- 100 * pca_146$sdev^2 / sum(pca_146$sdev^2)

col_data_874 <- data.frame(
  subject = factor(g874$subject),
  stage = factor(gsub("-", "_", g874$stage),
                 levels = c("pre_reperfusion", "post_reperfusion")),
  iri = factor(g874$iri, levels = c("negative", "positive")),
  row.names = g874$title
)
dds_874 <- DESeqDataSetFromMatrix(
  countData = g874_analysis_counts[keep_874, , drop = FALSE],
  colData = col_data_874,
  design = ~ subject + stage
)
vst_874 <- vst(dds_874, blind = TRUE)
pca_874 <- prcomp(t(assay(vst_874)))
pca_percent_874 <- 100 * pca_874$sdev^2 / sum(pca_874$sdev^2)

cleaning_summary <- data.frame(
  dataset = c("GSE146853", "GSE87487"),
  raw_features = c(g146_check$rows, g874_check$rows),
  annotation_columns = c(g146_check$annotation_cols, g874_check$annotation_cols),
  count_columns = c(g146_check$samples, g874_check$samples),
  metadata_samples = c(nrow(g146), nrow(g874)),
  sample_alignment = c(g146_aligned, g874_aligned),
  missing_key_metadata = c(sum(is.na(g146$cohort) | g146$cohort == "NA"),
                           sum(is.na(g874$stage) | is.na(g874$iri))),
  duplicate_gsm = c(g146_duplicate_gsm, g874_duplicate_gsm),
  duplicate_subject_stage = c(g146_duplicate_design, g874_duplicate_design),
  analysis_samples = c(ncol(g146_analysis_counts), ncol(g874_analysis_counts)),
  filtering_rule = c(
    sprintf("count >= 10 in at least %d first-visit samples", min_samples_146),
    sprintf("count >= 10 in at least %d samples", min_samples_874)
  ),
  features_before_filter = c(nrow(g146_analysis_counts), nrow(g874_analysis_counts)),
  features_after_filter = c(sum(keep_146), sum(keep_874)),
  deseq2_design = c("~ analysis_group", "~ subject + stage"),
  deseq2_dimensions = c(
    sprintf("%d features x %d samples", nrow(dds_146), ncol(dds_146)),
    sprintf("%d features x %d samples", nrow(dds_874), ncol(dds_874))
  ),
  vst_completed = c(TRUE, TRUE),
  stringsAsFactors = FALSE
)
write.table(cleaning_summary, file.path(output_dir, "cleaning_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

plot_colors <- c(Healthy = "#D97706", IBS = "#2563EB")
png(file.path(output_dir, "qc_vst_pca.png"),
    width = 2000, height = 1000, res = 180)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.8, 1))
plot(pca_146$x[, 1], pca_146$x[, 2],
     col = plot_colors[col_data_146$analysis_group], pch = 19, cex = 1.25,
     xlab = sprintf("PC1 (%.1f%%)", pca_percent_146[1]),
     ylab = sprintf("PC2 (%.1f%%)", pca_percent_146[2]),
     main = "GSE146853 first-visit PCA")
legend("topright", legend = names(plot_colors), col = plot_colors,
       pch = 19, bty = "n")

stage_colors <- c("pre_reperfusion" = "#6B7280",
                  "post_reperfusion" = "#7C3AED")
plot(pca_874$x[, 1], pca_874$x[, 2], type = "n",
     xlab = sprintf("PC1 (%.1f%%)", pca_percent_874[1]),
     ylab = sprintf("PC2 (%.1f%%)", pca_percent_874[2]),
     main = "GSE87487 paired-sample PCA")
for (subject_id in levels(col_data_874$subject)) {
  pair_index <- which(col_data_874$subject == subject_id)
  segments(pca_874$x[pair_index[1], 1], pca_874$x[pair_index[1], 2],
           pca_874$x[pair_index[2], 1], pca_874$x[pair_index[2], 2],
           col = "#D1D5DB", lwd = 1.2)
}
points(pca_874$x[, 1], pca_874$x[, 2],
       col = stage_colors[col_data_874$stage], pch = 19, cex = 1.25)
legend("topright", legend = c("pre-reperfusion", "post-reperfusion"), col = stage_colors,
       pch = 19, bty = "n")
dev.off()

pca_group <- as.character(col_data_146$analysis_group)
centroid_pc1 <- tapply(pca_146$x[, 1], pca_group, mean)
centroid_pc2 <- tapply(pca_146$x[, 2], pca_group, mean)
centroid_distance <- sqrt(diff(centroid_pc1)^2 + diff(centroid_pc2)^2)
within_distance <- mean(sqrt(
  (pca_146$x[, 1] - centroid_pc1[pca_group])^2 +
  (pca_146$x[, 2] - centroid_pc2[pca_group])^2
))

pca_stage <- as.character(col_data_874$stage)
centroid_pc1_874 <- tapply(pca_874$x[, 1], pca_stage, mean)
centroid_pc2_874 <- tapply(pca_874$x[, 2], pca_stage, mean)
centroid_distance_874 <- sqrt(diff(centroid_pc1_874)^2 + diff(centroid_pc2_874)^2)
within_distance_874 <- mean(sqrt(
  (pca_874$x[, 1] - centroid_pc1_874[pca_stage])^2 +
  (pca_874$x[, 2] - centroid_pc2_874[pca_stage])^2
))
qc_interpretation <- c(
  sprintf("After filtering and VST, PC1 and PC2 explain %.1f%% and %.1f%% of variance in GSE146853 and %.1f%% and %.1f%% in GSE87487.",
          pca_percent_146[1], pca_percent_146[2],
          pca_percent_874[1], pca_percent_874[2]),
  sprintf("The first two PCs show %s IBS-versus-healthy separation in GSE146853 and %s pre-versus-post separation in GSE87487; the gray lines identify paired biopsies, and neither panel is a differential-expression test.",
          if (centroid_distance <= within_distance) "limited" else "visible",
          if (centroid_distance_874 <= within_distance_874) "limited" else "visible")
)
writeLines(qc_interpretation, file.path(data_dir, "qc_interpretation.txt"), useBytes = TRUE)

retrieval_record <- data.frame(
  accession = c("GSE146853", "GSE146853", example_srx, "GSE87487", "GSE87487",
                example_srx_874, study_874),
  checked_on = as.character(Sys.Date()),
  file_name = basename(c(g146_soft, g146_counts, runinfo_path, g874_soft, g874_counts,
                         runinfo_path_874, study_runinfo_path_874)),
  file_type = c("GEO SOFT gzip", "raw gene-count TSV gzip", "SRA RunInfo CSV",
                "GEO SOFT gzip", "raw gene-count TSV gzip", "SRA RunInfo CSV",
                "SRA RunInfo CSV"),
  size_bytes = file.info(c(g146_soft, g146_counts, runinfo_path, g874_soft, g874_counts,
                           runinfo_path_874, study_runinfo_path_874))$size,
  stringsAsFactors = FALSE
)
write.table(retrieval_record, file.path(data_dir, "retrieval_record.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.csv(first_visit_146, file.path(data_dir, "GSE146853_first_visit_design.csv"), row.names = FALSE)
write.csv(paired_design_874, file.path(data_dir, "GSE87487_paired_design.csv"), row.names = FALSE)
write.csv(design_matrix_874, file.path(data_dir, "GSE87487_design_matrix.csv"), row.names = FALSE)

# 2. DISPLAY: print the values just obtained before checking expected answers.
summary_lines <- c(
  sprintf("GSE146853 count matrix: %d genes x %d samples; non-negative integers: %s.",
          g146_check$rows, g146_check$samples, g146_check$is_count),
  sprintf("GSE146853 metadata: %d samples from %d subjects.",
          nrow(g146), length(unique(g146$subject))),
  sprintf("GSE146853 repeat structure: %d subjects have two visits and %d have one visit.",
          visit_frequency_counts_146[["2"]], visit_frequency_counts_146[["1"]]),
  sprintf("GSE146853 cohorts: IBS-C %d, IBS-D %d, Healthy %d, missing cohort %d; T1 %d and T2 %d.",
          cohort_counts_146[["Constipation"]], cohort_counts_146[["Diarrhea"]],
          cohort_counts_146[["Healthy"]], cohort_counts_146[["NA"]],
          timepoint_counts_146[["T1"]], timepoint_counts_146[["T2"]]),
  sprintf("Primary first-visit contrast: IBS n=%d versus Healthy n=%d.",
          sum(first_visit_146$analysis_group == "IBS"), sum(first_visit_146$analysis_group == "Healthy")),
  sprintf("GSE146853 preprocessing: %d of %d features retained; DESeq2 object %d features x %d samples; VST and PCA completed without differential-expression testing.",
          sum(keep_146), nrow(g146_analysis_counts), nrow(dds_146), ncol(dds_146)),
  sprintf("Accession chain: GSE146853 -> %s -> %s -> %s; BioProject %s; SRA Study %s.",
          example_gsm, example_srx, example_srr, example_project, example_study),
  sprintf("GSE87487 count matrix: %d genes x %d samples; non-negative integers: %s.",
          g874_check$rows, g874_check$samples, g874_check$is_count),
  sprintf("GSE87487 matrix structure: %d annotation columns; compressed size %.1f MiB; file %s.",
          g874_check$annotation_cols, g874_count_size_mib,
          basename(g874_series$supplementary_file)),
  sprintf("GSE87487 GEO metadata: organism %s; series type %s; library strategy %s; tissue %s.",
          paste(g874_organisms, collapse = ", "), g874_series$type,
          paste(g874_library_strategies, collapse = ", "),
          paste(g874_tissues, collapse = ", ")),
  sprintf("GSE87487 material checks: liver allograft biopsies %s; bulk-tissue RNA-seq inference %s.",
          g874_is_liver_allograft_biopsy, g874_is_bulk_tissue),
  sprintf("GSE87487 design: %d samples from %d subjects, each sampled before and after reperfusion.",
          nrow(g874), length(unique(g874$subject))),
  sprintf("GSE87487 IRI status: %d IRI-positive subjects (%d samples) and %d IRI-negative subjects (%d samples).",
          iri_subject_counts_874[["positive"]], iri_counts_874[["positive"]],
          iri_subject_counts_874[["negative"]], iri_counts_874[["negative"]]),
  sprintf("GSE87487 preprocessing: %d of %d features retained; DESeq2 object %d features x %d samples with design ~ subject + stage; VST and PCA completed without differential-expression testing.",
          sum(keep_874), nrow(g874_analysis_counts), nrow(dds_874), ncol(dds_874)),
  sprintf("GSE87487 accession chain: GSE87487 -> %s -> %s -> %s; BioProject %s; SRA Study %s.",
          example_gsm_874, example_srx_874, example_srr_874,
          project_874, study_874),
  sprintf("The discovered study contains %d SRX experiments, %d SRR runs, and %d GSM samples; all GSE87487 samples have raw-read download paths: %s.",
          length(study_srx_874), length(study_srr_874), length(study_gsm_874),
          g874_raw_reprocessing_possible),
  comparison_answer_874
)

writeLines(summary_lines, file.path(data_dir, "verification_summary.txt"), useBytes = TRUE)
cat("DISCOVERED RESULTS\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\nGSE146853 cohort counts\n")
print(cohort_counts_146)
cat("\nGSE146853 timepoint counts\n")
print(timepoint_counts_146)
cat("\nGSE87487 stage counts\n")
print(stage_counts_874)
cat("\nGSE87487 IRI sample counts\n")
print(iri_counts_874)

# 3. VALIDATE: only now compare the discovered values with expected QC values.
cat("\nQC VALIDATION\n")
stopifnot(
  g146_check$rows == 64253L,
  g146_check$annotation_cols == 6L,
  g146_check$samples == 68L,
  g146_check$is_count,
  nrow(g146) == 68L,
  length(unique(g146$subject)) == 42L,
  unname(table(subject_frequency_146))[match(c("1", "2"), names(table(subject_frequency_146)))] == c(16L, 26L),
  unname(cohort_counts_146[c("Constipation", "Diarrhea", "Healthy", "NA")]) == c(23L, 23L, 21L, 1L),
  unname(timepoint_counts_146[c("T1", "T2")]) == c(40L, 28L),
  unname(table(first_visit_146$analysis_group)[c("IBS", "Healthy")]) == c(26L, 13L),
  example_gsm == "GSM4407890",
  example_srx == "SRX7899639",
  example_srr == "SRR11294080",
  example_project == "PRJNA612180",
  example_study == "SRP252532",
  g874_check$rows == 60498L,
  g874_check$annotation_cols == 6L,
  g874_check$samples == 20L,
  g874_check$is_count,
  nrow(g874) == 20L,
  length(unique(g874$subject)) == 10L,
  all(subject_frequency_874 == 2L),
  unname(stage_counts_874[c("pre-reperfusion", "post-reperfusion")]) == c(10L, 10L),
  unname(iri_counts_874[c("positive", "negative")]) == c(8L, 12L),
  unname(iri_subject_counts_874[c("positive", "negative")]) == c(4L, 6L),
  !g874_samples_independent,
  g874_complete_pairs,
  identical(g874_organisms, "Homo sapiens"),
  identical(g874_library_strategies, "RNA-Seq"),
  g874_is_liver_allograft_biopsy,
  g874_is_bulk_tissue,
  round(g874_count_size_mib, 1) == 3.7,
  basename(g874_series$supplementary_file) == "GSE87487_counts.20samples.txt.gz",
  example_gsm_874 == "GSM2332515",
  example_srx_874 == "SRX2199841",
  example_srr_874 == "SRR4305577",
  study_874 == "SRP090633",
  project_874 == "PRJNA344898",
  length(study_srx_874) == 20L,
  length(study_srr_874) == 76L,
  length(study_gsm_874) == 20L,
  g874_raw_reprocessing_possible
)
stopifnot(
  g146_aligned,
  g874_aligned,
  !anyNA(g146_analysis_index),
  !anyNA(g874_analysis_index),
  g146_duplicate_gsm == 0L,
  g874_duplicate_gsm == 0L,
  g146_duplicate_design == 0L,
  g874_duplicate_design == 0L,
  sum(keep_146) > 0L,
  sum(keep_874) > 0L,
  nrow(dds_146) == sum(keep_146),
  ncol(dds_146) == ncol(g146_analysis_counts),
  identical(dim(assay(vst_146)), dim(dds_146)),
  nrow(dds_874) == sum(keep_874),
  ncol(dds_874) == ncol(g874_analysis_counts),
  identical(dim(assay(vst_874)), dim(dds_874)),
  identical(all.vars(design(dds_146)), "analysis_group"),
  identical(all.vars(design(dds_874)), c("subject", "stage")),
  file.exists(file.path(output_dir, "cleaning_summary.tsv")),
  file.exists(file.path(output_dir, "qc_vst_pca.png"))
)
cat("All expected QC checks passed.\n")
