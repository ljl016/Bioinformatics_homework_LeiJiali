# Reproducible filtering for Week 4 Q4 synthetic teaching variants.
# Usage: Rscript q4_filter.R variants_q4.tsv results_directory

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript q4_filter.R <variants_q4.tsv> <results_dir>")

input_path <- args[1]
output_dir <- args[2]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

variants <- read.delim(input_path, comment.char = "#", stringsAsFactors = FALSE,
                       check.names = FALSE)
required <- c("CHROM", "POS", "REF", "ALT", "FILTER", "DP", "GQ", "AF",
              "GENE", "CONSEQUENCE", "CLINVAR_SIG", "CLINVAR_ID", "NOTE")
missing <- setdiff(required, names(variants))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))

variants$row_id <- seq_len(nrow(variants))
variants$pass_technical <- variants$FILTER == "PASS" & variants$DP >= 20 & variants$GQ >= 30
variants$pass_af_0.01 <- variants$AF <= 0.01
variants$pass_af_0.001 <- variants$AF <= 0.001

impact <- c("stop_gained", "frameshift_variant", "splice_acceptor_variant",
            "splice_donor_variant", "missense_variant")
variants$impact_consequence <- variants$CONSEQUENCE %in% impact

variants$consequence_points <- ifelse(
  variants$CONSEQUENCE %in% c("stop_gained", "frameshift_variant",
                              "splice_acceptor_variant", "splice_donor_variant"), 4,
  ifelse(variants$CONSEQUENCE == "missense_variant", 2, 0)
)
variants$rarity_points <- ifelse(variants$AF <= 0.0001, 2,
                                 ifelse(variants$AF <= 0.001, 1, 0))
variants$clinvar_points <- ifelse(
  variants$CLINVAR_SIG == "Pathogenic", 3,
  ifelse(variants$CLINVAR_SIG == "Conflicting_interpretations_of_pathogenicity", 1.5,
  ifelse(variants$CLINVAR_SIG == "Uncertain_significance", 1,
  ifelse(variants$CLINVAR_SIG == "Likely_benign", -2,
  ifelse(variants$CLINVAR_SIG == "Benign", -3, 0))))
)
variants$priority_score <- variants$consequence_points + variants$rarity_points +
                           variants$clinvar_points

variants$exclusion_reason <- ""
variants$exclusion_reason[variants$FILTER != "PASS"] <- "FILTER not PASS"
variants$exclusion_reason[variants$DP < 20] <- paste0(
  variants$exclusion_reason[variants$DP < 20],
  ifelse(variants$exclusion_reason[variants$DP < 20] == "", "", "; "), "DP < 20")
variants$exclusion_reason[variants$GQ < 30] <- paste0(
  variants$exclusion_reason[variants$GQ < 30],
  ifelse(variants$exclusion_reason[variants$GQ < 30] == "", "", "; "), "GQ < 30")

common <- variants$pass_technical & !variants$pass_af_0.01
variants$exclusion_reason[common] <- "AF > 0.01"
low_impact <- variants$pass_technical & variants$pass_af_0.01 & !variants$impact_consequence
variants$exclusion_reason[low_impact] <- "Lower-priority consequence without supporting evidence"

candidate <- variants$pass_technical & variants$pass_af_0.01 & variants$impact_consequence
shortlist <- variants[candidate, c(
  "row_id", "CHROM", "POS", "REF", "ALT", "GENE", "CONSEQUENCE", "DP", "GQ",
  "AF", "CLINVAR_SIG", "CLINVAR_ID", "consequence_points", "rarity_points",
  "clinvar_points", "priority_score", "NOTE"
)]
shortlist <- shortlist[order(-shortlist$priority_score, shortlist$AF, shortlist$row_id), ]
shortlist$rank <- seq_len(nrow(shortlist))
shortlist <- shortlist[, c("rank", setdiff(names(shortlist), "rank"))]

exclusions <- variants[!candidate, c(
  "row_id", "CHROM", "POS", "REF", "ALT", "GENE", "FILTER", "DP", "GQ", "AF",
  "CONSEQUENCE", "CLINVAR_SIG", "exclusion_reason"
)]

sensitivity <- data.frame(
  af_threshold = c(0.01, 0.001),
  candidate_count = c(
    sum(variants$pass_technical & variants$pass_af_0.01 & variants$impact_consequence),
    sum(variants$pass_technical & variants$pass_af_0.001 & variants$impact_consequence)
  ),
  candidate_genes = c(
    paste(variants$GENE[variants$pass_technical & variants$pass_af_0.01 &
                          variants$impact_consequence], collapse = ","),
    paste(variants$GENE[variants$pass_technical & variants$pass_af_0.001 &
                          variants$impact_consequence], collapse = ",")
  )
)

write.table(exclusions, file.path(output_dir, "q4_exclusion_log.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(shortlist, file.path(output_dir, "q4_shortlist.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(sensitivity, file.path(output_dir, "q4_sensitivity.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Small dataset-specific self-check: catches accidental changes to the locked rules.
stopifnot(nrow(variants) == 12, nrow(shortlist) == 4,
          shortlist$GENE[1] == "TP53", all(sensitivity$candidate_count == 4))

message("Wrote ", nrow(shortlist), " candidates and ", nrow(exclusions),
        " exclusions to ", normalizePath(output_dir))
