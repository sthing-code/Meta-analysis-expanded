# scripts/04_correlate.R
#
# Computes Pearson correlations between intratumoral bacterial abundance
# and gene expression for each cancer type.
# Applies BH FDR correction and filters by r threshold.
#
# Input:  data/processed/{PROJECT}_tpm_priority_genes.rds
#         data/processed/{PROJECT}_bacteria_processed.rds
# Output: data/processed/{PROJECT}_correlations.rds
#         (named list with full and filtered correlation matrices)

# ── Dependencies ──────────────────────────────────────────────────────────────

library(here)
library(dplyr)

source(here("R/utils.R"))

# ── Correlate one cancer type ─────────────────────────────────────────────────

correlate_cancer <- function(cancer_type) {

  cfg <- CANCER_CONFIG[[cancer_type]]
  project <- cfg$tcga_project

  message("\n========================================")
  message("  Correlating: ", project)
  message("========================================\n")

  # ── Load data ────────────────────────────────────────────────────────────────
  tcga_rds  <- here("data/processed", paste0(project, "_tpm_priority_genes.rds"))
  bact_rds  <- here("data/processed", paste0(project, "_bacteria_processed.rds"))

  if (!file.exists(tcga_rds)) stop("Run 02_preprocess_TCGA.R first.\n  Missing: ", tcga_rds)
  if (!file.exists(bact_rds)) stop("Run 03_preprocess_TCMbio.R first.\n  Missing: ", bact_rds)

  tcga_data <- readRDS(tcga_rds)
  bact_data <- readRDS(bact_rds)

  # Use log2(TPM+1) for gene expression
  # Use log10(abundance+1) for bacteria (already done in script 03)
  expr_mat <- tcga_data$log_tpm      # genes × samples
  bact_mat <- bact_data$abundance    # species × samples

  message("  Gene expression matrix: ", nrow(expr_mat), " genes × ", ncol(expr_mat), " samples")
  message("  Bacterial abundance matrix: ", nrow(bact_mat), " species × ", ncol(bact_mat), " samples")

  # ── Subset to priority genes only ────────────────────────────────────────────
  # Main heatmap uses priority genes; pathway genes are used in script 06.
  priority_present <- ALL_PRIORITY_GENES[ALL_PRIORITY_GENES %in% rownames(expr_mat)]
  expr_priority    <- expr_mat[priority_present, , drop = FALSE]
  message("  Priority genes for main heatmap: ", nrow(expr_priority))

  # ── Compute correlations ──────────────────────────────────────────────────────
  message("  Computing Pearson correlations (priority genes)...")
  cor_priority <- compute_pearson_matrix(bact_mat, expr_priority)

  message("  Computing Pearson correlations (pathway genes)...")
  pathway_present <- ALL_PATHWAY_GENES[ALL_PATHWAY_GENES %in% rownames(expr_mat)]
  expr_pathway    <- expr_mat[pathway_present, , drop = FALSE]
  cor_pathway     <- compute_pearson_matrix(bact_mat, expr_pathway)

  # ── Filter by r threshold ─────────────────────────────────────────────────────
  cor_r020 <- filter_by_r(cor_priority, r_threshold = 0.20)
  cor_r015 <- filter_by_r(cor_priority, r_threshold = 0.15)

  message("\n  Bacteria retained after |r| > 0.20 filter: ", nrow(cor_r020$r))
  message("  Bacteria retained after |r| > 0.15 filter: ", nrow(cor_r015$r))

  # ── Add literature benchmark annotation ──────────────────────────────────────
  # Uses flag_bacterium() from utils.R which normalises s__Genus_species format
  add_benchmark_flag <- function(cor_obj) {
    flags <- sapply(rownames(cor_obj$r),
                    function(b) flag_bacterium(b, cancer_type),
                    USE.NAMES = TRUE)
    cor_obj$benchmark_flag <- flags
    cor_obj
  }

  cor_priority <- add_benchmark_flag(cor_priority)
  cor_r020     <- add_benchmark_flag(cor_r020)
  cor_r015     <- add_benchmark_flag(cor_r015)

  # ── Summary statistics ────────────────────────────────────────────────────────
  sig_pairs <- sum(cor_r020$q < 0.05, na.rm = TRUE)
  message("\n  Significant pairs (|r| > 0.20, q < 0.05): ", sig_pairs)

  # Top positive correlations for IGFBP7 and Galectin-1
  for (gene in c("IGFBP7", "LGALS1")) {
    if (gene %in% colnames(cor_priority$r)) {
      top_pos <- sort(cor_priority$r[, gene], decreasing = TRUE)[1:5]
      message("\n  Top 5 positive correlations with ", gene, ":")
      print(round(top_pos, 3))
    }
  }

  # Top negative correlations for cathepsins
  cath_present <- PRIORITY_GENES$cathepsins[PRIORITY_GENES$cathepsins %in%
                                               colnames(cor_priority$r)]
  for (gene in cath_present) {
    top_neg <- sort(cor_priority$r[, gene], decreasing = FALSE)[1:5]
    message("\n  Top 5 negative correlations with ", gene, ":")
    print(round(top_neg, 3))
  }

  # ── Save ──────────────────────────────────────────────────────────────────────
  out <- list(
    full        = cor_priority,    # all bacteria × priority genes
    r020        = cor_r020,        # filtered |r| > 0.20
    r015        = cor_r015,        # filtered |r| > 0.15
    pathway     = cor_pathway,     # all bacteria × pathway genes (for script 06)
    cancer_type = cancer_type,
    project     = project,
    n_samples   = cor_priority$n_samples
  )

  rds_out <- here("data/processed", paste0(project, "_correlations.rds"))
  saveRDS(out, rds_out)
  message("\n  Saved: ", rds_out)

  return(invisible(out))
}

# ── Run ───────────────────────────────────────────────────────────────────────

# correlate_cancer("colon")
# correlate_cancer("breast")
# correlate_cancer("pancreatic")
correlate_cancer("prostate")

message("\nCorrelation analysis complete.")
