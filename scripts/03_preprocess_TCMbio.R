# scripts/03_preprocess_TCMbio.R
#
# Loads and preprocesses intratumoral bacterial abundance data downloaded
# from TCMbio (https://microbiomex.sdu.edu.cn/) for each cancer type.
#
# ── HOW TO DOWNLOAD TCMBIO DATA ───────────────────────────────────────────────
#
#  1. Go to https://microbiomex.sdu.edu.cn/
#  2. Select the cancer type (e.g. "COAD" for colon adenocarcinoma)
#  3. Navigate to the "Download" or "Data" section
#  4. Download the species-level bacterial abundance matrix as CSV
#     Expected format: rows = bacterial species, columns = TCGA sample barcodes
#     Abundance values are typically raw read counts or RPM from the
#     Kraken2/Bracken pipeline applied to TCGA RNA-seq non-human reads.
#  5. Save the file to:
#       data/raw/TCMbio_COAD_bacteria.csv   (for colon)
#       data/raw/TCMbio_BRCA_bacteria.csv   (for breast)
#       data/raw/TCMbio_PAAD_bacteria.csv   (for pancreatic)
#       data/raw/TCMbio_PRAD_bacteria.csv   (for prostate)
#
# ── ALTERNATIVE: cBioPortal microbiome data ───────────────────────────────────
#
#  If TCMbio download is unavailable, TCGA COAD microbiome data can be obtained
#  from cBioPortal (Poore et al. 2020 Nature pipeline):
#  https://www.cbioportal.org/
#  Select "Colon Adenocarcinoma (TCGA, PanCancer Atlas)" → download microbiome
#  abundance table. Column format and preprocessing steps are identical.
#
# ─────────────────────────────────────────────────────────────────────────────
#
# Input:  data/raw/TCMbio_{CANCER}_bacteria.csv
# Output: data/processed/{PROJECT}_bacteria_processed.rds
#         (named list: $abundance = species × samples log-transformed matrix,
#                       $raw       = species × samples raw matrix,
#                       $species   = character vector of retained species)

# ── Dependencies ──────────────────────────────────────────────────────────────

library(data.table)
library(dplyr)
library(here)

source(here("R/utils.R"))

# ── Load function ─────────────────────────────────────────────────────────────

#' Load and preprocess one TCMbio CSV file.
#'
#' cancer_type  : one of "colon", "breast", "pancreatic", "prostate"
#' min_prev     : minimum prevalence fraction for species retention (default 0.10)
#'
#' The function is permissive about CSV layout — it detects whether species
#' are in rows or columns and transposes if needed.

preprocess_tcmbio <- function(cancer_type, min_prev = 0.10) {

  cfg <- CANCER_CONFIG[[cancer_type]]
  if (is.null(cfg)) stop("Unknown cancer type: ", cancer_type)

  csv_path <- here(cfg$tcmbio_file)

  if (!file.exists(csv_path)) {
    stop(
      "TCMbio data file not found: ", csv_path,
      "\n\nPlease download from https://microbiomex.sdu.edu.cn/ and save to ",
      csv_path, "\nSee the instructions at the top of this script for details."
    )
  }

  message("\n========================================")
  message("  Loading TCMbio data: ", cancer_type)
  message("========================================\n")
  message("  File: ", csv_path)

  # ── Read CSV ────────────────────────────────────────────────────────────────
  raw_dt <- fread(csv_path, header = TRUE)
  message("  Raw dimensions: ", nrow(raw_dt), " rows × ", ncol(raw_dt), " cols")

  # Convert to matrix — first column expected to be species names or sample IDs
  first_col <- raw_dt[[1]]
  mat       <- as.matrix(raw_dt[, -1, with = FALSE])
  rownames(mat) <- first_col

  # ── Detect orientation ──────────────────────────────────────────────────────
  # TCMbio typically provides species as rows, TCGA barcodes as columns.
  # Detect by checking whether rownames or colnames contain "TCGA-".
  col_is_tcga <- grepl("^TCGA-", colnames(mat)[1])
  row_is_tcga <- grepl("^TCGA-", rownames(mat)[1])

  if (row_is_tcga && !col_is_tcga) {
    message("  Detected: samples in rows — transposing to species × samples")
    mat <- t(mat)
  } else if (!col_is_tcga && !row_is_tcga) {
    warning("  Could not determine orientation from barcode detection.",
            " Assuming species in rows.")
  } else {
    message("  Detected: species in rows, samples in columns — correct orientation")
  }

  message("  Matrix after orientation: ", nrow(mat), " species × ", ncol(mat), " samples")

  # Ensure numeric — apply() drops rownames so we save and restore them
  saved_rownames <- rownames(mat)
  mat            <- apply(mat, 2, as.numeric)
  rownames(mat)  <- saved_rownames

  # ── Standardise sample barcodes to 15 characters ────────────────────────────
  colnames(mat) <- sample_barcode(colnames(mat))

  # Remove any duplicate barcodes (keep first occurrence)
  if (any(duplicated(colnames(mat)))) {
    message("  Removing ", sum(duplicated(colnames(mat))), " duplicate sample barcodes")
    mat <- mat[, !duplicated(colnames(mat)), drop = FALSE]
  }

  # ── Filter to primary tumour samples ────────────────────────────────────────
  is_primary <- get_sample_type(colnames(mat)) == "01"
  mat        <- mat[, is_primary, drop = FALSE]
  message("  Primary tumour samples: ", ncol(mat))

  # ── Prevalence filter ────────────────────────────────────────────────────────
  n_species_pre <- nrow(mat)
  mat           <- filter_by_prevalence(mat, min_prev = min_prev)
  message("  Species after prevalence filter (>=", min_prev * 100, "%): ",
          nrow(mat), " (removed ", n_species_pre - nrow(mat), ")")

  # ── Log transformation ───────────────────────────────────────────────────────
  log_mat <- log_transform_abundance(mat)

  # ── Save ─────────────────────────────────────────────────────────────────────
  out <- list(
    raw         = mat,
    abundance   = log_mat,      # use this for all downstream correlations
    species     = rownames(mat),
    n_species   = nrow(mat),
    n_samples   = ncol(mat),
    cancer_type = cancer_type,
    project     = cfg$tcga_project
  )

  rds_out <- here("data/processed",
                  paste0(cfg$tcga_project, "_bacteria_processed.rds"))
  saveRDS(out, rds_out)
  message("  Saved: ", rds_out)

  # ── Report top species by mean abundance ─────────────────────────────────────
  mean_abund <- sort(rowMeans(log_mat), decreasing = TRUE)
  message("\n  Top 10 species by mean log-abundance:")
  print(head(mean_abund, 10))

  return(invisible(out))
}

# ── Run ───────────────────────────────────────────────────────────────────────

# v2.1: re-running original four cancer types on expanded gene panel.
# TCMbio CSVs must already be present in data/raw/ — see instructions above.

preprocess_tcmbio("colon")
preprocess_tcmbio("breast")
preprocess_tcmbio("pancreatic")
preprocess_tcmbio("prostate")

message("\nTCMbio preprocessing complete.")
