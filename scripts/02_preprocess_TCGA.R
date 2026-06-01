# scripts/02_preprocess_TCGA.R
#
# Takes the raw STAR counts from the downloaded SummarizedExperiment objects
# and produces clean, TPM-normalised expression matrices for each cancer type.
#
# Input:  data/processed/{PROJECT}_SE_priority_genes.rds
# Output: data/processed/{PROJECT}_tpm_priority_genes.rds
#         (named list: $tpm = genes × samples TPM matrix,
#                       $meta = sample metadata data.frame)

# ── Dependencies ──────────────────────────────────────────────────────────────

library(SummarizedExperiment)
library(here)
library(dplyr)

source(here("R/utils.R"))

# ── TPM computation from STAR counts ─────────────────────────────────────────
#
# STAR output in TCGA provides:
#  - unstranded counts (assay "unstranded")
#  - gene lengths can be derived from the rowRanges width
#
# TPM formula:
#   RPK_i = counts_i / (gene_length_i / 1000)
#   TPM_i = RPK_i / (sum(RPK) / 1e6)

compute_tpm_from_se <- function(se) {

  # STAR counts — use unstranded counts
  # Assay names in GDC STAR output: "unstranded", "stranded_first",
  #   "stranded_second", "tpm_unstrand" (some datasets include pre-computed TPM)

  assay_names <- assayNames(se)
  message("  Available assays: ", paste(assay_names, collapse = ", "))

  # Use pre-computed TPM if available (preferred — avoids length estimation)
  if ("tpm_unstrand" %in% assay_names) {
    message("  Using pre-computed TPM (tpm_unstrand assay)")
    tpm_mat <- assay(se, "tpm_unstrand")
    rownames(tpm_mat) <- rowData(se)$gene_name
    return(tpm_mat)
  }

  # Otherwise compute from unstranded counts + gene lengths
  message("  Computing TPM from unstranded counts + gene lengths")
  counts_mat  <- assay(se, "unstranded")
  rownames(counts_mat) <- rowData(se)$gene_name

  # Gene length: sum of exon widths from rowRanges
  # For genes with multiple rows (shouldn't happen post-GDCprepare), use first
  gene_lengths <- width(rowRanges(se))
  names(gene_lengths) <- rowData(se)$gene_name

  # Remove genes with zero length
  valid <- gene_lengths > 0
  counts_mat   <- counts_mat[valid, , drop = FALSE]
  gene_lengths <- gene_lengths[valid]

  tpm_mat <- counts_to_tpm(counts_mat, gene_lengths)
  return(tpm_mat)
}

# ── Preprocess one project ────────────────────────────────────────────────────

preprocess_project <- function(project) {

  message("\n========================================")
  message("  Preprocessing: ", project)
  message("========================================\n")

  # Load SE
  rds_in <- here("data/processed", paste0(project, "_SE_priority_genes.rds"))
  if (!file.exists(rds_in)) {
    stop("SE file not found — run 01_download_TCGA.R first.\n  Expected: ", rds_in)
  }
  se <- readRDS(rds_in)
  message("  Loaded SE: ", nrow(se), " genes × ", ncol(se), " samples")

  # Compute TPM
  tpm_mat <- compute_tpm_from_se(se)

  # Deduplicate gene symbols (keep row with highest mean expression if duplicates)
  if (any(duplicated(rownames(tpm_mat)))) {
    message("  Deduplicating gene symbols...")
    row_means      <- rowMeans(tpm_mat, na.rm = TRUE)
    tpm_mat        <- tpm_mat[order(row_means, decreasing = TRUE), ]
    tpm_mat        <- tpm_mat[!duplicated(rownames(tpm_mat)), ]
    message("  After dedup: ", nrow(tpm_mat), " genes")
  }

  # Keep only priority + pathway genes that are present
  genes_present <- rownames(tpm_mat)[rownames(tpm_mat) %in%
                                       unique(c(ALL_PRIORITY_GENES, ALL_PATHWAY_GENES))]
  tpm_mat <- tpm_mat[genes_present, , drop = FALSE]
  message("  Priority genes in TPM matrix: ", nrow(tpm_mat))

  # Standardise sample IDs to 15-character TCGA barcodes
  # (patient + sample type, e.g. "TCGA-AA-3977-01")
  colnames(tpm_mat) <- sample_barcode(colnames(tpm_mat))

  # Sample metadata
  meta <- as.data.frame(colData(se)) %>%
    mutate(
      sample_id_15 = sample_barcode(rownames(.)),
      patient_id   = shorten_barcode(rownames(.))
    )

  # Log-transform: log2(TPM + 1) for downstream correlation
  # (Pearson on raw TPM is dominated by highly expressed outliers)
  log_tpm <- log2(tpm_mat + 1)
  message("  Log2(TPM+1) transformation applied")

  # Save
  out <- list(
    tpm         = tpm_mat,
    log_tpm     = log_tpm,
    meta        = meta,
    project     = project,
    n_samples   = ncol(tpm_mat),
    genes_found = rownames(tpm_mat)
  )

  rds_out <- here("data/processed", paste0(project, "_tpm_priority_genes.rds"))
  saveRDS(out, rds_out)
  message("  Saved: ", rds_out)

  return(invisible(out))
}

# ── Run ───────────────────────────────────────────────────────────────────────

# v2.1: re-running original four cancer types on expanded gene panel.

projects_to_run <- c(
  "TCGA-COAD",   # colon adenocarcinoma
  "TCGA-BRCA",   # breast invasive carcinoma
  "TCGA-PAAD",   # pancreatic adenocarcinoma
  "TCGA-PRAD"    # prostate adenocarcinoma
)

for (proj in projects_to_run) {
  preprocess_project(proj)
}

message("\nPreprocessing complete.")
