# scripts/01_download_TCGA.R
#
# Downloads TCGA RNA-seq data (STAR — Counts) for each cancer type in this
# project using TCGAbiolinks. Run once; downloads are cached in data/raw/.
#
# Runtime: 15–40 minutes depending on connection speed and GDC server load.
# Output:  GDC query/download cache in data/raw/GDCdata/
#          Processed SummarizedExperiment RDS in data/processed/

# ── Dependencies ──────────────────────────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_pkgs <- c("TCGAbiolinks", "SummarizedExperiment")
missing_pkgs  <- required_pkgs[!sapply(required_pkgs, requireNamespace,
                                        quietly = TRUE)]
if (length(missing_pkgs) > 0)
  BiocManager::install(missing_pkgs)

library(TCGAbiolinks)
library(SummarizedExperiment)
library(here)

source(here("R/utils.R"))

# ── Config ────────────────────────────────────────────────────────────────────

# Set download cache inside the project data directory.
# GDCdownload will create data/raw/GDCdata/ automatically.
GDC_DIR <- here("data/raw/GDCdata")

# Genes to retain after download. We keep ALL priority + pathway genes.
GENES_TO_KEEP <- unique(c(ALL_PRIORITY_GENES, ALL_PATHWAY_GENES))

# ── Download function ─────────────────────────────────────────────────────────

#' Query, download and preprocess one TCGA project.
#'
#' project : e.g. "TCGA-COAD"
#' Returns a SummarizedExperiment with primary-tumour samples and all genes.
#' Saves a filtered RDS (priority genes only) to data/processed/.
download_tcga_project <- function(project) {

  message("\n========================================")
  message("  Processing: ", project)
  message("========================================\n")

  # ── 1. Query GDC ────────────────────────────────────────────────────────────
  query <- GDCquery(
    project           = project,
    data.category     = "Transcriptome Profiling",
    data.type         = "Gene Expression Quantification",
    workflow.type     = "STAR - Counts",        # GDC harmonised pipeline, hg38
    experimental.strategy = "RNA-Seq"
  )

  # ── 2. Download (skips files already on disk) ────────────────────────────────
  GDCdownload(
    query     = query,
    method    = "api",
    files.per.chunk = 10,
    directory = GDC_DIR
  )

  # ── 3. Prepare SummarizedExperiment ─────────────────────────────────────────
  se <- GDCprepare(
    query     = query,
    directory = GDC_DIR,
    save      = FALSE
  )

  # ── 4. Filter to primary tumour samples ─────────────────────────────────────
  # TCGA sample type "01" = Primary Solid Tumour
  sample_types <- colData(se)$sample_type
  is_primary   <- grepl("Primary", sample_types, ignore.case = TRUE)
  se           <- se[, is_primary]

  message("  Primary tumour samples retained: ", ncol(se))

  # ── 5. Save full SE (before gene filtering) ──────────────────────────────────
  rds_full <- here("data/processed", paste0(project, "_SE_full.rds"))
  saveRDS(se, rds_full)
  message("  Full SE saved: ", rds_full)

  # ── 6. Extract and save the priority-gene subset ─────────────────────────────
  # rowData(se)$gene_name contains HGNC symbols (from STAR/GDC)
  gene_names <- rowData(se)$gene_name
  keep_genes <- gene_names %in% GENES_TO_KEEP
  se_sub     <- se[keep_genes, ]

  message("  Priority genes found: ",
          sum(GENES_TO_KEEP %in% gene_names), " / ", length(GENES_TO_KEEP))

  missing_genes <- GENES_TO_KEEP[!GENES_TO_KEEP %in% gene_names]
  if (length(missing_genes) > 0)
    message("  Genes not found in dataset: ", paste(missing_genes, collapse = ", "))

  rds_sub <- here("data/processed", paste0(project, "_SE_priority_genes.rds"))
  saveRDS(se_sub, rds_sub)
  message("  Priority-gene SE saved: ", rds_sub)

  return(invisible(se_sub))
}

# ── Run downloads ─────────────────────────────────────────────────────────────

# v2.1: re-running original four cancer types on expanded gene panel.
# GBM was already downloaded in v2.0 — omit to avoid redundant re-download.
# To re-download GBM, add "TCGA-GBM" to the vector below.

projects_to_run <- c(
  "TCGA-COAD",   # colon adenocarcinoma
  "TCGA-BRCA",   # breast invasive carcinoma
  "TCGA-PAAD",   # pancreatic adenocarcinoma
  "TCGA-PRAD",   # prostate adenocarcinoma
  # "TCGA-GBM"     # glioblastoma multiforme — already downloaded in v2.0; GDCdownload skips cached files
)

for (proj in projects_to_run) {
  download_tcga_project(proj)
}

message("\nAll downloads complete.")
