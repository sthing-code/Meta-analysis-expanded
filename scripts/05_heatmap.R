# scripts/05_heatmap.R
#
# Generates the main Pearson correlation heatmaps for each cancer type.
# Produces two versions per cancer type:
#   (1) Main figure: |r| > 0.20 filter, significance annotation
#   (2) Supplementary: |r| > 0.15 filter
#
# Output: figures/{cancer}/heatmap_{PROJECT}_pearson_r0.2.pdf  + .png
#         figures/{cancer}/heatmap_{PROJECT}_pearson_r0.15.pdf + .png

# ── Dependencies ──────────────────────────────────────────────────────────────

if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
  BiocManager::install("ComplexHeatmap")
if (!requireNamespace("circlize", quietly = TRUE))
  install.packages("circlize")

library(ComplexHeatmap)
library(circlize)
library(grid)
library(RColorBrewer)
library(here)

source(here("R/utils.R"))
source(here("R/plot_themes.R"))

set_heatmap_options()   # apply global ComplexHeatmap settings

# ── Column (gene) order and annotation ───────────────────────────────────────

# Fixed column order: groups genes by biological function for readability
GENE_ORDER <- c(
  "IGFBP7", "LGALS1",                          # IGFBP7 / Galectin-1
  "VEGFA", "MMP9",                              # Angiogenesis
  "IL6", "STAT3",                               # JAK-STAT
  "CD274",                                      # Immune checkpoint
  "APC", "WNT3A", "TCF7", "LEF1",              # WNT
  "BRCA2",                                      # DNA repair
  "TP53", "ATM", "ATR", "MDM2",                # DNA damage
  "ERBB2", "ERBB3",                             # RTK signalling
  "HIF1A",                                      # Hypoxia
  "SOX2",                                       # Stemness
  "CTSB", "CTSD", "CTSL", "CTSS", "CTSV"       # Cathepsins
)

# Human-readable labels for column annotation bar
GENE_GROUPS_DISPLAY <- c(
  IGFBP7  = "IGFBP7/Galectin",
  LGALS1  = "IGFBP7/Galectin",
  VEGFA   = "Angiogenesis",
  MMP9    = "Angiogenesis",
  IL6     = "JAK-STAT",
  STAT3   = "JAK-STAT",
  CD274   = "Immune checkpoint",
  APC     = "WNT",
  WNT3A   = "WNT extended",
  TCF7    = "WNT extended",
  LEF1    = "WNT extended",
  BRCA2   = "DNA repair",
  TP53    = "DNA damage",
  ATM     = "DNA damage",
  ATR     = "DNA damage",
  MDM2    = "DNA damage",
  ERBB2   = "RTK signalling",
  ERBB3   = "RTK signalling",
  HIF1A   = "Hypoxia",
  SOX2    = "Stemness",
  CTSB    = "Cathepsins",
  CTSD    = "Cathepsins",
  CTSL    = "Cathepsins",
  CTSS    = "Cathepsins",
  CTSV    = "Cathepsins"
)

# ── Heatmap builder ───────────────────────────────────────────────────────────

#' Build and save a ComplexHeatmap correlation heatmap.
#'
#' cor_obj      : one element of the correlations list (e.g. $r020)
#' cancer_type  : "colon" | "breast" | "pancreatic" | "prostate"
#' r_threshold  : label for filename (0.20 or 0.15)
#' q_threshold  : FDR threshold for significance stars (default 0.05)

draw_correlation_heatmap <- function(cor_obj,
                                     cancer_type,
                                     r_threshold  = 0.20,
                                     q_threshold  = 0.05) {

  cfg     <- CANCER_CONFIG[[cancer_type]]
  project <- cfg$tcga_project
  label   <- cfg$label

  r_mat <- cor_obj$r
  q_mat <- cor_obj$q

  # Clean row names: strip s__ prefix, replace underscores with spaces
  # for manuscript-quality species labels
  clean_names <- function(x) {
    x <- sub("^s__", "", x)
    gsub("_", " ", x)
  }
  rownames(r_mat) <- clean_names(rownames(r_mat))
  rownames(q_mat) <- clean_names(rownames(q_mat))
  names(cor_obj$benchmark_flag) <- clean_names(names(cor_obj$benchmark_flag))

  if (nrow(r_mat) == 0) {
    message("  No bacteria pass the r = ", r_threshold, " threshold for ", project,
            " — skipping heatmap.")
    return(invisible(NULL))
  }

  # ── Reorder columns to fixed gene order ──────────────────────────────────────
  genes_present  <- intersect(GENE_ORDER, colnames(r_mat))
  r_mat          <- r_mat[, genes_present, drop = FALSE]
  q_mat          <- q_mat[, genes_present, drop = FALSE]

  message("\n  Building heatmap: ", project, " | r > ", r_threshold)
  message("  Bacteria: ", nrow(r_mat), "  |  Genes: ", ncol(r_mat))

  # ── Significance cell function ────────────────────────────────────────────────
  # Adds a star (*) to cells where q < q_threshold AND the correlation direction
  # matches the expected flag for that gene.

  cell_fun <- function(j, i, x, y, width, height, fill) {
    gene_name <- colnames(r_mat)[j]
    r_val     <- r_mat[i, j]
    q_val     <- q_mat[i, j]
    sig       <- !is.na(q_val) && q_val < q_threshold

    # Positive flag genes: mark strong positive correlations
    if (sig && is_positive_flag_gene(gene_name) && !is.na(r_val) && r_val > r_threshold) {
      grid.text("★", x, y,
                gp = gpar(fontsize = 8, col = FLAG_COLOURS["positive_flag"], fontface = "bold"))

    # Negative flag genes (cathepsins): mark strong negative correlations
    } else if (sig && is_negative_flag_gene(gene_name) && !is.na(r_val) && r_val < -r_threshold) {
      grid.text("★", x, y,
                gp = gpar(fontsize = 8, col = FLAG_COLOURS["negative_flag"], fontface = "bold"))

    # General significance marker for all other genes
    } else if (sig) {
      grid.text("·", x, y,
                gp = gpar(fontsize = 10, col = "grey30"))
    }
  }

  # ── Row (bacterium) annotation: benchmark flag ────────────────────────────────
  benchmark_flags  <- cor_obj$benchmark_flag
  # Align to current row order (may differ after filter)
  flags_aligned    <- benchmark_flags[rownames(r_mat)]

  row_anno <- rowAnnotation(
    Evidence = flags_aligned,
    col = list(
      Evidence = c(
        "literature-supported" = unname(FLAG_COLOURS["literature"]),
        "novel"                = unname(FLAG_COLOURS["novel"])
      )
    ),
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size   = unit(3, "mm"),
    show_legend        = TRUE
  )

  # ── Column (gene) annotation: functional group ────────────────────────────────
  gene_groups <- GENE_GROUPS_DISPLAY[genes_present]

  col_anno <- HeatmapAnnotation(
    Group = gene_groups,
    col   = list(Group = GENE_GROUP_COLS),
    annotation_name_gp = gpar(fontsize = 8),
    simple_anno_size   = unit(3, "mm"),
    show_legend        = TRUE
  )

  # ── Row clustering ────────────────────────────────────────────────────────────
  # Hierarchical clustering on Pearson r values, Ward linkage
  row_dend <- as.dendrogram(
    hclust(dist(r_mat, method = "euclidean"), method = "ward.D2")
  )

  # ── Heatmap title ─────────────────────────────────────────────────────────────
  ht_title <- paste0(
    label,
    "\nPearson r | Bacteria × Gene Expression | |r| > ", r_threshold,
    " | n = ", cor_obj$n_samples, " tumours"
  )

  # ── Draw heatmap ─────────────────────────────────────────────────────────────
  ht <- Heatmap(
    r_mat,
    name                   = "Pearson r",
    col                    = HEATMAP_COL_FUN,
    cluster_rows           = row_dend,
    cluster_columns        = FALSE,       # columns fixed to biological order
    show_row_names         = TRUE,
    show_column_names      = TRUE,
    row_names_gp           = gpar(fontsize = 6.5, fontface = "italic"),
    column_names_gp        = gpar(fontsize = 8, fontface = "italic"),
    column_names_rot       = 45,
    column_names_max_height = unit(4, "cm"),
    row_names_side         = "left",
    cell_fun               = cell_fun,
    top_annotation         = col_anno,
    right_annotation       = row_anno,
    rect_gp                = gpar(col = "white", lwd = 0.5),
    column_title           = ht_title,
    column_title_gp        = gpar(fontsize = 9, fontface = "bold"),
    heatmap_legend_param   = list(
      title     = "Pearson r",
      at        = c(-0.5, -0.25, 0, 0.25, 0.5),
      labels    = c("−0.5", "−0.25", "0", "0.25", "0.5"),
      direction = "vertical"
    ),
    border                 = TRUE,
    border_gp              = gpar(col = "grey60", lwd = 0.5)
  )

  # ── Legend for significance stars ─────────────────────────────────────────────
  star_legend <- Legend(
    labels     = c(
      paste0("★ IGFBP7/Galectin: r > ", r_threshold, ", q < ", q_threshold),
      paste0("★ Cathepsin: r < −", r_threshold, ", q < ", q_threshold),
      paste0("· Other: q < ", q_threshold)
    ),
    type       = "points",
    pch        = c("★", "★", "·"),
    legend_gp  = gpar(col = c(FLAG_COLOURS["positive_flag"],
                               FLAG_COLOURS["negative_flag"],
                               "grey30"),
                      fontsize = 8),
    title      = "Significance",
    title_gp   = gpar(fontsize = 8, fontface = "bold")
  )

  # ── Save to file ──────────────────────────────────────────────────────────────
  r_str    <- gsub("\\.", "", as.character(r_threshold))
  base_fn  <- file.path(here(cfg$figures_dir),
                         paste0("heatmap_", project, "_pearson_r", r_str))

  # Compute figure dimensions: scale with number of bacteria
  fig_h <- max(7, 2.5 + nrow(r_mat) * 0.18)
  fig_w <- 10

  # PDF (vector — for manuscript submission)
  pdf(paste0(base_fn, ".pdf"), width = fig_w, height = fig_h)
  draw(ht, annotation_legend_list = list(star_legend), merge_legend = FALSE)
  dev.off()

  # PNG (raster — for quick review)
  png(paste0(base_fn, ".png"), width = fig_w, height = fig_h,
      units = "in", res = 300)
  draw(ht, annotation_legend_list = list(star_legend), merge_legend = FALSE)
  dev.off()

  message("  Saved: ", base_fn, ".pdf / .png")
  return(invisible(ht))
}

# ── Run ───────────────────────────────────────────────────────────────────────

make_heatmaps_for_cancer <- function(cancer_type) {

  cfg     <- CANCER_CONFIG[[cancer_type]]
  project <- cfg$tcga_project

  rds_in <- here("data/processed", paste0(project, "_correlations.rds"))
  if (!file.exists(rds_in)) {
    stop("Run 04_correlate.R first.\n  Missing: ", rds_in)
  }

  cor_data <- readRDS(rds_in)

  # Use cancer-type-specific thresholds from utils.R
  thresholds <- HEATMAP_THRESHOLDS[[cancer_type]]
  primary_r  <- thresholds$primary
  suppl_r    <- thresholds$supplementary

  # Select the correct pre-filtered correlation object.
  # For most cancer types supplementary < primary (looser threshold).
  # PAAD is the exception: primary = 0.15, supplementary = 0.20 (stricter view).
  primary_cor <- if (primary_r == 0.20) cor_data$r020 else cor_data$r015
  suppl_cor   <- if (suppl_r  == 0.20) cor_data$r020 else cor_data$r015

  message("  Using primary threshold: |r| > ", primary_r,
          " | supplementary: |r| > ", suppl_r)

  # Main figure
  draw_correlation_heatmap(primary_cor, cancer_type, r_threshold = primary_r)

  # Supplementary (only draw if different from primary)
  if (primary_r != suppl_r) {
    draw_correlation_heatmap(suppl_cor, cancer_type, r_threshold = suppl_r)
  }
}

# make_heatmaps_for_cancer("colon")
# make_heatmaps_for_cancer("breast")
# make_heatmaps_for_cancer("pancreatic")
# make_heatmaps_for_cancer("prostate")
make_heatmaps_for_cancer("glioblastoma")

message("\nHeatmap generation complete.")
