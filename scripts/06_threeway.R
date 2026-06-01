# scripts/06_threeway.R
#
# Three-way correlation analysis: bacterium abundance × gene × downstream
# pathway gene. Targets the three axes defined in the study brief:
#
#   IL-6  → STAT3
#   VEGFA → MAPK signalling (MAPK1, MAPK3, KRAS)
#   APC   → WNT/β-catenin (CTNNB1, MYC, CCND1)
#
# For each bacterium that shows |r| > 0.20 with the primary gene (IL6, VEGFA,
# or APC), we test whether that bacterium also correlates with the downstream
# genes. The output is a bubble plot per axis, one per cancer type.
#
# Input:  data/processed/{PROJECT}_correlations.rds
#         data/processed/{PROJECT}_tpm_priority_genes.rds
#         data/processed/{PROJECT}_bacteria_processed.rds
# Output: figures/{cancer}/threeway_{axis}_{PROJECT}.pdf + .png

# ── Dependencies ──────────────────────────────────────────────────────────────

library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(patchwork)
library(here)

source(here("R/utils.R"))
source(here("R/plot_themes.R"))

# ── Axis definitions ──────────────────────────────────────────────────────────

THREEWAY_AXES <- list(

  # ── Original three axes (all cancer types) ───────────────────────────────
  il6_stat3 = list(
    primary    = "IL6",
    downstream = c("STAT3"),
    label      = "IL-6 → STAT3 axis"
  ),
  vegfa_mapk = list(
    primary    = "VEGFA",
    downstream = c("MAPK1", "MAPK3", "KRAS"),
    label      = "VEGFA → MAPK axis"
  ),
  apc_wnt = list(
    primary    = "APC",
    downstream = c("CTNNB1", "MYC", "CCND1"),
    label      = "APC → WNT axis"
  ),

  # ── NEW v2.0 — expanded axes ──────────────────────────────────────────────

  # Glycolytic / metabolic reprogramming axis
  # LDHA is the canonical Warburg endpoint gene; downstream genes cover
  # glucose transport and lactate export.
  metabolic = list(
    primary    = "LDHA",
    downstream = c("SLC2A1", "SLC16A1", "SLC16A3", "PFKP", "TIGAR", "SESN2"),
    label      = "LDHA → metabolic reprogramming axis"
  ),

  # RTK signalling axis
  # ERBB2 as primary; FGFR2/3 as co-receptor downstream comparators.
  rtk_extended = list(
    primary    = "ERBB2",
    downstream = c("FGFR2", "FGFR3"),
    label      = "ERBB2 → RTK signalling axis"
  ),

  # TCA / IDH axis — particularly relevant for GBM (IDH mutation status)
  # IDH1 as primary (most commonly mutated in glioma); IDH2 and GLUD1
  # as downstream TCA cycle nodes.
  tca_idh = list(
    primary    = "IDH1",
    downstream = c("IDH2", "GLUD1"),
    label      = "IDH1 → TCA cycle axis"
  ),

  # NF-kB → IL-6 → STAT3 axis
  # NFKB1 as primary; IL6 and STAT3 as downstream effectors.
  # NF-kB drives IL-6 transcription, which activates STAT3 — a well-established
  # inflammatory signalling cascade in GBM (Brennan et al. 2013; TCGA GBM).
  nfkb_il6_stat3 = list(
    primary    = "NFKB1",
    downstream = c("IL6", "STAT3"),
    label      = "NF-kB → IL-6 → STAT3 axis"
  ),

  # p53 DNA damage response axis
  # TP53 as primary; MDM2 (negative regulator), MDM4 (co-regulator),
  # ATM and ATR (upstream kinases) as downstream context.
  p53_axis = list(
    primary    = "TP53",
    downstream = c("MDM2", "MDM4", "ATM", "ATR"),
    label      = "TP53 → DNA damage response axis"
  ),

  # Stemness axis
  # SOX2 as primary stem cell transcription factor;
  # KLF4 and POU5F1 (OCT4) as co-expressed pluripotency factors.
  stemness_ext = list(
    primary    = "SOX2",
    downstream = c("KLF4", "POU5F1"),
    label      = "SOX2 → stemness axis"
  )
)

# ── Three-way analysis for one axis ──────────────────────────────────────────

run_threeway_axis <- function(cancer_type, axis_name,
                               r_primary    = 0.20,
                               q_threshold  = 0.05) {

  axis <- THREEWAY_AXES[[axis_name]]
  cfg  <- CANCER_CONFIG[[cancer_type]]
  project <- cfg$tcga_project

  # ── Load data ────────────────────────────────────────────────────────────────
  cor_rds  <- here("data/processed", paste0(project, "_correlations.rds"))
  tcga_rds <- here("data/processed", paste0(project, "_tpm_priority_genes.rds"))
  bact_rds <- here("data/processed", paste0(project, "_bacteria_processed.rds"))

  cor_data  <- readRDS(cor_rds)
  tcga_data <- readRDS(tcga_rds)
  bact_data <- readRDS(bact_rds)

  # Full correlation matrix (all bacteria × all pathway genes)
  pathway_cor <- cor_data$pathway

  # Check primary gene exists
  if (!axis$primary %in% colnames(pathway_cor$r)) {
    message("  ", axis$primary, " not in correlation matrix for ", project, " — skipping")
    return(invisible(NULL))
  }

  # ── Select bacteria with |r| > r_primary for the primary gene ────────────────
  r_primary_vec <- pathway_cor$r[, axis$primary]
  selected_bact <- names(r_primary_vec)[abs(r_primary_vec) >= r_primary &
                                           !is.na(r_primary_vec)]

  if (length(selected_bact) == 0) {
    message("  No bacteria meet |r| > ", r_primary, " threshold with ",
            axis$primary, " in ", project)
    return(invisible(NULL))
  }

  message("\n  Three-way: ", axis$label, " | ", project)
  message("  Bacteria selected (|r_primary| > ", r_primary, "): ", length(selected_bact))

  # ── Gather r values for primary + downstream genes ────────────────────────────
  all_genes <- c(axis$primary, axis$downstream)
  genes_avail <- all_genes[all_genes %in% colnames(pathway_cor$r)]

  r_sub <- pathway_cor$r[selected_bact, genes_avail, drop = FALSE]
  q_sub <- pathway_cor$q[selected_bact, genes_avail, drop = FALSE]

  # ── Build tidy data frame for plotting ───────────────────────────────────────
  df <- as.data.frame(r_sub) %>%
    tibble::rownames_to_column("bacterium") %>%
    pivot_longer(-bacterium,
                 names_to  = "gene",
                 values_to = "r") %>%
    left_join(
      as.data.frame(q_sub) %>%
        tibble::rownames_to_column("bacterium") %>%
        pivot_longer(-bacterium, names_to = "gene", values_to = "q"),
      by = c("bacterium", "gene")
    ) %>%
    mutate(
      gene_type   = ifelse(gene == axis$primary, "Primary", "Downstream"),
      sig         = !is.na(q) & q < q_threshold,
      evidence    = sapply(bacterium, function(b) flag_bacterium(b, cancer_type)),
      label_bact  = ifelse(abs(r_primary) >= 0.25 & abs(r) >= 0.20 & sig,
                           sub("^s__", "", gsub("_", " ", bacterium)),
                           NA_character_)
    )

  # ── Bubble plot ───────────────────────────────────────────────────────────────
  # x = r with primary gene, y = r with each downstream gene, size = |r|
  # One panel per downstream gene, plus the primary gene on x-axis

  primary_r  <- df %>% filter(gene == axis$primary) %>%
    select(bacterium, r_primary = r, evidence)

  downstream_df <- df %>%
    filter(gene != axis$primary) %>%
    left_join(primary_r, by = c("bacterium", "evidence")) %>%
    mutate(
      bacterium_clean = sub("^s__", "", gsub("_", " ", bacterium))
    )

  if (nrow(downstream_df) == 0 || length(unique(downstream_df$gene)) == 0) {
    message("  No downstream gene data available for ", axis$label)
    return(invisible(NULL))
  }

  p <- ggplot(downstream_df,
              aes(x = r_primary, y = r, size = abs(r),
                  colour = evidence)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.4) +
    geom_hline(yintercept =  c(-r_primary, r_primary),
               linetype = "dotted", colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = c(-r_primary, r_primary),
               linetype = "dotted", colour = "grey80", linewidth = 0.3) +
    geom_point(alpha = 0.75) +
    geom_text_repel(
      aes(label = label_bact),
      size = 2.5,
      max.overlaps = 15,
      segment.colour = "grey50",
      segment.size   = 0.3,
      fontface       = "italic"
    ) +
    scale_size_continuous(range = c(1.5, 6), guide = "none") +
    scale_colour_manual(
      values = c("literature-supported" = unname(FLAG_COLOURS["literature"]),
                 "novel"                = unname(FLAG_COLOURS["novel"])),
      name   = "Evidence",
      labels = c("literature-supported" = "Literature-supported",
                 "novel"                = "Novel")
    ) +
    guides(colour = guide_legend(override.aes = list(size = 4))) +
    facet_wrap(~ gene, scales = "free_y") +
    labs(
      x        = paste0("r (bacterium × ", axis$primary, ")"),
      y        = "r (bacterium × downstream gene)",
      title    = paste0(cfg$label, " — ", axis$label),
      subtitle = paste0("Bacteria with |r| > ", r_primary, " with ", axis$primary,
                        " | Stars = q < ", q_threshold,
                        " | n = ", cor_data$n_samples, " tumours"),
      caption  = paste0("★ = significant (BH q < ", q_threshold, "). ",
                        "Dashed lines at r = ±", r_primary, ".")
    ) +
    theme_manuscript()

  # ── Save ──────────────────────────────────────────────────────────────────────
  base_fn <- file.path(here(cfg$figures_dir),
                        paste0("threeway_", axis_name, "_", project))

  fig_h <- 5
  fig_w <- max(8, 4 * length(unique(downstream_df$gene)))

  ggsave(paste0(base_fn, ".pdf"), p, width = fig_w, height = fig_h)
  ggsave(paste0(base_fn, ".png"), p, width = fig_w, height = fig_h, dpi = 300)

  message("  Saved: ", base_fn, ".pdf / .png")
  return(invisible(p))
}

# ── Run ───────────────────────────────────────────────────────────────────────

run_threeway_for_cancer <- function(cancer_type) {
  # Use cancer-type-specific primary threshold from utils.R
  threshold <- HEATMAP_THRESHOLDS[[cancer_type]]$primary
  message("  Three-way threshold for ", cancer_type, ": |r| > ", threshold)
  for (axis_name in names(THREEWAY_AXES)) {
    tryCatch(
      run_threeway_axis(cancer_type, axis_name, r_primary = threshold),
      error = function(e) message("  ERROR in ", axis_name, ": ", e$message)
    )
  }
}

# ── Run ───────────────────────────────────────────────────────────────────────

# v2.1: re-running original four cancer types on expanded gene panel.

run_threeway_for_cancer("colon")
run_threeway_for_cancer("breast")
run_threeway_for_cancer("pancreatic")
run_threeway_for_cancer("prostate")

message("\nThree-way correlation analysis complete.")
