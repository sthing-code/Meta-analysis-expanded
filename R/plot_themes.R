# R/plot_themes.R
# Consistent visual theme for all figures in this project.
# Source alongside R/utils.R at the top of plotting scripts.

library(ggplot2)
library(RColorBrewer)
library(circlize)        # for ComplexHeatmap colour scale

# ── Colour palettes ───────────────────────────────────────────────────────────

# Diverging palette for correlation heatmap (blue = negative, red = positive)
# Centred at 0, clipped at ±0.5
HEATMAP_COL_FUN <- colorRamp2(
  breaks = c(-0.5, -0.25, 0, 0.25, 0.5),
  colors = c("#2166AC", "#92C5DE", "white", "#F4A582", "#D6604D")
)

# Annotation colours
FLAG_COLOURS <- c(
  "positive_flag"   = "#C0392B",    # red star — strong positive r for IGFBP7/Galectin-1
  "negative_flag"   = "#2471A3",    # blue star — strong negative r for cathepsins
  "literature"      = "#1A936F",    # green — benchmark bacterium
  "novel"           = "#8E44AD"     # purple — novel observation
)

# Column group colours for gene annotation bar
GENE_GROUP_COLS <- c(
  "IGFBP7/Galectin"  = "#E74C3C",
  "Angiogenesis"     = "#E67E22",
  "JAK-STAT"         = "#8E44AD",
  "WNT"              = "#1ABC9C",
  "DNA repair"       = "#2C3E50",
  "Cathepsins"       = "#2980B9"
)

# ── ggplot2 base theme ────────────────────────────────────────────────────────

theme_manuscript <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      text               = element_text(family = "sans", colour = "black"),
      axis.text          = element_text(size = base_size - 1, colour = "black"),
      axis.title         = element_text(size = base_size,     colour = "black"),
      legend.text        = element_text(size = base_size - 2),
      legend.title       = element_text(size = base_size - 1, face = "bold"),
      plot.title         = element_text(size = base_size + 1, face = "bold"),
      plot.subtitle      = element_text(size = base_size - 1, colour = "grey40"),
      strip.text         = element_text(size = base_size - 1, face = "bold"),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.major.x = element_blank()
    )
}

# Set globally
theme_set(theme_manuscript())

# ── ComplexHeatmap global options ─────────────────────────────────────────────

# Call this once before drawing any heatmap
set_heatmap_options <- function() {
  ht_opt(
    heatmap_row_names_gp = gpar(fontsize = 7),
    heatmap_column_names_gp = gpar(fontsize = 8, fontface = "italic"),
    legend_title_gp = gpar(fontsize = 8, fontface = "bold"),
    legend_labels_gp = gpar(fontsize = 7)
  )
}
