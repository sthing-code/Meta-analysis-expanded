# R/utils.R
# Shared helper functions for the microbiome-cancer correlation pipeline.
# Source this file at the top of every analysis script.
#
# VERSION HISTORY
# v1.0  — Original four-cancer analysis (COAD, BRCA, PAAD, PRAD)
# v2.0  — Expanded gene panel; GBM added; new gene groups; GBM benchmarks added
#          GBM heatmap thresholds: TBD — to be set after script 04 runs on TCGA-GBM

# ── Gene lists ────────────────────────────────────────────────────────────────

# Priority genes, grouped by biological function.
# These are the column variables in every heatmap.
#
# NEW IN v2.0 (expanded panel, confirmed by supervisor):
#   immune_checkpoint : CD274 (PD-L1)
#   hypoxia           : HIF1A
#   rtk_signalling    : ERBB2, ERBB3
#   dna_damage        : TP53, ATM, ATR, MDM2
#   wnt_extended      : WNT3A, TCF7, LEF1
#   stemness          : SOX2
#
# Original groups retained unchanged.

PRIORITY_GENES <- list(
  igfbp_galectin    = c("IGFBP7", "LGALS1"),
  angiogenesis      = c("VEGFA", "MMP9"),
  jak_stat          = c("IL6", "STAT3"),
  wnt               = c("APC"),
  dna_repair        = c("BRCA2"),
  cathepsins        = c("CTSB", "CTSD", "CTSL", "CTSS", "CTSV"),
  # ── NEW v2.0 ──────────────────────────────────────────────────────────────
  immune_checkpoint = c("CD274"),
  hypoxia           = c("HIF1A"),
  rtk_signalling    = c("ERBB2", "ERBB3"),
  dna_damage        = c("TP53", "ATM", "ATR", "MDM2"),
  wnt_extended      = c("WNT3A", "TCF7", "LEF1"),
  stemness          = c("SOX2")
)

# Flat vector of all priority genes (for TCGA query)
ALL_PRIORITY_GENES <- unlist(PRIORITY_GENES, use.names = FALSE)

# Downstream pathway genes for three-way correlation analysis.
#
# NEW IN v2.0 (metabolic, RTK, stemness, and other axes):
#   metabolic    : LDHA, SLC2A1 (GLUT1), SLC16A1 (MCT1), SLC16A3 (MCT4),
#                  PFKP, TIGAR, SESN2
#   rtk_extended : FGFR2, FGFR3
#   tca_idh      : GLUD1, IDH1, IDH2
#   nfkb         : NFKB1
#   p53_axis     : MDM4
#   stemness_ext : KLF4, POU5F1 (OCT4)
#   glioma_marker: GFAP, RBFOX3
#
# MYC, KRAS, CTNNB1 already present — no changes required.

PATHWAY_GENES <- list(
  il6_stat3    = c("IL6",    "STAT3"),
  vegfa_mapk   = c("VEGFA",  "MAPK1", "MAPK3", "KRAS"),
  apc_wnt      = c("APC",    "CTNNB1", "MYC",  "CCND1"),
  # ── NEW v2.0 ──────────────────────────────────────────────────────────────
  metabolic    = c("LDHA",  "SLC2A1", "SLC16A1", "SLC16A3", "PFKP",
                   "TIGAR", "SESN2"),
  rtk_extended = c("FGFR2", "FGFR3"),
  tca_idh      = c("GLUD1", "IDH1",  "IDH2"),
  nfkb         = c("NFKB1"),
  p53_axis     = c("MDM4"),
  stemness_ext = c("KLF4",  "POU5F1"),
  glioma_marker = c("GFAP", "RBFOX3")
)

ALL_PATHWAY_GENES <- unique(unlist(PATHWAY_GENES, use.names = FALSE))

# ── TCGA cancer-type config ───────────────────────────────────────────────────

CANCER_CONFIG <- list(
  colon = list(
    tcga_project = "TCGA-COAD",
    tcmbio_file  = "data/raw/TCMbio_COAD_bacteria.csv",
    figures_dir  = "figures/colon",
    label        = "Colon adenocarcinoma (COAD)"
  ),
  breast = list(
    tcga_project = "TCGA-BRCA",
    tcmbio_file  = "data/raw/TCMbio_BRCA_bacteria.csv",
    figures_dir  = "figures/breast",
    label        = "Breast invasive carcinoma (BRCA)"
  ),
  pancreatic = list(
    tcga_project = "TCGA-PAAD",
    tcmbio_file  = "data/raw/TCMbio_PAAD_bacteria.csv",
    figures_dir  = "figures/pancreatic",
    label        = "Pancreatic adenocarcinoma (PAAD)"
  ),
  prostate = list(
    tcga_project = "TCGA-PRAD",
    tcmbio_file  = "data/raw/TCMbio_PRAD_bacteria.csv",
    figures_dir  = "figures/prostate",
    label        = "Prostate adenocarcinoma (PRAD)"
  ),
  # ── NEW v2.0 ──────────────────────────────────────────────────────────────
  glioblastoma = list(
    tcga_project = "TCGA-GBM",
    tcmbio_file  = "data/raw/TCMbio_GBM_bacteria.csv",
    figures_dir  = "figures/glioblastoma",
    label        = "Glioblastoma multiforme (GBM)"
  )
)

# ── Cancer-type-specific heatmap thresholds ──────────────────────────────────
# COAD: strong correlations (104 bacteria at |r| > 0.20); strict primary threshold.
# BRCA: weaker correlations overall (9 bacteria at |r| > 0.20); |r| > 0.15 primary.
# PAAD: small dataset (n=178, 23 species post-prevalence filter); |r| > 0.15 primary
#       consistent with BRCA. Supplementary |r| > 0.20 (17 species).
# PRAD: n=497, 32 species post-prevalence filter. |r| > 0.15 primary (22 species).
#       Supplementary |r| > 0.20 (18 species). Confirmed after script 04.
# GBM:  TBD — thresholds to be determined after script 04 runs on TCGA-GBM.
#       Placeholder values set to NA; update once script 04 output is reviewed.

HEATMAP_THRESHOLDS <- list(
  colon        = list(primary = 0.20, supplementary = 0.15),
  breast       = list(primary = 0.15, supplementary = 0.10),
  pancreatic   = list(primary = 0.15, supplementary = 0.20),
  prostate     = list(primary = 0.15, supplementary = 0.20),
  # ── NEW v2.0 — PLACEHOLDER ────────────────────────────────────────────────
  glioblastoma = list(primary = NA_real_,  supplementary = NA_real_)
)

# ── Literature-supported bacteria per cancer type ─────────────────────────────
# Used to colour-flag correlations as literature-supported (green) vs novel (grey).
# Sources cited inline; confidence levels in comments (HIGH / MODERATE).

BENCHMARK_BACTERIA <- list(

  colon = c(
    # ── Enriched in CRC (pro-tumorigenic) ──────────────────────────────────
    # Kostic et al. Genome Res 2012; Rubinstein et al. Cell Host Microbe 2013
    "Fusobacterium nucleatum",
    # Yu et al. Gut 2017; Proctor et al. NAR Cancer 2021
    "Parvimonas micra",
    "Peptostreptococcus stomatis",
    "Solobacterium moorei",
    # Dalal et al. 2021; Wu et al. Nat Med 2004
    "Bacteroides fragilis",
    # Dalal et al. 2021
    "Streptococcus gallolyticus",
    "Enterococcus faecalis",
    "Peptostreptococcus anaerobius",
    "Helicobacter pylori",
    "Clostridium septicum",
    "Salmonella enterica",
    # Bautista et al. 2026; Sheng et al. 2024
    "Porphyromonas gingivalis",
    # ── Depleted in CRC (protective commensals) ────────────────────────────
    # Bautista et al. 2026; Cao et al. 2025
    "Faecalibacterium prausnitzii",
    # Bautista et al. 2026
    "Roseburia intestinalis",
    "Roseburia hominis"
  ),

  breast = c(
    # ── Enriched in breast cancer / estrobolome-active ─────────────────────
    # Larnder et al. 2025; Mahno et al. 2024; Sheng et al. 2024
    "Escherichia coli",
    # Larnder et al. 2025 — broadest estrobolome enzyme activity
    "Bacteroides fragilis",
    # Sheng et al. 2024
    "Fusobacterium nucleatum",
    "Lactobacillus iners",
    "Staphylococcus aureus",       # contamination caveat — flag in figure
    "Staphylococcus epidermidis",  # contamination caveat — flag in figure
    # Mahno et al. 2024; Larnder et al. 2025 (conflicting evidence)
    "Ruminococcus gnavus",
    # Larnder et al. 2025
    "Klebsiella pneumoniae",
    "Citrobacter freundii",
    "Enterobacter hormaechei",
    # Sheng et al. 2024
    "Sphingomonas panacisoli",
    "Sphingomonas paucimobilis",
    # ── Depleted in breast cancer / protective ──────────────────────────────
    # Ruo et al. 2021; Larnder et al. 2025
    "Faecalibacterium prausnitzii",
    # Larnder et al. 2025
    "Roseburia inulinivorans",
    "Roseburia hominis",
    "Roseburia intestinalis",
    "Bifidobacterium longum",
    "Bifidobacterium adolescentis",
    "Collinsella aerofaciens",
    # ── TMAO-producing Clostridiales — antitumour in TNBC ──────────────────
    # Wang et al. 2022 Cell Metabolism — TMAO → pyroptosis → CD8+ T cells
    "Blautia obeum",
    "Dorea longicatena",
    "Ruminococcus torques"
  ),

  pancreatic = c(
    # ── Pro-tumorigenic: oral-to-pancreas translocation ───────────────────
    # Wang et al. 2024 (K-rasG12D mouse model); Sheng et al. 2024
    "Porphyromonas gingivalis",
    # Sheng et al. 2024 — oral microbiota enriched in early PDAC precursors
    "Fusobacterium nucleatum",
    "Prevotella oris",
    "Prevotella melaninogenica",
    "Veillonella parvula",
    # ── Pro-tumorigenic: gemcitabine inactivation ─────────────────────────
    # Geller et al. Science 2017, cited in Sheng et al. 2024
    "Mycoplasma hyorhinis",
    # Bautista et al. 2026 — Gammaproteobacteria inactivate gemcitabine
    # via cytidine deaminase; class-level benchmark
    "Pseudomonas aeruginosa",     # Gammaproteobacteria representative
    "Klebsiella pneumoniae",      # Gammaproteobacteria representative
    # ── TMAO-producing bacteria — antitumour in PAAD ──────────────────────
    # Liu et al. 2023; W. Zhang et al. 2024
    "Blautia obeum",
    "Dorea longicatena",
    "Faecalibacterium prausnitzii",
    # ── Protective ────────────────────────────────────────────────────────
    # Wang et al. 2024
    "Lactobacillus acidophilus",
    "Lactobacillus rhamnosus"
  ),

  prostate = c(
    # ── Pro-tumorigenic: intratumoral PRAD pathogens ──────────────────────
    # HIGH confidence:
    # Ashida et al. 2024 — detected in 19/20 PC tissues by IHC; invades prostate
    # epithelial cells; significantly downregulates BRCA2; induces BRCAness.
    # Davidsson et al. 2021 — stimulates M2 polarisation; upregulates PD-L1,
    # CCL17, CCL18; associated with increased Tregs in stroma (P=0.0004) and
    # epithelia (P=0.046).
    "Cutibacterium acnes",
    # MODERATE confidence:
    # Ashida et al. 2024 — detected by PCR in 20/20 PC tissues. Detection only.
    "Moraxella osloensis",
    # Ashida et al. 2024 — detected by PCR in 16/20 PC tissues. Detection only.
    "Micrococcus luteus"
    # Note: Uncultured Chroococcidiopsis excluded — absent from TCMbio namespace.
  ),

  # ── NEW v2.0 — GBM benchmark ────────────────────────────────────────────────
  # Only two confirmed entries. C. acnes explicitly excluded — no direct GBM
  # tissue evidence in any reviewed paper (see handoff Section 6).
  glioblastoma = c(
    # HIGH confidence:
    # Li et al. 2024 (mSystems, DOI: 10.1128/msystems.00457-24)
    # Genus Fusobacterium significantly enriched in glioma tissue vs adjacent
    # normal brain (n=50 patients); bacterial RNA and LPS confirmed by IHC,
    # immunofluorescence, and FISH; F. nucleatum promotes glioma cell proliferation
    # and upregulates CCL2, CXCL1, CXCL2 in xenograft and organoid models.
    # Note: study includes mixed glioma grades, not exclusively GBM.
    "Fusobacterium nucleatum",
    # MODERATE confidence:
    # Mehelleb et al. 2025 (Int J Mol Sci, n=9 GBM patients, 16S rRNA fresh tissue)
    # Sphingomonas ranked as second most abundant genus (12.90% relative abundance).
    # No genus reached statistical significance in differential-abundance testing.
    # Treat as abundance ranking only — not validated enrichment.
    "Sphingomonas melonis"
    # Note: Burkholderia-Caballeronia-Paraburkholderia (27.43%), Helicobacter (4.16%),
    # and Leifsonia (8.01%) also reported by Mehelleb et al. 2025 but excluded here
    # as they are unlikely to match TCMbio species namespace at species level.
    # Revisit if these genera appear post-prevalence filter in script 03.
  )
)

# ── Barcode utilities ─────────────────────────────────────────────────────────

#' Truncate a TCGA barcode to 12-character patient ID
#' e.g. "TCGA-AA-3977-01A-01R-1635-07" → "TCGA-AA-3977"
shorten_barcode <- function(x) substr(x, 1, 12)

#' Truncate to 15 characters (patient + sample type)
#' e.g. "TCGA-AA-3977-01A-01R-1635-07" → "TCGA-AA-3977-01"
sample_barcode <- function(x) substr(x, 1, 15)

#' Extract sample type code from TCGA barcode (positions 14–15)
#' "01" = Primary tumour, "11" = solid normal tissue
get_sample_type <- function(x) substr(x, 14, 15)

#' Filter a vector of barcodes to primary tumour samples only
filter_primary_tumour <- function(barcodes) {
  barcodes[get_sample_type(barcodes) == "01"]
}

# ── TPM computation ───────────────────────────────────────────────────────────

#' Compute TPM from a raw counts matrix and a gene-length vector.
#' counts : genes × samples integer matrix (STAR unstranded counts)
#' gene_lengths : named numeric vector (gene name → length in bp)
#' Returns a genes × samples TPM matrix.
counts_to_tpm <- function(counts, gene_lengths) {
  stopifnot(all(rownames(counts) %in% names(gene_lengths)))
  lengths  <- gene_lengths[rownames(counts)]
  rpk      <- counts / (lengths / 1e3)
  col_sums <- colSums(rpk, na.rm = TRUE)
  tpm      <- t(t(rpk) / (col_sums / 1e6))
  return(tpm)
}

# ── Bacterial abundance preprocessing ────────────────────────────────────────

#' Log-transform bacterial abundance: log10(x + 1)
log_transform_abundance <- function(mat) {
  log10(mat + 1)
}

#' Filter bacteria by prevalence across samples.
#' Keeps only species present (abundance > 0) in >= min_prev fraction of samples.
filter_by_prevalence <- function(mat, min_prev = 0.10) {
  prev <- rowMeans(mat > 0)
  mat[prev >= min_prev, , drop = FALSE]
}

# ── Correlation utilities ─────────────────────────────────────────────────────

#' Compute Pearson r and two-sided p-value between two numeric vectors.
#' Returns NA silently if fewer than 10 finite paired observations exist.
pearson_pair <- function(x, y) {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 10) return(c(r = NA_real_, p = NA_real_))
  ct <- cor.test(x[valid], y[valid], method = "pearson")
  c(r = unname(ct$estimate), p = ct$p.value)
}

#' Compute full Pearson correlation matrix (bacteria × genes) with p-values.
#'
#' bact_mat : bacteria × samples (rows = species, cols = samples)
#' expr_mat : genes × samples   (rows = genes, cols = samples)
#'
#' Returns a list:
#'   $r  — bacteria × genes correlation matrix
#'   $p  — bacteria × genes raw p-value matrix
#'   $q  — bacteria × genes BH-adjusted q-value matrix
compute_pearson_matrix <- function(bact_mat, expr_mat) {

  common_samples <- intersect(colnames(bact_mat), colnames(expr_mat))
  if (length(common_samples) < 10) {
    stop("Fewer than 10 matched samples — check barcode formatting.")
  }
  message("  Matched samples: ", length(common_samples))

  bact_mat <- bact_mat[, common_samples, drop = FALSE]
  expr_mat <- expr_mat[, common_samples, drop = FALSE]

  n_bact <- nrow(bact_mat)
  n_gene <- nrow(expr_mat)

  r_mat <- matrix(NA_real_, nrow = n_bact, ncol = n_gene,
                  dimnames = list(rownames(bact_mat), rownames(expr_mat)))
  p_mat <- r_mat

  for (i in seq_len(n_bact)) {
    for (j in seq_len(n_gene)) {
      res          <- pearson_pair(bact_mat[i, ], expr_mat[j, ])
      r_mat[i, j] <- res["r"]
      p_mat[i, j] <- res["p"]
    }
  }

  q_mat <- matrix(p.adjust(as.vector(p_mat), method = "BH"),
                  nrow = n_bact, ncol = n_gene,
                  dimnames = dimnames(p_mat))

  list(r = r_mat, p = p_mat, q = q_mat,
       n_samples = length(common_samples))
}

#' Filter correlation results to retain only |r| >= threshold.
#' Returns the same list structure with rows (bacteria) that have at least
#' one gene meeting the threshold.
filter_by_r <- function(cor_list, r_threshold = 0.20) {
  keep <- rowSums(abs(cor_list$r) >= r_threshold, na.rm = TRUE) > 0
  list(
    r = cor_list$r[keep, , drop = FALSE],
    p = cor_list$p[keep, , drop = FALSE],
    q = cor_list$q[keep, , drop = FALSE],
    n_samples = cor_list$n_samples
  )
}

# ── Flagging helpers ──────────────────────────────────────────────────────────

#' Given a gene name, return whether it should be flagged for strong POSITIVE r.
#' Covers original and expanded panel positive-axis genes.
is_positive_flag_gene <- function(gene) {
  gene %in% c("IGFBP7", "LGALS1", "CD274", "HIF1A",
              "ERBB2",  "ERBB3",  "SOX2",  "WNT3A")
}

#' Given a gene name, return whether it should be flagged for strong NEGATIVE r.
#' Covers cathepsins and expanded DNA-damage / tumour-suppressor genes.
is_negative_flag_gene <- function(gene) {
  gene %in% c(PRIORITY_GENES$cathepsins,
              PRIORITY_GENES$dna_damage,   # TP53, ATM, ATR, MDM2
              "IDH1", "IDH2")              # glioma tumour suppressors
}

#' Given a bacterium name and a cancer type, classify as benchmark or novel.
#' Handles TCMbio-format names (s__Genus_species) automatically.
flag_bacterium <- function(bacterium, cancer_type = "colon") {
  benchmarks <- BENCHMARK_BACTERIA[[cancer_type]]
  if (is.null(benchmarks)) return("novel")
  plain <- normalise_species_name(bacterium)
  genus <- strsplit(plain, " ")[[1]][1]
  if (plain %in% benchmarks ||
      any(grepl(paste0("^", genus, " "), benchmarks))) {
    "literature-supported"
  } else {
    "novel"
  }
}

#' Normalise a TCMbio species name to plain text for benchmark matching.
#' Strips "s__" prefix and converts underscores to spaces.
#' e.g. "s__Bacteroides_fragilis" -> "Bacteroides fragilis"
normalise_species_name <- function(x) {
  x <- sub("^s__", "", x)
  x <- gsub("_", " ", x)
  trimws(x)
}
