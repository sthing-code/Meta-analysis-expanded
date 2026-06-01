# Intratumoral Microbiome × Cancer Gene Expression (Expanded Correlation Analysis)

## Overview

This repository contains the full analysis pipeline correlating intratumoral bacterial abundance
(from TCMbio) with the expression of cancer-related genes (from TCGA) across multiple cancer types.
The primary outputs are Pearson correlation heatmaps and three-way correlation plots per cancer type.

This is the **expanded panel version (v2.0)** of the analysis. The original four-cancer analysis
(COAD, BRCA, PAAD, PRAD) used a 13-gene panel. This repository extends to a 49-gene panel and
adds Glioblastoma multiforme (GBM) as a fifth cancer type.

**Cancer types covered:**
- Colon adenocarcinoma (COAD) — original panel
- Breast invasive carcinoma (BRCA) — original panel
- Pancreatic adenocarcinoma (PAAD) — original panel
- Prostate adenocarcinoma (PRAD) — original panel
- Glioblastoma multiforme (GBM) — expanded panel v2.0

---

## Priority Genes (v2.0 — 49 genes)

### Original panel

| Gene symbol | Protein | Group | Notes |
|---|---|---|---|
| IGFBP7 | IGF-binding protein 7 | IGFBP7/Galectin | Strong positive correlations flagged |
| LGALS1 | Galectin-1 | IGFBP7/Galectin | Strong positive correlations flagged |
| VEGFA | VEGF-A | Angiogenesis | Downstream: MAPK pathway |
| MMP9 | MMP-9 | Angiogenesis | |
| IL6 | IL-6 | JAK-STAT | Downstream: STAT3 |
| STAT3 | STAT3 | JAK-STAT | |
| APC | APC | WNT | Downstream: WNT/β-catenin |
| BRCA2 | BRCA2 | DNA repair | Strong negative correlations flagged |
| CTSB | Cathepsin B | Cathepsins | Strong negative correlations flagged |
| CTSD | Cathepsin D | Cathepsins | Strong negative correlations flagged |
| CTSL | Cathepsin L | Cathepsins | Strong negative correlations flagged |
| CTSS | Cathepsin S | Cathepsins | Strong negative correlations flagged |
| CTSV | Cathepsin V | Cathepsins | Strong negative correlations flagged |

### Expanded panel (v2.0 — new genes)

| Gene symbol | Protein | Group | Notes |
|---|---|---|---|
| CD274 | PD-L1 | Immune checkpoint | Immune evasion marker |
| HIF1A | HIF-1α | Hypoxia | Master hypoxia regulator |
| ERBB2 | HER2 | RTK signalling | Downstream: FGFR2, FGFR3 |
| ERBB3 | HER3 | RTK signalling | |
| TP53 | p53 | DNA damage | Downstream: MDM2, MDM4, ATM, ATR |
| ATM | ATM | DNA damage | DNA damage kinase |
| ATR | ATR | DNA damage | Replication stress kinase |
| MDM2 | MDM2 | DNA damage | Negative regulator of p53 |
| WNT3A | WNT3A | WNT extended | Canonical WNT ligand |
| TCF7 | TCF1 | WNT extended | WNT transcription factor |
| LEF1 | LEF1 | WNT extended | WNT transcription factor |
| SOX2 | SOX2 | Stemness | Pluripotency transcription factor |

---

## Pathway Genes (three-way axes)

| Group | Genes | Three-way axis |
|---|---|---|
| il6_stat3 | IL6, STAT3 | IL-6 → STAT3 |
| vegfa_mapk | VEGFA, MAPK1, MAPK3, KRAS | VEGFA → MAPK |
| apc_wnt | APC, CTNNB1, MYC, CCND1 | APC → WNT |
| metabolic | LDHA, SLC2A1, SLC16A1, SLC16A3, PFKP, TIGAR, SESN2 | LDHA → metabolic reprogramming |
| rtk_extended | ERBB2, FGFR2, FGFR3 | ERBB2 → RTK signalling |
| tca_idh | GLUD1, IDH1, IDH2 | IDH1 → TCA cycle |
| nfkb | NFKB1 | NF-κB → IL-6 → STAT3 |
| p53_axis | TP53, MDM4 | TP53 → DNA damage response |
| stemness_ext | SOX2, KLF4, POU5F1 | SOX2 → stemness |
| glioma_marker | GFAP, RBFOX3 | Heatmap only (tissue markers) |

---

## Data Sources

### Gene expression
- **Source:** TCGA via Bioconductor `TCGAbiolinks`
- **Workflow:** STAR — Counts (GDC harmonised, hg38)
- **Normalisation:** TPM (transcripts per million), computed from STAR raw counts + gene lengths
- **Sample filter:** Primary tumour samples only (TCGA barcode suffix `-01`)

### Bacterial abundance
- **Source:** TCMbio (https://microbiomex.sdu.edu.cn/) — cancer-type-specific download
- **Format:** Species × sample abundance matrix (CSV), downloaded manually per cancer type
- **Preprocessing:** log₁₀(x + 1) transformation; prevalence filter ≥ 10% of samples

---

## Project Structure

```
.
├── README.md
├── .gitignore
│
├── R/                          # Shared helper functions
│   ├── utils.R                 # Gene lists, cancer configs, thresholds, benchmarks,
│   │                           # barcode matching, TPM computation, Pearson helpers
│   └── plot_themes.R           # ggplot2 / ComplexHeatmap theme + gene group colours
│
├── scripts/
│   ├── 01_download_TCGA.R      # Download TCGA data for target cancer type
│   ├── 02_preprocess_TCGA.R    # Normalise, filter, extract priority genes
│   ├── 03_preprocess_TCMbio.R  # Load and clean TCMbio bacterial abundance data
│   ├── 04_correlate.R          # Pearson correlations + BH FDR correction
│   ├── 05_heatmap.R            # Correlation heatmaps (primary + supplementary)
│   └── 06_threeway.R           # Three-way correlation plots (all axes)
│
├── data/
│   ├── raw/                    # Downloaded files (not tracked by git — see .gitignore)
│   └── processed/              # RDS files output by preprocessing scripts
│                               # (SE full objects excluded from git — too large)
│
└── figures/
    ├── colon/
    ├── breast/
    ├── pancreatic/
    ├── prostate/
    └── glioblastoma/
```

---

## Reproduction Steps

Run scripts in order. Change the active cancer type string at the bottom of each script
before sourcing (e.g. `"glioblastoma"`, `"colon"`, `"breast"`, `"pancreatic"`, `"prostate"`).

```r
source("scripts/01_download_TCGA.R")
source("scripts/02_preprocess_TCGA.R")
# Manually download TCMbio CSV — see instructions in script 03
source("scripts/03_preprocess_TCMbio.R")
source("scripts/04_correlate.R")
# Review script 04 output to set heatmap thresholds in R/utils.R before continuing
source("scripts/05_heatmap.R")
source("scripts/06_threeway.R")
```

---

## Dependencies

```r
# Bioconductor
BiocManager::install(c("TCGAbiolinks", "SummarizedExperiment", "DESeq2",
                       "ComplexHeatmap", "circlize"))

# CRAN
install.packages(c("tidyverse", "data.table", "corrplot",
                   "RColorBrewer", "ggrepel", "patchwork", "here"))
```

R version: ≥ 4.3.0 recommended.

---

## Correlation Thresholds and Flagging

Thresholds are set per cancer type based on observed signal strength and sample size, and are
defined in `R/utils.R` (`HEATMAP_THRESHOLDS`). The rationale for each cancer type is as follows.

**COAD** yields the strongest correlations in the dataset (104 bacteria at |r| > 0.20),
justifying the strictest primary threshold.

**BRCA** shows substantially weaker correlations overall (9 bacteria at |r| > 0.20). A relaxed
primary threshold of |r| > 0.15 is used and documented transparently.

**PAAD** is constrained by small dataset size (n = 178, 23 species post-filter). The |r| > 0.15
primary threshold is used, consistent with BRCA. The |r| > 0.20 view is retained as supplementary.

**PRAD** follows the same rationale as PAAD (n = 497, 32 species post-filter). Primary |r| > 0.15;
supplementary |r| > 0.20.

**GBM** has the weakest overall signal in the dataset (n = 154, 45 species post-filter, 4
FDR-significant pairs). Primary threshold |r| > 0.15 (34 bacteria); supplementary |r| > 0.20
(19 bacteria). The weak signal is consistent with the small matched sample size and a substantial
contamination burden from skin/environmental species introduced during craniotomy-based tissue
access.

---

## GBM-Specific Notes

### Contamination caveat
GBM tissue is accessed via craniotomy. Skin and environmental bacteria are a documented
contamination source in brain tumour microbiome studies (Mehelleb et al. 2025). The following
species dominate the top-abundance ranking and should be interpreted with high scepticism:
*Pseudomonas* sp. CIP-10 (rank 1), *Cutibacterium acnes* (rank 2), *Staphylococcus aureus*
(rank 4), *Cutibacterium granulosum* (rank 5), *Staphylococcus epidermidis* (rank 9),
*Corynebacterium kroppenstedtii* (rank 8).

### Key GBM findings
- Single FDR-significant hit: *Priestia megaterium* × CD274 (r = 0.155, q < 0.05) — novel,
  no GBM literature, requires independent validation
- Strongest three-way signal: NF-κB → STAT3 negative co-suppression cluster (13 bacteria)
- Second strongest: SOX2 → stemness co-suppression (4 bacteria; KLF4 and POU5F1 concordant)
- ERBB2 × FGFR2/FGFR3 divergence: positive ERBB2 tracks FGFR2; negative ERBB2 tracks FGFR3
- No benchmark bacteria (F. nucleatum, Sphingomonas, Burkholderia, Helicobacter) confirmed
  above |r| > 0.15 threshold

---

## Literature Benchmarks (Colon — COAD)

### Enriched in CRC (pro-tumorigenic)

| Species | Mechanism | Citation |
|---|---|---|
| *Fusobacterium nucleatum* | E-cadherin/β-catenin invasion; NF-κB/IL-6 activation; FadA adhesin; pan-cancer prognostic risk | Kostic et al. 2012; Rubinstein et al. 2013 |
| *Parvimonas micra* | Enriched in CRC tissue; pan-cancer prognostic risk | Yu et al. 2017; Proctor et al. 2021 |
| *Peptostreptococcus stomatis* | Enriched in CRC tissue | Yu et al. 2017 |
| *Solobacterium moorei* | Enriched in CRC tissue | Yu et al. 2017 |
| *Bacteroides fragilis* | Enterotoxin BFT activates STAT3; promotes colitis-associated CRC | Dalal et al. 2021; Wu et al. 2004 |
| *Streptococcus gallolyticus* | Associated with colonic polyps and CRC | Dalal et al. 2021 |
| *Enterococcus faecalis* | Extracellular superoxide production; DNA damage | Dalal et al. 2021 |
| *Porphyromonas gingivalis* | Oral pathobiont; enriched in CRC; pan-cancer prognostic risk | Bautista et al. 2026; Sheng et al. 2024 |
| *Helicobacter pylori* | CagA oncoprotein; carcinogenic in upper GI; associated with CRC | Dalal et al. 2021 |

### Depleted in CRC (protective)

| Species | Mechanism | Citation |
|---|---|---|
| *Faecalibacterium prausnitzii* | Major butyrate producer; anti-inflammatory; consistently depleted in CRC | Bautista et al. 2026; Cao et al. 2025 |
| *Roseburia intestinalis* | Butyrate producer; depleted in CRC-associated dysbiosis | Bautista et al. 2026 |
| *Roseburia hominis* | Butyrate producer; depleted in CRC-associated dysbiosis | Bautista et al. 2026 |

---

## Literature Benchmarks (Breast — BRCA)

The dominant mechanism is **oestrogen metabolism via the estrobolome** — bacteria producing
β-glucuronidase, sulfatase, and hydroxysteroid dehydrogenase enzymes that reactivate conjugated
oestrogens, elevating circulating hormone levels and driving HR+ tumour growth. A second
mechanism operates in TNBC via **TMAO-producing Clostridiales**.

### Enriched in breast cancer (pro-tumorigenic or risk-associated)

| Species/Genus | Mechanism | Notes | Citation |
|---|---|---|---|
| *Escherichia coli* | β-glucuronidase; oestrogen deconjugation → elevated HR+ tumour growth | Strong, consistent | Larnder et al. 2025; Mahno et al. 2024; Sheng et al. 2024 |
| *Bacteroides fragilis* | Broadest estrobolome enzyme activity (β-glucuronidase + sulfatase + 3β-HSD + 17β-HSD) | Strong | Larnder et al. 2025 |
| *Fusobacterium nucleatum* | Enriched in breast tumour tissue; pro-inflammatory | Moderate | Sheng et al. 2024 |
| *Ruminococcus* spp. | β-glucuronidase activity; **conflicting evidence** across cohorts | Mixed | Larnder et al. 2025; Mahno et al. 2024 |
| *Klebsiella* spp. | β-glucuronidase activity; enriched in cases | Moderate | Mahno et al. 2024 |
| *Sphingomonas* spp. | Enriched in breast tumour tissue; **cross-cancer signal — also top IGFBP7 correlator in COAD** | Moderate | Sheng et al. 2024 |

### Depleted in breast cancer (protective)

| Species/Genus | Mechanism | Citation |
|---|---|---|
| *Faecalibacterium prausnitzii* | Suppresses breast cancer cell growth via IL-6/STAT3 inhibition; butyrate producer | Ruo et al. 2021; Larnder et al. 2025 |
| *Roseburia* spp. | Butyrate producers; β-glucosidase activates protective phytoestrogens | Larnder et al. 2025 |
| *Bifidobacterium* spp. | Anti-inflammatory; consistently protective | Larnder et al. 2025; Mahno et al. 2024 |
| *Collinsella aerofaciens* | Depleted in cases; protective gut commensal | Larnder et al. 2025 |

### TMAO-producing Clostridiales — antitumour in TNBC

| Genera | Evidence | Citation |
|---|---|---|
| *Blautia*, *Dorea*, *Ruminococcus*, *Tyzzerella*, *Roseburia* | TMAO → PERK ER stress → GSDME pyroptosis → CD8+ T cell immunity; n = 360 TNBC cohort | Wang et al. 2022 |

---

## Literature Benchmarks (Pancreatic — PAAD)

Three defining features: (1) oral-to-pancreas bacterial translocation; (2) intratumoral
gemcitabine inactivation; (3) TMAO anti-tumour immunity (cross-cancer bridge with BRCA).

### Pro-tumorigenic: oral-to-pancreas translocation

| Species/Genus | Mechanism | Citation |
|---|---|---|
| *Porphyromonas gingivalis* | Promotes PDAC in K-rasG12D mouse models; induces immunosuppressive TME | Sheng et al. 2024; Wang et al. 2024 |
| *Fusobacterium nucleatum* | Oral pathobiont enriched in PAAD tissue; pro-inflammatory | Sheng et al. 2024 |
| *Prevotella oris* | Enriched in early cystic PDAC precursors | Sheng et al. 2024 |
| *Veillonella parvula* | Oral commensal; pan-cancer prognostic risk | Sheng et al. 2024 |

### Pro-tumorigenic: gemcitabine inactivation

| Species/Group | Mechanism | Citation |
|---|---|---|
| *Mycoplasma* spp. | Inactivates gemcitabine intracellularly; key chemotherapy resistance determinant | Geller et al. 2017, cited in Sheng et al. 2024 |
| Gammaproteobacteria (class) | Cytidine deaminase inactivates gemcitabine in TME | Bautista et al. 2026 |

### Antitumour: TMAO-producing bacteria

| Genera | Mechanism | Citation |
|---|---|---|
| Clostridiales (*Blautia*, *Dorea*, *Ruminococcus*, *Roseburia*) | TMAO → type I IFN in macrophages → cytotoxic T cell activation → ICI sensitisation | Mirji et al. 2022, cited in Liu et al. 2023; W. Zhang et al. 2024 |

> **Cross-cancer note:** This is the same Clostridiales–TMAO mechanism documented in BRCA
> (Wang et al. *Cell Metabolism* 2022). Consistent positive gene correlations in both BRCA
> and PAAD heatmaps would constitute a strong cross-cancer mechanistic story.

---

## Literature Benchmarks (Prostate — PRAD)

PRAD benchmarks are drawn from **direct intratumoral detection** in prostate cancer tissue, not
gut microbiome dysbiosis studies. Two mechanistic axes are relevant to the priority gene panel:
*C. acnes* directly downregulates BRCA2 and drives immunosuppression via PD-L1/CCL17/CCL18.

| # | Species | Mechanism | Confidence | Citation |
|---|---|---|---|---|
| 1 | *Cutibacterium acnes* | Invades prostate epithelial cells; significantly downregulates BRCA2; impairs homologous recombination (BRCAness); promotes DNA double-strand breaks | High | Ashida et al. 2024 |
| 2 | *Cutibacterium acnes* | Stimulates macrophage M2 polarisation; upregulates PD-L1, CCL17, CCL18; increases Tregs in stroma (P = 0.0004) and epithelia (P = 0.046) | High | Davidsson et al. 2021 |
| 3 | *Moraxella osloensis* | De novo RNA-seq detection; validated by PCR in 20/20 PC tissue samples; no mechanistic data | Moderate | Ashida et al. 2024 |
| 4 | *Micrococcus luteus* | De novo RNA-seq detection; validated by PCR in 16/20 PC tissue samples; no mechanistic data | Moderate | Ashida et al. 2024 |

> **Note:** *C. acnes* is excluded from the GBM benchmark — it is a well-established skin
> commensal and a documented craniotomy contamination source. No direct GBM tissue evidence
> exists for this species.

---

## Literature Benchmarks (GBM)

| Species (representative) | Confidence | Source | Notes |
|---|---|---|---|
| *Fusobacterium nucleatum* | HIGH | Li et al. 2024 (mSystems, DOI: 10.1128/msystems.00457-24) | Enrichment confirmed by IHC, immunofluorescence, FISH in glioma tissue (n = 50); functional evidence: upregulates CCL2, CXCL1, CXCL2 in xenograft and organoid models. Mixed glioma grades — not exclusively GBM. |
| *Burkholderia cepacia* | MODERATE | Mehelleb et al. 2025 (Int J Mol Sci) | Burkholderia-Caballeronia-Paraburkholderia most abundant genus (27.43%); no genus reached statistical significance; n = 9 GBM patients; 16S rRNA fresh tissue |
| *Sphingomonas melonis* | MODERATE | Mehelleb et al. 2025 | Second most abundant genus (12.90%); abundance ranking only |
| *Leifsonia xyli* | MODERATE | Mehelleb et al. 2025 | 8.01% abundance; abundance ranking only |
| *Helicobacter pylori* | MODERATE | Mehelleb et al. 2025 | 4.16% abundance; abundance ranking only |
| *Cutibacterium acnes* | EXCLUDED | — | Skin contamination source; no direct GBM tissue evidence |

All MODERATE entries use one representative species per genus to trigger genus-level matching
in the pipeline flagging function. TCMbio operates at species level; genus matching is applied.

---

## Key References

- Ashida et al. 2024 — *Cutibacterium acnes* in prostate cancer; BRCAness mechanism
- Bautista et al. 2026 — Pan-cancer microbiome; gemcitabine resistance; protective commensals
- Cao et al. 2025 — *F. prausnitzii* in CRC
- Davidsson et al. 2021 — *C. acnes* immunosuppression in prostate cancer
- Kostic et al. 2012 / Rubinstein et al. 2013 — *F. nucleatum* in CRC
- Larnder et al. 2025 — Estrobolome and breast cancer
- Li et al. 2024 (mSystems) — *Fusobacterium* in glioma
- Liu et al. 2023 — TMAO in pancreatic cancer
- Mahno et al. 2024 — Gut microbiome and breast cancer systematic review
- Mehelleb et al. 2025 — GBM intratumoral microbiome (16S rRNA, fresh tissue)
- Mirji et al. 2022 — TMAO and anti-tumour immunity
- Ruo et al. 2021 — *F. prausnitzii* in breast cancer
- Sheng et al. 2024 — Pan-cancer tumour-resident microbiome atlas
- Wang et al. 2022 — TMAO in TNBC
- Wang et al. 2024 — Oral microbiota in pancreatic cancer
- Yu et al. 2017 — CRC microbiome
