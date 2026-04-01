#===========================
# Figure 3 Honours
#=============================

# notes....

#===========================
# loading libraries
#===========================

library(Seurat)
library(dplyr)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(patchwork)
library(stringr)
library(forcats)
library(scales)
library(scales)
library(readxl)
library(babelgene)


setwd("D:/Finn/Honours/figures/Figure_3")
outdir <- getwd()
set.seed(1)

# ==========================
# laoding Dataset 
# ==========================

seu <- readRDS("D:/Finn/Honours/figures/Figure_1/MS_CTRL_Whole_OBJECT.RDS")

# ==========================
# Figure 3A umap of AEMs in MS and CTRL
# =========================

# --------------------------------------------------
# Create grouped labels in seu
# --------------------------------------------------
seu$AEMS_Astro_group <- dplyr::case_when(
  seu$assigned_label == "Astrocyte_5" ~ "AEMS",
  seu$assigned_label %in% c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "Astrocyte_4", "Astrocyte_6") ~ "Astrocytes",
  TRUE ~ NA_character_
)

seu$AEMS_Astro_group <- factor(
  seu$AEMS_Astro_group,
  levels = c("Astrocytes", "AEMS")
)

# --------------------------------------------------
# Subset to only Astrocytes and AEMS
# --------------------------------------------------
seu_astro_aems <- subset(
  seu,
  subset = AEMS_Astro_group %in% c("Astrocytes", "AEMS")
)

# --------------------------------------------------
# Set identities
# --------------------------------------------------
Idents(seu_astro_aems) <- seu_astro_aems$AEMS_Astro_group

# --------------------------------------------------
# Plot UMAP with zoomed coordinates
# --------------------------------------------------
library(ggplot2)

AEMS_Astro_umap <- DimPlot(
  seu_astro_aems,
  reduction = "umap_MP",
  group.by = "AEMS_Astro_group",
  cols = c(
    "Astrocytes" = "#9ACD32",
    "AEMS" = "forestgreen"
  ),
  pt.size = 1
) +
  coord_cartesian(xlim = c(-15, 5), ylim = c(-6, 5)) +
  theme_classic(base_size = 16) +
  labs(title = "Astrocytes and AEMS (Zoomed UMAP)")


ggsave(AEMS_Astro_umap, file = "AEMS_Astro_umap.png", width = 10, height = 8)


ggsave(AEMS_Astro_umap, file = "AEMS_Astro_legend_umap.png", width = 10, height = 8)

# ==========================
# Figure 3B analysis and plotting code to add
# =========================


# ==========================
# Figure 3C analysis and plotting code to add
# ==========================

# ==========================
# Figure 3D analysis and plotting code to add
# ==========================

# ==========================
# Figure 3E analysis and plotting code to add
# ==========================

# ==========================
# Figure 3F PAMs MS vs CTRL
# ==========================

# ==========================================================
# PAMs signature UCell scoring and sample-level normalization
# ==========================================================

# -------------------------------
# Libraries
# -------------------------------
library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)
library(UCell)
library(tidyr)

# -------------------------------
# PAMs gene set
# -------------------------------
pams_genes <- c(
  "KHDRBS3","SIPA1L2","GAB3","TNFRSF1B","SAMSN1","ETV6","PIK3AP1","ITPR2",
  "TGFBR2","SPP1","BACH1","RAPGEF1","SCIN","GALNT2","MTHFD1L","ST6GAL1",
  "ADAM28","RBM47","A2M","SLC11A1","ALOX5","TGFBR1","RHBDF2","MSR1","RUNX2",
  "RIN3","AKAP13","SLC4A7","DOCK2","IKZF1","PELI1","SLCO2B1","PREX1","DOCK4",
  "IPCEF1","FAM102B","RAP1A","CD83","LYN","SRGAP2","SLC1A3","RASSF3","BASP1",
  "CTSB","USP15","WDFY4","TMEM156","SP100","IL6ST","CYFIP1","PAG1","PALD1",
  "OXR1","SFMBT2","LRRK1","SRGAP2B","SLC2A5","MGAT4A","KCNQ1","PLXDC2","PTPRJ",
  "CHST11","ATP8B4","TBC1D22A","LPAR6","GRB2","SH3KBP1","HIF1A","DISC1",
  "LDLRAD4","SYNDIG1","SKAP2","PDE3B","ST8SIA4","SRGAP2C","INPP5D","MEF2A",
  "SOCS6","TBXAS1","ACSL1","EPB41L2","RNASET2","MCF2L","CPVL","KYNU","FMN1",
  "FAM49B","FLI1","SSH2","SLA","FAM149A","FGD2","PICALM","ARHGAP25","CSF1R",
  "CACNA1D","CMTM7","KLHL6","FOXN3","BMP2K","DOCK8","PLA2G4A","MANBA","LY86",
  "SMAP2","KCNQ3","ARHGAP26","TRPM2","CTSC","PIK3R5","CPED1","KCNIP1","ENTPD1",
  "CD86","IGSF21","IRAK3","FGD4","HS3ST4","CACNA1A","FKBP5","TFEC"
)

pams_genes <- unique(pams_genes)
pams_genes <- pams_genes[!is.na(pams_genes)]
pams_genes <- pams_genes[pams_genes != ""]

# -------------------------------
# Keep only PAMs genes present in RNA assay
# -------------------------------
DefaultAssay(seu) <- "RNA"

pams_present <- intersect(pams_genes, rownames(seu[["RNA"]]))
pams_missing <- setdiff(pams_genes, pams_present)

cat("Total PAMs genes:", length(pams_genes), "\n")
cat("PAMs genes present in RNA assay:", length(pams_present), "\n")
cat("PAMs genes missing from RNA assay:", length(pams_missing), "\n\n")

if (length(pams_present) < 3) {
  stop("Too few PAMs genes are present in the RNA assay to score reliably.")
}

cat("First present PAMs genes:\n")
print(head(pams_present, 25))

# -------------------------------
# Score with UCell
# -------------------------------
sig_name <- "PAMs_Signature"
sig_list <- list(PAMs_Signature = pams_present)

seu <- AddModuleScore_UCell(
  seu,
  features = sig_list,
  assay = "RNA",
  name = NULL
)

score_col <- sig_name

if (!score_col %in% colnames(seu@meta.data)) {
  stop(paste0("Could not find score column: ", score_col))
}

# -------------------------------
# Define AEMS vs regular astrocytes
# -------------------------------
seu$AEMS_Astro_group <- dplyr::case_when(
  seu$assigned_label == "Astrocyte_5" ~ "AEMS",
  seu$assigned_label %in% c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "Astrocyte_4", "Astrocyte_6") ~ "Astrocytes",
  TRUE ~ NA_character_
)

seu$AEMS_Astro_group <- factor(
  seu$AEMS_Astro_group,
  levels = c("Astrocytes", "AEMS")
)

# -------------------------------
# Build metadata table
# -------------------------------
meta_df <- seu@meta.data %>%
  dplyr::select(orig.ident, assigned_label, AEMS_Astro_group, all_of(score_col)) %>%
  dplyr::filter(!is.na(AEMS_Astro_group)) %>%
  dplyr::mutate(
    Condition = stringr::str_extract(orig.ident, "CTRL|MS"),
    Condition = factor(Condition, levels = c("CTRL", "MS"))
  ) %>%
  dplyr::filter(!is.na(Condition))

colnames(meta_df)[colnames(meta_df) == score_col] <- "PAMs_UCell_score"

# -------------------------------
# Compute mean score per sample per group
# -------------------------------
sample_group_scores <- meta_df %>%
  dplyr::group_by(orig.ident, Condition, AEMS_Astro_group) %>%
  dplyr::summarise(
    mean_score = mean(PAMs_UCell_score, na.rm = TRUE),
    median_score = median(PAMs_UCell_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

print(sample_group_scores)

# -------------------------------
# Normalize AEMS score to Astrocyte score within sample
# -------------------------------
plot_df <- sample_group_scores %>%
  dplyr::select(orig.ident, Condition, AEMS_Astro_group, mean_score, median_score, n_cells) %>%
  tidyr::pivot_wider(
    names_from = AEMS_Astro_group,
    values_from = c(mean_score, median_score, n_cells),
    values_fill = NA
  ) %>%
  dplyr::mutate(
    normalized_score = mean_score_AEMS / mean_score_Astrocytes
  ) %>%
  dplyr::filter(
    !is.na(normalized_score),
    is.finite(normalized_score)
  )

cat("\nSample-level normalized PAMs scores:\n")
print(plot_df)

# -------------------------------
# Statistical test
# -------------------------------
wilcox_res <- wilcox.test(
  normalized_score ~ Condition,
  data = plot_df,
  exact = FALSE
)

cat("\nWilcoxon p-value:", wilcox_res$p.value, "\n")

# -------------------------------
# Plot
# -------------------------------
p <- ggplot(plot_df, aes(x = Condition, y = normalized_score, color = Condition)) +
  
  # Bigger points
  geom_jitter(width = 0.12, size = 5.5, alpha = 0.9) +
  
  # Median (thin black line)
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    color = "black",
    fatten = 0,
    linewidth = 0.8
  ) +
  
  # Mean (thicker + more prominent)
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "black",
    linewidth = 0.7
  ) +
  
  scale_color_manual(values = c("CTRL" = "blue2", "MS" = "red3")) +
  
  theme_classic(base_size = 18) +
  
  labs(
    x = "",
    y = "Mean PAMs UCell score in AEMS / Mean PAMs UCell score in Astrocytes",
    subtitle = paste0("Wilcoxon p = ", signif(wilcox_res$p.value, 3))
  ) +
  
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 25),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"
  )

print(p)

# -------------------------------
# Save plot
# -------------------------------
ggsave(
  "PAMs_signature_normalized_to_astrocytes_by_sample.png",
  plot = p,
  width = 7.5,
  height = 8,
  dpi = 300
)

# ==========================
# Figure 3G LEE ET AL IN ms VS ctrl
# ==========================

# ==========================================================
# Load mouse gene profile, convert to human orthologs,
# score with UCell on RNA, plot violin in MS vs CTRL,
# and print a p-value
# ==========================================================

# -------------------------------
# Libraries
# -------------------------------
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)
library(Seurat)
library(UCell)
library(babelgene)

# -------------------------------
# Input file
# -------------------------------
path <- "D:/Finn/Honours/figures/Figure_2/41586_2024_7187_MOESM4_ESM.xlsx"

# -------------------------------
# Load gene list from first column
# -------------------------------
genes <- read_excel(path, col_names = FALSE)[[1]]
genes <- genes[!is.na(genes)]
genes <- unique(as.character(genes))
genes <- genes[genes != ""]

cat("First 20 input genes:\n")
print(head(genes, 20))
cat("Total input genes:", length(genes), "\n\n")

# -------------------------------
# Convert mouse -> human
# -------------------------------
genes_mouse <- unique(genes)
genes_mouse <- genes_mouse[!is.na(genes_mouse)]
genes_mouse <- genes_mouse[genes_mouse != ""]

mapping <- babelgene::orthologs(
  genes = genes_mouse,
  species = "mouse",
  human = FALSE
)

mapping_clean <- mapping %>%
  dplyr::select(mouse_gene = symbol, human_gene = human_symbol, support) %>%
  dplyr::filter(!is.na(human_gene), human_gene != "") %>%
  dplyr::distinct()

mapping_best <- mapping_clean %>%
  dplyr::mutate(n_support = lengths(strsplit(support, "\\|"))) %>%
  dplyr::arrange(mouse_gene, dplyr::desc(n_support), human_gene) %>%
  dplyr::group_by(mouse_gene) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

genes_human <- unique(mapping_best$human_gene)

mapped_mouse <- unique(mapping_best$mouse_gene)
unmapped_mouse <- setdiff(genes_mouse, mapped_mouse)

cat("Input mouse genes:", length(genes_mouse), "\n")
cat("Mapped mouse genes:", length(mapped_mouse), "\n")
cat("Unmapped mouse genes:", length(unmapped_mouse), "\n\n")

cat("Top mappings:\n")
print(head(mapping_best, 25))

cat("\nUnmapped genes:\n")
print(unmapped_mouse)

cat("\nTop human genes:\n")
print(head(genes_human, 25))
cat("\nTotal mapped human genes:", length(genes_human), "\n\n")

# -------------------------------
# Keep only genes present in seu RNA assay
# -------------------------------
DefaultAssay(seu) <- "RNA"

genes_present <- intersect(genes_human, rownames(seu[["RNA"]]))
genes_missing <- setdiff(genes_human, genes_present)

cat("Human genes present in seu RNA assay:", length(genes_present), "\n")
cat("Human genes missing from seu RNA assay:", length(genes_missing), "\n\n")

if (length(genes_present) < 3) {
  stop("Too few mapped genes are present in the RNA assay to score reliably.")
}

# -------------------------------
# Build UCell signature
# -------------------------------
sig_name <- "ConvertedSignature"
sig_list <- list()
sig_list[[sig_name]] <- genes_present

# -------------------------------
# Score with UCell
# -------------------------------
seu <- AddModuleScore_UCell(
  seu,
  features = sig_list,
  assay = "RNA",
  name = NULL
)

# UCell adds the score as metadata with the list name
score_col <- sig_name
if (!score_col %in% colnames(seu@meta.data)) {
  stop(paste0("Could not find UCell score column: ", score_col))
}


# -------------------------------
# Make sure the AEMS UCell score exists
# Replace this if your score column has a different name
# -------------------------------
score_col <- "ConvertedSignature"

if (!score_col %in% colnames(seu@meta.data)) {
  stop(paste0("Could not find score column: ", score_col))
}

# -------------------------------
# Define AEMS vs regular astrocytes
# -------------------------------
seu$AEMS_Astro_group <- dplyr::case_when(
  seu$assigned_label == "Astrocyte_5" ~ "AEMS",
  seu$assigned_label %in% c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "Astrocyte_4", "Astrocyte_6") ~ "Astrocytes",
  TRUE ~ NA_character_
)

seu$AEMS_Astro_group <- factor(
  seu$AEMS_Astro_group,
  levels = c("Astrocytes", "AEMS")
)

# -------------------------------
# Build metadata table
# -------------------------------
meta_df <- seu@meta.data %>%
  dplyr::select(orig.ident, assigned_label, AEMS_Astro_group, all_of(score_col)) %>%
  dplyr::filter(!is.na(AEMS_Astro_group)) %>%
  dplyr::mutate(
    Condition = stringr::str_extract(orig.ident, "CTRL|MS"),
    Condition = factor(Condition, levels = c("CTRL", "MS"))
  ) %>%
  dplyr::filter(!is.na(Condition))

colnames(meta_df)[colnames(meta_df) == score_col] <- "UCell_score"

# -------------------------------
# Compute mean score per sample per group
# -------------------------------
sample_group_scores <- meta_df %>%
  dplyr::group_by(orig.ident, Condition, AEMS_Astro_group) %>%
  dplyr::summarise(
    mean_score = mean(UCell_score, na.rm = TRUE),
    median_score = median(UCell_score, na.rm = TRUE),
    n_cells = dplyr::n(),
    .groups = "drop"
  )

print(sample_group_scores)

# -------------------------------
# Convert to wide format
# Normalize AEMS score to Astrocyte score within sample
# -------------------------------
plot_df <- sample_group_scores %>%
  dplyr::select(orig.ident, Condition, AEMS_Astro_group, mean_score, n_cells) %>%
  tidyr::pivot_wider(
    names_from = AEMS_Astro_group,
    values_from = c(mean_score, n_cells),
    values_fill = NA
  ) %>%
  dplyr::mutate(
    normalized_score = mean_score_AEMS / mean_score_Astrocytes
  ) %>%
  dplyr::filter(
    !is.na(normalized_score),
    is.finite(normalized_score)
  )

# -------------------------------
# Inspect per-sample normalized values
# -------------------------------
print(plot_df)

# -------------------------------
# Statistical test at sample level
# -------------------------------
wilcox_res <- wilcox.test(
  normalized_score ~ Condition,
  data = plot_df,
  exact = FALSE
)

print(wilcox_res)

cat("Wilcoxon p-value:", wilcox_res$p.value, "\n")

# -------------------------------
# Plot
# -------------------------------
p <- ggplot(plot_df, aes(x = Condition, y = normalized_score, color = Condition)) +
  
  # Bigger points
  geom_jitter(width = 0.12, size = 5.5, alpha = 0.9) +
  
  # Median (thin black line)
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    color = "black",
    fatten = 0,
    linewidth = 0.8
  ) +
  
  # Mean (thicker + more prominent)
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "black",
    linewidth = 0.7
  ) +
  
  scale_color_manual(values = c("CTRL" = "blue2", "MS" = "red3")) +
  
  theme_classic(base_size = 18) +
  
  labs(
    x = "",
    y = "AEMS Mean /Astrocyte Mean, Lee et al. Memory",
  ) +
  
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 25),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"
  )

print(p)

# -------------------------------
# Save plot
# -------------------------------
ggsave(
  "AEMS_signature_normalized_to_astrocytes_by_sample.png",
  plot = p,
  width = 7.5,
  height = 8,
  dpi = 300
)

# -------------------------------
# Save plot
# -------------------------------
ggsave(
  filename = "UCell_converted_signature_MS_vs_CTRL_violin.png",
  plot = p,
  width = 7,
  height = 7,
  dpi = 300
)

# -------------------------------
# Optional: print medians
# -------------------------------
cat("Median UCell score by condition:\n")
print(
  plot_df %>%
    dplyr::group_by(Condition) %>%
    dplyr::summarise(
      n_cells = dplyr::n(),
      median_score = median(UCell_score, na.rm = TRUE),
      mean_score = mean(UCell_score, na.rm = TRUE),
      .groups = "drop"
    )
)

# ==========================
# Figure 3H aems AMOUNT IN ms vs CTRL
# ==========================

# --------------------------------------------------
# Libraries
# --------------------------------------------------
library(dplyr)
library(ggplot2)
library(stringr)

# --------------------------------------------------
# Create grouped labels from assigned_label
# --------------------------------------------------
seu$AEMS_Astro_group <- dplyr::case_when(
  seu$assigned_label == "Astrocyte_5" ~ "AEMS",
  seu$assigned_label %in% c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "Astrocyte_4", "Astrocyte_6") ~ "Astrocytes",
  TRUE ~ NA_character_
)

# --------------------------------------------------
# Pull metadata
# --------------------------------------------------
meta_df <- seu@meta.data %>%
  dplyr::select(orig.ident, AEMS_Astro_group) %>%
  dplyr::filter(!is.na(AEMS_Astro_group))

# --------------------------------------------------
# Count AEMS and Astrocytes per sample
# --------------------------------------------------
count_df <- meta_df %>%
  dplyr::group_by(orig.ident, AEMS_Astro_group) %>%
  dplyr::summarise(n_cells = n(), .groups = "drop")

# --------------------------------------------------
# Convert to wide format
# --------------------------------------------------
count_wide <- count_df %>%
  tidyr::pivot_wider(
    names_from = AEMS_Astro_group,
    values_from = n_cells,
    values_fill = 0
  )

# --------------------------------------------------
# Compute normalized AEMS abundance
# AEMS cells / Astrocyte cells within each sample
# --------------------------------------------------
count_wide <- count_wide %>%
  dplyr::mutate(
    AEMS_per_Astrocyte = AEMS / Astrocytes
  )

# --------------------------------------------------
# Extract condition from orig.ident
# Assumes names like sample1_CTRL or sample3_MS
# --------------------------------------------------
count_wide <- count_wide %>%
  dplyr::mutate(
    Condition = stringr::str_extract(orig.ident, "CTRL|MS"),
    Condition = factor(Condition, levels = c("CTRL", "MS"))
  ) %>%
  dplyr::filter(!is.na(Condition))

# --------------------------------------------------
# Print table so you can inspect sample-level values
# --------------------------------------------------
print(count_wide)

# --------------------------------------------------
# Wilcoxon test
# --------------------------------------------------
wilcox_res <- wilcox.test(
  AEMS_per_Astrocyte ~ Condition,
  data = count_wide,
  exact = FALSE
)

print(wilcox_res)

# Optional: also run t-test
ttest_res <- t.test(
  AEMS_per_Astrocyte ~ Condition,
  data = count_wide
)

print(ttest_res)

# --------------------------------------------------
# Make plotting dataframe
# --------------------------------------------------
plot_df <- count_wide

# --------------------------------------------------
# Plot
# --------------------------------------------------
p <- ggplot(plot_df, aes(x = Condition, y = AEMS_per_Astrocyte, color = Condition)) +
  
  # Bigger points
  geom_jitter(width = 0.12, size = 5.5, alpha = 0.9) +
  
  # Median (thin black line)
  stat_summary(
    fun = median,
    geom = "crossbar",
    width = 0.45,
    color = "black",
    fatten = 0,
    linewidth = 0.8
  ) +
  
  # Mean (thicker + more prominent)
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.6,
    color = "black",
    linewidth = 0.7
  ) +
  
  scale_color_manual(values = c("CTRL" = "blue2", "MS" = "red3")) +
  
  theme_classic(base_size = 18) +
  
  labs(
    x = "",
    y = "AEMS cells / Astrocyte cells"
  ) +
  
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    plot.subtitle = element_text(size = 25),
    axis.title = element_text(size = 25),
    axis.text = element_text(size = 25),
    legend.position = "none"
  )
print(p)

# --------------------------------------------------
# Save plot
# --------------------------------------------------
ggsave("AEMS_per_Astrocyte_by_condition.png", plot = p, width = 7.5, height = 8, dpi = 300)














