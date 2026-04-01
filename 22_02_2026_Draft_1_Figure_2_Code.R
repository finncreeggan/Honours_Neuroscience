#===========================
# Figure 2 Honours
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
library(readxl)
library(babelgene)



setwd("D:/Finn/Honours/figures/Figure_2")
outdir <- getwd()
set.seed(1)


# ==========================
# laoding Data and listing properties
# ==========================

seu <- readRDS("D:/Finn/Honours/figures/Figure_1/MS_CTRL_Whole_OBJECT.RDS")

head(seu@meta.data)

path <- "D:/Finn/Honours/figures/Figure_2/41586_2024_7187_MOESM4_ESM.xlsx"

genes <- read_excel(path, col_names = FALSE)[[1]]
genes <- genes[!is.na(genes)]
genes <- unique(as.character(genes))

print(head(genes, 20))
print(length(genes))

# =======================================
# convert mouse genes to human genes
# =======================================

genes_mouse <- genes
genes_mouse <- unique(genes_mouse)
genes_mouse <- genes_mouse[!is.na(genes_mouse)]



# 2. Map mouse -> human
#    species can be common name or taxon

?orthologs

mapping <- babelgene::orthologs(
  genes = genes_mouse,
  species = "mouse",
  human = FALSE
)

# Look at columns
print(colnames(mapping))
head(mapping)
View(mapping)


# 3. Keep the main useful columns
#    In babelgene output:
#    - symbol is the input-species gene (mouse here)
#    - human_symbol is the mapped human ortholog

mapping_clean <- mapping %>%
  dplyr::select(mouse_gene = symbol, human_gene = human_symbol, support) %>%
  dplyr::filter(!is.na(human_gene), human_gene != "") %>%
  dplyr::distinct()


# 4. For genes with multiple human orthologs,
#    keep the one with the most support sources

mapping_best <- mapping_clean %>%
  dplyr::mutate(n_support = lengths(strsplit(support, "\\|"))) %>%
  dplyr::arrange(mouse_gene, dplyr::desc(n_support), human_gene) %>%
  dplyr::group_by(mouse_gene) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()


# 5. Final human gene vector

genes_human <- unique(mapping_best$human_gene)


# 6. Report what mapped and what did not

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


# 7. Save both outputs to working directory

write.csv(mapping_best, "mouse_to_human_ortholog_mapping.csv", row.names = FALSE)
write.csv(data.frame(human_gene = genes_human), "human_gene_list.csv", row.names = FALSE)

# ==========================
# scoring RNA and GA 
# ==========================
# 
# genes <- genes_human
# 
# 
# # 1. Score gene profile at RNA level
# 
# DefaultAssay(seu) <- "RNA"
# 
# rna_genes_present <- intersect(genes, rownames(seu[["RNA"]]))
# cat("RNA genes present:", length(rna_genes_present), "of", length(genes), "\n")
# 
# seu <- AddModuleScore(
#   object = seu,
#   features = list(rna_genes_present),
#   assay = "RNA",
#   name = "GeneProfile_RNA"
# )
# 
# # Rename the default added column to something cleaner
# seu$GeneProfile_RNA_Score <- seu$GeneProfile_RNA1
# 
# 
# # 2. Score gene profile at GeneActivity level
# 
# DefaultAssay(seu) <- "GeneActivity"
# 
# ga_genes_present <- intersect(genes, rownames(seu[["GeneActivity"]]))
# cat("GeneActivity genes present:", length(ga_genes_present), "of", length(genes), "\n")
# 
# seu <- AddModuleScore(
#   object = seu,
#   features = list(ga_genes_present),
#   assay = "GeneActivity",
#   name = "GeneProfile_GA"
# )
# 
# # Rename the default added column to something cleaner
# seu$GeneProfile_GA_Score <- seu$GeneProfile_GA1
# 
# 
# # 3. Subset to CTRL cells only
# 
# seu_CTRL <- subset(seu, subset = condition == "CTRL")
# 
# 
# # 4. Quick checks
# 
# table(seu_CTRL$condition)
# 
# summary(seu$GeneProfile_RNA_Score)
# summary(seu$GeneProfile_GA_Score)
# 
# summary(seu_CTRL$GeneProfile_RNA_Score)
# summary(seu_CTRL$GeneProfile_GA_Score)
# 
# 
# # 5. Optional: visualize scores in CTRL only
# 
# DefaultAssay(seu_CTRL) <- "RNA"
# 
# p1 <- FeaturePlot(
#   seu_CTRL,
#   reduction = "umap_MP",
#   features = "GeneProfile_RNA_Score"
# )
# 
# DefaultAssay(seu_CTRL) <- "GeneActivity"
# 
# p2 <- FeaturePlot(
#   seu_CTRL,
#   reduction = "umap_MP",
#   features = "GeneProfile_GA_Score"
# )
# 
# p1
# p2


# ==========================
# Score gene profile with UCell
# ==========================


library(UCell)

genes <- genes_human

# Clean gene list
genes <- unique(as.character(genes))
genes <- genes[!is.na(genes) & genes != ""]

# =========================================================
# 1. Score gene profile at RNA level with UCell
# =========================================================
DefaultAssay(seu) <- "RNA"

rna_genes_present <- intersect(genes, rownames(seu[["RNA"]]))
cat("RNA genes present:", length(rna_genes_present), "of", length(genes), "\n")

seu <- AddModuleScore_UCell(
  seu,
  features = list(GeneProfile_RNA = rna_genes_present),
  assay = "RNA"
)

# =========================================================
# 2. Score gene profile at GeneActivity level with UCell
# =========================================================
DefaultAssay(seu) <- "GeneActivity"

ga_genes_present <- intersect(genes, rownames(seu[["GeneActivity"]]))
cat("GeneActivity genes present:", length(ga_genes_present), "of", length(genes), "\n")

seu <- AddModuleScore_UCell(
  seu,
  features = list(GeneProfile_GA = ga_genes_present),
  assay = "GeneActivity"
)

# =========================================================
# 3. Subset to CTRL cells only
# =========================================================
seu_CTRL <- subset(seu, subset = condition == "CTRL")

# =========================================================
# 4. Quick checks
# =========================================================
table(seu_CTRL$condition)

summary(seu$GeneProfile_RNA_UCell)
summary(seu$GeneProfile_GA_UCell)

summary(seu_CTRL$GeneProfile_RNA_UCell)
summary(seu_CTRL$GeneProfile_GA_UCell)

# =========================================================
# 5. Optional: visualize scores in CTRL only
# =========================================================
DefaultAssay(seu_CTRL) <- "RNA"

p1 <- FeaturePlot(
  seu_CTRL,
  reduction = "umap_MP",
  features = "GeneProfile_RNA_UCell"
)

DefaultAssay(seu_CTRL) <- "GeneActivity"

p2 <- FeaturePlot(
  seu_CTRL,
  reduction = "umap_MP",
  features = "GeneProfile_GA_UCell"
)

p1
p2


# ==========================
# Figure 2B conocial Astrocyte and lee et al memory profile RNA
# =========================

table(seu_CTRL@meta.data$assigned_label)


# --------------------------------------------------
# 1. Subset to the 4 groups of interest
# --------------------------------------------------
table(seu_CTRL@meta.data$assigned_label)

seu_fig1b <- subset(
  seu_CTRL,
  subset = assigned_label %in% c(
    "Astrocyte_1",
    "Astrocyte_2",
    "Astrocyte_3",
    "Astrocyte_4",
    "Astrocyte_5",
    "Astrocyte_6",
    "Microglia",
    "Neuron"
  )
)

# --------------------------------------------------
# 2. Create grouped labels
# --------------------------------------------------
seu_fig1b$Figure1B_group <- dplyr::case_when(
  seu_fig1b$assigned_label == "Astrocyte_5" ~ "AEMS",
  seu_fig1b$assigned_label %in% c("Astrocyte_1", "Astrocyte_2", "Astrocyte_3", "Astrocyte_4", "Astrocyte_6") ~ "Astrocytes",
  seu_fig1b$assigned_label == "Microglia" ~ "Microglia",
  seu_fig1b$assigned_label == "Neuron" ~ "Neuron",
  TRUE ~ NA_character_
)

seu_fig1b$Figure1B_group <- factor(
  seu_fig1b$Figure1B_group,
  levels = c("Astrocytes", "AEMS", "Microglia", "Neuron")
)

table(seu_fig1b$Figure1B_group, useNA = "ifany")

# --------------------------------------------------
# 3. Define colors
# --------------------------------------------------
fig1b_cols <- c(
  "Astrocytes" = "#9ACD32",
  "AEMS"       = "forestgreen",
  "Microglia"  = "orange3",
  "Neuron"     = "pink"
)

# --------------------------------------------------
# 4. Score canonical astrocyte RNA markers with UCell
# --------------------------------------------------
DefaultAssay(seu_fig1b) <- "RNA"

astro_markers <- c("GFAP", "SLC1A2", "AQP4", "ALDH1A1")
astro_markers_present <- intersect(astro_markers, rownames(seu_fig1b[["RNA"]]))

cat("Astrocyte RNA markers present:", length(astro_markers_present), "of", length(astro_markers), "\n")
print(astro_markers_present)

seu_fig1b <- AddModuleScore_UCell(
  seu_fig1b,
  features = list(AstrocyteMarkers_RNA = astro_markers_present),
  assay = "RNA"
)

summary(seu_fig1b$AstrocyteMarkers_RNA_UCell)

# --------------------------------------------------
# 5. Check that GeneProfile_GA_UCell already exists
# --------------------------------------------------
if (!"GeneProfile_GA_UCell" %in% colnames(seu_fig1b@meta.data)) {
  stop("GeneProfile_GA_UCell is not in seu_CTRL/seu_fig1b metadata. Score it first before running this block.")
}

summary(seu_fig1b$GeneProfile_GA_UCell)

# --------------------------------------------------
# 6. Make violin plot 1: canonical astrocyte RNA UCell
# --------------------------------------------------
p_astro_rna <- ggplot(
  seu_fig1b@meta.data,
  aes(x = Figure1B_group, y = AstrocyteMarkers_RNA_UCell, fill = Figure1B_group)
) +
  geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = fig1b_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Astrocyte Score RNA"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

p_astro_rna

ggsave(
  filename = "Figure_2B_canonical_astrocyte_RNA_UCell_violin.png",
  plot = p_astro_rna,
  width = 4.5,
  height = 4,
  dpi = 600,
  bg = "white"
)

# --------------------------------------------------
# 7. Make violin plot 2: Lee et al / human rna UCell
# --------------------------------------------------
P_rna <- ggplot(
  seu_fig1b@meta.data,
  aes(x = Figure1B_group, y = GeneProfile_RNA_UCell, fill = Figure1B_group)
) +
  geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = fig1b_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Lee et al. Memory Profile RNA "
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

P_rna

ggsave(
  filename = "Figure_2B_memory_profile_RNA_UCell_violin.png",
  plot = P_rna,
  width = 4.5,
  height = 4,
  dpi = 600,
  bg = "white"
)

# --------------------------------------------------
# 8. Optional quick checks
# --------------------------------------------------
table(seu_fig1b$Figure1B_group)

aggregate(
  AstrocyteMarkers_RNA_UCell ~ Figure1B_group,
  data = seu_fig1b@meta.data,
  summary
)

aggregate(
  GeneProfile_GA_UCell ~ Figure1B_group,
  data = seu_fig1b@meta.data,
  summary
)


# ==========================
# Figure 2C conocial Astrocyte and lee et al memory profile ATAC
# =========================


# --------------------------------------------------
# 4. Score canonical astrocyte GeneActivity markers with UCell
# --------------------------------------------------
DefaultAssay(seu_fig1b) <- "GeneActivity"

astro_markers <- c("GFAP", "SLC1A2", "AQP4", "ALDH1A1")
astro_markers_present_GA <- intersect(astro_markers, rownames(seu_fig1b[["GeneActivity"]]))

cat("Astrocyte GeneActivity markers present:", length(astro_markers_present_GA), "of", length(astro_markers), "\n")
print(astro_markers_present_GA)

seu_fig1b <- AddModuleScore_UCell(
  seu_fig1b,
  features = list(AstrocyteMarkers_GA = astro_markers_present_GA),
  assay = "GeneActivity"
)

summary(seu_fig1b$AstrocyteMarkers_GA_UCell)

# --------------------------------------------------
# 5. Check that Lee et al memory profile GeneActivity UCell already exists
# --------------------------------------------------
if (!"GeneProfile_GA_UCell" %in% colnames(seu_fig1b@meta.data)) {
  stop("GeneProfile_GA_UCell is not in seu_fig1b metadata. Score it first before running this block.")
}

summary(seu_fig1b$GeneProfile_GA_UCell)

# --------------------------------------------------
# 6. Make violin plot 1: canonical astrocyte GeneActivity UCell
# --------------------------------------------------
p_astro_ga <- ggplot(
  seu_fig1b@meta.data,
  aes(x = Figure1B_group, y = AstrocyteMarkers_GA_UCell, fill = Figure1B_group)
) +
  geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = fig1b_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Astrocyte Score GA"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

p_astro_ga

ggsave(
  filename = "Figure_2c_canonical_astrocyte_GA_UCell_violin.png",
  plot = p_astro_ga,
  width = 4.5,
  height = 4,
  dpi = 600,
  bg = "white"
)

# --------------------------------------------------
# 7. Make violin plot 2: Lee et al. memory profile GeneActivity UCell
# --------------------------------------------------
p_memory_ga <- ggplot(
  seu_fig1b@meta.data,
  aes(x = Figure1B_group, y = GeneProfile_GA_UCell, fill = Figure1B_group)
) +
  geom_violin(trim = TRUE, scale = "width", color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.25) +
  scale_fill_manual(values = fig1b_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Lee et al. Memory Profile GA"
  ) +
  theme_classic(base_size = 18) +
  theme(
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

p_memory_ga

ggsave(
  filename = "Figure_2c_memory_profile_GA_UCell_violin.png",
  plot = p_memory_ga,
  width = 4.5,
  height = 4,
  dpi = 600,
  bg = "white"
)

# --------------------------------------------------
# 8. Optional quick checks
# --------------------------------------------------
table(seu_fig1b$Figure1B_group)

aggregate(
  AstrocyteMarkers_GA_UCell ~ Figure1B_group,
  data = seu_fig1b@meta.data,
  summary
)

aggregate(
  GeneProfile_GA_UCell ~ Figure1B_group,
  data = seu_fig1b@meta.data,
  summary
)


# =================================
# # Figure 2G Python code to add to github
# ===============================

# =================================
# # Figure 2H heatmap of significant DEGs between AEMs and Astrocyte in the assgined label column 
# ===============================


# --------------------------------------------------
# Setup identities
# --------------------------------------------------
seu_sub$Figure1B_group <- factor(
  seu_sub$Figure1B_group,
  levels = c("Astrocytes", "AEMS")
)
Idents(seu_sub) <- seu_sub$Figure1B_group

# --------------------------------------------------
# Differential expression
# --------------------------------------------------
deg <- FindMarkers(
  seu_sub,
  ident.1 = "AEMS",
  ident.2 = "Astrocytes",
  logfc.threshold = 0,
  min.pct = 0,
  test.use = "wilcox"
)

deg$gene <- rownames(deg)

# Avoid Inf values in plotting
deg$plot_p <- ifelse(deg$p_val_adj == 0, .Machine$double.xmin, deg$p_val_adj)

# --------------------------------------------------
# Define coloring
# --------------------------------------------------
deg$highlight <- "Not significant"
deg$highlight[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0] <- "Up in AEMS"
deg$highlight[deg$p_val_adj < 0.05 & deg$avg_log2FC < 0] <- "Down in AEMS"

deg$highlight <- factor(
  deg$highlight,
  levels = c("Not significant", "Down in AEMS", "Up in AEMS")
)

# --------------------------------------------------
# Volcano plot (NO LABELS + BIG TEXT)
# --------------------------------------------------
library(ggplot2)

volcano_plot <- ggplot(deg, aes(x = avg_log2FC, y = -log10(plot_p))) +
  
  geom_point(aes(color = highlight), alpha = 0.8, size = 2) +
  
  scale_color_manual(values = c(
    "Up in AEMS" = "#d62728",
    "Down in AEMS" = "#1f77b4",
    "Not significant" = "grey85"
  )) +
  
  theme_classic(base_size = 18) +  # 🔥 increases all base text
  
  labs(
    title = "Volcano Plot: AEMS vs Astrocytes",
    subtitle = "Red = up in AEMS | Blue = down in AEMS",
    x = "Log2 Fold Change (left ← Astrocyte | AEMS → right)",
    y = "-log10(adj p-value)"
  ) +
  
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    plot.subtitle = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 16),
    legend.title = element_blank()
  )

ggsave("volcano_plot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)



# =================================
# # Figure 2I Python code to add to github
# ===============================








