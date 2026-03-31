#===========================
# Figure 1 Honours
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

setwd("D:/Finn/Honours/figures/Figure_1")
outdir <- getwd()
set.seed(1)


# ==========================
# laoding Data and listing properties
# ==========================

seu <- readRDS("D:/CNS_immune_memory/archive/2025_13_2025/atac_on_NMF/ATAC_info_addded_NMF_QC_PV_glia_reactivity.rds")

seu

head(seu@meta.data)

meta <- seu@meta.data

# 2. Identify non-numeric columns
cat_cols <- colnames(meta)[!sapply(meta, is.numeric) & !sapply(meta, is.logical)]

# 3. Print the factors/unique values for those columns
for (col in cat_cols) {
  cat("\n---", col, "---\n")
  # Get unique values, sorted
  unique_vals <- sort(unique(as.character(meta[[col]])))
  print(unique_vals)
}

# Get Metadata Summary
summary(seu[[]])

# List Assay names and data types
Assays(seu)

# Get Dimensional Reductions information
Reductions(seu)

# If clusters exist, get cluster sizes
table(Idents(seu))

# ==========================
# helper functions
# ==========================

# Clean/save theme (BIG TEXT, no titles)
theme_qc <- function() {
  theme_classic(base_size = 22) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 18),
      axis.text.y = element_text(size = 18),
      axis.title.x = element_blank(),  # REMOVE X LABEL
      axis.title.y = element_text(size = 20, face = "bold"),
      plot.title = element_blank(),    # REMOVE TITLE
      legend.title = element_text(size = 18, face = "bold"),
      legend.text = element_text(size = 16)
    )
}

save_plot <- function(p, filename, width = 6.48, height = 6, dpi = 400) {
  ggsave(
    filename = file.path(outdir, filename),
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
}

# For nice patient ordering: all CTRL first, then MS
get_patient_order <- function(meta) {
  patient_df <- meta %>%
    distinct(orig.ident, condition) %>%
    mutate(
      sample_num = suppressWarnings(as.numeric(str_extract(orig.ident, "\\d+"))),
      condition = factor(condition, levels = c("CTRL", "MS"))
    ) %>%
    arrange(condition, sample_num, orig.ident)
  
  patient_df$orig.ident
}

# Manual color palette: cool controls, warm MS
make_patient_palette <- function(meta) {
  patient_df <- meta %>%
    distinct(orig.ident, condition) %>%
    mutate(
      sample_num = suppressWarnings(as.numeric(str_extract(orig.ident, "\\d+"))),
      condition = factor(condition, levels = c("CTRL", "MS"))
    ) %>%
    arrange(condition, sample_num, orig.ident)
  
  ctrl_ids <- patient_df$orig.ident[patient_df$condition == "CTRL"]
  ms_ids   <- patient_df$orig.ident[patient_df$condition == "MS"]
  
  ctrl_cols <- c("#4C78A8", "#72B7B2", "#54A24B", "#A0CBE8", "#8CD17D")
  ms_cols   <- c("#E45756", "#F58518", "#FF9DA6", "#B279A2", "#D37295")
  
  ctrl_cols <- ctrl_cols[seq_along(ctrl_ids)]
  ms_cols   <- ms_cols[seq_along(ms_ids)]
  
  pal <- c(ctrl_cols, ms_cols)
  names(pal) <- c(ctrl_ids, ms_ids)
  pal
}

# Make one violin plot per metric
make_qc_violin <- function(seu, metric, ylab = NULL, palette = NULL) {
  meta <- seu@meta.data
  
  if (!metric %in% colnames(meta)) {
    stop(paste0("Metric '", metric, "' not found in seu@meta.data"))
  }
  
  meta$orig.ident <- factor(meta$orig.ident, levels = get_patient_order(meta))
  
  p <- ggplot(meta, aes(x = orig.ident, y = .data[[metric]], fill = orig.ident)) +
    geom_violin(scale = "width", trim = TRUE, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = palette, drop = FALSE) +
    labs(
      y = ifelse(is.null(ylab), metric, ylab)
    ) +
    theme_qc() +
    theme(legend.position = "none")
  
  p
}

# Saving Metadata 
meta <- seu@meta.data

# Patient ordering + colors
meta <- seu@meta.data
patient_order <- get_patient_order(meta)
patient_palette <- make_patient_palette(meta)

seu$orig.ident <- factor(seu$orig.ident, levels = patient_order)
seu$condition  <- factor(seu$condition, levels = c("CTRL", "MS"))

# ==========================
# Figure 1C nCountRNA
# =========================
p_1C <- make_qc_violin(
  seu, metric = "nCount_RNA",
  ylab = "nCount RNA",
  palette = patient_palette
)
save_plot(p_1C, "Figure_1C_nCountRNA_violin.png")

# ==========================
# Figure 1D nFeatureRNA
# =========================
p_1D <- make_qc_violin(
  seu, metric = "nFeature_RNA",
  ylab = "nFeature RNA",
  palette = patient_palette
)
save_plot(p_1D, "Figure_1D_nFeatureRNA_violin.png")

# ==========================
# Figure 1E Percent.mt
# =========================
p_1E <- make_qc_violin(
  seu, metric = "percent.mt",
  ylab = "Percent mitochondrial",
  palette = patient_palette
)
save_plot(p_1E, "Figure_1E_percent_mt_violin.png")

# ==========================
# Figure 1F Percent.ribosomal
# =========================
p_1F <- make_qc_violin(
  seu, metric = "percent.ribosomal",
  ylab = "Percent ribosomal",
  palette = patient_palette
)
save_plot(p_1F, "Figure_1F_percent_ribosomal_violin.png")

# ==========================
# Figure 1G nCountATAC
# =========================
p_1G <- make_qc_violin(
  seu, metric = "nCount_ATAC",
  ylab = "nCount ATAC",
  palette = patient_palette
)
save_plot(p_1G, "Figure_1G_nCountATAC_violin.png")

# ==========================
# Figure 1H nFeatureATAC
# =========================
p_1H <- make_qc_violin(
  seu, metric = "nFeature_ATAC",
  ylab = "nFeature ATAC",
  palette = patient_palette
)
save_plot(p_1H, "Figure_1H_nFeatureATAC_violin.png")

# ==========================
# Figure 1I Nucleosome signal
# =========================
p_1I <- make_qc_violin(
  seu, metric = "nucleosome_signal",
  ylab = "Nucleosome signal",
  palette = patient_palette
)
save_plot(p_1I, "Figure_1I_nucleosome_signal_violin.png")

# ==========================
# Figure 1J TSS.enrichment
# =========================
p_1J <- make_qc_violin(
  seu, metric = "TSS.enrichment",
  ylab = "TSS enrichment",
  palette = patient_palette
)
save_plot(p_1J, "Figure_1J_TSS_enrichment_violin.png")

# ==========================
# Figure 1K NMF Reduced UMAP
# =========================


# 1. Check the reduction

seu@reductions[["MPsignatures"]]


# 2. Build graph from MPsignatures dimensions

seu <- FindNeighbors(
  object = seu,
  reduction = "MPsignatures",
  dims = 1:10,
  graph.name = "MP_snn"
)


# 3. Cluster cells based on MPsignatures

seu <- FindClusters(
  object = seu,
  graph.name = "MP_snn",
  resolution = 0.6,
  algorithm = 1
)

# By default this writes to seurat_clusters
table(seu$seurat_clusters)


# 4. Save these as a separate metadata column

seu$MP_clusters <- seu$seurat_clusters

# optional: make them factors in numeric order
seu$MP_clusters <- factor(
  seu$MP_clusters,
  levels = sort(unique(as.numeric(as.character(seu$MP_clusters))))
)


# 5. Plot clusters on your existing UMAP

p_MP_clusters <- DimPlot(
  object = seu,
  reduction = "umap_MP",
  group.by = "MP_clusters",
  label = TRUE,
  repel = TRUE,
  pt.size = 2,
  raster = FALSE
) +
  labs(title = "Clusters computed from MPsignatures") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )

p_MP_clusters



# 6. Save

ggsave(
  filename = "Figure_1K_MPsignature_clusters.png",
  plot = p_MP_clusters,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)

# ------------------------
# 1. Define new cluster -> cell type mapping

# for 0.6 Res  
cluster2celltype <- c(
  "0"  = "Oligo",
  "1"  = "Oligo",
  "2"  = "Astrocyte",
  "3"  = "Oligo",
  "4"  = "Oligo",
  "5"  = "Oligo",
  "6"  = "Oligo",
  "7"  = "Astrocyte",
  "8"  = "Oligo",
  "9"  = "Microglia",
  "10" = "Neuron",
  "11" = "Oligo",
  "12" = "Astrocyte",
  "13" = "Unknown",
  "14" = "Ependymal",
  "15" = "OPC",
  "16" = "Astrocyte"
)


# 2. Filter out the OPC cluster at the start


# # remove cluster 15 before the rest of the analysis
# seu_noOPC <- subset(seu, subset = MP_clusters != "15")
# 
# seu <- seu_noOPC

# check
table(seu$MP_clusters)

# Pull MP cluster IDs as character
cluster_ids <- as.character(seu$MP_clusters)

# Map clusters to broad cell types
mapped_celltypes <- cluster2celltype[cluster_ids]

# Name by cell barcode so it aligns properly
mapped_celltypes <- unname(mapped_celltypes)
names(mapped_celltypes) <- colnames(seu)

# Add as metadata
seu$MP_cell_type <- mapped_celltypes

# Check mapping
table(seu$MP_clusters, seu$MP_cell_type, useNA = "ifany")


# 3. Define cell-type colors


celltype_cols <- c(
  "Oligo" = "#440154FF",  # dark purple
  "Astrocyte"       = "#9ACD32",
  "Microglia"       = "orange3",
  "Neuron"          = "pink",
  "Ependymal"       = "purple",
  "Unknown"         = "grey"
)

# Put in a nice legend order
seu$MP_cell_type <- factor(
  seu$MP_cell_type,
  levels = c(
    "Oligo",
    "Astrocyte",
    "Ependymal",
    "Microglia",
    "Neuron",
    "Unknown"
  )
)


# 4. Plot cell type on UMAP


p_MP_celltype <- DimPlot(
  object = seu,
  reduction = "umap_MP",
  group.by = "MP_cell_type",
  cols = celltype_cols,
  label = FALSE,
  repel = TRUE,
  pt.size = 0.8,
  raster = FALSE
) +
  labs(title = "MP-signature clusters mapped to broad cell types") +
  theme_classic(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 12)
  )

p_MP_celltype

ggsave(
  filename = "Figure_1K_MPsignature_broad_celltypes.png",
  plot = p_MP_celltype,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)

# ==========================
# Figure 1L Cell type markers dotplot
# ==========================

DefaultAssay(seu) <- "RNA"

# ordered vector instead of named list
marker_vec <- c(
  "MOBP", "PLP1", "MOG",        # Oligo
  "GFAP", "AQP4", "ALDH1L1",    # Astrocyte
  "FOXJ1", "PIFO", "TPPP3",     # Ependymal
  "P2RY12", "CX3CR1", "CSF1R",  # Microglia
  "RBFOX3", "SNAP25", "SYT1",   # Neuron
  "SOX2", "VIM", "HES1"         # Unknown
)

# make sure cell types are ordered nicely
seu$MP_cell_type <- factor(
  seu$MP_cell_type,
  levels = c(
    "Oligo",
    "Astrocyte",
    "Ependymal",
    "Microglia",
    "Neuron",
    "Unknown"
  )
)

p_1L_dot <- DotPlot(
  object = seu,
  features = marker_vec,
  group.by = "MP_cell_type",
  assay = "RNA",
  dot.scale = 22
) + RotatedAxis() +
  scale_color_gradient2(low = "black", mid = "purple3", high = "red2") + 
  labs(
    x = "Marker genes",
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 30 , face = "bold"),
    axis.text.y = element_text(size = 30 , face = "bold"),
    axis.title = element_text(size = 0),
    legend.title = element_text(size = 30, face = "bold"),
    legend.text = element_text(size = 30),
    legend.key.size = unit(1.5, "cm")
  )

p_1L_dot

ggsave(
  filename = "Figure_1L_celltype_marker_dotplot.png",
  plot = p_1L_dot,
  width = 25.92,
  height = 12,
  dpi = 600,
  bg = "white"
)

# ==========================
# Figure 1M per patient stack bar
# ==========================


# 1. Define colors

celltype_cols <- c(
  "Oligo" = "#440154FF",
  "Astrocyte"       = "#9ACD32",
  "Microglia"       = "orange3",
  "Neuron"          = "pink",
  "Ependymal"       = "purple",
  "Unknown"         = "grey"
)


# 2. Make sure labels exist

table(seu$MP_cell_type, useNA = "ifany")
table(seu$orig.ident, useNA = "ifany")

# optional: replace NA cell types with Unknown
seu$MP_cell_type <- as.character(seu$MP_cell_type)
seu$MP_cell_type[is.na(seu$MP_cell_type)] <- "Unknown"

# put in desired legend/order
seu$MP_cell_type <- factor(
  seu$MP_cell_type,
  levels = c(
    "Oligo",
    "Astrocyte",
    "Ependymal",
    "Microglia",
    "Neuron",
    "Unknown"
  )
)


# 3. Compute per-patient proportions

stack_df <- seu@meta.data %>%
  as.data.frame() %>%
  count(orig.ident, MP_cell_type, name = "n_cells") %>%
  group_by(orig.ident) %>%
  mutate(prop = n_cells / sum(n_cells)) %>%
  ungroup()


# 4. Optional patient ordering

# If you already have condition in metadata and want CTRL left, MS right:
if ("condition" %in% colnames(seu@meta.data)) {
  patient_order <- seu@meta.data %>%
    as.data.frame() %>%
    distinct(orig.ident, condition) %>%
    mutate(condition = factor(condition, levels = c("CTRL", "MS"))) %>%
    arrange(condition, orig.ident) %>%
    pull(orig.ident)
  
  stack_df$orig.ident <- factor(stack_df$orig.ident, levels = patient_order)
} else {
  stack_df$orig.ident <- factor(
    stack_df$orig.ident,
    levels = unique(stack_df$orig.ident)
  )
}


# 5. Plot stacked bar plot

p_1M <- ggplot(stack_df, aes(x = orig.ident, y = prop, fill = MP_cell_type)) +
  geom_col(width = 0.8, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = celltype_cols, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "Figure 1M. Cell-type composition per patient",
    x = "Patient",
    y = "Proportion of cells",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11)
  )

p_1M


# 6. Save

ggsave(
  filename = "Figure_1M_per_patient_stacked_barplot.png",
  plot = p_1M,
  width = 11,
  height = 7,
  dpi = 600,
  bg = "white"
)

saveRDS(seu, file = "MS_CTRL_Whole_OBJECT.RDS")
