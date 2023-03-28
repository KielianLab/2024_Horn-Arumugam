## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure 1 in Horn et al., 2023
## Created: 2023-01-23
## Updated: 2023-03-28

## Setting up environment ----
# Load in packages
packages <- c(
  "tidyverse",
  "Seurat",
  "patchwork",
  "plotly",
  "ggsci",
  "readxl",
  "data.table",
  "msigdbr",
  "here",
  "ComplexHeatmap",
  "grid",
  "circlize",
  "reshape2",
  "ggpubr",
  "ggrepel",
  "ggvenn"
)

# Check overlap between existing & new packages; install new (if any)
new.packages <- packages[!(packages %in% installed.packages()[, "Package"])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(packages, library, character.only = TRUE))

rm(
  packages,
  new.packages
)

gc()

# Set options
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(12345)

# Save session info
writeLines(capture.output(sessionInfo()), here("sessionInfo.txt"))

## Figure 1B ----
# Load in data
integrated <- readRDS(here("Figure 1", "Figure 1.rds"))
dimred <- readRDS(here("Figure 1", "Figure 1 dim red.rds")) # Was calculated outside of the normal Seurat workflow & used for downstream analysis

# Add custom dimensional reduction
integrated[["umap"]] <- CreateDimReducObject(embeddings = dimred, key = "UMAP_", assay = DefaultAssay(integrated))

# Create plot
p <- DimPlot(
  integrated,
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 1.5,
  label.color = "black") +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank()
)

p <- p + theme(legend.position = "bottom")

# Save plot
ggsave(
  "Figure 1B.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure 1")
)

# Clear the environment
rm(list = ls())
gc()

## Figure 1C ----
# Load in data
integrated <- readRDS(here("Figure 1", "Figure 1.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

# Calculate cluster-averaged normalized counts
avg <- AverageExpression(integrated, return.seurat = TRUE)
avg.norm_counts <- avg@assays$RNA@data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = "gene")

# Filter averaged normalized counts for genes of interest & convert back into a matrix
gene_list <- c(
  "Il1b",
  "Clec4e",
  "Junb",
  "Ctsd",
  "Wfdc17",
  "Il1f9",
  "Pla2g7",
  "Arg2",
  "Cd84",
  "Lcn2",
  "Prdx5",
  "Ngp",
  "Camp",
  "Ltf",
  "Arhgdib",
  "Anxa1",
  "Plbd1",
  "Tkt",
  "Aldh2",
  "Ly6c2",
  "Adpgk",
  "Cd177"
)

heatmap.data <- avg.norm_counts %>%
  filter(gene %in% gene_list)
heatmap.data <- column_to_rownames(heatmap.data, var = "gene")
heatmap.data <- as.matrix(heatmap.data)

## Create vectors for column/row annotations
col_anno <- c(
  "G-MDSC",
  "G-MDSC",
  "G-MDSC",
  "G-MDSC",
  "G-MDSC",
  "G-MDSC",
  "G-MDSC",
  
  "PMN",
  "PMN",
  "PMN"
)

row_order <- gene_list
heatmap.data <- heatmap.data[row_order, ]

row_anno <- c(
  "G-MDSC genes",
  "G-MDSC genes",
  "G-MDSC genes",
  "G-MDSC genes",
  "G-MDSC genes",
  "G-MDSC genes",
  "G-MDSC genes",
  "G-MDSC genes",
  "G-MDSC genes",
  
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes",
  "PMN genes"
)

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c(0, htmp_range), c("white", "red"))

# Create & save heatmap
png(
  here("Figure 1", "Figure 1C.png"),
  width = 400,
  height = 550
)

ht <- Heatmap(
  heatmap.data,
  name = "Average\nNormalized\nCounts",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  column_split = col_anno,
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  column_gap = unit(5, "mm"),
  row_split = row_anno,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  row_gap = unit(5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  row_names_gp = gpar(fontface = "bold.italic"),
  row_title_gp = gpar(fontface = "bold"),
  column_names_gp = gpar(fontface = "bold"),
  column_title_gp = gpar(fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", heatmap.data[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)

# Adjust heatmap paddings
draw(ht, padding = unit(c(5, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())

gc()
