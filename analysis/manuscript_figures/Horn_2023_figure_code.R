## Made by Christopher M. Horn, PhD
## Kielian Lab data
## Code for figures in Horn et al., 2024
## Created: 2023-01-23
## Updated: 2024-02-02

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
  "ggvenn",
  "svglite"
)

# Check overlap between existing & new packages; install new (if any)
installed <- suppressWarnings(sapply(packages, require, character.only = TRUE))
BiocManager::install(names(installed)[!installed])
install.packages(names(installed)[!installed])
invisible(lapply(packages, library, character.only = TRUE))

rm(
  packages,
  installed
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
  "svg__Figure 1B.svg",
  plot = p,
  device = "svg",
  width = 3587,
  height = 2799,
  unit = "px",
  dpi = 600,
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
avg.norm_counts <- avg[["RNA"]]$data
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
svglite::svglite(
  here("Figure 1", "svg__Figure 1C.svg"),
  width = 6,
  height = 6
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

## Figure 2A ----
# Read in pathway data & select only granulocytic data
integrated <- readRDS(here("Figure 2", "Figure 2.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

avg <- AverageExpression(integrated, return.seurat = TRUE)
avg.norm_counts <- avg[["RNA"]]$data
avg.norm_counts <- as.data.frame(avg.norm_counts)
avg.norm_counts <- rownames_to_column(avg.norm_counts, var = "gene")

gene_list <- c(
  "Gapdh",
  "Pkm",
  "Pgam1",
  "Eno1",
  "Pgk1",
  "Tpi1",
  "Gpi1",
  "Taldo1",
  "Pgd",
  "Tkt",
  "G6pdx",
  "Pgls"
)

# Create matrix of NES values, fill with 0 if not found
heatmap.data <- avg.norm_counts %>%
  filter(gene %in% gene_list)
heatmap.data <- column_to_rownames(heatmap.data, var = "gene")
heatmap.data <- as.matrix(heatmap.data)

# Set custom row order
col_order <- c(
  "G1",
  "G2",
  "G3",
  "G4",
  "G5",
  "G6",
  "G7",
  "G8",
  "G9",
  "G10"
)

row_order <- gene_list
heatmap.data <- heatmap.data[row_order, col_order]

# Make the second heatmap match the row & column order of the first
col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

# Create a vector to split the heatmap according to cluster
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

row_anno <- c(
  "Glycolysis genes",
  "Glycolysis genes",
  "Glycolysis genes",
  "Glycolysis genes",
  "Glycolysis genes",
  "Glycolysis genes",
  "Glycolysis genes",
  
  "PPP genes",
  "PPP genes",
  "PPP genes",
  "PPP genes",
  "PPP genes"
)

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c( 0, htmp_range), c("white", "red"))

# Create heatmap
svglite::svglite(
  here("Figure 2", "svg__Figure 2A.svg"),
  width = 6,
  height = 5.5
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

# Adjust heatmap padding
draw(ht, padding = unit(c(5, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure 2B ----
# Read in pathway data & combine
h <- read_xlsx(here("Figure 2", "Figure 2B.xlsx"), sheet = "H pathways")
c2 <- read_xlsx(here("Figure 2", "Figure 2B.xlsx"), sheet = "C2 pathways")

pathways <- rbind(
  h,
  c2
)

pathway_list <- c(
  "HALLMARK_GLYCOLYSIS",
  "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
  "HALLMARK_HYPOXIA",
  "SEMENZA_HIF1_TARGETS",
  
  "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
  "REACTOME_PENTOSE_PHOSPHATE_PATHWAY"
)

# Create matrix of NES values, fill with 0 if not found
heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- dcast(heatmap.data, pathway ~ cluster, value.var = "NES", fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = "pathway")
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- t(heatmap.data)

# Create another matrix of adjusted p-values, fill with 1 if not found
heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- dcast(heatmap2.data, pathway ~ cluster, value.var = "padj", fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = "pathway")
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- t(heatmap2.data)

supp_mat <- matrix(
  0,
  nrow = 4,
  ncol = 6,
  dimnames = list(
    c("G5", "G6", "G7", "G8"),
    c(
      "HALLMARK_GLYCOLYSIS",
      "HALLMARK_HYPOXIA",
      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
      "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
      "REACTOME_PENTOSE_PHOSPHATE_PATHWAY",
      "SEMENZA_HIF1_TARGETS"
      )
  )
)

supp_mat2 <- matrix(
  1,
  nrow = 4,
  ncol = 6,
  dimnames = list(
    c("G5", "G6", "G7", "G8"),
    c(
      "HALLMARK_GLYCOLYSIS",
      "HALLMARK_HYPOXIA",
      "KEGG_GLYCOLYSIS_GLUCONEOGENESIS",
      "KEGG_PENTOSE_PHOSPHATE_PATHWAY",
      "REACTOME_PENTOSE_PHOSPHATE_PATHWAY",
      "SEMENZA_HIF1_TARGETS"
    )
  )
)

heatmap.data <- rbind(heatmap.data, supp_mat)
heatmap2.data <- rbind(heatmap2.data, supp_mat2)

# Set custom row order
row_order <- c(
  "G1",
  "G2",
  "G3",
  "G4",
  "G5",
  "G6",
  "G7",
  
  "G8",
  "G9",
  "G10"
)

heatmap.data <- heatmap.data[row_order, ]

# Make the second heatmap match the row & column order of the first
col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

# Create a vector to split the heatmap according to cluster
row_anno <- c(
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

# Change formatting of pathway names from the defaults
col_names <- c(
  "Hallmark: Glycolysis",
  "Hallmark: Hypoxia",
  "KEGG: Glycolysis & Gluconeogenesis",
  "KEGG: Pentose Phosphate Pathway",
  "Reactome: Pentose Phosphate Pathway",
  "Semenza HIF-1a Targets"
)

colnames(heatmap.data) <- col_names

heatmap.data <- heatmap.data[, colnames(heatmap.data) != "KEGG: Pentose Phosphate Pathway"]
heatmap2.data <- heatmap2.data[, colnames(heatmap2.data) != "KEGG_PENTOSE_PHOSPHATE_PATHWAY"]

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c(0, htmp_range), c("white", "red"))

# Create heatmap
svglite::svglite(
  here("Figure 2", "svg__Figure 2B.svg"),
  width = 3,
  height = 5.5
)

ht <- Heatmap(
  heatmap.data,
  name = "NES",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_gap = unit(5, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_anno,
  row_gap = unit(5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  row_names_gp = gpar(fontface = "bold"),
  row_title_gp = gpar(fontface = "bold"),
  column_names_gp = gpar(fontface = "bold", fontsize = 10),
  column_title_gp = gpar(fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(heatmap2.data[i, j] < 0.0001) {
      grid.text("****", x, y)
    } else if(heatmap2.data[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(heatmap2.data[i, j] < 0.01) {
      grid.text("**", x, y)
    } else if(heatmap2.data[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
)

# Adjust heatmap padding
draw(ht, padding = unit(c(10, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure 2C ----
# Load in data
compass_results <- read.csv(here("Figure 2", "Figure 2C.csv"))

aam <- c(
  "Alanine and aspartate metabolism",
  "Arginine and Proline Metabolism",
  "beta-Alanine metabolism",
  "Cysteine Metabolism",
  "D-alanine metabolism",
  "Folate metabolism",
  "Glutamate metabolism",
  "Glycine, serine, alanine and threonine metabolism",
  "Histidine metabolism",
  "Lysine metabolism",
  "Methionine and cysteine metabolism",
  "Taurine and hypotaurine metabolism",
  "Tryptophan metabolism",
  "Tyrosine metabolism",
  "Urea cycle",
  "Valine, leucine, and isoleucine metabolism"
)

p <- ggplot(compass_results, aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "gray", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  annotate("text", x = 2.75, y = 2, label = "G-MDSC") +
  annotate("text", x = -2.75, y = 2, label = "PMN") +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(bold(-log[10](BH-adjusted~Wilcoxon~rank~sum~p)))) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"))

p1 <- ggplot(compass_results %>% filter(subsystem =="Glycolysis/gluconeogenesis"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#BC3C29FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.5, y = 2.5, label = "G-MDSC") +
  annotate("text", x = -2.5, y = 2.5, label = "PMN") +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Glycolysis") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p2 <- ggplot(compass_results %>% filter(subsystem =="Citric acid cycle"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#0072B5FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.5, y = 2.5, label = "G-MDSC") +
  annotate("text", x = -2.5, y = 2.5, label = "PMN") +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("TCA cycle") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p3 <- ggplot(compass_results %>% filter(subsystem %in% aam), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#E18727FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.5, y = 2.5, label = "G-MDSC") +
  annotate("text", x = -2.5, y = 2.5, label = "PMN") +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Amino acid metabolism") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p4 <- ggplot(compass_results %>% filter(subsystem =="Fatty acid oxidation"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#20854EFF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.5, y = 2.5, label = "G-MDSC") +
  annotate("text", x = -2.5, y = 2.5, label = "PMN") +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Fatty acid oxidation") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

figure <- ggarrange(
  p1 + rremove("ylab") + rremove("xlab"), 
  p2 + rremove("ylab") + rremove("xlab"), 
  p3 + rremove("ylab") + rremove("xlab"), 
  p4 + rremove("ylab") + rremove("xlab"),
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
)

p5 <- annotate_figure(
  figure,
  left = textGrob(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), rot = 90, vjust = 0.5, gp = gpar(cex = 1.3)),
  bottom = textGrob("Cohen's D", gp = gpar(cex = 1.3)))

final <- p | figure

ggsave(
  "svg__Figure 2C.svg",
  plot = final,
  device = "svg",
  path = here("Figure 2"),
  width = 14,
  height = 7,
  units = "in",
  dpi = 600
)

# Clear the environment
rm(list = ls())
gc()

## Figure 2D ----
# Read in data
compass.cor <- read.csv(here("Figure 2", "Figure 2D.csv")) %>%
  select(-c(X, label_point))

compass.cor <- compass.cor %>%
  mutate(label_repel = case_when(
    (neg_log_p >= -log10(0.05)) &
      (abs(patho_spearman) >= 0.5) &
      ((subsystem == "Glycolysis/gluconeogenesis") | (subsystem == "Pentose phosphate pathway"))
    ~ metadata_r_id,
    TRUE ~ "")) %>%
  filter(!row_number() %in% c(547, 551))

compass.cor <- compass.cor %>%
  mutate(label_point = case_when(
    subsystem == "Glycolysis/gluconeogenesis" ~ "red",
    subsystem == "Pentose phosphate pathway" ~ "blue",
    subsystem != "Glycolysis/gluconeogenesis" | subsystem != "Pentose phosphate pathway" ~ "gray"
  ))

compass.cor$label_point <- factor(x = compass.cor$label_point, levels = c("gray", "blue", "red"))

p1 <- ggplot(compass.cor %>% arrange(label_point), aes(x = signed_neg_log_p, y = patho_spearman)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(
    aes(fill = label_point),
    shape = 21,
    color = "black",
    size = 3
  ) +
  geom_label_repel(
    aes(label = label_repel),
    max.overlaps = Inf
  ) +
  annotate("label", x = 5.0, y = 1.0, label = "G-MDSC", fontface = "bold") +
  annotate("label", x = -5.0, y = -1.0, label = "PMN", fontface = "bold") +
  theme_classic() +
  labs(x = bquote(bold(Signed~-log[10](BH-adjusted~Wilcoxon~rank~sum~p))), y = "Spearman correlation with pathogenicity score") +
  ylim(c(-1, 1)) +
  scale_fill_manual(
    name = "",
    labels = c("Other", "Pentose Phosphate Pathway", "Glycolysis/Gluconeogenesis"),
    values = c("gray", "blue", "red")
  ) +
  theme(axis.title = element_text(face = "bold"), legend.position = "bottom")

p2_d <- compass.cor %>%
  filter(label_repel != "") %>%
  select(c(metadata_r_id, cohens_d, adjusted_pval, patho_spearman, subsystem)) %>%
  mutate(cohens_d = round(cohens_d, 4)) %>%
  mutate(patho_spearman = round(patho_spearman, 4)) %>%
  arrange(desc(cohens_d))

colnames(p2_d) <- c("", "Cohen's D", "Adjusted p", "Spearman Rho w/Pathogenicity", "Metabolic Subsystem")

p2 <- ggtexttable(p2_d, rows = NULL, theme = ttheme("light"))

final <- p1 / p2

ggsave(
  "svg__Figure 2D.svg",
  plot = final,
  device = "svg",
  units = "in",
  width = 8,
  height = 8,
  path = here("Figure 2"),
  dpi = 600
)

# Clear the environment
rm(list = ls())
gc()

## Figure 4A ----
#Load in data, normalize, scale, & subset
integrated <- readRDS(here("Figure 4", "Figure 4.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

new.cluster.ids <- c(
  "G1",
  "G2",
  "G3",
  "G4",
  "G5",
  "G6",
  "G7",
  "G8",
  "G9",
  "G10"
)

names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)

# Create plot
vln1 <- VlnPlot(
  integrated,
  features = "Sod2",
  idents = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10")
) +
  NoLegend() +
  scale_fill_npg() +
  theme(plot.title = element_text(face = "bold.italic"))

vln2 <- VlnPlot(
  integrated,
  features = "Nfe2l2",
  idents = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10")
) +
  NoLegend() +
  scale_fill_npg() +
  theme(plot.title = element_text(face = "bold.italic"))

vln3 <- VlnPlot(
  integrated,
  features = "Hmox1",
  idents = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9", "G10")
) +
  NoLegend() +
  scale_fill_npg() +
  theme(plot.title = element_text(face = "bold.italic"))

vln <- vln1 | vln2 | vln3

ggsave(
  "svg__Figure 4A.svg",
  plot = vln,
  device = "svg",
  path = here("Figure 4"),
  width = 14,
  height = 7,
  units = "in",
  dpi = 600
)

# Clear the environment
rm(list = ls())
gc()

## Figure 9A ----
# Load in data
integrated <- readRDS(here("Figure 9", "Figure 9.rds"))

# Create plot
integrated_umap.clus <- DimPlot(
  integrated,
  label = TRUE,
  label.box = TRUE,
  label.size = 6,
  repel = TRUE,
  pt.size = 1.5) +
  theme_classic() +
  NoLegend() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(face = "bold", size = 12)
)

ggsave(
  "svg__Figure 9A.svg",
  plot = integrated_umap.clus,
  device = "svg",
  units = "in",
  width = 8,
  height = 8,
  path = here("Figure 9"),
  dpi = 600
)

# Clear the environment
rm(list = ls())
gc()

## Figure 9B ----
# Load in data, normalize, scale, & subset
integrated <- readRDS(here("Figure 9", "Figure 9.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

new.cluster.ids <- c(
  "G1",
  "G2",
  "G3",
  "G4",
  "G5",
  "G6",
  "G7",
  "G8",
  "GMM",
  "L",
  "G9"
)

names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)

# Create plot
vln1 <- VlnPlot(
  integrated,
  features = "Ldha",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2")
) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

vln2 <- VlnPlot(
  integrated,
  features = "Eno1",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2")
) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

vln3 <- VlnPlot(
  integrated,
  features = "Gapdh",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2")
) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

vln4 <- VlnPlot(
  integrated,
  features = "Aldoa",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2")
) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

vln5 <- VlnPlot(
  integrated,
  features = "Pfkp",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2")
) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

vln6 <- VlnPlot(
  integrated,
  features = "Hk2",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2")
) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

vln1$data$split <- factor(x = vln1$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))
vln2$data$split <- factor(x = vln2$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))
vln3$data$split <- factor(x = vln3$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))
vln4$data$split <- factor(x = vln4$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))
vln5$data$split <- factor(x = vln5$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))
vln6$data$split <- factor(x = vln6$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))

vln1 <- vln1 + scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6)
vln2 <- vln2 + scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6)
vln3 <- vln3 + scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6)
vln4 <- vln4 + scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6)
vln5 <- vln5 + scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6)
vln6 <- vln6 + scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6)

vln <- (vln1 | vln2 | vln3) / (vln4 | vln5 | vln6) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(face = "bold", size = 12))

ggsave(
  "svg__Figure 9B.svg",
  plot = vln,
  device = "svg",
  units = "in",
  width = 10,
  height = 8,
  path = here("Figure 9"),
  dpi = 600
)

# Clear the environment
rm(list = ls())
gc()

## Figure 9C ----
# Read in pathway data & combine
h <- read_xlsx(here("Figure 9", "Figure 9C.xlsx"), sheet = "HM pathways")
c2 <- read_xlsx(here("Figure 9", "Figure 9C.xlsx"), sheet = "C2 pathways")

pathways <- rbind(
  h,
  c2
)

# Top 10 pathways by NES
pathways <- pathways %>%
  filter(str_detect(cluster_name, "Granulocytes [:digit:]")) %>%
  filter(cluster_name != "Granulocytes 9")

top_pathways <- pathways %>%
  arrange(desc(NES)) %>%
  group_by(cluster_num, cluster_subset) %>%
  do(head(., n = 5))

pathway_list <- unique(top_pathways$pathway)

heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- dcast(heatmap.data, pathway ~ cluster_name + cluster_subset, value.var = "NES", fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = "pathway")
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- dcast(heatmap2.data, pathway ~ cluster_name + cluster_subset, value.var = "padj", fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = "pathway")
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

row_order <- c(
  "Granulocytes 1_d3_null",
  "Granulocytes 1_d3_cre",
  "Granulocytes 1_d14_null",
  "Granulocytes 1_d14_cre",
  
  "Granulocytes 2_d3_null",
  "Granulocytes 2_d3_cre",
  "Granulocytes 2_d14_null",
  "Granulocytes 2_d14_cre",
  
  "Granulocytes 3_d3_null",
  "Granulocytes 3_d3_cre",
  "Granulocytes 3_d14_null",
  "Granulocytes 3_d14_cre",
  
  "Granulocytes 4_d3_null",
  "Granulocytes 4_d3_cre",
  "Granulocytes 4_d14_null",
  "Granulocytes 4_d14_cre",
  
  "Granulocytes 5_d3_null",
  "Granulocytes 5_d3_cre",
  "Granulocytes 5_d14_null",
  "Granulocytes 5_d14_cre",
  
  "Granulocytes 6_d3_null",
  "Granulocytes 6_d3_cre",
  "Granulocytes 6_d14_null",
  "Granulocytes 6_d14_cre",
  
  "Granulocytes 7_d3_null",
  "Granulocytes 7_d3_cre",
  "Granulocytes 7_d14_null",
  "Granulocytes 7_d14_cre",
  
  "Granulocytes 8_d3_null",
  "Granulocytes 8_d3_cre",
  "Granulocytes 8_d14_null",
  "Granulocytes 8_d14_cre"
)

heatmap.data <- heatmap.data[row_order, ]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

row_anno <- c(
  "G1",
  "G1",
  "G1",
  "G1",
  
  "G2",
  "G2",
  "G2",
  "G2",
  
  "G3",
  "G3",
  "G3",
  "G3",
  
  "G4",
  "G4",
  "G4",
  "G4",
  
  "G5",
  "G5",
  "G5",
  "G5",
  
  "G6",
  "G6",
  "G6",
  "G6",
  
  "G7",
  "G7",
  "G7",
  "G7",
  
  "G8",
  "G8",
  "G8",
  "G8"
)

row_names <- c(
  "G1_D3 WT",
  "G1_D3 HIF-1a cKO",
  "G1_D14 WT",
  "G1_D14 HIF-1a cKO",
  
  "G2_D3 WT",
  "G2_D3 HIF-1a cKO",
  "G2_D14 WT",
  "G2_D14 HIF-1a cKO",
  
  "G3_D3 WT",
  "G3_D3 HIF-1a cKO",
  "G3_D14 WT",
  "G3_D14 HIF-1a cKO",
  
  "G4_D3 WT",
  "G4_D3 HIF-1a cKO",
  "G4_D14 WT",
  "G4_D14 HIF-1a cKO",
  
  "G5_D3 WT",
  "G5_D3 HIF-1a cKO",
  "G5_D14 WT",
  "G5_D14 HIF-1a cKO",
  
  "G6_D3 WT",
  "G6_D3 HIF-1a cKO",
  "G6_D14 WT",
  "G6_D14 HIF-1a cKO",
  
  "G7_D3 WT",
  "G7_D3 HIF-1a cKO",
  "G7_D14 WT",
  "G7_D14 HIF-1a cKO",
  
  "G8_D3 WT",
  "G8_D3 HIF-1a cKO",
  "G8_D14 WT",
  "G8_D14 HIF-1a cKO"
)

col_names <- c(
  "Hallmark: Apical Junction",                                                                                                 
  "Hallmark: Complement",                                                                                                      
  "Hallmark: Glycolysis",                                                                                                      
  "Hallmark: Hypoxia",                                                                                                       
  "Hallmark: Inflammatory Response",                                                                                           
  "Hallmark: IFNa Response",                                                                                       
  "Hallmark: IFNg Response",                                                                                       
  "Hallmark: MTORC1 Signaling",                                                                                                
  "Hallmark: TNFA Signalling via NFKB",                                                                                         
  "KEGG: Antigen Processing and Presentation",                                                                                 
  "KEGG: Axon Guidance",                                                                                                       
  "KEGG: Cytokine-Cytokine Receptor Interaction",                                                                              
  "KEGG: FC Epsilon RI Signaling Pathway",                                                                                     
  "KEGG: Glycolysis/Gluconeogenesis",                                                                                          
  "KEGG: Leukocyte Transendotheliall Migration",                                                                                
  "KEGG: Lysosome",                                                                                                            
  "KEGG: NK Cell Mediated Cytotoxicity",                                                                           
  "KEGG: Oxidative Phosphorylation",                                                                                           
  "KEGG: Phosphatidylinositol Signaling System",                                                                               
  "KEGG: Regulation of Actin Cytoskeleton",                                                                                    
  "NABA: Secreted Factors",                                                                                                    
  "PID: HIF-1 TF Pathway",                                                                                                       
  "Reactome: Adaptive Immune System",                                                                                          
  "Reactome: Death Receptor Signaling",                                                                                       
  "Reactome: Diseases of Immune System",                                                                                       
  "Reactome: EPH Ephrin Signaling",                                                                                            
  "Reactome: Glucose Metabolism",                                                                                              
  "Reactome: Glycolysis",                                                                                                      
  "Reactome: Immunoregulatory Interactions Between a Lymphoid and a Non-Lymphoid Cell",                                        
  "Reactome: Innate Immune System",                                                                                            
  "Reactome: Interferon Alpha Beta Signaling",                                                                                 
  "Reactome: IFNg Signaling",                                                                                      
  "Reactome: IL-10 Signaling",                                                                                        
  "Reactome: Iron Uptake and Transport",                                                                                       
  "Reactome: Metabolism of Carbohydrates",                                                                                     
  "Reactome: Muscle Contraction",                                                                                              
  "Reactome: Neutrophil Degranulation",                                                                                        
  "Reactome: Post-Translational Protein Modification",                                                                         
  "Reactome: RAC1 GTPase Cycle",                                                                                               
  "Reactome: ATP Synthesis by Chemiosmotic Coupling in the ETC & Heat Production by Uncoupling Proteins",
  "Reactome: Rho GTPase Effectors",                                                                                            
  "Reactome: Rho GTPases Activate NADPH Oxidases",                                                                             
  "Reactome: RhoA GTPase Cycle",                                                                                               
  "Reactome: Signaling by GPCR",                                                                                               
  "Reactome: Signaling by Interleukins",                                                                                       
  "WP: Clear Cell Renal Cell Carcinoma Pathways",                                                                              
  "WP: Electron Transport Chain OXPhos System in Mitochondria",                                                                
  "WP: Ferroptosis",                                                                                                           
  "WP: Glycolysis and Gluconeogenesis",                                                                                        
  "WP: Hippo-Merlin Signaling Dysregulation",                                                                                   
  "WP: IL-18 Signaling Pathway",                                                                                                
  "WP: Metabolic Reprogramming in Colon Cancer",                                                                               
  "WP: Microglia Pathogen Phagocytosis Pathway",                                                                               
  "WP: miRNAs Involvement in the Immune Response in Sepsis",                                                                   
  "WP: Overview of Proinflammatory and Profibrotic Mediators",                                                                 
  "WP: Proteasome Degradation",                                                                                                
  "WP: Type II Interferon Signaling IFNg"
)

rownames(heatmap.data) <- row_names
colnames(heatmap.data) <- col_names

row_labels <- c(
  "G1 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G1 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G1 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G1 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  
  "G2 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G2 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G2 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G2 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  
  "G3 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G3 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G3 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G3 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  
  "G4 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G4 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G4 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G4 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  
  "G5 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G5 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G5 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G5 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  
  "G6 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G6 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G6 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G6 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  
  "G7 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G7 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G7 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G7 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  
  "G8 D3 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G8 D3 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G8 D14 *Mrp8*<sup>*Null*</sup>/*Hif1a*<sup>*fl/fl*</sup>",
  "G8 D14 *Mrp8*<sup>*Cre*</sup>/*Hif1a*<sup>*fl/fl*</sup>"
)

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c(-htmp_range, 0, htmp_range), c("blue", "white", "red"))

# Create heatmap
svglite::svglite(
  here("Figure 9", "svg__Figure 9C.svg"),
  width = 20,
  height = 15
)

ht <- Heatmap(
  heatmap.data,
  row_labels = gt_render(row_labels),
  name = "NES",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_gap = unit(5, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_anno,
  row_gap = unit(2.5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  row_names_gp = gpar(fontface = "bold"),
  row_title_gp = gpar(fontface = "bold"),
  column_names_gp = gpar(fontface = "bold", fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(heatmap2.data[i, j] < 0.0001) {
      grid.text("****", x, y)
    } else if(heatmap2.data[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(heatmap2.data[i, j] < 0.01) {
      grid.text("**", x, y)
    } else if(heatmap2.data[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
)

draw(ht, padding = unit(c(125, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure 10A ----
# Read in pathway data & select only granulocytic data
de <- read.csv(here("Figure 10", "Figure 10A.csv")) %>%
  filter(str_detect(cluster, "Granulocytes [:digit:]")) %>%
  select(-c(1, 2))

# Find overlap between hallmark_glycolysis & hallmark_hypoxia
HM.set <- msigdbr(species = "Homo sapiens", category = "H")

hypoxia <- HM.set %>%
  filter(gs_name == "HALLMARK_HYPOXIA")

glycolysis <- HM.set %>%
  filter(gs_name == "HALLMARK_GLYCOLYSIS")

hypoxia <- unique(hypoxia$human_gene_symbol)
glycolysis <- unique(glycolysis$human_gene_symbol)
gene_list <- intersect(hypoxia, glycolysis)
# Drop some non-informative genes
to_drop <- c(
  "GYS1",
  "LDHC",
  "ISG20",
  "HDLBP",
  "ANKZF1",
  "CHST2",
  "UGP2",
  "SDC2",
  "B3GALT6",
  "TGFBI",
  "PPFIA4",
  "COL5A1",
  "DCN",
  "PAM",
  "GALK1"
)

gene_list <- setdiff(gene_list, to_drop)

# Create matrix of NES values, fill with 0 if not found
heatmap.data <- de %>%
  filter(gene %in% gene_list)
heatmap.data <- dcast(heatmap.data, gene ~ cluster, value.var = "avg_log2FC", fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = "gene")
heatmap.data <- as.matrix(heatmap.data)

# Create another matrix of adjusted p-values, fill with 1 if not found
heatmap2.data <- de %>%
  filter(gene %in% gene_list)
heatmap2.data <- dcast(heatmap2.data, gene ~ cluster, value.var = "p_val_adj", fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = "gene")
heatmap2.data <- as.matrix(heatmap2.data)

# Set custom row order
col_order <- c(
  "Granulocytes 1_blood",
  "Granulocytes 1_tissue",
  
  "Granulocytes 2_blood",
  "Granulocytes 2_tissue",
  
  "Granulocytes 3_blood",
  "Granulocytes 3_tissue",
  
  "Granulocytes 4_blood",
  "Granulocytes 4_tissue",
  
  "Granulocytes 5_blood",
  "Granulocytes 5_tissue",
  
  "Granulocytes 6_blood",
  "Granulocytes 6_tissue",
  
  "Granulocytes 7_blood",
  "Granulocytes 7_tissue",
  
  "Granulocytes 8_blood",
  "Granulocytes 8_tissue",
  
  "Granulocytes 9_blood",
  "Granulocytes 9_tissue",
  
  "Granulocytes 10_blood",
  "Granulocytes 10_tissue",
  
  "Granulocytes 11_blood",
  "Granulocytes 11_tissue",
  
  "Granulocytes 12_blood",
  "Granulocytes 12_tissue"
)

heatmap.data <- heatmap.data[, col_order]

# Make the second heatmap match the row & column order of the first
col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

# Create a vector to split the heatmap according to cluster
col_anno <- c(
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue",
  
  "Blood",
  "Tissue"
)

# Create new shorter rownames
col_names <- c(
  "b_G1",
  "t_G1",
  
  "b_G2",
  "t_G2",
  
  "b_G3",
  "t_G3",
  
  "b_G4",
  "t_G4",
  
  "b_G5",
  "t_G5",
  
  "b_G6",
  "t_G6",
  
  "b_G7",
  "t_G7",
  
  "b_G8",
  "t_G8",
  
  "b_G9",
  "t_G9",
  
  "b_G10",
  "t_G10",
  
  "b_G11",
  "t_G11",
  
  "b_G12",
  "t_G12"
)

colnames(heatmap.data) <- col_names

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c(-htmp_range, 0, htmp_range), c("blue", "white", "red"))

# Create heatmap
svglite::svglite(
  here("Figure 10", "svg__Figure 10A.svg"),
  width = 10,
  height = 10
)

ht <- Heatmap(
  heatmap.data,
  name = "Log2\nFold Change",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  column_gap = unit(5, "mm"),
  column_split = col_anno,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  row_gap = unit(5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  row_names_gp = gpar(fontface = "bold.italic"),
  row_title_gp = gpar(fontface = "bold"),
  column_names_gp = gpar(fontface = "bold", fontsize = 10),
  column_title_gp = gpar(fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(heatmap2.data[i, j] < 0.0001) {
      grid.text("****", x, y)
    } else if(heatmap2.data[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(heatmap2.data[i, j] < 0.01) {
      grid.text("**", x, y)
    } else if(heatmap2.data[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
)

# Adjust heatmap padding
draw(ht, padding = unit(c(5, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure 10B ----
# Read in pathway data & combine
h_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "H_tissue")
c2_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "C2_tissue")

h_blood <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "H_blood")
c2_blood <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "C2_blood")

pathways <- rbind(
  h_tissue,
  c2_tissue,
  h_blood,
  c2_blood
)

tissue_up <- pathways %>%
  filter(sample_subset == "tissue") %>%
  filter(NES > 0) %>%
  filter(padj <= 0.05)

blood_up <- pathways %>%
  filter(sample_subset == "blood") %>%
  filter(NES > 0) %>%
  filter(padj <= 0.05)

up <- list(Tissue = tissue_up$pathway, Blood = blood_up$pathway)

names(up) <- c("Tissue - Up", "Blood - Up")

up <- ggvenn(
  up,
  fill_color = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  fill_alpha = 0.6,
  stroke_size = 0.5,
  set_name_size = 4
)

ggsave(
  "svg__Figure 10B_1.svg",
  plot = up,
  device = "svg",
  units = "in",
  width = 8,
  height = 8,
  path = here("Figure 10"),
  dpi = 600
)

tissue_dn <- pathways %>%
  filter(sample_subset == "tissue") %>%
  filter(NES < 0) %>%
  filter(padj <= 0.05)

blood_dn <- pathways %>%
  filter(sample_subset == "blood") %>%
  filter(NES < 0) %>%
  filter(padj <= 0.05)

dn <- list(Tissue = tissue_dn$pathway, Blood = blood_dn$pathway)

names(dn) <- c("Tissue - Down", "Blood - Down")

dn <- ggvenn(
  dn,
  fill_color = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  fill_alpha = 0.6,
  stroke_size = 0.5,
  set_name_size = 4
)

ggsave(
  "svg__Figure 10B_2.svg",
  plot = dn,
  device = "svg",
  units = "in",
  width = 8,
  height = 8,
  path = here("Figure 10"),
  dpi = 600
)

# Clear the environment
rm(list = ls())
gc()

## Figure 10C ----
# Read in pathway data & combine
h_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "H_tissue")
c2_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "C2_tissue")

h_blood <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "H_blood")
c2_blood <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "C2_blood")

pathways <- rbind(
  h_tissue,
  c2_tissue,
  h_blood,
  c2_blood
)

# Select only granulocytic data
pathways <- pathways %>%
  filter(str_detect(cluster_name, "Granulocytes [:digit:]"))

# Find top 5 enriched pathways per cluster
top_pathways <- pathways %>%
  arrange(desc(NES)) %>%
  group_by(cluster_num, sample_subset) %>%
  do(head(., n = 5))

# Find unique pathways from top enriched
pathway_list <- unique(top_pathways$pathway)

# Create matrix of NES values, fill with 0 if not found
heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- dcast(heatmap.data, pathway ~ cluster_name + sample_subset, value.var = "NES", fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = "pathway")
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

# Create another matrix of adjusted p-values, fill with 1 if not found
heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- dcast(heatmap2.data, pathway ~ cluster_name + sample_subset, value.var = "padj", fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = "pathway")
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

# Set custom row order
row_order <- c(
  "Granulocytes 1_blood",
  "Granulocytes 1_tissue",
  
  "Granulocytes 2_blood",
  "Granulocytes 2_tissue",
  
  "Granulocytes 3_blood",
  "Granulocytes 3_tissue",
  
  "Granulocytes 4_blood",
  "Granulocytes 4_tissue",
  
  "Granulocytes 5_blood",
  "Granulocytes 5_tissue",
  
  "Granulocytes 6_blood",
  "Granulocytes 6_tissue",
  
  "Granulocytes 7_blood",
  "Granulocytes 7_tissue",
  
  "Granulocytes 8_blood",
  "Granulocytes 8_tissue",
  
  "Granulocytes 9_blood",
  "Granulocytes 9_tissue",
  
  "Granulocytes 10_blood",
  "Granulocytes 10_tissue",
  
  "Granulocytes 11_blood",
  "Granulocytes 11_tissue",
  
  "Granulocytes 12_blood",
  "Granulocytes 12_tissue"
)

heatmap.data <- heatmap.data[row_order, ]

# Make the second heatmap match the row & column order of the first
col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

# Create a vector to split the heatmap according to cluster
row_anno <- c(
  "G01",
  "G01",
  
  "G02",
  "G02",
  
  "G03",
  "G03",
  
  "G04",
  "G04",
  
  "G05",
  "G05",
  
  "G06",
  "G06",
  
  "G07",
  "G07",
  
  "G08",
  "G08",
  
  "G09",
  "G09",
  
  "G10",
  "G10",
  
  "G11",
  "G11",
  
  "G12",
  "G12"
)

# Create new shorter rownames
row_names <- c(
  "G1_Blood",
  "G1_Tissue",
  
  "G2_Blood",
  "G2_Tissue",
  
  "G3_Blood",
  "G3_Tissue",
  
  "G4_Blood",
  "G4_Tissue",
  
  "G5_Blood",
  "G5_Tissue",
  
  "G6_Blood",
  "G6_Tissue",
  
  "G7_Blood",
  "G7_Tissue",
  
  "G8_Blood",
  "G8_Tissue",
  
  "G9_Blood",
  "G9_Tissue",
  
  "G10_Blood",
  "G10_Tissue",
  
  "G11_Blood",
  "G11_Tissue",
  
  "G12_Blood",
  "G12_Tissue"
)

# Change formatting of pathway names from the defaults
col_names <- c(
  "Biocarta: Biopeptides Pathway",                                                     
  "Hallmark: Epithelial Mesenchymal Transition",                                       
  "Hallmark: Hypoxia",                                                                 
  "Hallmark: Inflammatory Response",                                                   
  "Hallmark: IFNa Response",                                               
  "Hallmark: TNFa Signaling via NFKB",                                                 
  "KEGG: Antigen Processing and Presentation",                                         
  "KEGG: Cell Adhesion Molecules CAMs",                                                
  "KEGG: Gap Junction",                                                                
  "KEGG: Glioma",                                                                      
  "KEGG: GnRH Signaling Pathway",                                                      
  "KEGG: Long Term Potentiation",                                                      
  "KEGG: Melanogenesis",                                                               
  "KEGG: Neurotrophin Signaling Pathway",                                              
  "KEGG: Pathogenic E. coli Infection",                                       
  "KEGG: TGF beta Signaling Pathway",                                                  
  "KEGG: Vascular Smooth Muscle Contraction",                                          
  "KEGG: Wnt Signaling Pathway",                                                       
  "PID: AP1 Pathway",                                                                  
  "PID: NFAT 3 Pathway",                                                                
  "PID: P38 Alpha Beta Pathway",                                                       
  "PID: TRAIL Pathway",                                                                
  "Reactome: Activation of NMDA Receptors and Postsynaptic Events",                    
  "Reactome: Antigen Processing Cross Presentation",                                   
  "Reactome: Biological Oxidations",                                                   
  "Reactome: CDC42 GTPase Cycle",                                                      
  "Reactome: Cellular Response to Heat Stress",                                        
  "Reactome: G Protein Mediated Events",                                               
  "Reactome: HSF1 Dependent Transactivation",                                          
  "Reactome: Immunoregulatory Interactions between a Lymphoid and a Non-lymphoid Cell",
  "Reactome: Integration of Energy Metabolism",                                        
  "Reactome: IL-10 Signaling",                                                
  "Reactome: Neurotransmitter Receptors and Postsynaptic Signal Transmission",         
  "Reactome: Neutrophil Degranulation",                                                
  "Reactome: Opioid Signalling",                                                       
  "Reactome: Rho GTPase Effectors",                                                    
  "Reactome: Signaling by Wnt",                                                        
  "Reactome: TCF Dependent Signaling in Response to WNT",                              
  "WP: Angiopoietinlike Protein 8 Regulatory Pathway",                                 
  "WP: Glycogen Synthesis and Degradation",                                            
  "WP: Immune Response to Tuberculosis",                                               
  "WP: Integrated Cancer Pathway",                                                     
  "WP: Lung Fibrosis",                                                                 
  "WP: miRNAs in Cardiomyocyte Hypertrophy",                                        
  "WP: Network Map of SARS-CoV2 Signaling Pathway",                                     
  "WP: Orexin Receptor Pathway",                                                       
  "WP: Overview of Proinflammatory and Profibrotic Mediators",                         
  "WP: Pathogenic E. coli Infection",                                         
  "WP: Pathways Affected in Adenoid Cystic Carcinoma",                                
  "WP: Regucalcin in Proximal Tubule Epithelial Kidney Cells",                         
  "WP: Renin-Angiotensin-Aldosterone System (RAAS)",                                       
  "WP: Spinal Cord Injury"    
)

rownames(heatmap.data) <- row_names
colnames(heatmap.data) <- col_names

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c(-htmp_range, 0, htmp_range), c("blue", "white", "red"))

# Create heatmap
svglite::svglite(
  here("Figure 10", "svg__Figure 10C.svg"),
  width = 20,
  height = 15
)

ht <- Heatmap(
  heatmap.data,
  name = "NES",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_gap = unit(5, "mm"),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_anno,
  row_gap = unit(2.5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  row_names_gp = gpar(fontface = "bold"),
  row_title_gp = gpar(fontface = "bold"),
  column_names_gp = gpar(fontface = "bold", fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(heatmap2.data[i, j] < 0.0001) {
      grid.text("****", x, y)
    } else if(heatmap2.data[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(heatmap2.data[i, j] < 0.01) {
      grid.text("**", x, y)
    } else if(heatmap2.data[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
)

# Adjust heatmap padding
draw(ht, padding = unit(c(90, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure 11A ----
# Read in pathway data & select only granulocytic data
c00 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_0")
c01 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_1")
c02 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_2")
c03 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_3")
c04 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_4")
c05 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_5")
c06 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_6")
c07 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_7")
c08 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_8")
c09 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_9")
c10 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_10")
c11 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_11")
c12 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_12")
c13 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_13")
c14 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_14")
c15 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_15")
c16 <- read_xlsx(here("Figure 11", "Figure 11A.xlsx"), sheet = "DE_16")

de <- rbind(
  c00,
  c01,
  c02,
  c03,
  c04,
  c05,
  c06,
  c07,
  c08,
  c09,
  c10,
  c11,
  c12,
  c13,
  c14,
  c15,
  c16
)

de <- de %>%
  filter(str_detect(cluster, "Granulocytes [:digit:]")) %>%
  filter(str_detect(cluster, "tissue"))

# Find overlap between hallmark_glycolysis & hallmark_hypoxia
HM.set <- msigdbr(species = "Homo sapiens", category = "H")

hypoxia <- HM.set %>%
  filter(gs_name == "HALLMARK_HYPOXIA")

glycolysis <- HM.set %>%
  filter(gs_name == "HALLMARK_GLYCOLYSIS")

hypoxia <- unique(hypoxia$human_gene_symbol)
glycolysis <- unique(glycolysis$human_gene_symbol)
gene_list <- intersect(hypoxia, glycolysis)

# Drop some non-informative genes
to_drop <- c(
  "HDLBP",
  "TGFBI",
  "DCN",
  "PPFIA4",
  "PFKP",
  "PGAM2",
  "ANKZF1",
  "CITED2",
  "GALK1",
  "GYS1"
)

gene_list <- setdiff(gene_list, to_drop)

# Create matrix of NES values, fill with 0 if not found
heatmap.data <- de %>%
  filter(gene %in% gene_list)
heatmap.data <- dcast(heatmap.data, gene ~ cluster, value.var = "avg_log2FC", fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = "gene")
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]

# Create another matrix of adjusted p-values, fill with 1 if not found
heatmap2.data <- de %>%
  filter(gene %in% gene_list)
heatmap2.data <- dcast(heatmap2.data, gene ~ cluster, value.var = "p_val_adj", fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = "gene")
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]

# Set custom row order
col_order <- c(
  "Granulocytes 1_s1 tissue",
  "Granulocytes 1_s2 tissue",
  "Granulocytes 1_s5 tissue",
  
  "Granulocytes 2_s1 tissue",
  "Granulocytes 2_s2 tissue",
  "Granulocytes 2_s5 tissue",
  
  "Granulocytes 3_s1 tissue",
  "Granulocytes 3_s2 tissue",
  "Granulocytes 3_s5 tissue",
  
  "Granulocytes 4_s1 tissue",
  "Granulocytes 4_s2 tissue",
  "Granulocytes 4_s5 tissue",
  
  "Granulocytes 5_s1 tissue",
  "Granulocytes 5_s2 tissue",
  "Granulocytes 5_s5 tissue",
  
  "Granulocytes 6_s1 tissue",
  "Granulocytes 6_s2 tissue",
  "Granulocytes 6_s5 tissue",
  
  "Granulocytes 7_s1 tissue",
  "Granulocytes 7_s2 tissue",
  "Granulocytes 7_s5 tissue",
  
  "Granulocytes 8_s1 tissue",
  "Granulocytes 8_s2 tissue",
  "Granulocytes 8_s5 tissue",
  
  "Granulocytes 9_s1 tissue",
  "Granulocytes 9_s2 tissue",
  "Granulocytes 9_s5 tissue",
  
  "Granulocytes 10_s2 tissue",
  "Granulocytes 10_s5 tissue",
  
  "Granulocytes 11_s1 tissue",
  "Granulocytes 11_s2 tissue",
  "Granulocytes 11_s5 tissue",
  
  "Granulocytes 12_s1 tissue",
  "Granulocytes 12_s2 tissue",
  "Granulocytes 12_s5 tissue"
)

# Make the second heatmap match the row & column order of the first
heatmap.data <- heatmap.data[, col_order]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

# Create a vector to split the heatmap according to cluster
col_anno <- c(
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)"
)

# Create new shorter rownames
col_names <- c(
  "s1_G1",
  "s2_G1",
  "s5_G1",
  
  "s1_G2",
  "s2_G2",
  "s5_G2",
  
  "s1_G3",
  "s2_G3",
  "s5_G3",
  
  "s1_G4",
  "s2_G4",
  "s5_G4",
  
  "s1_G5",
  "s2_G5",
  "s5_G5",
  
  "s1_G6",
  "s2_G6",
  "s5_G6",
  
  "s1_G7",
  "s2_G7",
  "s5_G7",
  
  "s1_G8",
  "s2_G8",
  "s5_G8",
  
  "s1_G9",
  "s2_G9",
  "s5_G9",

  "s1_G10",
  "s5_G10",
  
  "s1_G11",
  "s2_G11",
  "s5_G11",
  
  "s1_G12",
  "s2_G12",
  "s5_G12"
)

colnames(heatmap.data) <- col_names

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c(-htmp_range, 0, htmp_range), c("blue", "white", "red"))

# Create heatmap
svglite::svglite(
  here("Figure 11", "svg__Figure 11A.svg"),
  width = 15,
  height = 10
)

ht <- Heatmap(
  heatmap.data,
  name = "Log2\nFold Change",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  cluster_columns = TRUE,
  cluster_column_slices = TRUE,
  column_gap = unit(2.5, "mm"),
  column_split = col_anno,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  row_gap = unit(2.5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  row_names_gp = gpar(fontface = "bold.italic"),
  row_title_gp = gpar(fontface = "bold"),
  column_names_gp = gpar(fontface = "bold", fontsize = 10),
  column_title_gp = gpar(fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(heatmap2.data[i, j] < 0.0001) {
      grid.text("****", x, y)
    } else if(heatmap2.data[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(heatmap2.data[i, j] < 0.01) {
      grid.text("**", x, y)
    } else if(heatmap2.data[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
)

# Adjust heatmap padding
draw(ht, padding = unit(c(5, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure 11B ----
# Read in pathway data & combine
h_s1_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_H_s1 tissue")
h_s2_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_H_s2 tissue")
h_s5_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_H_s5 tissue")

c2_s1_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_C2_s1 tissue")
c2_s2_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_C2_s2 tissue")
c2_s5_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_C2_s5 tissue")

pathways <- rbind(
  h_s1_tissue,
  h_s2_tissue,
  h_s5_tissue,
  
  c2_s1_tissue,
  c2_s2_tissue,
  c2_s5_tissue
)

s1_up <- pathways %>%
  filter(sample_subset == "s1 tissue") %>%
  filter(NES > 0) %>%
  filter(padj <= 0.05)

s2_up <- pathways %>%
  filter(sample_subset == "s2 tissue") %>%
  filter(NES > 0) %>%
  filter(padj <= 0.05)

s5_up <- pathways %>%
  filter(sample_subset == "s5 tissue") %>%
  filter(NES > 0) %>%
  filter(padj <= 0.05)

up <- list(S1 = s1_up$pathway, S2 = s2_up$pathway, S5 = s5_up$pathway)

names(up) <- c("S1 - Up", "S2 - Up", "S5 - Up")

up <- ggvenn(
  up,
  fill_color = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  fill_alpha = 0.6,
  stroke_size = 0.5,
  set_name_size = 4
)

ggsave(
  "svg__Figure 11B_1.svg",
  plot = up,
  device = "svg",
  units = "in",
  width = 8,
  height = 8,
  path = here("Figure 11"),
  dpi = 600
)

s1_dn <- pathways %>%
  filter(sample_subset == "s1 tissue") %>%
  filter(NES < 0) %>%
  filter(padj <= 0.05)

s2_dn <- pathways %>%
  filter(sample_subset == "s2 tissue") %>%
  filter(NES < 0) %>%
  filter(padj <= 0.05)

s5_dn <- pathways %>%
  filter(sample_subset == "s5 tissue") %>%
  filter(NES < 0) %>%
  filter(padj <= 0.05)

dn <- list(S1 = s1_dn$pathway, S2 = s2_dn$pathway, S5 = s5_dn$pathway)

names(dn) <- c("S1 - Down", "S2 - Down", "S5 - Down")

dn <- ggvenn(
  dn,
  fill_color = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  fill_alpha = 0.6,
  stroke_size = 0.5,
  set_name_size = 4
)

ggsave(
  "svg__Figure 11B_2.svg",
  plot = dn,
  device = "svg",
  units = "in",
  width = 8,
  height = 8,
  path = here("Figure 11"),
  dpi = 600
)

# Clear the environment
rm(list = ls())
gc()

## Figure 11C ----
# Read in pathway data & combine
h_s1_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_H_s1 tissue")
h_s2_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_H_s2 tissue")
h_s5_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_H_s5 tissue")

c2_s1_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_C2_s1 tissue")
c2_s2_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_C2_s2 tissue")
c2_s5_tissue <- read_xlsx(here("Figure 11", "Figure 11C.xlsx"), sheet = "full_C2_s5 tissue")

pathways <- rbind(
  h_s1_tissue,
  h_s2_tissue,
  h_s5_tissue,
  
  c2_s1_tissue,
  c2_s2_tissue,
  c2_s5_tissue
)

# Select only granulocytic data
pathways <- pathways %>%
  filter(str_detect(cluster_name, "Granulocytes [:digit:]"))

# Find top 5 enriched pathways per cluster
top_pathways <- pathways %>%
  arrange(desc(NES)) %>%
  group_by(cluster_num, sample_subset) %>%
  do(head(., n = 3))

# Find unique pathways from top enriched
pathway_list <- unique(top_pathways$pathway)

# Create matrix of NES values, fill with 0 if not found
heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- dcast(heatmap.data, pathway ~ cluster_name + sample_subset, value.var = "NES", fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = "pathway")
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

# Create another matrix of adjusted p-values, fill with 1 if not found
heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- dcast(heatmap2.data, pathway ~ cluster_name + sample_subset, value.var = "padj", fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = "pathway")
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

# Set custom row order
row_order <- c(
  "Granulocytes 1_s1 tissue",
  "Granulocytes 1_s2 tissue",
  "Granulocytes 1_s5 tissue",
  
  "Granulocytes 2_s1 tissue",
  "Granulocytes 2_s2 tissue",
  "Granulocytes 2_s5 tissue",
  
  "Granulocytes 3_s1 tissue",
  "Granulocytes 3_s2 tissue",
  "Granulocytes 3_s5 tissue",
  
  "Granulocytes 4_s1 tissue",
  "Granulocytes 4_s2 tissue",
  "Granulocytes 4_s5 tissue",
  
  "Granulocytes 5_s1 tissue",
  "Granulocytes 5_s2 tissue",
  "Granulocytes 5_s5 tissue",
  
  "Granulocytes 6_s1 tissue",
  "Granulocytes 6_s2 tissue",
  "Granulocytes 6_s5 tissue",
  
  "Granulocytes 7_s1 tissue",
  "Granulocytes 7_s2 tissue",
  "Granulocytes 7_s5 tissue",
  
  "Granulocytes 8_s1 tissue",
  "Granulocytes 8_s2 tissue",
  "Granulocytes 8_s5 tissue",
  
  "Granulocytes 9_s1 tissue",
  "Granulocytes 9_s2 tissue",
  "Granulocytes 9_s5 tissue",
  
  "Granulocytes 10_s2 tissue",
  "Granulocytes 10_s5 tissue",
  
  "Granulocytes 11_s1 tissue",
  "Granulocytes 11_s2 tissue",
  "Granulocytes 11_s5 tissue",
  
  "Granulocytes 12_s1 tissue",
  "Granulocytes 12_s2 tissue",
  "Granulocytes 12_s5 tissue"
)

heatmap.data <- heatmap.data[row_order, ]

# Make the second heatmap match the row & column order of the first
col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

# Create a vector to split the heatmap according to cluster
row_anno <- c(
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)",
  
  "Subject 1 (G+)",
  "Subject 2 (G-)",
  "Subject 5 (G+)"
)

# Create new shorter rownames
row_names <- c(
  "G1_s1 tissue",
  "G1_s2 tissue",
  "G1_s5 tissue",
  
  "G2_s1 tissue",
  "G2_s2 tissue",
  "G2_s5 tissue",
  
  "G3_s1 tissue",
  "G3_s2 tissue",
  "G3_s5 tissue",
  
  "G4_s1 tissue",
  "G4_s2 tissue",
  "G4_s5 tissue",
  
  "G5_s1 tissue",
  "G5_s2 tissue",
  "G5_s5 tissue",
  
  "G6_s1 tissue",
  "G6_s2 tissue",
  "G6_s5 tissue",
  
  "G7_s1 tissue",
  "G7_s2 tissue",
  "G7_s5 tissue",
  
  "G8_s1 tissue",
  "G8_s2 tissue",
  "G8_s5 tissue",
  
  "G9_s1 tissue",
  "G9_s2 tissue",
  "G9_s5 tissue",
  
  "G10_s2 tissue",
  "G10_s5 tissue",
  
  "G11_s1 tissue",
  "G11_s2 tissue",
  "G11_s5 tissue",
  
  "G12_s1 tissue",
  "G12_s2 tissue",
  "G12_s5 tissue"
)

# Change formatting of pathway names from the defaults
col_names <- c(
 "Biocarta: IL-1R pathway",                                                                    
 "Hallmark: Epithelial Mesenchymal Transition",                                               
 "Hallmark: Hypoxia",                                                                         
 "Hallmark: IFNa Response",                                                       
 "Hallmark: IFNg Response",                                                       
 "Hallmark: Oxidative Phosphorylation",                                                       
 "Hallmark: TNFa Signaling via NFKB",                                                         
 "Hallmark: UV Response Down",                                                                  
 "KEGG: Neurotrophin Signaling Pathway",                                                      
 "KEGG: Oxidative Phosphorylation",                                                           
 "KEGG: Ribosome",                                                                            
 "Naba: Matrisome",                                                                           
 "PID: AP1 Pathway",                                                                          
 "PID: Fra Pathway",                                                                          
 "PID: IL-1 Pathway",                                                                          
 "PID: IL-12 2 Pathway",                                                                        
 "PID: NFAT TF Pathway",                                                                       
 "Reactome: Attenuation Phase",                                                               
 "Reactome: Cellular Responses to Stimuli",                                                   
 "Reactome: Chemokine Receptors Bind Chemokines",                                             
 "Reactome: Diseases of Signal Transduction by GF Receptors & Second Messengers",
 "Reactome: Eukaryotic Translation Elongation",                                               
 "Reactome: HSF1 Activation",                                                                 
 "Reactome: HSF1 Dependent Transactivation",                                                  
 "Reactome: Innate Immune System",                                                            
 "Reactome: IL-1 Family Signaling",                                                  
 "Reactome: IL-10 Signaling",                                                        
 "Reactome: Metabolism of Amino Acids & Derivatives",                                       
 "Reactome: Neutrophil Degranulation",                                                        
 "Reactome: NGF Stimulated Transcription",                                                    
 "Reactome: Nuclear Events Kinase & Transcription Factor Activation",                       
 "Reactome: Peptide Ligand Binding Receptors",                                                
 "Reactome: Protein Localization",                                                            
 "Reactome: RAB GEFs Exchange GTP for GDP on RABs",                                           
 "Reactome: RAB Regulation of Trafficking",                                                   
 "Reactome: Response of EIF2AK4 (GCN2) to Amino Acid Deficiency",                               
 "Reactome: Rho GTPase CCycle",                                                                
 "Reactome: rRNA Processing",                                                                 
 "Reactome: SARS-CoV Infections",                                                             
 "Reactome: Signaling by Rho GTPases, Miro GTPases, & RHOBTB3",                               
 "Reactome: SRP Dependent Cotranslational Protein Targeting to Membrane",                     
 "WP: Ectoderm Differentiation",                                                              
 "WP: ETC OxPhos System in Mitochondria",                               
 "WP: Fatty Acid Beta Oxidation",                                                              
 "WP: IL-1 Signaling Pathway",                                                                 
 "WP: IL-18 Signaling Pathway",                                                                
 "WP: Network Map of SARS-CoV2 Signaling Pathway",                                             
 "WP: Notch Signaling Pathway",                                                               
 "WP: Orexin Receptor Pathway",                                                               
 "WP: Overview of Proinflammatory & Profibrotic Mediators",                                 
 "WP: Proximal Tubule Transport",                                                             
 "WP: Signal Transduction Through IL-1R"   
)

rownames(heatmap.data) <- row_names
colnames(heatmap.data) <- col_names

# Find the NES value that is furthest from zero & use it to set heatmap scale range
htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = colorRamp2(c(-htmp_range, 0, htmp_range), c("blue", "white", "red"))

# Create heatmap
svglite::svglite(
  here("Figure 11", "svg__Figure 11C.svg"),
  width = 20,
  height = 15
)

ht <- Heatmap(
  heatmap.data,
  name = "NES",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_gap = unit(2.5, "mm"),
  cluster_rows = TRUE,
  cluster_row_slices = TRUE,
  row_split = row_anno,
  row_gap = unit(2.5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  row_names_gp = gpar(fontface = "bold"),
  row_title_gp = gpar(fontface = "bold"),
  column_names_gp = gpar(fontface = "bold", fontsize = 10),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(heatmap2.data[i, j] < 0.0001) {
      grid.text("****", x, y)
    } else if(heatmap2.data[i, j] < 0.001) {
      grid.text("***", x, y)
    } else if(heatmap2.data[i, j] < 0.01) {
      grid.text("**", x, y)
    } else if(heatmap2.data[i, j] < 0.05) {
      grid.text("*", x, y)
    }
  }
)

# Adjust heatmap padding
draw(ht, padding = unit(c(85, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure S1A ----
# Load in data
integrated <- readRDS(here("Figure S1", "Figure S1.rds"))
dimred <- readRDS(here("Figure S1", "Figure S1 dim red.rds")) # Was calculated outside of the normal Seurat workflow & used for downstream analysis

# Add custom dimensional reduction
integrated[["umap"]] <- CreateDimReducObject(embeddings = dimred, key = "UMAP_", assay = DefaultAssay(integrated))

# Create plot
p <- DimPlot(
  integrated,
  pt.size = 1.5,
  label.color = "white",
  group.by = "Sample.ID") +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  ggtitle("") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(face = "bold"))

p$data$Sample.ID <- factor(x = p$data$Sample.ID, levels = c("D03", "D07", "D14"))
p <- p + scale_color_nejm(labels = c("D03", "D07", "D14"), alpha = 0.6)
p <- p + theme(legend.position = "bottom")

ggsave(
  "Figure S1A.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S1")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S1B ----
# Read in data
integrated <- readRDS(here("Figure S1", "Figure S1.rds"))

ggData = data.frame(prop.table(table(integrated$Sample.ID, Idents(integrated)), margin = 2))
colnames(ggData) = c("Sample", "cluster", "value")
ggData$Sample <- factor(x = ggData$Sample, levels = c("D03", "D07", "D14"))

p1 <- ggplot(ggData, aes(cluster, value, fill = Sample)) +
  geom_col() +
  xlab("Cluster") +
  ylab("Proportion of Cells (%)") +
  scale_fill_nejm(labels = c("D03", "D07", "D14"), alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

ggData = data.frame(table(integrated$Sample.ID, Idents(integrated)))
colnames(ggData) = c("Sample", "cluster", "value")
ggData$Sample <- factor(x = ggData$Sample, levels = c("D03", "D07", "D14"))

p2 <- ggplot(ggData, aes(cluster, value, fill = Sample)) +
  geom_col() +
  xlab("Cluster") +
  ylab("Cell Number") +
  scale_fill_nejm(labels = c("D03", "D07", "D14"), alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

p <- p1 / p2 + plot_layout(guides = "collect")

ggsave(
  "Figure S1B.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S1")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S2 ----
# Load in data
compass_results <- read.csv(here("Figure S2", "Figure S2.csv")) %>%
  group_by(subsystem) %>%
  mutate(alpha = case_when(adjusted_pval > 0.1 ~ "N", adjusted_pval <= 0.1 ~ "Y"))

p <- ggplot(compass_results, aes(x = cohens_d, y = subsystem)) +
  geom_point(data = subset(compass_results, cohens_d > 0), aes(alpha = alpha), shape = 21, fill = "red", size = 2) +
  geom_point(data = subset(compass_results, cohens_d < 0), aes(alpha = alpha), shape = 21, fill = "blue", size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(c(-5, 5)) +
  labs(x = "Cohen's d", y = "Metabolic Subsystem") +
  theme_bw() +
  scale_alpha_manual(name = "", labels = c("NS", bquote(p[adj]<=~0.1)), values = c(0.25, 1.00)) +
  guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
  theme(legend.position = "bottom")

ggsave(
  "Figure S2.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S2")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S3A ----
# Read in data
d <- read.csv(here("Figure S3", "Figure S3A.csv")) %>%
  select(-1) %>%
  mutate(cell_type = str_replace(cell_type, "MDSC", "G-MDSC"))

p <- ggplot(d, aes(x = Seurat_patho_score, y = Seurat_maturity_score, fill = cell_type)) +
  geom_point(shape = 21, color = "black", size = 3) +
  scale_fill_manual(name = "Cell type", values = c("tomato", "slateblue1")) +
  labs(x = "Pathogenicity score", y = "Maturity score") +
  theme_classic() +
  theme(legend.position = "bottom", axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

ggsave(
  "Figure S3A.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S3")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S3B ----
# Read in data
cor.table <- read.csv(here("Figure S3", "Figure S3B_1.csv")) %>%
  select(-1)

compass.result <- read.csv(here("Figure S3", "Figure S3B_2.csv")) %>%
  filter(adjusted_pval < 0.1)

cor.table <- column_to_rownames(cor.table, var = "reaction")
cor.table <- cor.table %>%
  filter(rownames(.) %in% compass.result$X)

heatmap.data <- as.matrix(cor.table)

ha <- HeatmapAnnotation(
  "Gene markers" = c(rep("G-MDSC", 9), rep("PMN", 13)),
  col = list("Gene markers" = c("G-MDSC" = "tomato", "PMN" = "slateblue1")),
  gp = gpar(col = "black"),
  show_annotation_name = FALSE
)

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

tiff(
  here("Figure S3", "Figure S3B.tiff"),
  width = 860,
  height = 670
)

Heatmap(
  heatmap.data,
  name = "Spearman correlation",
  col = col_fun,
  column_names_side = "top",
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_parent_dend_line = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  column_names_rot = 45,
  border = TRUE,
  heatmap_legend_param = list(border = "black"),
  top_annotation = ha
)

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure S6A ----
# Load in data
integrated <- readRDS(here("Figure S6", "Figure S6.rds"))

# Create plot
p <- DimPlot(
  integrated,
  pt.size = 1.5,
  label.color = "white",
  group.by = "sample_origin",
  order = c("d3_null", "d3_cre", "d14_null", "d14_cre")) +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  ggtitle("") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(face = "bold"))

p$data$sample_origin <- factor(x = p$data$sample_origin, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))
p <- p + scale_color_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6)
p <- p + theme(legend.position = "bottom")

ggsave(
  "Figure S6A.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S6")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S6B ----
# Read in data
integrated <- readRDS(here("Figure S6", "Figure S6.rds"))

ggData = data.frame(prop.table(table(integrated$sample_origin, Idents(integrated)), margin = 2))
colnames(ggData) = c("Sample", "cluster", "value")
ggData$Sample <- factor(x = ggData$Sample, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))

p1 <- ggplot(ggData, aes(cluster, value, fill = Sample)) +
  geom_col() +
  xlab("Cluster") +
  ylab("Proportion of Cells (%)") +
  scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

ggData = data.frame(table(integrated$sample_origin, Idents(integrated)))
colnames(ggData) = c("Sample", "cluster", "value")
ggData$Sample <- factor(x = ggData$Sample, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))

p2 <- ggplot(ggData, aes(cluster, value, fill = Sample)) +
  geom_col() +
  xlab("Cluster") +
  ylab("Cell Number") +
  scale_fill_nejm(labels = c(bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})), bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))), alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

p <- p1 / p2 + plot_layout(guides = "collect")

ggsave(
  "Figure S6B.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S6")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S7A ----
# Load in data, normalize, scale, & subset
integrated <- readRDS(here("Figure S7", "Figure S7.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

new.cluster.ids <- c(
  "G1",
  "G2",
  "G3",
  "G4",
  "G5",
  "G6",
  "G7",
  "G8",
  "GMM",
  "L",
  "G9"
)

names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)

# Create plot
vln1 <- VlnPlot(
  integrated,
  features = "Hif1a",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9")
)

vln2 <- VlnPlot(
  integrated,
  features = "S100a8",
  split.by = "sample_origin",
  cols = c("#BC3C2999", "#0072B599", "#E1872799", "#20854E99"),
  idents = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8", "G9")
)

vln1$data$split <- factor(x = vln1$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))
vln2$data$split <- factor(x = vln2$data$split, levels = c("d3_null", "d3_cre", "d14_null", "d14_cre"))

vln1 <- vln1 + scale_fill_nejm(labels = c("D3 WT", "D3 HIF-1a cKO", "D14 WT", "D14 HIF-1a cKO"), alpha = 0.6)
vln2 <- vln2 + scale_fill_nejm(labels = c("D3 WT", "D3 HIF-1a cKO", "D14 WT", "D14 HIF-1a cKO"), alpha = 0.6)

vln <- (vln1 / vln2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(
  "Figure S7A.tiff",
  plot = vln,
  device = "tiff",
  path = here("Figure S7")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S7B ----
# Load in data
integrated <- readRDS(here("Figure S7", "Figure S7.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

integrated <- subset(integrated, idents = c(
  "Granulocytes 1",
  "Granulocytes 2",
  "Granulocytes 3",
  "Granulocytes 4",
  "Granulocytes 5",
  "Granulocytes 6",
  "Granulocytes 7",
  "Granulocytes 8",
  "Granulocytes 9"
))

DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

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

# Create heatmap
png(
  here("Figure S6", "Figure S6B.png"),
  width = 500,
  height = 670
)

ht <- Heatmap(
  heatmap.data,
  name = "Average\nNormalized\nCounts",
  col = col_fun,
  rect_gp = gpar(col = "black", lwd = 2),
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
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", heatmap.data[i, j]), x, y, gp = gpar(fontsize = 10))
  }
)

# Adjust heatmap padding
draw(ht, padding = unit(c(5, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()

# Clear the environment
rm(list = ls())
gc()

## Figure S8 ----
# Load in data
compass_results <- read.csv(here("Figure S8", "Figure S8_1.csv"))

aam <- c(
  "Alanine and aspartate metabolism",
  "Arginine and Proline Metabolism",
  "beta-Alanine metabolism",
  "Cysteine Metabolism",
  "D-alanine metabolism",
  "Folate metabolism",
  "Glutamate metabolism",
  "Glycine, serine, alanine and threonine metabolism",
  "Histidine metabolism",
  "Lysine metabolism",
  "Methionine and cysteine metabolism",
  "Taurine and hypotaurine metabolism",
  "Tryptophan metabolism",
  "Tyrosine metabolism",
  "Urea cycle",
  "Valine, leucine, and isoleucine metabolism"
)

p <- ggplot(compass_results, aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "gray", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  annotate("text", x = 2.75, y = 3.50, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.75, y = 3.50, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(bold(-log[10](BH-adjusted~Wilcoxon~rank~sum~p)))) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"))

p1 <- ggplot(compass_results %>% filter(subsystem =="Glycolysis/gluconeogenesis"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#BC3C29FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Glycolysis") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p2 <- ggplot(compass_results %>% filter(subsystem =="Citric acid cycle"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#0072B5FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("TCA cycle") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p3 <- ggplot(compass_results %>% filter(subsystem %in% aam), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#E18727FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Amino acid metabolism") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p4 <- ggplot(compass_results %>% filter(subsystem =="Fatty acid oxidation"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#20854EFF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's d") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Fatty acid oxidation") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

figure <- ggarrange(
  p1 + rremove("ylab") + rremove("xlab"), 
  p2 + rremove("ylab") + rremove("xlab"), 
  p3 + rremove("ylab") + rremove("xlab"), 
  p4 + rremove("ylab") + rremove("xlab"),
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
)

p5 <- annotate_figure(
  figure,
  left = textGrob(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), rot = 90, vjust = 0.5, gp = gpar(cex = 1.3)),
  bottom = textGrob("Cohen's D", gp = gpar(cex = 1.3)))

d3 <- p | figure

# Load in data
compass_results <- read.csv(here("Figure S8", "Figure S8_2.csv"))

aam <- c(
  "Alanine and aspartate metabolism",
  "Arginine and Proline Metabolism",
  "beta-Alanine metabolism",
  "Cysteine Metabolism",
  "D-alanine metabolism",
  "Folate metabolism",
  "Glutamate metabolism",
  "Glycine, serine, alanine and threonine metabolism",
  "Histidine metabolism",
  "Lysine metabolism",
  "Methionine and cysteine metabolism",
  "Taurine and hypotaurine metabolism",
  "Tryptophan metabolism",
  "Tyrosine metabolism",
  "Urea cycle",
  "Valine, leucine, and isoleucine metabolism"
)

p <- ggplot(compass_results, aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "gray", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  annotate("text", x = 2.75, y = 3.50, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.75, y = 3.50, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(bold(-log[10](BH-adjusted~Wilcoxon~rank~sum~p)))) +
  theme_classic() +
  theme(axis.title = element_text(face = "bold"))

p1 <- ggplot(compass_results %>% filter(subsystem =="Glycolysis/gluconeogenesis"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#BC3C29FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Glycolysis") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p2 <- ggplot(compass_results %>% filter(subsystem =="Citric acid cycle"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#0072B5FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("TCA cycle") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p3 <- ggplot(compass_results %>% filter(subsystem %in% aam), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#E18727FF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Amino acid metabolism") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

p4 <- ggplot(compass_results %>% filter(subsystem =="Fatty acid oxidation"), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = "#20854EFF", size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 1.3, linetype = "dashed") +
  annotate("text", x = 2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Null), bolditalic(Hif1a^{fl/fl})))) +
  annotate("text", x = -2.25, y = 6.00, label = bquote(atop(bolditalic(Mrp8^Cre), bolditalic(Hif1a^{fl/fl})))) +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's D") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle("Fatty acid oxidation") +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA)
  )

figure <- ggarrange(
  p1 + rremove("ylab") + rremove("xlab"), 
  p2 + rremove("ylab") + rremove("xlab"), 
  p3 + rremove("ylab") + rremove("xlab"), 
  p4 + rremove("ylab") + rremove("xlab"),
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = "bottom",
  align = "hv", 
  font.label = list(size = 10, color = "black", face = "bold", family = NULL, position = "top")
)

p5 <- annotate_figure(
  figure,
  left = textGrob(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), rot = 90, vjust = 0.5, gp = gpar(cex = 1.3)),
  bottom = textGrob("Cohen's D", gp = gpar(cex = 1.3)))

d14 <- p | figure

final <- d3 / d14



ggsave(
  "Figure S8.tiff",
  plot = final,
  device = "tiff",
  path = here("Figure S8"),
  width = 14,
  height = 10,
  units = "in"
)

# Clear the environment
rm(list = ls())
gc()

## Figure S9 ----
# Load in data
compass_results <- read.csv(here("Figure S9", "Figure S9.csv")) %>%
  group_by(subsystem) %>%
  mutate(alpha = case_when(adjusted_pval > 0.1 ~ "N", adjusted_pval <= 0.1 ~ "Y"))

p <- ggplot(compass_results, aes(x = cohens_d, y = subsystem)) +
  geom_point(data = subset(compass_results, cohens_d > 0), aes(alpha = alpha), shape = 21, fill = "red", size = 2) +
  geom_point(data = subset(compass_results, cohens_d < 0), aes(alpha = alpha), shape = 21, fill = "blue", size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(c(-2.5, 2.5)) +
  labs(x = "Cohen's D", y = "Metabolic Subsystem") +
  theme_bw() +
  scale_alpha_manual(name = "", labels = c("NS", bquote(p[adj]<=~0.1)), values = c(0.25, 1.00)) +
  guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
  theme(legend.position = "bottom")

ggsave(
  "Figure S9.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S9")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S10 ----
# Load in data
compass_results <- read.csv(here("Figure S10", "Figure S10.csv")) %>%
  group_by(subsystem) %>%
  mutate(alpha = case_when(adjusted_pval > 0.1 ~ "N", adjusted_pval <= 0.1 ~ "Y"))

p <- ggplot(compass_results, aes(x = cohens_d, y = subsystem)) +
  geom_point(data = subset(compass_results, cohens_d > 0), aes(alpha = alpha), shape = 21, fill = "red", size = 2) +
  geom_point(data = subset(compass_results, cohens_d < 0), aes(alpha = alpha), shape = 21, fill = "blue", size = 2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  xlim(c(-2.5, 2.5)) +
  labs(x = "Cohen's D", y = "Metabolic Subsystem") +
  theme_bw() +
  scale_alpha_manual(name = "", labels = c("NS", bquote(p[adj]<=~0.1)), values = c(0.25, 1.00)) +
  guides(alpha = guide_legend(override.aes = list(fill = "black"))) +
  theme(legend.position = "bottom")

ggsave(
  "Figure S10.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S10")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S11A ----
# Load in data
integrated <- readRDS(here("Figure S10", "Figure S11.rds"))

# Create plot
p <- DimPlot(
  integrated,
  label = TRUE,
  label.box = TRUE,
  label.size = 6,
  repel = TRUE,
  pt.size = 1.5,
  label.color = "black") +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(face = "bold"))

p <- p + theme(legend.position = "bottom")

ggsave(
  "Figure S10A.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S10")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S11B ----
# Load in data
integrated <- readRDS(here("Figure S10", "Figure S11.rds"))

# Create plot
p <- DimPlot(
  integrated,
  pt.size = 1.5,
  label.color = "white",
  group.by = "Sample_origin") +
  theme_classic() +
  labs(x = "UMAP 1", y = "UMAP 2") +
  ggtitle("") +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_text(size = 12, face = "bold"), legend.text = element_text(face = "bold"))

p$data$sample_origin <- factor(x = p$data$Sample_origin, levels = c("s1 blood", "s1 tissue", "s2 blood", "s2 tissue", "s5 blood", "s5 tissue"))
p <- p + scale_color_nejm(labels = c(
  "Subject 1 - blood",
  "Subject 1 - tissue",
  "Subject 2 - blood",
  "Subject 2 - tissue",
  "Subject 5 - blood",
  "Subject 5 - tissue"
),
alpha = 0.6
)

p <- p + theme(legend.position = "bottom")

ggsave(
  "Figure S10B.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S10")
)

# Clear the environment
rm(list = ls())
gc()

## Figure S11C ----
# Read in data
integrated <- readRDS(here("Figure S10", "Figure S11.rds"))

ggData = data.frame(prop.table(table(integrated$Sample_origin, Idents(integrated)), margin = 2))
colnames(ggData) = c("Sample", "cluster", "value")
ggData$Sample <- factor(x = ggData$Sample, levels = c("s1 blood", "s1 tissue", "s2 blood", "s2 tissue", "s5 blood", "s5 tissue"))

p1 <- ggplot(ggData, aes(cluster, value, fill = Sample)) +
  geom_col() +
  xlab("Cluster") +
  ylab("Proportion of Cells (%)") +
  scale_fill_nejm(labels = c("Subject 1 - blood", "Subject 1 - tissue", "Subject 2 - blood", "Subject 2 - tissue", "Subject 5 - blood", "Subject 5 - tissue"), alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

ggData = data.frame(table(integrated$Sample_origin, Idents(integrated)))
colnames(ggData) = c("Sample", "cluster", "value")
ggData$Sample <- factor(x = ggData$Sample, levels = c("s1 blood", "s1 tissue", "s2 blood", "s2 tissue", "s5 blood", "s5 tissue"))

p2 <- ggplot(ggData, aes(cluster, value, fill = Sample)) +
  geom_col() +
  xlab("Cluster") +
  ylab("Cell Number") +
  scale_fill_nejm(labels = c("Subject 1 - blood", "Subject 1 - tissue", "Subject 2 - blood", "Subject 2 - tissue", "Subject 5 - blood", "Subject 5 - tissue"), alpha = 0.6) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  coord_flip() +
  theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(face = "bold"), legend.title = element_text(face = "bold"), legend.text = element_text(face = "bold"))

p <- p1 / p2 + plot_layout(guides = "collect")

ggsave(
  "Figure S10C.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S10")
)

# Clear the environment
rm(list = ls())
gc()