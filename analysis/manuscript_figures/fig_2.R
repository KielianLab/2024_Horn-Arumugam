## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure 2 in Horn et al., 2023
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

## Figure 2A ----
# Read in pathway data & select only granulocytic data
integrated <- readRDS(here("Figure 2", "Figure 2.rds"))
DefaultAssay(integrated) <- "RNA"
integrated <- NormalizeData(integrated)
integrated <- ScaleData(integrated)

avg <- AverageExpression(integrated, return.seurat = TRUE)
avg.norm_counts <- avg@assays$RNA@data
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

heatmap2.data <- heatmap2.data[row_order, col_order]

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
png(
  here("Figure 2", "Figure 2A.png"),
  width = 400,
  height = 300
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
png(
  here("Figure 2", "Figure 2B.png"),
  width = 200,
  height = 400
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
  "Figure 2C.png",
  plot = final,
  device = "png",
  path = here("Figure 2"),
  width = 14,
  height = 7,
  units = "in"
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
  "Figure 2D.tiff",
  plot = final,
  device = "tiff",
  path = here("Figure 2")
)

# Clear the environment
rm(list = ls())
gc()
