## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure S2 in Horn et al., 2023
## Created: 2023-01-23
## Updated: 2023-04-11

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

## Figure S2A ----
#Load in data
integrated <- readRDS(here("Figure S2", "Figure S2.rds"))
dimred <- readRDS(here("Figure S2", "Figure S2 dim red.rds")) # Was calculated outside of the normal Seurat workflow & used for downstream analysis

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
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(face = "bold")
  )

p$data$Sample.ID <- factor(x = p$data$Sample.ID, levels = c("D03", "D07", "D14"))
p <- p + scale_color_nejm(labels = c("D03", "D07", "D14"), alpha = 0.6)
p <- p + theme(legend.position = "bottom")



ggsave(
  "Figure S2A.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S2")
)

# Clear the environment
rm(list = ls())

gc()

## Figure S2B ----
# Read in data
integrated <- readRDS(here("Figure S2", "Figure S2.rds"))

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
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

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
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

p <- p1 / p2 + plot_layout(guides = "collect")

ggsave(
  "Figure S2B.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S2")
)

# Clear the environment
rm(list = ls())

gc()
