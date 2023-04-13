## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure 4A in Horn et al., 2023
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
  "Figure 4A.tiff",
  plot = vln,
  device = "tiff",
  path = here("Figure 4")
)

# Clear the environment
rm(list = ls())
gc()

