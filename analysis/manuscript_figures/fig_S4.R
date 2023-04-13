## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure S4 in Horn et al., 2023
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

## Figure S4A ----
# Read in data
d <- read.csv(here("Figure S4", "Figure S4A.csv")) %>%
  select(-1) %>%
  mutate(cell_type = str_replace(cell_type, "MDSC", "G-MDSC"))

p <- ggplot(d, aes(x = Seurat_patho_score, y = Seurat_maturity_score, fill = cell_type)) +
  geom_point(shape = 21, color = "black", size = 3) +
  scale_fill_manual(name = "Cell type", values = c("tomato", "slateblue1")) +
  labs(x = "Pathogenicity score", y = "Maturity score") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"),
    legend.text = element_text(face = "bold")
  )

ggsave(
  "Figure S4A.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S4")
)

# Clear the environment
rm(list = ls())

gc()

## Figure S4B ----
# Read in data
cor.table <- read.csv(here("Figure S4", "Figure S4B_1.csv")) %>%
  select(-1)

compass.result <- read.csv(here("Figure S4", "Figure S4B_2.csv")) %>%
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
  here("Figure S4", "Figure S4B.tiff"),
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
