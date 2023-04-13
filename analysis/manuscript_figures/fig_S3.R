## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure S3 in Horn et al., 2023
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

## Figure S3 ----
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
  "Figure S3.tiff",
  plot = p,
  device = "tiff",
  path = here("Figure S3")
)

# Clear the environment
rm(list = ls())

gc()
