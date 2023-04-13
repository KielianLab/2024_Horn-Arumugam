## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure S10 in Horn et al., 2023
## Created: 2023-01-23
## Updated: 2023-04-13

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

## Figure S10A ----
# Load in data
compass_results <- read.csv(here("Figure S10", "Figure S10A.csv"))

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
compass_results <- read.csv(here("Figure S10", "Figure S10B.csv"))

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
  "Figure S10.tiff",
  plot = final,
  device = "tiff",
  path = here("Figure S10"),
  width = 14,
  height = 10,
  units = "in"
)

# Clear the environment
rm(list = ls())

gc()
