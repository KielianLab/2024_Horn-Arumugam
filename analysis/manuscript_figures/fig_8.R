## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure 8 in Horn et al., 2023
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

## Figure 8A ----
#Load in data
integrated <- readRDS(here("Figure 8", "Figure 8.rds"))

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
  "Figure 8A.tiff",
  plot = integrated_umap.clus,
  device = "tiff",
  path = here("Figure 8")
)

# Clear the environment
rm(list = ls())
gc()

## Figure 8B ----
#Load in data, normalize, scale, & subset
integrated <- readRDS(here("Figure 8", "Figure 8.rds"))
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

vln1 <- vln1 + scale_fill_nejm(labels = c(
  bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))),
  alpha = 0.6
)

vln2 <- vln2 + scale_fill_nejm(labels = c(
  bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))),
  alpha = 0.6
)

vln3 <- vln3 + scale_fill_nejm(labels = c(
  bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))),
  alpha = 0.6
)

vln4 <- vln4 + scale_fill_nejm(labels = c(
  bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))),
  alpha = 0.6
)

vln5 <- vln5 + scale_fill_nejm(labels = c(
  bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))),
  alpha = 0.6
)

vln6 <- vln6 + scale_fill_nejm(labels = c(
  bquote(bolditalic(D3~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D3~Mrp8^Cre/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Null/Hif1a^{fl/fl})),
  bquote(bolditalic(D14~Mrp8^Cre/Hif1a^{fl/fl}))),
  alpha = 0.6
)



vln <- (vln1 | vln2 | vln3) / (vln4 | vln5 | vln6) + plot_layout(guides = "collect") & theme(legend.position = "bottom", legend.text = element_text(face = "bold", size = 12))



ggsave(
  "Figure 8B.tiff",
  plot = vln,
  device = "tiff",
  path = here("Figure 8")
)



# Clear the environment
rm(list = ls())

gc()

## Figure 8C ----
# Read in pathway data & combine
h <- read_xlsx(here("Figure 8", "Figure 8C.xlsx"), sheet = "HM pathways")
c2 <- read_xlsx(here("Figure 8", "Figure 8C.xlsx"), sheet = "C2 pathways")

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
png(
  here("Figure 8", "Figure 8C.png"),
  width = 1600,
  height = 1400
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

draw(ht, padding = unit(c(135, 5, 2, 5), "mm")) # bottom, left, top, right paddings

# Close the device
dev.off()


# Clear the environment
rm(list = ls())

gc()
