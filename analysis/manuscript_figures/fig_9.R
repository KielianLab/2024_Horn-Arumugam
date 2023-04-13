## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure 9 in Horn et al., 2023
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

## Figure 9A ----
# Read in pathway data & combine
h_tissue <- read_xlsx(here("Figure 9", "Figure 9A.xlsx"), sheet = "H_tissue")
c2_tissue <- read_xlsx(here("Figure 9", "Figure 9A.xlsx"), sheet = "C2_tissue")

h_blood <- read_xlsx(here("Figure 9", "Figure 9A.xlsx"), sheet = "H_blood")
c2_blood <- read_xlsx(here("Figure 9", "Figure 9A.xlsx"), sheet = "C2_blood")

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
  fill_color = c("#BC3C2999", "#0072B599"),
  fill_alpha = 0.6,
  stroke_size = 0.5,
  set_name_size = 4
)

# Clear the environment
rm(list = ls())

gc()

## Figure 9B ----
# Read in pathway data & select only granulocytic data
de <- read.csv(here("Figure 9", "Figure 9B.csv")) %>%
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
png(
  here("Figure 9", "Figure 9B.png"),
  width = 800,
  height = 600
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

## Figure 9C ----
# Read in pathway data & combine
h_tissue <- read_xlsx(here("Figure 9", "Figure 9C.xlsx"), sheet = "H_tissue")
c2_tissue <- read_xlsx(here("Figure 9", "Figure 9C.xlsx"), sheet = "C2_tissue")

h_blood <- read_xlsx(here("Figure 9", "Figure 9C.xlsx"), sheet = "H_blood")
c2_blood <- read_xlsx(here("Figure 9", "Figure 9C.xlsx"), sheet = "C2_blood")

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
png(
  here("Figure 9", "Figure 9C.png"),
  width = 1600,
  height = 1000
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
