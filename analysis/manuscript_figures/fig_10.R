## Made by Christopher M. Horn, MS
## Kielian Lab data
## Code for Figure 10 in Horn et al., 2023
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

## Figure 10A ----
# Read in pathway data & combine
h_s1_tissue <- read_xlsx(here("Figure 10", "Figure 10A.xlsx"), sheet = "full_H_s1 tissue")
h_s2_tissue <- read_xlsx(here("Figure 10", "Figure 10A.xlsx"), sheet = "full_H_s2 tissue")
h_s5_tissue <- read_xlsx(here("Figure 10", "Figure 10A.xlsx"), sheet = "full_H_s5 tissue")

c2_s1_tissue <- read_xlsx(here("Figure 10", "Figure 10A.xlsx"), sheet = "full_C2_s1 tissue")
c2_s2_tissue <- read_xlsx(here("Figure 10", "Figure 10A.xlsx"), sheet = "full_C2_s2 tissue")
c2_s5_tissue <- read_xlsx(here("Figure 10", "Figure 10A.xlsx"), sheet = "full_C2_s5 tissue")

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

# Clear the environment
rm(list = ls())

gc()

## Figure 10B ----
# Read in pathway data & select only granulocytic data
c00 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_0")
c01 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_1")
c02 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_2")
c03 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_3")
c04 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_4")
c05 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_5")
c06 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_6")
c07 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_7")
c08 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_8")
c09 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_9")
c10 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_10")
c11 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_11")
c12 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_12")
c13 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_13")
c14 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_14")
c15 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_15")
c16 <- read_xlsx(here("Figure 10", "Figure 10B.xlsx"), sheet = "DE_16")

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
png(
  here("Figure 10", "Figure 10B.png"),
  width = 1200,
  height = 600
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

## Figure 10C ----
# Read in pathway data & combine
h_s1_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "full_H_s1 tissue")
h_s2_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "full_H_s2 tissue")
h_s5_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "full_H_s5 tissue")

c2_s1_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "full_C2_s1 tissue")
c2_s2_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "full_C2_s2 tissue")
c2_s5_tissue <- read_xlsx(here("Figure 10", "Figure 10C.xlsx"), sheet = "full_C2_s5 tissue")

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
png(
  here("Figure 10", "Figure 10C.png"),
  width = 1650,
  height = 1100
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
