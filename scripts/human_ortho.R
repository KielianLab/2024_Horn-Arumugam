## Made by Christopher M. Horn, MS
## Kielian Lab data
## Created: 2022-11-27
## Updated: 2022-11-27

## Notes: Blood & tissue data from a patient with a PJI




# Setting up environment -----

## Load in packages
packages <- c(
  'tidyverse',
  'Seurat',
  'patchwork',
  'SingleR',
  'SingleCellExperiment',
  'sctransform',
  'MAST',
  'plotly',
  'ggsci',
  'readxl',
  'fgsea',
  'data.table',
  'clusterExperiment',
  'pheatmap',
  'mgcv',
  'slingshot',
  'msigdbr',
  'celldex',
  'dittoSeq',
  'here'
)

new.packages <- packages[!(packages %in% installed.packages()[,'Package'])]
if(length(new.packages)) install.packages(new.packages)
invisible(lapply(packages, library, character.only = TRUE))

rm(packages)

## Set options
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(12345)



# Initializing objects -----

## Setting up Seurat objects
## Load in the data
blood.data <- Read10X(data.dir = here('raw data', 'blood'))
blood.data <- blood.data$`Gene Expression`
tissue.data <- Read10X(data.dir = here('raw data', 'tissue'))
tissue.data <- tissue.data$`Gene Expression`

## Initialize Seurat object w/raw data
blood <- CreateSeuratObject(counts = blood.data, project = 'blood', min.cells = 3, min.features = 200)
tissue <- CreateSeuratObject(counts = tissue.data, project = 'tissue', min.cells = 3, min.features = 200)

## Pre-processing and QC
blood[['percent.mt']] <- PercentageFeatureSet(blood, pattern = '^MT-')
blood <- subset(blood, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

tissue[['percent.mt']] <- PercentageFeatureSet(tissue, pattern = '^MT-')
tissue <- subset(tissue, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Create timing metadata
sample.id_1 <- rep('blood', length(blood@meta.data$orig.ident))
sample.id_2 <- rep('tissue', length(tissue@meta.data$orig.ident))
names(sample.id_1) <- rownames(blood@meta.data)
names(sample.id_2) <- rownames(tissue@meta.data)
blood <- AddMetaData(blood, sample.id_1, col.name = 'Sample_origin')
tissue <- AddMetaData(tissue, sample.id_2, col.name = 'Sample_origin')

## Write individual object metadata to file
write.csv(blood@meta.data, here('output', 'blood_metadata.csv'))
write.csv(tissue@meta.data, here('output', 'tissue_metadata.csv'))

## Write QC metrics
blood.QC_mets <- dim(blood)
tissue.QC_mets <- dim(tissue)

write.csv(blood.QC_mets, here('output', 'QC', 'blood QC_mets.csv'))
write.csv(tissue.QC_mets, here('output', 'QC', 'tissue QC_mets.csv'))

blood.QC_mets.plot  <- VlnPlot(blood, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
tissue.QC_mets.plot  <- VlnPlot(tissue, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

ggsave('blood QC_mets.plot.png', plot = blood.QC_mets.plot, device = 'png', path = here('output', 'QC'))
ggsave('tissue QC_mets.plot.png', plot = tissue.QC_mets.plot, device = 'png', path = here('output', 'QC'))

## Remove temp objects
rm(blood.data,
   tissue.data,
   sample.id_1,
   sample.id_2,
   blood.QC_mets,
   blood.QC_mets.plot,
   tissue.QC_mets,
   tissue.QC_mets.plot
)

gc()



# Annotating cell types -----

## Annotate the cells w/SingleR
blood <- NormalizeData(blood)
tissue <- NormalizeData(tissue)

ref.se <- HumanPrimaryCellAtlasData() # snapshot date: 2020-10-27

blood_sce <- as.SingleCellExperiment(blood)
tissue_sce <- as.SingleCellExperiment(tissue)

commonGenes.1 <- intersect(rownames(blood_sce), rownames(ref.se))
commonGenes.2 <- intersect(rownames(tissue_sce), rownames(ref.se))

ref.se_1 <- ref.se[commonGenes.1,]
ref.se_2 <- ref.se[commonGenes.2,]

blood_sce <- blood_sce[commonGenes.1,]
tissue_sce <- tissue_sce[commonGenes.2,]

pred.blood <- SingleR(test = blood_sce, ref = ref.se_1, labels = ref.se_1$label.main)
pred.tissue <- SingleR(test = tissue_sce, ref = ref.se_2, labels = ref.se_2$label.main)

blood[['celltype']] <- pred.blood$pruned.labels
tissue[['celltype']] <- pred.tissue$pruned.labels



## Remove temp objects
rm(ref.se,
   blood_sce,
   tissue_sce,
   commonGenes.1,
   commonGenes.2,
   ref.se_1,
   ref.se_2,
   pred.blood,
   pred.tissue
)

gc()



# Object integration -----

## Integrating all cells from all days
sample.list <- c(blood, tissue)
names(sample.list) <- c('blood', 'tissue')
for (i in 1:length(sample.list))
  {sample.list[[i]] <- SCTransform(sample.list[[i]], verbose = TRUE)}
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- PrepSCTIntegration(object.list = sample.list, anchor.features = sample.features, verbose = TRUE)
sample.anchors <- FindIntegrationAnchors(object.list = sample.list, normalization.method = 'SCT', anchor.features = sample.features, verbose = TRUE)
ortho.integrated <- IntegrateData(anchorset = sample.anchors, normalization.method = 'SCT', verbose = TRUE)
ortho.integrated <- RunPCA(ortho.integrated, verbose = TRUE)
ortho.integrated <- RunUMAP(ortho.integrated, dims = 1:30)
ortho.integrated <- FindNeighbors(ortho.integrated, dims = 1:30)
ortho.integrated <- FindClusters(ortho.integrated, resolution = 0.5)
saveRDS(ortho.integrated, here('output', 'ortho.integrated.rds'))



## Remove temp objects
rm(sample.list,
   sample.features,
   sample.anchors,
   blood,
   tissue,
   i
)

gc()



## Output celltype composition of each cluster
sample.comp_origin <- table(ortho.integrated$Sample_origin)
clust.comp_origin <- table(Idents(ortho.integrated), ortho.integrated$Sample_origin)
clust.comp_celltype <- prop.table(table(Idents(ortho.integrated), ortho.integrated$celltype), margin = 1)



write.csv(sample.comp_origin, here('output', 'sample.comp_origin.csv'))
write.csv(clust.comp_origin, here('output', 'clust.comp_origin.csv'))
write.csv(clust.comp_celltype, here('output', 'clust.comp_celltype.csv'))



rm(sample.comp_origin,
   clust.comp_origin,
   clust.comp_celltype)

gc()



# Renaming/plotting -----

## Renaming clusters
new.cluster.ids <- c('Granulocytes 1',
                     'Granulocytes 2',
                     'Granulocytes 3',
                     'Granulocytes 4',
                     'Granulocytes 5',
                     'T Cells',
                     'Granulocytes 6',
                     'Granulocytes 7',
                     'Granulocytes/Monocytes'
)

names(new.cluster.ids) <- levels(ortho.integrated)
ortho.integrated <- RenameIdents(ortho.integrated, new.cluster.ids)



rm(new.cluster.ids)

gc()



## Save seurat objects
saveRDS(ortho.integrated, here('output', 'ortho.integrated.rds'))



## Plotting w/new labels
ortho_pca.origin <- DimPlot(
  ortho.integrated,
  group.by = 'Sample_origin',
  reduction = 'pca',
  cols = c(
    '#BC3C2999',
    '#0072B599',
    '#E1872799',
    '#20854E99'
  ),
  pt.size = 1.5
) +
theme_classic() +
ggtitle('Split by sample origin') +
labs(
  x = 'PC 1',
  y = 'PC 2'
) +
theme(
  axis.text = element_blank(),
  axis.ticks = element_blank()
)



ortho_pca.clus <- DimPlot(
  ortho.integrated,
  reduction = 'pca',
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 1.5
) +
theme_classic() +
NoLegend() +
labs(
  x = 'PC 1',
  y = 'PC 2'
) +
theme(
  axis.text = element_blank(),
  axis.ticks = element_blank()
)



ortho_umap.origin <- DimPlot(
  ortho.integrated,
  group.by = 'Sample_origin',
  cols = c(
    '#BC3C2999',
    '#0072B599',
    '#E1872799',
    '#20854E99'
  ),
  pt.size = 1.5
) +
theme_classic() +
ggtitle('Split by sample origin') +
labs(
  x = 'UMAP 1',
  y = 'UMAP 2'
) +
theme(
  axis.text = element_blank(),
  axis.ticks = element_blank()
)



ortho_umap.clus <- DimPlot(
  ortho.integrated,
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 1.5
) +
theme_classic() +
NoLegend() +
labs(
  x = 'UMAP 1',
  y = 'UMAP 2'
) +
theme(
  axis.text = element_blank(),
  axis.ticks = element_blank()
)



ggsave('ortho_pca origin.png', plot = ortho_pca.origin, device = 'png', path = here('output'))
ggsave('ortho_pca cluster.png', plot = ortho_pca.clus, device = 'png', path = here('output'))
ggsave('ortho_umap origin.png', plot = ortho_umap.origin, device = 'png', path = here('output'))
ggsave('ortho_umap cluster.png', plot = ortho_umap.clus, device = 'png', path = here('output'))



rm(
  ortho_pca.origin,
  ortho_pca.clus,
  ortho_umap.origin,
  ortho_umap.clus
)

gc()



## Output cluster number to name cheat sheet
num2name <- data.frame(
  'Cluster_name' = levels(ortho.integrated),
  'Cluster_num' = levels(ortho.integrated$seurat_clusters)
)
write.csv(num2name, here('output', 'num 2 name cheat sheet.csv'))

## Write the metadata to a file
write.csv(ortho.integrated@meta.data, here('output', 'ortho integrated_metadata.csv'))



rm(num2name)

gc()



## Output m vs p for each cluster
sample.comp_origin <- table(ortho.integrated$Sample_origin)
clust.comp_origin <- table(Idents(ortho.integrated), ortho.integrated$Sample_origin)
clust.comp_celltype <- table(ortho.integrated$celltype, Idents(ortho.integrated))

write.csv(sample.comp_origin, here('output', 'sample.comp_origin.csv'))
write.csv(clust.comp_origin, here('output', 'clust.comp_origin.csv'))
write.csv(clust.comp_celltype, here('output', 'clust.comp_celltype.csv'))



rm(sample.comp_origin,
   clust.comp_origin,
   clust.comp_celltype)

gc()



## Subsetting out
brain.microglia <- subset(brain.integrated, idents = c('Microglia 1',
                                                       'Microglia 2',
                                                       'Microglia 3'))
brain.granulocyte <- subset(brain.integrated, idents = c('Granulocytes 1',
                                                         'Granulocytes 2',
                                                         'Granulocytes 3'))
brain.monomac <- subset(brain.integrated, idents = c('Mono/Mac 1',
                                                     'Mono/Mac 2'))




# Cluster-level DE -----
## New (post 2022-10-05)
DefaultAssay(ortho.integrated) <- 'RNA'
ortho.integrated <- NormalizeData(ortho.integrated, verbose = TRUE)
ortho.integrated <- ScaleData(ortho.integrated, verbose = TRUE)

DE <- FindAllMarkers(ortho.integrated, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')

write.csv(DE, here('output', 'DE', 'full DE.csv'))

for (i in seq_along(levels(DE$cluster))) {
  if(i != length(DE$cluster)) {
    single.DE <- DE %>% filter(DE$cluster %in% c(levels(DE$cluster)[i]))
    write.csv(single.DE, here('output', 'DE', paste0('DE_', i-1, '.csv')))
  }
}



rm(
  i,
  single.DE,
  DE
)

gc()




## Within cluster DE Analysis
c01 <- subset(ortho.integrated, idents = c('Granulocytes 1'))
c02 <- subset(ortho.integrated, idents = c('Granulocytes 2'))
c03 <- subset(ortho.integrated, idents = c('Granulocytes 3'))
c04 <- subset(ortho.integrated, idents = c('Granulocytes 4'))
c05 <- subset(ortho.integrated, idents = c('Granulocytes 5'))
c06 <- subset(ortho.integrated, idents = c('T Cells'))
c07 <- subset(ortho.integrated, idents = c('Granulocytes 6'))
c08 <- subset(ortho.integrated, idents = c('Granulocytes 7'))
c09 <- subset(ortho.integrated, idents = c('Granulocytes/Monocytes'))


c.list <- list(
  c01,
  c02,
  c03,
  c04,
  c05,
  c06,
  c07,
  c08,
  c09
)



rm(
  c01,
  c02,
  c03,
  c04,
  c05,
  c06,
  c07,
  c08,
  c09
)

gc()



for(i in 1:length(c.list)){
  DefaultAssay(c.list[[i]]) <- 'RNA'
  c.list[[i]] <- NormalizeData(c.list[[i]], verbose = TRUE)
  c.list[[i]] <- ScaleData(c.list[[i]], verbose = TRUE)
  c.list[[i]]$cell_origin <- paste(Idents(c.list[[i]]), c.list[[i]]$Sample_origin, sep = '_')
  c.list[[i]]$cell <- Idents(c.list[[i]])
  Idents(c.list[[i]]) <- 'cell_origin'
}

for(i in 1:length(c.list)){
  DE <- FindAllMarkers(c.list[[i]], min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
  DE <- DE %>%
    arrange(desc(avg_log2FC))
  write.csv(DE, here('output', 'DE', 'within cluster', paste0('DE_', i-1, '.csv')))
}



rm(
  i,
  c.list,
  DE
)

gc()



GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # Canonical Pathways
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

geneSet_list <- list(GO.set, HM.set, CP.set)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:length(levels(ortho.integrated))) {
    glist <- DE %>% filter(cluster == levels(ortho.integrated)[ii])
    glist <- column_to_rownames(glist, var = 'gene')
    stats <- glist$avg_log2FC
    names(stats) <- toupper(rownames(glist))
    stats <- sort(stats, decreasing = TRUE)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(eaRes, file = here('output', 'GSEA', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '.tsv')), sep="\t", sep2=c("", " ", ""))
    temp <- read.table(file = here('output', 'GSEA', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '.tsv')), sep = '\t', header = TRUE)
    temp <- temp %>% mutate(cluster_num = rep(ii-1, length(nrow(temp))), cluster_name = rep(levels(ortho.integrated)[ii], length(nrow(temp))), pathway_db = rep(unique(geneSet_list[[i]]$gs_cat), length(nrow(temp))))
    write.csv(temp, here('output', 'GSEA', 'csv', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '.csv')))
  }
}



rm(geneSet_list,
   m_list,
   stats,
   glist,
   eaRes,
   DE,
   GO.set,
   HM.set,
   CP.set,
   i,
   ii)

gc()



GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

geneSet_list <- list(GO.set, HM.set, CP.set)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:length(DE_list)) {
    for(iii in 1:length(levels(as.factor(DE_list[[ii]]$cluster)))) {
      glist <- DE_list[[ii]] %>% filter(cluster == levels(as.factor(DE_list[[ii]]$cluster))[iii])
      glist <- column_to_rownames(glist, var = 'gene')
      stats <- glist$avg_log2FC
      names(stats) <- toupper(rownames(glist))
      stats <- sort(stats, decreasing = TRUE)
      eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
      eaRes <- arrange(eaRes, desc(NES))
      fwrite(eaRes, file = here('output', 'GSEA', 'within cluster', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '_', str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'), '.tsv')), sep="\t", sep2=c("", " ", ""))
      temp <- read.table(file = here('output', 'GSEA', 'within cluster', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '_', str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'), '.tsv')), sep = '\t', header = TRUE)
      temp <- temp %>% mutate(cluster_num = rep(ii-1, length(nrow(temp))), cluster_name = rep(levels(ortho.integrated)[ii], length(nrow(temp))), pathway_db = rep(unique(geneSet_list[[i]]$gs_cat), length(nrow(temp))), sample_subset = rep(str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'), length(nrow(temp))))
      write.csv(temp, here('output', 'GSEA', 'within cluster', 'csv', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '_', str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'), '.csv')))
    }
  }
}



rm(geneSet_list,
   m_list,
   stats,
   glist,
   eaRes,
   DE_list,
   GO.set,
   HM.set,
   CP.set,
   i,
   ii,
   iii)

gc()



c2_00 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_0.csv'))
c2_01 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_1.csv'))
c2_02 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_2.csv'))
c2_03 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_3.csv'))
c2_04 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_4.csv'))
c2_05 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_5.csv'))
c2_06 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_6.csv'))
c2_07 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_7.csv'))

h_00 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_0.csv'))
h_01 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_1.csv'))
h_02 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_2.csv'))
h_03 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_3.csv'))
h_04 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_4.csv'))
h_05 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_5.csv'))
h_06 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_6.csv'))
h_07 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_7.csv'))

c5_00 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_0.csv'))
c5_01 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_1.csv'))
c5_02 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_2.csv'))
c5_03 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_3.csv'))
c5_04 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_4.csv'))
c5_05 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_5.csv'))
c5_06 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_6.csv'))
c5_07 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_7.csv'))


c2 <- rbind(c2_00, c2_01, c2_02, c2_03, c2_04, c2_05, c2_06, c2_07)
h <- rbind(h_00, h_01, h_02, h_03, h_04, h_05, h_06, h_07)
c5 <- rbind(c5_00, c5_01, c5_02, c5_03, c5_04, c5_05, c5_06, c5_07)



write.csv(c2, here('output', 'GSEA', 'full_C2.csv'))
write.csv(h, here('output', 'GSEA', 'full_H.csv'))
write.csv(c5, here('output', 'GSEA', 'full_C5.csv'))




c2_00_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_blood.csv'))
c2_00_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_tissue.csv'))

c2_01_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_blood.csv'))
c2_01_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_tissue.csv'))

c2_02_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_blood.csv'))
c2_02_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_tissue.csv'))

c2_03_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_blood.csv'))
c2_03_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_tissue.csv'))

c2_04_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_blood.csv'))
c2_04_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_tissue.csv'))

c2_05_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_blood.csv'))
c2_05_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_tissue.csv'))

c2_06_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_blood.csv'))
c2_06_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_tissue.csv'))

c2_07_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_blood.csv'))
c2_07_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_tissue.csv'))



h_00_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_blood.csv'))
h_00_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_tissue.csv'))

h_01_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_blood.csv'))
h_01_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_tissue.csv'))

h_02_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_blood.csv'))
h_02_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_tissue.csv'))

h_03_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_blood.csv'))
h_03_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_tissue.csv'))

h_04_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_blood.csv'))
h_04_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_tissue.csv'))

h_05_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_blood.csv'))
h_05_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_tissue.csv'))

h_06_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_blood.csv'))
h_06_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_tissue.csv'))

h_07_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_blood.csv'))
h_07_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_tissue.csv'))



c5_00_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_blood.csv'))
c5_00_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_tissue.csv'))

c5_01_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_blood.csv'))
c5_01_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_tissue.csv'))

c5_02_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_blood.csv'))
c5_02_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_tissue.csv'))

c5_03_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_blood.csv'))
c5_03_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_tissue.csv'))

c5_04_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_blood.csv'))
c5_04_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_tissue.csv'))

c5_05_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_blood.csv'))
c5_05_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_tissue.csv'))

c5_06_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_blood.csv'))
c5_06_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_tissue.csv'))

c5_07_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_blood.csv'))
c5_07_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_tissue.csv'))



c2_blood <- rbind(c2_00_blood, c2_01_blood, c2_02_blood, c2_03_blood, c2_04_blood, c2_05_blood, c2_06_blood, c2_07_blood)
c2_tissue <- rbind(c2_00_tissue, c2_01_tissue, c2_02_tissue, c2_03_tissue, c2_04_tissue, c2_05_tissue, c2_06_tissue, c2_07_tissue)

h_blood <- rbind(h_00_blood, h_01_blood, h_02_blood, h_03_blood, h_04_blood, h_05_blood, h_06_blood, h_07_blood)
h_tissue <- rbind(h_00_tissue, h_01_tissue, h_02_tissue, h_03_tissue, h_04_tissue, h_05_tissue, h_06_tissue, h_07_tissue)

c5_blood <- rbind(c5_00_blood, c5_01_blood, c5_02_blood, c5_03_blood, c5_04_blood, c5_05_blood, c5_06_blood, c5_07_blood)
c5_tissue <- rbind(c5_00_tissue, c5_01_tissue, c5_02_tissue, c5_03_tissue, c5_04_tissue, c5_05_tissue, c5_06_tissue, c5_07_tissue)



write.csv(c2_blood, here('output', 'GSEA', 'within cluster', 'full_C2_blood.csv'))
write.csv(c2_tissue, here('output', 'GSEA', 'within cluster', 'full_C2_tissue.csv'))

write.csv(h_blood, here('output', 'GSEA', 'within cluster', 'full_H_blood.csv'))
write.csv(h_tissue, here('output', 'GSEA', 'within cluster', 'full_H_tissue.csv'))

write.csv(c5_blood, here('output', 'GSEA', 'within cluster', 'full_C5_blood.csv'))
write.csv(c5_tissue, here('output', 'GSEA', 'within cluster', 'full_C5_tissue.csv'))



top_pathways <- pathways %>%
  arrange(desc(NES)) %>%
  group_by(cluster_num, sample_subset) %>%
  do(head(., n = 10))

clust_num <- 19 # Choose your cluster of interest here

path_clust <- top_pathways %>% filter(cluster_num == clust_num)
clust_name <- unique(path_clust$cluster_name)
pathway_list <- unique(path_clust$pathway)

heatmap.data <- path_clust %>%
  filter(pathway %in% pathway_list)
heatmap.data <- reshape2::dcast(heatmap.data, pathway ~ sample_subset, value.var = 'NES', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'pathway')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

heatmap2.data <- path_clust %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- reshape2::dcast(heatmap2.data, pathway ~ sample_subset, value.var = 'padj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'pathway')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

row_order <- c(
  'wt',
  'rag',
  'rag.th1',
  'rag.th17'
)

heatmap.data <- heatmap.data[row_order, ]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = circlize::colorRamp2(c(0, 0, htmp_range), c('white', 'white', 'red')) # Check to make sure that the upper bound of this range isn't cutting off the expression of some genes in the matrix
ht <- ComplexHeatmap::Heatmap(heatmap.data,
                              name = 'NES',
                              column_title = clust_name,
                              col = col_fun,
                              rect_gp = grid::gpar(col = 'black', lwd = 2),
                              cluster_columns = TRUE,
                              cluster_column_slices = FALSE,
                              column_gap = grid::unit(5, 'mm'),
                              column_names_gp = grid::gpar(fontsize = 8),
                              cluster_rows = FALSE,
                              cluster_row_slices = FALSE,
                              row_gap = grid::unit(5, 'mm'),
                              show_parent_dend_line = F,
                              heatmap_legend_param = list(border = 'black'),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                if(heatmap2.data[i, j] < 0.0001) {
                                  grid::grid.text('****', x, y)
                                } else if(heatmap2.data[i, j] < 0.001) {
                                  grid::grid.text('***', x, y)
                                } else if(heatmap2.data[i, j] < 0.01) {
                                  grid::grid.text('**', x, y)
                                } else if(heatmap2.data[i, j] < 0.05) {
                                  grid::grid.text('*', x, y)
                                }
                              })

ComplexHeatmap::draw(ht, padding = unit(c(70, 5, 2, 5), 'mm')) # bottom, left, top, right paddings



top_genes_up <- de %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  do(head(., n = 10))

top_genes_down <- de %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  do(tail(., n = 10))

top_genes <- rbind(top_genes_up, top_genes_down)

clust <- 19 # Choose your cluster of interest here

gene_clust <- top_genes %>% filter(cluster_num == clust)
clust_name <- unique(str_extract(gene_clust$cluster, '[^_]+'))
gene_list <- unique(gene_clust$gene)

heatmap.data <- gene_clust %>%
  filter(gene %in% gene_list)
heatmap.data <- reshape2::dcast(heatmap.data, gene ~ cluster, value.var = 'avg_log2FC', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)

heatmap2.data <- gene_clust %>%
  filter(gene %in% gene_list)
heatmap2.data$p_val_adj <- as.numeric(heatmap2.data$p_val_adj)
heatmap2.data <- reshape2::dcast(heatmap2.data, gene ~ cluster, value.var = 'p_val_adj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'gene')
heatmap2.data <- as.matrix(heatmap2.data)

col_order <- c(
  paste0(clust_name, '_wt'),
  paste0(clust_name, '_rag'),
  paste0(clust_name, '_rag.th1'),
  paste0(clust_name, '_rag.th17')
)

heatmap.data <- heatmap.data[, col_order]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)



col_fun = circlize::colorRamp2(c(-htmp_range, 0, htmp_range), c('blue', 'white', 'red'))

ht <- ComplexHeatmap::Heatmap(heatmap.data,
                              name = 'Average Log2\nFold Change',
                              column_title = clust_name,
                              col = col_fun,
                              rect_gp = grid::gpar(col = 'black', lwd = 2),
                              cluster_columns = FALSE,
                              cluster_column_slices = FALSE,
                              column_gap = grid::unit(5, 'mm'),
                              column_names_gp = grid::gpar(fontsize = 8),
                              cluster_rows = TRUE,
                              cluster_row_slices = FALSE,
                              row_gap = grid::unit(5, 'mm'),
                              row_names_gp = grid::gpar(fontsize = 8),
                              show_parent_dend_line = FALSE,
                              heatmap_legend_param = list(border = 'black'),
                              cell_fun = function(j, i, x, y, width, height, fill) {
                                if(heatmap2.data[i, j] < 0.0001) {
                                  grid::grid.text('****', x, y)
                                } else if(heatmap2.data[i, j] < 0.001) {
                                  grid::grid.text('***', x, y)
                                } else if(heatmap2.data[i, j] < 0.01) {
                                  grid::grid.text('**', x, y) 
                                } else if(heatmap2.data[i, j] < 0.05) {
                                  grid::grid.text('*', x, y)
                                }
                              }
)

ComplexHeatmap::draw(ht, padding = unit(c(2, 5, 2, 5), 'mm')) # bottom, left, top, right paddings



## Old (pre 2022-10-06)
## Between cluster DE Analysis
DefaultAssay(cellplex.integrated) <- 'RNA'
cellplex.integrated <- NormalizeData(cellplex.integrated, verbose = T)
cellplex.integrated <- ScaleData(cellplex.integrated, verbose = T)

DE <- FindAllMarkers(cellplex.integrated, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')

write.csv(DE, here('output', 'DE', 'full DE.csv'))

for (i in seq_along(levels(DE$cluster))) {
  if(i != length(DE$cluster)) {
    single.DE <- DE %>% filter(DE$cluster %in% c(levels(DE$cluster)[i]))
    write.csv(single.DE, here('output', 'DE', paste0('DE_', i-1, '.csv')))
  }
}



rm(i,
   single.DE,
   DE)

gc()



## Within cluster DE Analysis
c01 <- subset(cellplex.integrated, idents = c('Microglia 1'))
c02 <- subset(cellplex.integrated, idents = c('Microglia 2'))
c03 <- subset(cellplex.integrated, idents = c('Granulocytes 1'))
c04 <- subset(cellplex.integrated, idents = c('Granulocytes 2'))
c05 <- subset(cellplex.integrated, idents = c('Microglia/Macrophages'))
c06 <- subset(cellplex.integrated, idents = c('Granulocytes 3'))
c07 <- subset(cellplex.integrated, idents = c('Microglia 3'))
c08 <- subset(cellplex.integrated, idents = c('T/NKT Cells'))
c09 <- subset(cellplex.integrated, idents = c('Monocytes/DC/B Cells'))
c10 <- subset(cellplex.integrated, idents = c('NK/NKT/ILC'))
c11 <- subset(cellplex.integrated, idents = c('Granulocytes 4'))
c12 <- subset(cellplex.integrated, idents = c('Granulocytes 5'))
c13 <- subset(cellplex.integrated, idents = c('Microglia 4'))
c14 <- subset(cellplex.integrated, idents = c('Macrophages'))
c15 <- subset(cellplex.integrated, idents = c('DC'))
c16 <- subset(cellplex.integrated, idents = c('Monocytes 1'))
c17 <- subset(cellplex.integrated, idents = c('Microglia 5'))
c18 <- subset(cellplex.integrated, idents = c('Granulocytes/Monocytes/Macrophages'))
c19 <- subset(cellplex.integrated, idents = c('DC/pro B Cells'))
c20 <- subset(cellplex.integrated, idents = c('Monocytes 2'))

c.list <- list(c01,
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
               c16,
               c17,
               c18,
               c19,
               c20)



rm(c01,
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
   c16,
   c17,
   c18,
   c19,
   c20)

gc()



for(i in 1:length(c.list)){
  DefaultAssay(c.list[[i]]) <- 'RNA'
  c.list[[i]] <- NormalizeData(c.list[[i]], verbose = T)
  c.list[[i]] <- ScaleData(c.list[[i]], verbose = T)
  c.list[[i]]$cell_origin <- paste(Idents(c.list[[i]]), c.list[[i]]$Sample_origin, sep = '_')
  c.list[[i]]$cell <- Idents(c.list[[i]])
  Idents(c.list[[i]]) <- 'cell_origin'
}

for(i in 1:length(c.list)){
  DE <- FindAllMarkers(c.list[[i]], min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
  DE <- DE %>%
    arrange(desc(avg_log2FC))
  write.csv(DE, here('output', 'DE', 'within cluster', paste0('DE_', i-1, '.csv')))
}



rm(i,
   c.list,
   DE)

gc()



## fgsea
GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # Canonical Pathways
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

geneSet_list <- list(GO.set, HM.set, CP.set)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:length(levels(cellplex.integrated))) {
    glist <- DE %>% filter(cluster == levels(cellplex.integrated)[ii])
    glist <- column_to_rownames(glist, var = 'gene')
    stats <- glist$avg_log2FC
    names(stats) <- toupper(rownames(glist))
    stats <- sort(stats, decreasing = T)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(eaRes, file = here('output', 'GSEA', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '.tsv')), sep="\t", sep2=c("", " ", ""))
  }
}



rm(geneSet_list,
   m_list,
   stats,
   glist,
   eaRes,
   DE,
   GO.set,
   HM.set,
   CP.set,
   i,
   ii)

gc()



GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

geneSet_list <- list(GO.set, HM.set, CP.set)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:length(DE_list)) {
    for(iii in 1:length(levels(as.factor(DE_list[[ii]]$cluster)))) {
      glist <- DE_list[[ii]] %>% filter(cluster == levels(as.factor(DE_list[[ii]]$cluster))[iii])
      glist <- column_to_rownames(glist, var = 'gene')
      stats <- glist$avg_log2FC
      names(stats) <- toupper(rownames(glist))
      stats <- sort(stats, decreasing = T)
      eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
      eaRes <- arrange(eaRes, desc(NES))
      fwrite(eaRes, file = here('output', 'GSEA', 'within cluster', unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '_', iii, '.tsv')), sep="\t", sep2=c("", " ", ""))
    }
  }
}



rm(geneSet_list,
   m_list,
   stats,
   glist,
   eaRes,
   DE_list,
   GO.set,
   HM.set,
   CP.set,
   i,
   ii,
   iii)

gc()



# Single-cell GSEA ------

GS <- getGeneSets(library = 'H', species = 'Mus musculus')
# GS <- getGeneSets(library = 'C5', species = 'Mus musculus')
ES <- enrichIt(neutro.integrated, gene.sets = GS, groups = 1000, cores = 2)
neutro.integrated <- AddMetaData(neutro.integrated, ES)

colors <- colorRampPalette(c("#FF4B20", "#FFB433", "#C6FDEC", "#7AC5FF", "#0348A6"))

multi_dittoPlot(neutro.integrated, vars = c('HALLMARK_APOPTOSIS', 'HALLMARK_FATTY_ACID_METABOLISM', 'HALLMARK_GLYCOLYSIS', 'HALLMARK_HEME_METABOLISM', 'HALLMARK_HYPOXIA', 'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_MTORC1_SIGNALING', 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', 'HALLMARK_PEROXISOME', 'HALLMARK_PI3K_AKT_MTOR_SIGNALING', 'HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY'), 
                group.by = 'seurat_clusters', plots = c('jitter', 'vlnplot', 'boxplot'), 
                ylab = 'Enrichment Scores', 
                theme = theme_classic() + theme(plot.title = element_text(size = 10)))



ES2 <- data.frame(neutro.integrated[[]], Idents(neutro.integrated))
colnames(ES2)[ncol(ES2)] <- 'cluster'
ridgeEnrichment(ES2, gene.set = 'HALLMARK_HYPOXIA', group = 'seurat_clusters', add.rug = TRUE) # Add 'facet = 'Sample.ID'' to split by day



PCA <- performPCA(enriched = ES2, groups = c('seurat_clusters', 'Sample.ID'))
pcaEnrichment(PCA, PCx = 'PC1', PCy = 'PC2', contours = TRUE)
pcaEnrichment(PCA, PCx = 'PC1', PCy = 'PC2', contours = FALSE, facet = 'seurat_clusters')

output <- getSignificance(ES2, group = 'seurat_clusters', fit = 'ANOVA') # Can use linear.model, T.test, or ANOVA



# Complex Heatmap ------

## Average single cell data & pull out normalized counts
## Counts are normalized by dividing the counts for a given feature by the total counts per cell, multiplying by a scale factor (default == 10,000), and then taking the natural log using log1p()
cellplex.wt <- subset(cellplex.integrated, Sample_origin =='wt')
cellplex.rag <- subset(cellplex.integrated, Sample_origin =='rag')
cellplex.ragTh1 <- subset(cellplex.integrated, Sample_origin =='rag.th1')
cellplex.ragTh17 <- subset(cellplex.integrated, Sample_origin =='rag.th17')

cellplex.avg <- AverageExpression(cellplex.integrated, return.seurat = T)
cellplex.wt.avg <- AverageExpression(cellplex.wt, return.seurat = T)
cellplex.rag.avg <- AverageExpression(cellplex.rag, return.seurat = T)
cellplex.ragTh1.avg <- AverageExpression(cellplex.ragTh1, return.seurat = T)
cellplex.ragTh17.avg <- AverageExpression(cellplex.ragTh17, return.seurat = T)

avg.norm_counts <- cellplex.avg@assays$RNA@data
wt.avg.norm_counts <- cellplex.wt.avg@assays$RNA@data
rag.avg.norm_counts <- cellplex.rag.avg@assays$RNA@data
ragTh1.avg.norm_counts <- cellplex.ragTh1.avg@assays$RNA@data
ragTh17.avg.norm_counts <- cellplex.ragTh17.avg@assays$RNA@data

avg.norm_counts <- as.data.frame(avg.norm_counts)
wt.avg.norm_counts <- as.data.frame(wt.avg.norm_counts)
rag.avg.norm_counts <- as.data.frame(rag.avg.norm_counts)
ragTh1.avg.norm_counts <- as.data.frame(ragTh1.avg.norm_counts)
ragTh17.avg.norm_counts <- as.data.frame(ragTh17.avg.norm_counts)

avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')
wt.avg.norm_counts <- rownames_to_column(wt.avg.norm_counts, var = 'gene')
rag.avg.norm_counts <- rownames_to_column(rag.avg.norm_counts, var = 'gene')
ragTh1.avg.norm_counts <- rownames_to_column(ragTh1.avg.norm_counts, var = 'gene')
ragTh17.avg.norm_counts <- rownames_to_column(ragTh17.avg.norm_counts, var = 'gene')

## Grabbing pathway sets
# GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
# CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
# HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

## Grabbing gene lists for desired pathways
# pathways <- CP.set %>%
#   filter(gs_name == 'KEGG_GLYCOLYSIS_GLUCONEOGENESIS') # Set this to be whatever pathway you'd like

# gene_list <- pathways$human_gene_symbol
gene_list <- c('Ilb1',
               'Clec4e',
               'Junb',
               'Ctsd',
               'Wfdc17',
               'Il1f9',
               'Pla2g7',
               'Arg2',
               'Cd84',
               'Lcn2',
               'Prdx5',
               'Ngp',
               'Camp',
               'Ltf',
               'Arhgdib',
               'Anxa1',
               'Plbd1',
               'Tkt',
               'Aldh2',
               'Ly6c2',
               'Adpgk',
               'Cd177')

## Filtering average Seurat object for genes in pathway
heatmap.data <- avg.norm_counts %>%
  filter(gene %in% gene_list)
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]

## Create vectors for column/row annotations
# mdsc_pmn <- c('')
gene_anno <- c('MDSC',
               'MDSC',
               'PMN',
               'PMN',
               'MDSC',
               'PMN',
               'PMN',
               'PMN',
               'MDSC',
               'MDSC',
               'PMN',
               'PMN',
               'PMN',
               'PMN',
               'MDSC',
               'MDSC',
               'PMN',
               'PMN',
               'MDSC',
               'PMN',
               'PMN') # Check to make sure annotation order matches how the genes appear in the matrix

which(heatmap.data == max(heatmap.data), arr.ind = TRUE) # Use this to find the max value within the matrix so that you can set your upper bound
which(heatmap.data == min(heatmap.data), arr.ind = TRUE) # Use this to find the min value within the matrix so that you can set your lower bound

col_fun = circlize::colorRamp2(c(0, 5), c("white", "blue4")) # Check to make sure that the upper bound of this range isn't cutting off the expression of some genes in the matrix
ComplexHeatmap::Heatmap(heatmap.data,
                        name = 'Average\nExpression',
                        col = col_fun,
                        rect_gp = grid::gpar(col = 'black', lwd = 2),
                        cluster_columns = T,
                        cluster_column_slices = F,
                        column_gap = grid::unit(5, 'mm'),
                        row_split = gene_anno,
                        cluster_rows = T,
                        cluster_row_slices = F,
                        row_gap = grid::unit(5, 'mm'),
                        show_parent_dend_line = F,
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid::grid.text(sprintf("%.1f", heatmap.data[i, j]), x, y, gp = grid::gpar(fontsize = 10))
                        })



# Top pathways by NES
pathways <- pathways %>%
  filter(str_detect(cluster_name, 'Granulocytes [:digit:]'))

top_pathways <- pathways %>%
  arrange(desc(NES)) %>%
  group_by(cluster_num, sample_subset) %>%
  do(head(., n = 5))

pathway_list <- unique(top_pathways$pathway)

heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- reshape2::dcast(heatmap.data, pathway ~ cluster_name + sample_subset, value.var = 'NES', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'pathway')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- reshape2::dcast(heatmap2.data, pathway ~ cluster_name + sample_subset, value.var = 'padj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'pathway')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)



row_order <- c(
  'Granulocytes 1_blood',
  'Granulocytes 1_tissue',
  
  'Granulocytes 2_blood',
  'Granulocytes 2_tissue',
  
  'Granulocytes 3_blood',
  'Granulocytes 3_tissue',
  
  'Granulocytes 4_blood',
  'Granulocytes 4_tissue',
  
  'Granulocytes 5_blood',
  'Granulocytes 5_tissue',
  
  'Granulocytes 6_blood',
  'Granulocytes 6_tissue',
  
  'Granulocytes 7_blood',
  'Granulocytes 7_tissue',
  
  'Granulocytes 8_blood',
  'Granulocytes 8_tissue',
  
  'Granulocytes 9_blood',
  'Granulocytes 9_tissue',
  
  'Granulocytes 10_blood',
  'Granulocytes 10_tissue',
  
  'Granulocytes 11_blood',
  'Granulocytes 11_tissue',
  
  'Granulocytes 12_blood',
  'Granulocytes 12_tissue'
)



heatmap.data <- heatmap.data[row_order, ]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]



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



row_names <- c(
  'G1_Blood',
  'G1_Tissue',
  
  'G2_Blood',
  'G2_Tissue',
  
  'G3_Blood',
  'G3_Tissue',
  
  'G4_Blood',
  'G4_Tissue',
  
  'G5_Blood',
  'G5_Tissue',
  
  'G6_Blood',
  'G6_Tissue',
  
  'G7_Blood',
  'G7_Tissue',
  
  'G8_Blood',
  'G8_Tissue',
  
  'G9_Blood',
  'G9_Tissue',
  
  'G10_Blood',
  'G10_Tissue',
  
  'G11_Blood',
  'G11_Tissue',
  
  'G12_Blood',
  'G12_Tissue'
)



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

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = circlize::colorRamp2(c(-htmp_range, 0, htmp_range), c("blue", "white", "red"))

ht <- ComplexHeatmap::Heatmap(
  heatmap.data,
  name = "NES",
  column_title = "Top enriched pathways per granulocyte cluster",
  col = col_fun,
  rect_gp = grid::gpar(col = "black", lwd = 2),
  cluster_columns = TRUE,
  cluster_column_slices = FALSE,
  column_gap = grid::unit(5, "mm"),
  column_names_gp = grid::gpar(fontsize = 8),
  cluster_rows = FALSE,
  cluster_row_slices = FALSE,
  row_split = row_anno,
  row_gap = grid::unit(5, "mm"),
  show_parent_dend_line = FALSE,
  heatmap_legend_param = list(border = "black"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(heatmap2.data[i, j] < 0.0001) {
      grid::grid.text("****", x, y)
    } else if(heatmap2.data[i, j] < 0.001) {
      grid::grid.text("***", x, y)
    } else if(heatmap2.data[i, j] < 0.01) {
      grid::grid.text("**", x, y)
    } else if(heatmap2.data[i, j] < 0.05) {
      grid::grid.text("*", x, y)
    }
  }
)

ComplexHeatmap::draw(ht, padding = unit(c(55, 5, 2, 5), "mm")) # bottom, left, top, right paddings



# Trajectory Analysis -----

# Finding trajectories =====

## Inferring trajectories
object_counts <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@data), 'sparseMatrix'))
neutro.integrated_dyn <- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

rm(object_counts, object_expression)

## Add a dimensionality reduction
neutro.integrated_dimred <- dyndimred::dimred_umap(neutro.integrated_dyn$expression)

## Infer the trajectory
neutro.integrated_model <- infer_trajectory(neutro.integrated_dyn, ti_slingshot(), verbose = T)

## Plot trajectory & pseudotime
neutro.integrated_milestone.umap <- plot_dimred(neutro.integrated_model, label_milestones = T, dimred = neutro.integrated_dimred, hex_cells = F) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'UMAP 1', y = 'UMAP 2') ## Check this prior to rooting
neutro.integrated_milestone.pca <- plot_dimred(neutro.integrated_model, label_milestones = T, dimred = 'pca', hex_cells = F) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'PC 1', y = 'PC 2') ## Check this prior to rooting
neutro.integrated_traj.umap <- plot_dimred(neutro.integrated_model, dimred = neutro.integrated_dimred, grouping = neutro.integrated@active.ident, color_density = 'grouping', hex_cells = F) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'UMAP 1', y = 'UMAP 2')
neutro.integrated_traj.pca <- plot_dimred(neutro.integrated_model, dimred = 'pca', grouping = neutro.integrated@active.ident, color_density = 'grouping', hex_cells = F) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'PC 1', y = 'PC 2')
neutro.integrated_traj2.umap <- plot_dimred(neutro.integrated_model, dimred = neutro.integrated_dimred, grouping = neutro.integrated@meta.data$Sample.ID, color_density = 'grouping', hex_cells = F) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'UMAP 1', y = 'UMAP 2')
neutro.integrated_traj2.pca <- plot_dimred(neutro.integrated_model, dimred = 'pca', grouping = neutro.integrated@meta.data$Sample.ID, color_density = 'grouping', hex_cells = F) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'PC 1', y = 'PC 2')
neutro.integrated_pseudo.umap <- plot_dimred(neutro.integrated_model, "pseudotime", pseudotime = calculate_pseudotime(neutro.integrated_model), dimred = neutro.integrated_dimred, hex_cells = F) + theme_classic() + labs(x = 'UMAP 1', y = 'UMAP 2') + theme(axis.text = element_blank(), axis.ticks = element_blank())
neutro.integrated_pseudo.pca <- plot_dimred(neutro.integrated_model, "pseudotime", pseudotime = calculate_pseudotime(neutro.integrated_model), dimred = 'pca', hex_cells = F) + theme_classic() + labs(x = 'PC 1', y = 'PC 2') + theme(axis.text = element_blank(), axis.ticks = element_blank())

## Root trajectory if necessary
neutro.integrated_model <- add_root(neutro.integrated_model, root_milestone_id = "4")

## Simplify trajectory
simp <- simplify_trajectory(neutro.integrated_model)
simp_lab <- simp %>% label_milestones(c('1' = 'Stuck', '3' = 'Bifurcation', '4' = 'Mature', '5' = 'Immature'))



# Plotting gene expression over pseudotime =====

# Slingshot #####

## Get trajectory & clustering information for lineage acquisition (pseudotime)
expression <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@data), 'sparseMatrix'))
ndim <- 20L
max_clusters <- min(nrow(expression)-1, 10)
pca <- irlba::prcomp_irlba(expression, n = ndim)

# Select optimal number of dimensions if ndim is large enough
if (ndim > 3) {
  # This code is adapted from the expermclust() function in TSCAN
  # The only difference is in how PCA is performed
  # (they specify scale. = TRUE and we leave it as FALSE)
  x <- 1:ndim
  optpoint1 <- which.min(sapply(2:10, function(i) {
    x2 <- pmax(0, x - i)
    sum(lm(pca$sdev[1:ndim] ~ x + x2)$residuals^2 * rep(1:2,each = 10))
  }))
  
  # This is a simple method for finding the "elbow" of a curve, from
  # https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
  x <- cbind(1:ndim, pca$sdev[1:ndim])
  line <- x[c(1, nrow(x)),]
  proj <- princurve::project_to_curve(x, line)
  optpoint2 <- which.max(proj$dist_ind)-1
  
  # We will take more than 3 PCs only if both methods recommend it
  optpoint <- max(c(min(c(optpoint1, optpoint2)), 3))
}

dimred <- pca$x[, seq_len(optpoint)]
rownames(dimred) <- rownames(expression)

clusterings <- lapply(3:max_clusters, function(K){
  cluster::pam(dimred, K) # we generally prefer PAM as a more robust alternative to k-means
})
wh.cl <- which.max(sapply(clusterings, function(x){ x$silinfo$avg.width })) + 1
labels <- clusterings[[min(c(wh.cl, 8))]]$clustering

rm(ndim, optpoint, optpoint1, optpoint2, x, pca, proj, line, expression, clusterings, wh.cl, max_clusters)



## Find trajectory
lineages <- getLineages(dimred, labels, start.clus = '4')
sling <- getCurves(lineages, shrink = 1L, reweight = T, reassign = T, thresh = 0.001, maxit = 10L, stretch = 2L, smoother = 'smooth.spline', shrink.method = 'cosine') ## Parameters taken from dynverse ti_slingshot code

rm(dimred, labels, lineages)



# TradeSeq #####

## Find optimal number of knots to use for GAM fitting
counts <- neutro.integrated@assays$RNA@counts
filt <- rowSums(counts > 1) >= 120 # Filtering transcripts (genes w/count of at least two in at least 120 different cells)
counts.filt <- counts[filt, ]
rm(filt)

icMat <- evaluateK(counts = as.matrix(counts.filt), sds = sling, k = 3:15, nGenes = 300, verbose = T, plot = T)

## Fit GAM
sce <- fitGAM(as.matrix(counts.filt), sling, nknots = 10) # SCE method
converge <- table(rowData(sce)$tradeSeq$converged) # Check whether genes converged

# List method
# control <- gam.control()
# control$maxit <- 1000 # Set maximum number of iterations to 1,000
# gamList <- fitGAM(counts = as.matrix(counts.filt), pseudotime = slingPseudotime(sling, na = F), cellWeights = slingCurveWeights(sling), control = control, sce = F)
# pvalLineage <- getSmootherPvalues(gamList)
# statLineage <- getSmootherTestStats(gamList)

## Testing
assoRes <- associationTest(sce, lineages = T) # Testing whether genes are significantly changed along pseudotime (independent lineages)
startRes <- startVsEndTest(sce, lineages = T) # Testing whether genes are significantly changed between the start & end of pseudotime (independent lineages)
endRes <- diffEndTest(sce) # Testing whether genes are significantly changed at end of pseudotime
patternRes <- patternTest(sce) # Testing whether genes have significantly different expression patterns throughout pseudotime
earlyDE.plot <- plotGeneCount(curve = sling, counts = as.matrix(counts.filt), clusters = apply(slingClusterLabels(sling), 1, which.max), models = sce) # Visualize where knots are
earlyDERes <- earlyDETest(sce, knots = c(1, 3))

## Identifying significant changes in assoRes
assoRes.sig.lin1 <- rownames(assoRes)[which(p.adjust(assoRes$pvalue_1, "fdr") <= 0.05)]
assoRes.sig.lin2 <- rownames(assoRes)[which(p.adjust(assoRes$pvalue_2, "fdr") <= 0.05)]

assoRes.sig.upset <- upset(fromList(list('Lineage 1' = assoRes.sig.lin1 , 'Lineage 2' = assoRes.sig.lin2)))

yhatSmooth.1 <- predictSmooth(sce, gene = assoRes.sig.lin1, nPoints = 50)
heatSmooth.1 <- pheatmap(t(scale(t(yhatSmooth.1[, 1:50]))), cluster_cols = F, show_rownames = F, show_colnames = F)

yhatSmooth.2 <- predictSmooth(sce, gene = assoRes.sig.lin2, nPoints = 50)
heatSmooth.2 <- pheatmap(t(scale(t(yhatSmooth.2[, 1:50]))), cluster_cols = F, show_rownames = F, show_colnames = F)

## Combining End & Pattern testing
patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes.comp <- patternRes[, c('Gene', 'pattern')]

endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
endRes.comp <- endRes[, c('Gene', 'end')]

compare <- merge(patternRes.comp, endRes.comp, by = 'Gene', all = F)
compare$transientScore <- rank(-compare$end, ties.method = 'min')^2 + rank(compare$pattern, ties.method = 'random')^2

rm(patternRes.comp, endRes.comp)

## Plotting
## Expression vs Pseudotime
plotSmoothers(sce, assays(sce)$counts, gene = 'Ltf')

## End vs Pattern
ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = 'patternTest Wald Statistic (log scale)',
       y = 'diffEndTest Wald Statistic (log scale)') +
  scale_color_continuous(low = 'grey', high = 'blue') +
  theme_classic()

## Clustering genes with similar expression along pseudotime
nPointsClus <- 50
clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus, genes = rownames(as.matrix(counts.filt)))
clusterLabels <- primaryCluster(clusPat$rsec)
cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes

gene <- rownames(as.matrix(counts.filt))
Gene.Cluster <- as.data.frame(gene)
Gene.Cluster$cluster <- primaryCluster(clusPat$rsec)
Gene.Cluster <- arrange(Gene.Cluster, cluster)
write.csv(Gene.Cluster, here('Integrated', 'TradeSeq', 'Real time', 'gene cluster', 'full_cluster.csv'))

for (i in seq_along(unique(Gene.Cluster$cluster))) {
  if(i != length(Gene.Cluster$cluster)) {
    clust <- Gene.Cluster %>% filter(Gene.Cluster$cluster %in% c(unique(Gene.Cluster$cluster)[i]))
    write.csv(clust, here('Integrated', 'TradeSeq', 'Real time', 'gene cluster', paste0('cluster_', i-1, '.csv')))
  }
}

rm(clust, i)

for (xx in sort(cUniq)[1:6]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0('Cluster ', xx),  x = 'Pseudotime', y = 'Normalized expression') +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 2),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:1, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  p <- p + guides(color = FALSE) +
    scale_color_manual(values = c('#440154FF', '#FDE725FF'),
                       breaks = c('0', '1'))
  assign(paste0('c', xx), p)
}

rm(geneId, ii, xx, p, cId, nPointsClus, cUniq, clusterLabels)



# TradeSeq GSEA #####

## fgsea
geneSets <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
# geneSets <- msigdbr(species = 'Mus musculus', category = 'C2', subcategory = 'CP:KEGG') # KEGG
# geneSets <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark
geneSets$gene_symbol <- toupper(geneSets$gene_symbol)
m_list <- geneSets %>% split(x = .$gene_symbol, f = .$gs_name)
stats.lin1 <- assoRes$waldStat_1
stats.lin2 <- assoRes$waldStat_2
names(stats.lin1) <- toupper(rownames(assoRes))
names(stats.lin2) <- toupper(rownames(assoRes))
eaRes.lin1 <- fgsea(pathways = m_list, stats = stats.lin1, eps = 0.0, scoreType = 'pos', minSize = 10)
eaRes.lin2 <- fgsea(pathways = m_list, stats = stats.lin2, eps = 0.0, scoreType = 'pos', minSize = 10)
eaRes.lin1 <- filter(eaRes.lin1, padj <= 0.05)
eaRes.lin2 <- filter(eaRes.lin2, padj <= 0.05)
eaRes.lin1_unique <- anti_join(eaRes.lin1, eaRes.lin2, by = 'pathway')
eaRes.lin2_unique <- anti_join(eaRes.lin2, eaRes.lin1, by = 'pathway')
fwrite(eaRes.lin1, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_lineage 1 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))
fwrite(eaRes.lin2, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_lineage 2 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))
fwrite(eaRes.lin1_unique, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_unique lineage 1 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))
fwrite(eaRes.lin2_unique, file = here('Integrated', 'Pseudotime', 'Real time', 'TradeSeq', 'fgsea', 'GO_unique lineage 2 enrichment analysis.tsv'), sep="\t", sep2=c("", " ", ""))

rm(geneSets, m_list, stats.lin1, stats.lin2)



# Notes -----

# If you'd like to do everything in three dimensions, use RumUMAP(object, dims = 1:n, n.components = 3L). Bringing this into
# Slingshot requires you to create objects that contain the cell embeddings (object@reductions$umap@cell.embeddings) & the
# clustering info (object@active.iden (or any grouping you desire)). You can then supply these to slingshot as you normally
# would. If you would like to visualize, you could do the following (requires rgl package):

# gg_color_hue <- function(n) {
#   hues = seq(15, 375, length = n + 1)
#   hcl(h = hues, l = 65, c = 100)[1:n]
# }
# plot3d.SlingshotDataSet(sds)
# plot3d(rd, col = gg_color_hue(10)[cl], aspect = 'iso', add = T)


# I'm not sure if using three dimensions has any clear benefit over using two, but it may.



# Test area -----

files_to_read <- list.files(path = 'output/brain/GSEA/monomac', pattern = '\\.tsv$', full.names = T)
all_files <- lapply(files_to_read, function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = T)
})

for(i in 1:length(all_files)) {
  write.csv(all_files[i], here('output', 'brain', 'GSEA', 'monomac', 'converted', paste0(list.files('output/brain/GSEA/monomac')[i], '.csv')))
}



