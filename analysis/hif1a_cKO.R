## Made by Christopher M. Horn, MS
## Kielian Lab data
## Created: 2022-07-01
## Updated: 2022-08-30

## Notes: WT & Hif1a cKO mice were infected & tissue samples were prepared for scRNA-seq at D3 & D14 post-infection; 2-Mar == Marc2 (row #892) & 2-Mar == March2 (row #13629)

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
  'here',
  'sqldf'
  )

invisible(lapply(packages, library, character.only = TRUE))

rm(packages)

## Set options
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(12345)

# Initializing objects -----
## Setting up Seurat objects
## Load in the data
d3_null.data <- Read10X(data.dir = here('raw_data', 'd3_null'))
d3_null.data <- d3_null.data$`Gene Expression`

d3_cre.data <- Read10X(data.dir = here('raw_data', 'd3_cre'))
d3_cre.data <- d3_cre.data$`Gene Expression`

d14_null.data <- Read10X(data.dir = here('raw_data', 'd14_null'))
d14_null.data <- d14_null.data$`Gene Expression`

d14_cre.data <- Read10X(data.dir = here('raw_data', 'd14_cre'))
d14_cre.data <- d14_cre.data$`Gene Expression`

## Initialize Seurat object w/raw data
d3_null <- CreateSeuratObject(
  counts = d3_null.data,
  project = 'd3_null',
  min.cells = 3,
  min.features = 200
  )

d3_cre <- CreateSeuratObject(
  counts = d3_cre.data,
  project = 'd3_cre',
  min.cells = 3,
  min.features = 200
  )

d14_null <- CreateSeuratObject(
  counts = d14_null.data,
  project = 'd14_null',
  min.cells = 3,
  min.features = 200
  )

d14_cre <- CreateSeuratObject(
  counts = d14_cre.data,
  project = 'd14_cre',
  min.cells = 3,
  min.features = 200
  )

## Pre-processing and QC
d3_null[['percent.mt']] <- PercentageFeatureSet(d3_null, pattern = '^mt-')
d3_null <- subset(d3_null, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

d3_cre[['percent.mt']] <- PercentageFeatureSet(d3_cre, pattern = '^mt-')
d3_cre <- subset(d3_cre, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

d14_null[['percent.mt']] <- PercentageFeatureSet(d14_null, pattern = '^mt-')
d14_null <- subset(d14_null, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

d14_cre[['percent.mt']] <- PercentageFeatureSet(d14_cre, pattern = '^mt-')
d14_cre <- subset(d14_cre, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Create timing metadata
sample.id_1 <- rep('d3_null', length(d3_null@meta.data$orig.ident))
sample.id_2 <- rep('d3_cre', length(d3_cre@meta.data$orig.ident))
sample.id_3 <- rep('d14_null', length(d14_null@meta.data$orig.ident))
sample.id_4 <- rep('d14_cre', length(d14_cre@meta.data$orig.ident))

names(sample.id_1) <- rownames(d3_null@meta.data)
names(sample.id_2) <- rownames(d3_cre@meta.data)
names(sample.id_3) <- rownames(d14_null@meta.data)
names(sample.id_4) <- rownames(d14_cre@meta.data)

d3_null <- AddMetaData(
  d3_null,
  sample.id_1,
  col.name = 'sample_origin'
  )

d3_cre <- AddMetaData(
  d3_cre,
  sample.id_2,
  col.name = 'sample_origin'
  )

d14_null <- AddMetaData(
  d14_null,
  sample.id_3,
  col.name = 'sample_origin'
  )

d14_cre <- AddMetaData(
  d14_cre,
  sample.id_4,
  col.name = 'sample_origin'
  )

## Write individual object metadata to file
write.csv(d3_null@meta.data, here('output', 'd3_null metadata.csv'))
write.csv(d3_cre@meta.data, here('output', 'd3_cre metadata.csv'))
write.csv(d14_null@meta.data, here('output', 'd14_null metadata.csv'))
write.csv(d14_cre@meta.data, here('output', 'd14_cre metadata.csv'))

## Write QC metrics
d3_null.QC_mets <- dim(d3_null)
names(d3_null.QC_mets) <- c('features', 'cells')

d3_cre.QC_mets <- dim(d3_cre)
names(d3_cre.QC_mets) <- c('features', 'cells')

d14_null.QC_mets <- dim(d14_null)
names(d14_null.QC_mets) <- c('features', 'cells')

d14_cre.QC_mets <- dim(d14_cre)
names(d14_cre.QC_mets) <- c('features', 'cells')

write.csv(d3_null.QC_mets, here('output', 'QC', 'd3_null QC_mets.csv'))
write.csv(d3_cre.QC_mets, here('output', 'QC', 'd3_cre QC_mets.csv'))
write.csv(d14_null.QC_mets, here('output', 'QC', 'd14_null QC_mets.csv'))
write.csv(d14_cre.QC_mets, here('output', 'QC', 'd14_cre QC_mets.csv'))

d3_null.QC_mets.plot  <- VlnPlot(
  d3_null,
  features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
  ncol = 3
  )

d3_cre.QC_mets.plot  <- VlnPlot(
  d3_cre,
  features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
  ncol = 3
  )

d14_null.QC_mets.plot  <- VlnPlot(
  d14_null,
  features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
  ncol = 3
  )

d14_cre.QC_mets.plot  <- VlnPlot(
  d14_cre,
  features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'),
  ncol = 3
  )

ggsave(
  'd3_null QC_mets.plot.png',
  plot = d3_null.QC_mets.plot,
  device = 'png',
  path = here('output', 'QC')
  )

ggsave(
  'd3_cre QC_mets.plot.png',
  plot = d3_cre.QC_mets.plot,
  device = 'png',
  path = here('output', 'QC')
  )

ggsave(
  'd14_null QC_mets.plot.png',
  plot = d14_null.QC_mets.plot,
  device = 'png',
  path = here('output', 'QC')
  )

ggsave(
  'd14_cre QC_mets.plot.png',
  plot = d14_cre.QC_mets.plot,
  device = 'png',
  path = here('output', 'QC')
  )

## Remove temp objects
rm(
  d3_null.data,
  d3_cre.data,
  d14_null.data,
  d14_cre.data,
  sample.id_1,
  sample.id_2,
  sample.id_3,
  sample.id_4,
  d3_null.QC_mets,
  d3_null.QC_mets.plot,
  d3_cre.QC_mets,
  d3_cre.QC_mets.plot,
  d14_null.QC_mets,
  d14_null.QC_mets.plot,
  d14_cre.QC_mets,
  d14_cre.QC_mets.plot
  )

gc()

# Annotating cell types -----
## Annotate the cells w/SingleR
d3_null <- NormalizeData(d3_null)
d3_cre <- NormalizeData(d3_cre)
d14_null <- NormalizeData(d14_null)
d14_cre <- NormalizeData(d14_cre)

ref.se <- ImmGenData() # snapshot date: 2020-10-27

d3_null_sce <- as.SingleCellExperiment(d3_null)
d3_cre_sce <- as.SingleCellExperiment(d3_cre)
d14_null_sce <- as.SingleCellExperiment(d14_null)
d14_cre_sce <- as.SingleCellExperiment(d14_cre)

commonGenes.1 <- intersect(rownames(d3_null_sce), rownames(ref.se))
commonGenes.2 <- intersect(rownames(d3_cre_sce), rownames(ref.se))
commonGenes.3 <- intersect(rownames(d14_null_sce), rownames(ref.se))
commonGenes.4 <- intersect(rownames(d14_cre_sce), rownames(ref.se))

ref.se_1 <- ref.se[commonGenes.1,]
ref.se_2 <- ref.se[commonGenes.2,]
ref.se_3 <- ref.se[commonGenes.3,]
ref.se_4 <- ref.se[commonGenes.4,]

d3_null_sce <- d3_null_sce[commonGenes.1,]
d3_cre_sce <- d3_cre_sce[commonGenes.2,]
d14_null_sce <- d14_null_sce[commonGenes.3,]
d14_cre_sce <- d14_cre_sce[commonGenes.4,]

pred.d3_null <- SingleR(
  test = d3_null_sce,
  ref = ref.se_1,
  labels = ref.se_1$label.main)

pred.d3_cre <- SingleR(
  test = d3_cre_sce,
  ref = ref.se_2,
  labels = ref.se_2$label.main)

pred.d14_null <- SingleR(
  test = d14_null_sce,
  ref = ref.se_3,
  labels = ref.se_3$label.main)

pred.d14_cre <- SingleR(
  test = d14_cre_sce,
  ref = ref.se_4,
  labels = ref.se_4$label.main)

d3_null[['celltype']] <- pred.d3_null$pruned.labels
d3_cre[['celltype']] <- pred.d3_cre$pruned.labels
d14_null[['celltype']] <- pred.d14_null$pruned.labels
d14_cre[['celltype']] <- pred.d14_cre$pruned.labels

## Remove temp objects
rm(
  ref.se,
  d3_null_sce,
  d3_cre_sce,
  d14_null_sce,
  d14_cre_sce,
  commonGenes.1,
  commonGenes.2,
  commonGenes.3,
  commonGenes.4,
  ref.se_1,
  ref.se_2,
  ref.se_3,
  ref.se_4,
  pred.d3_null,
  pred.d3_cre,
  pred.d14_null,
  pred.d14_cre
  )

gc()

# Object integration -----
## Integrating all cells from all days
sample.list <- c(
  d3_null,
  d3_cre,
  d14_null,
  d14_cre
  )
names(sample.list) <- c(
  'd3_null',
  'd3_cre',
  'd14_null',
  'd14_cre'
  )
for (i in 1:length(sample.list))
  {sample.list[[i]] <- SCTransform(sample.list[[i]], verbose = TRUE)}
sample.features <- SelectIntegrationFeatures(object.list = sample.list, nfeatures = 3000)
sample.list <- PrepSCTIntegration(
  object.list = sample.list,
  anchor.features = sample.features,
  verbose = TRUE
  )
sample.anchors <- FindIntegrationAnchors(
  object.list = sample.list,
  normalization.method = 'SCT',
  anchor.features = sample.features,
  verbose = TRUE
  )
integrated <- IntegrateData(
  anchorset = sample.anchors,
  normalization.method = 'SCT',
  verbose = TRUE
  )
integrated <- RunPCA(integrated, verbose = TRUE)
integrated <- RunUMAP(integrated, dims = 1:30)
integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5)

saveRDS(integrated, here('output', 'integrated.rds'))

## Remove temp objects
rm(
  sample.list,
  sample.features,
  sample.anchors,
  d3_null,
  d3_cre,
  d14_null,
  d14_cre,
  i
  )

gc()

## Output celltype composition of each cluster
sample.comp_origin <- table(integrated$sample_origin)
clust.comp_origin <- table(Idents(integrated), integrated$sample_origin)
clust.comp_celltype <- t(prop.table(table(Idents(integrated), integrated$celltype), margin = 1))

write.csv(sample.comp_origin, here('output', 'sample.comp_origin.csv'))
write.csv(clust.comp_origin, here('output', 'clust.comp_origin.csv'))
write.csv(clust.comp_celltype, here('output', 'clust.comp_celltype.csv'))

rm(
  sample.comp_origin,
  clust.comp_origin,
  clust.comp_celltype
  )

gc()

# Renaming/plotting -----
## Renaming clusters
new.cluster.ids <- c(
  'Granulocytes 1',
  'Granulocytes 2',
  'Granulocytes 3',
  'Granulocytes 4',
  'Granulocytes 5',
  'Granulocytes 6',
  'Granulocytes 7',
  'Granulocytes 8',
  'Granulocyte/Monocyte/Macrophage',
  'Lymphocytes',
  'Granulocytes 9'
  )

names(new.cluster.ids) <- levels(integrated)
integrated <- RenameIdents(integrated, new.cluster.ids)

rm(new.cluster.ids)

gc()

## Save seurat objects
saveRDS(integrated, here('output', 'integrated.rds'))

## Plotting w/new labels
integrated_pca.origin <- DimPlot(
  integrated,
  group.by = 'sample_origin',
  reduction = 'pca',
  cols = c('#BC3C2999', '#0072B599', '#E1872799', '#20854E99'),
  pt.size = 1.5) +
  theme_classic() +
  ggtitle('Split by sample origin') +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank()
  )

integrated_pca.clus <- DimPlot(
  integrated,
  reduction = 'pca',
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 1.5) +
  theme_classic() +
  NoLegend() +
  labs(x = 'PC 1', y = 'PC 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank()
  )

integrated_umap.origin <- DimPlot(
  integrated,
  group.by = 'sample_origin',
  cols = c('#BC3C2999', '#0072B599', '#E1872799', '#20854E99'),
  pt.size = 1.5) +
  theme_classic() +
  ggtitle('Split by sample origin') +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank()
  )

integrated_umap.clus <- DimPlot(
  integrated,
  label = TRUE,
  label.box = TRUE,
  repel = TRUE,
  pt.size = 1.5) +
  theme_classic() +
  NoLegend() +
  labs(x = 'UMAP 1', y = 'UMAP 2') +
  theme(axis.text = element_blank(), axis.ticks = element_blank()
  )

ggsave(
  'integrated_pca origin.png',
  plot = integrated_pca.origin,
  device = 'png',
  path = here('output')
  )

ggsave(
  'integrated_pca cluster.png',
  plot = integrated_pca.clus,
  device = 'png',
  path = here('output')
  )

ggsave(
  'integrated_umap origin.png',
  plot = integrated_umap.origin,
  device = 'png',
  path = here('output')
  )

ggsave(
  'integrated_umap cluster.png',
  plot = integrated_umap.clus,
  device = 'png',
  path = here('output')
  )

rm(
  integrated_pca.origin,
  integrated_pca.clus,
  integrated_umap.origin,
  integrated_umap.clus
  )

gc()

## Output cluster number to name cheat sheet
num2name <- data.frame(
  'Cluster_name' = levels(integrated),
  'Cluster_num' = levels(integrated$seurat_clusters)
  )

write.csv(num2name, here('output', 'num 2 name cheat sheet.csv'))

## Write the metadata to a file
write.csv(integrated@meta.data, here('output', 'integrated_metadata.csv'))

rm(num2name)

gc()

## Output origin for each cluster
sample.comp_origin <- t(table(integrated$Sample_origin))
clust.comp_origin <- table(Idents(integrated), integrated$Sample_origin)
clust.comp_celltype <- t(table(integrated$celltype, Idents(integrated)))

write.csv(sample.comp_origin, here('output', 'sample.comp_origin.csv'))
write.csv(clust.comp_origin, here('output', 'clust.comp_origin.csv'))
write.csv(clust.comp_celltype, here('output', 'clust.comp_celltype.csv'))

rm(
  sample.comp_origin,
  clust.comp_origin,
  clust.comp_celltype
  )

gc()

## Subsetting out
integrated.granulocyte <- subset(brain.integrated, idents = c(
  'Microglia 1',
  'Microglia 2',
  'Microglia 3')
  )

# Cluster-level DE -----
## Between cluster DE Analysis
DefaultAssay(integrated) <- 'RNA'
integrated <- NormalizeData(integrated, verbose = TRUE)
integrated <- ScaleData(integrated, verbose = TRUE)

DE <- FindAllMarkers(integrated, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')

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
c01 <- subset(integrated, idents = c('Granulocytes 1'))
c02 <- subset(integrated, idents = c('Granulocytes 2'))
c03 <- subset(integrated, idents = c('Granulocytes 3'))
c04 <- subset(integrated, idents = c('Granulocytes 4'))
c05 <- subset(integrated, idents = c('Granulocytes 5'))
c06 <- subset(integrated, idents = c('Granulocytes 6'))
c07 <- subset(integrated, idents = c('Granulocytes 7'))
c08 <- subset(integrated, idents = c('Granulocytes 8'))
c09 <- subset(integrated, idents = c('Granulocyte/Monocyte/Macrophage'))
c10 <- subset(integrated, idents = c('Lymphocytes'))
c11 <- subset(integrated, idents = c('Granulocytes 9'))

c.list <- list(
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
  c11
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
  c09,
  c10,
  c11
  )

gc()

for(i in 1:length(c.list)){
  DefaultAssay(c.list[[i]]) <- 'RNA'
  c.list[[i]] <- NormalizeData(c.list[[i]], verbose = TRUE)
  c.list[[i]] <- ScaleData(c.list[[i]], verbose = TRUE)
  c.list[[i]]$cell_origin <- paste(Idents(c.list[[i]]), c.list[[i]]$sample_origin, sep = '_')
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

## fgsea
GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # Canonical Pathways
HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

geneSet_list <- list(GO.set, HM.set, CP.set)

for (i in 1:length(geneSet_list)) {
  geneSet_list[[i]]$gene_symbol <- toupper(geneSet_list[[i]]$gene_symbol)
  m_list <- split(x = geneSet_list[[i]]$gene_symbol, f = geneSet_list[[i]]$gs_name)
  for (ii in 1:length(levels(integrated))) {
    glist <- DE %>% filter(cluster == levels(integrated)[ii])
    glist <- column_to_rownames(glist, var = 'gene')
    stats <- glist$avg_log2FC
    names(stats) <- toupper(rownames(glist))
    stats <- sort(stats, decreasing = TRUE)
    eaRes <- fgsea(pathways = m_list, stats = stats, eps = 0.0, minSize = 10, maxSize = 500)
    eaRes <- arrange(eaRes, desc(NES))
    fwrite(
      eaRes,
      file = here(
        'output',
        'GSEA',
        unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '.tsv')), sep="\t", sep2=c("", " ", "")
      )
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
      fwrite(
        eaRes,
        file = here(
          'output',
          'GSEA',
          'within cluster',
          unique(geneSet_list[[i]]$gs_cat), paste0(unique(geneSet_list[[i]]$gs_cat), '_eaRes_', ii-1, '_', iii, '.tsv')), sep="\t", sep2=c("", " ", "")
        )
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

HM <- read_xlsx(here('output', 'GSEA', 'within cluster', 'Hif1a cKO split within cluster pathways.xlsx'), sheet = 'HM pathways')
C2 <- read_xlsx(here('output', 'GSEA', 'within cluster', 'Hif1a cKO split within cluster pathways.xlsx'), sheet = 'C2 pathways')
C5 <- read_xlsx(here('output', 'GSEA', 'within cluster', 'Hif1a cKO split within cluster pathways.xlsx'), sheet = 'C5 pathways')

for(i in 1:length(unique(C5$cluster_num))) {
  d3_null <- C5 %>% filter(cluster_num == i-1) %>% filter(str_detect(cluster_subset, 'd3_null'))
  d3_cre <- C5 %>% filter(cluster_num == i-1) %>% filter(str_detect(cluster_subset, 'd3_cre'))
  d14_null <- C5 %>% filter(cluster_num == i-1) %>% filter(str_detect(cluster_subset, 'd14_null'))
  d14_cre <- C5 %>% filter(cluster_num == i-1) %>% filter(str_detect(cluster_subset, 'd14_cre'))
  
  d3_sql <- sqldf(
    '
    SELECT d3_null.pathway, d3_null.NES AS d3_null_NES, d3_cre.NES AS d3_cre_NES, ABS(d3_null.NES-d3_cre.NES) AS NES_abs_diff
    FROM d3_null
    INNER JOIN d3_cre ON d3_null.pathway = d3_cre.pathway
    ORDER BY NES_abs_diff DESC
    '
  )
  
  d14_sql <- sqldf(
    '
    SELECT d14_null.pathway, d14_null.NES AS d14_null_NES, d14_cre.NES AS d14_cre_NES, ABS(d14_null.NES-d14_cre.NES) AS NES_abs_diff
    FROM d14_null
    INNER JOIN d14_cre ON d14_null.pathway = d14_cre.pathway
    ORDER BY NES_abs_diff DESC
    '
  )
  
  write.csv(
    d3_sql,
    file = here(
      'output',
      'GSEA',
      'within cluster',
      'NES diff',
      'C5',
      paste0('Hif1a cKO d3 NES diff_G', unique(C5$cluster_num)[i], '.csv')
    )
  )

  write.csv(
    d14_sql,
    file = here(
      'output',
      'GSEA',
      'within cluster',
      'NES diff',
      'C5',
      paste0('Hif1a cKO d14 NES diff_G', unique(C5$cluster_num)[i], '.csv')
    )
  )
}

# Compass ------
DefaultAssay(integrated) <- 'RNA'
integrated <- NormalizeData(integrated, verbose = TRUE)
integrated <- ScaleData(integrated, verbose = TRUE)

first3 <- subset(integrated, idents = c('Granulocytes 1', 'Granulocytes 2', 'Granulocytes 3'))

first3_d3_null <- subset(first3, subset = sample_origin == 'd3_null')
first3_d3_cre <- subset(first3, subset = sample_origin == 'd3_cre')
first3_d14_null <- subset(first3, subset = sample_origin == 'd14_null')
first3_d14_cre <- subset(first3, subset = sample_origin == 'd14_cre')

first3_d3_null <- first3_d3_null@assays$RNA@data
first3_d3_cre <- first3_d3_cre@assays$RNA@data
first3_d14_null <- first3_d14_null@assays$RNA@data
first3_d14_cre <- first3_d14_cre@assays$RNA@data

first3_d3_null <- as.data.frame(first3_d3_null)
first3_d3_cre <- as.data.frame(first3_d3_cre)
first3_d14_null <- as.data.frame(first3_d14_null)
first3_d14_cre <- as.data.frame(first3_d14_cre)

fwrite(first3_d3_null, file = '/Users/christopherhorn/Documents/Research/Compass/hif1a/d3_null/expression.tsv', sep="\t", sep2=c("", " ", ""), row.names = TRUE)
fwrite(first3_d3_cre, file = '/Users/christopherhorn/Documents/Research/Compass/hif1a/d3_cre/expression.tsv', sep="\t", sep2=c("", " ", ""), row.names = TRUE)
fwrite(first3_d14_null, file = '/Users/christopherhorn/Documents/Research/Compass/hif1a/d14_null/expression.tsv', sep="\t", sep2=c("", " ", ""), row.names = TRUE)
fwrite(first3_d14_cre, file = '/Users/christopherhorn/Documents/Research/Compass/hif1a/d14_cre/expression.tsv', sep="\t", sep2=c("", " ", ""), row.names = TRUE)

# In terminal
# compass \
# --data path/to/data.tsv
# --num-processes 10 \
# --species mus_musculus \
# --output-dir path/to/output/dir \
# --temp-dir path/to/temp/dir \
# --microcluster-size x (set x so that data is reduced to ~100 microclusters)

# Comparison volcano plots
compass_results <- read.csv(here('compass', 'd14 null vs d14 cre', 'd14 null vs cre compass analysis results.csv'))

aam <- c(
  'Alanine and aspartate metabolism',
  'Arginine and Proline Metabolism',
  'beta-Alanine metabolism',
  'Cysteine Metabolism',
  'D-alanine metabolism',
  'Folate metabolism',
  'Glutamate metabolism',
  'Glycine, serine, alanine and threonine metabolism',
  'Histidine metabolism',
  'Lysine metabolism',
  'Methionine and cysteine metabolism',
  'Taurine and hypotaurine metabolism',
  'Tryptophan metabolism',
  'Tyrosine metabolism',
  'Urea cycle',
  'Valine, leucine, and isoleucine metabolism'
)

ggplot(compass_results, aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = 'gray', size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 1.3, linetype = 'dashed') +
  annotate('text', x = 2.75, y = 2, label = 'Null') +
  annotate('text', x = -2.75, y = 2, label = 'Cre') +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's d") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  theme_classic()

p1 <- ggplot(compass_results %>% filter(subsystem =='Glycolysis/gluconeogenesis'), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = '#BC3C29FF', size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 1.3, linetype = 'dashed') +
  annotate('text', x = 2.75, y = 2, label = 'Null') +
  annotate('text', x = -2.75, y = 2, label = 'Cre') +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's d") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle('Glycolysis') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA)
    )

p2 <- ggplot(compass_results %>% filter(subsystem =='Citric acid cycle'), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = '#0072B5FF', size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 1.3, linetype = 'dashed') +
  annotate('text', x = 2.75, y = 2, label = 'Null') +
  annotate('text', x = -2.75, y = 2, label = 'Cre') +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's d") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle('TCA cycle') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA)
    )

p3 <- ggplot(compass_results %>% filter(subsystem %in% aam), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = '#E18727FF', size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 1.3, linetype = 'dashed') +
  annotate('text', x = 2.75, y = 2, label = 'Null') +
  annotate('text', x = -2.75, y = 2, label = 'Cre') +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's d") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle('Amino acid metabolism') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA)
    )

p4 <- ggplot(compass_results %>% filter(subsystem =='Fatty acid oxidation'), aes(x = cohens_d, y = -log10(adjusted_pval))) +
  geom_point(shape = 21, fill = '#20854EFF', size = 3) +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  geom_hline(yintercept = 1.3, linetype = 'dashed') +
  annotate('text', x = 2.75, y = 2, label = 'Null') +
  annotate('text', x = -2.75, y = 2, label = 'Cre') +
  xlim(c(-3, 3)) +
  ylim(c(0, 30)) +
  xlab("Cohen's d") +
  ylab(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p))) +
  ggtitle('Fatty acid oxidation') +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = 'bold'),
    panel.border = element_rect(color = 'black', fill = NA)
    )

figure <- ggpubr::ggarrange(
  p1 + ggpubr::rremove('ylab') + ggpubr::rremove('xlab'), 
  p2 + ggpubr::rremove('ylab') + ggpubr::rremove('xlab'), 
  p3 + ggpubr::rremove('ylab') + ggpubr::rremove('xlab'), 
  p4 + ggpubr::rremove('ylab') + ggpubr::rremove('xlab'),
  labels = NULL,
  ncol = 2, nrow = 2,
  common.legend = TRUE, legend = 'bottom',
  align = 'hv', 
  font.label = list(size = 10, color = 'black', face = 'bold', family = NULL, position = 'top')
  )

ggpubr::annotate_figure(
  figure,
  left = grid::textGrob(bquote(-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), rot = 90, vjust = 0.5, gp = grid::gpar(cex = 1.3)),
  bottom = grid::textGrob("Cohen's d", gp = grid::gpar(cex = 1.3)))

# Calculating metabolic activity
null_micropool <- read.table(here('compass', 'd14_null', 'micropooled_data.tsv'), sep = '\t', header = TRUE, row.names = 1)
cre_micropool <- read.table(here('compass', 'd14_cre', 'micropooled_data.tsv'), sep = '\t', header = TRUE, row.names = 1)

null_micropool <- t(null_micropool)
cre_micropool <- t(cre_micropool)

null_micropool <- as.data.frame(null_micropool)
cre_micropool <- as.data.frame(cre_micropool)

null_micropool <- null_micropool %>%
  rowwise() %>%
  mutate(protein_sum = sum(c_across(everything())))

cre_micropool <- cre_micropool %>%
  rowwise() %>%
  mutate(protein_sum = sum(c_across(everything())))

metabolic_genes <- read.csv('/Users/christopherhorn/Library/Python/3.8/lib/python/site-packages/compass/Resources/Recon2_export/gene_md.csv')

metabolic_genes <- metabolic_genes %>%
  filter(MGI.symbol != 'N/A')

metabolic_genes <- str_to_title(metabolic_genes$MGI.symbol)

null_micropool <- null_micropool %>%
  rowwise() %>%
  mutate(metabolic_sum = sum(c_across(any_of((metabolic_genes)))))

cre_micropool <- cre_micropool %>%
  rowwise() %>%
  mutate(metabolic_sum = sum(c_across(any_of((metabolic_genes)))))

null_micropool <- null_micropool %>%
  rowwise() %>%
  mutate(metabolic_activity = metabolic_sum / protein_sum)

cre_micropool <- cre_micropool %>%
  rowwise() %>%
  mutate(metabolic_activity = metabolic_sum / protein_sum)

null_micropool <- null_micropool %>%
  select(c('protein_sum', 'metabolic_sum', 'metabolic_activity'))

cre_micropool <- cre_micropool %>%
  select(c('protein_sum', 'metabolic_sum', 'metabolic_activity'))

write.csv(null_micropool, here('compass', 'd14 null vs d14 cre', 'd14 null metabolic activity.csv'))
write.csv(cre_micropool,here('compass', 'd14 null vs d14 cre', 'd14 cre metabolic activity.csv'))

rm(
  null_micropool,
  cre_micropool,
  metabolic_genes
)

# PCA of compass meta rxn consistencies
meta_rxn_consist <- read.csv(here('compass', 'd14 null vs d14 cre', 'd14 null vs cre compass metareaction consistencies.csv'), header = TRUE)
meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'X')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))
meta_rxn_consist.pca <- prcomp(meta_rxn_consist, scale. = TRUE)
meta_rxn_consist.pca_coord <- meta_rxn_consist.pca$x[ , 1:3]
meta_rxn_consist.pca_load <- meta_rxn_consist.pca$rotation[ , 1:3]

write.csv(meta_rxn_consist.pca_coord, here('compass', 'd14 null vs d14 cre', 'meta_rxn_pca_coord.csv'))
write.csv(meta_rxn_consist.pca_load, here('compass', 'd14 null vs d14 cre', 'meta_rxn_pca_loadings.csv'))

met_activity <- read.csv(here('compass', 'd14 null vs d14 cre', 'metabolic activity calculations.csv'), header = TRUE) # PCA coords have already been added as columns
met_activity <- select(met_activity, -1)

# PCA plots
p1 <- ggplot(met_activity, aes(x = PC1, y = PC2)) +
  geom_point(shape = 21, color = 'black', size = 3, aes(fill = cell_genotype)) +
  theme_classic() +
  labs(x = 'PC1 (62.66%)', y = 'PC2 (6.44%)') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(fill = guide_legend(title = 'Cell genotype'))

p2 <- ggplot(met_activity, aes(x = PC1, y = PC3)) +
  geom_point(shape = 21, color = 'black', size = 3, aes(fill = cell_genotype)) +
  theme_classic() +
  labs(x = 'PC1 (62.66%)', y = 'PC3 (4.52%)') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(fill = guide_legend(title = 'Cell genotype'))

p3 <- ggplot(met_activity, aes(x = PC2, y = PC3)) +
  geom_point(shape = 21, color = 'black', size = 3, aes(fill = cell_genotype)) +
  theme_classic() +
  labs(x = 'PC2 (6.44%)', y = 'PC3 (4.52%)') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(fill = guide_legend(title = 'Cell genotype'))

meta_rxn_consist.pca_plot <- ggpubr::ggarrange(p1, p2, p3, common.legend = TRUE, legend = 'bottom', nrow = 1, ncol = 3)

rm(
  p1,
  p2,
  p3
  )

# Computing maturity & pathogenicity scores based on Seurat's AddModuleScore function
combined.express <- read.csv(here::here('compass', 'd14 null vs d14 cre', 'combined d14 micropooled_data.csv'))
combined.express[892, 1] <- 'Marc2'
combined.express[13629, 1] <- 'March2'
combined.express <- tibble::column_to_rownames(combined.express, var = 'X')

path_list <- c('Il1b', 'Clec4e', 'Junb', 'Ctsd', 'Wfdc17', 'Il1f9', 'Pla2g7', 'Arg2')
mature_list <- c('Lcn2', 'Ngp', 'Camp', 'Ltf', 'Arhgdib', 'Anxa1', 'Plbd1', 'Tkt', 'Aldh2', 'Ly6c2', 'Adpgk', 'Cd177')

object = combined.express
features = list(path_list) # Make sure to change to the gene list of interest
pool = rownames(combined.express)
nbin = 24
ctrl = 100
k = FALSE
name = 'path_score' # Make sure the name matches the selected gene list
seed = 1

cluster.length <- length(x = features)
assay.data <- object
data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
data.avg <- data.avg[order(data.avg)]
data.cut <- ggplot2::cut_number(
  x = data.avg + rnorm(n = length(data.avg)) / 1e+30,
  n = nbin,
  labels = FALSE,
  right = FALSE
  )

names(x = data.cut) <- names(x = data.avg)
ctrl.use <- vector(mode = 'list', length = cluster.length)

for (i in 1:cluster.length) {
  features.use <- features[[i]]
  for (j in 1:length(x = features.use)) {
    ctrl.use[[i]] <- c(
      ctrl.use[[i]],
      names(x = sample(x = data.cut[which(x = data.cut == data.cut[features.use[j]])], size = ctrl, replace = FALSE))
      )
  }
}

ctrl.use <- lapply(X = ctrl.use, FUN = unique)
ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), ncol = ncol(x = object))

for (i in 1:length(ctrl.use)) {
  features.use <- ctrl.use[[i]]
  ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use,])
}

features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, ncol = ncol(x = object))

for (i in 1:cluster.length) {
  features.use <- features[[i]]
  data.use <- assay.data[features.use, , drop = FALSE]
  features.scores[i, ] <- Matrix::colMeans(x = data.use)
}

features.scores.use <- features.scores - ctrl.scores
rownames(x = features.scores.use) <- paste0(name, 1:cluster.length)
features.scores.use <- as.data.frame(x = t(x = features.scores.use))
rownames(x = features.scores.use) <- colnames(x = object)

rm(
  ctrl.scores,
  ctrl.use,
  data.use,
  features,
  features.scores,
  cluster.length,
  ctrl,
  data.avg,
  data.cut,
  features.use,
  i,
  j,
  k,
  name,
  nbin,
  pool,
  seed,
  assay.data,
  object,
  combined.express
  )

# PCs vs metabolic activity plots
met_activity <- read.csv(here('compass', 'd14 null vs d14 cre', 'metabolic activity calculations.csv'), header = TRUE) # PCA coords have already been added as columns

p1 <- ggplot(met_activity, aes(x = PC1, y = metabolic_activity)) +
  geom_point(aes(shape = cell_genotype, color = maturity_score), size = 3) +
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic() +
  labs(x = 'PC1 (62.66%)', y = 'Metabolic activity') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(shape = guide_legend(title = 'Cell genotype')) +
  viridis::scale_color_viridis(option = 'inferno', 'Maturity score') +
  annotate('text', x = 0, y = 0.12, label = paste('Pearson correlation:', round(cor(met_activity$metabolic_activity, met_activity$PC1), 3)))

p2 <- ggplot(met_activity, aes(x = PC2, y = metabolic_activity)) +
  geom_point(aes(shape = cell_genotype, color = maturity_score), size = 3) +
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic() +
  labs(x = 'PC2 (6.44%)', y = 'Metabolic activity') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(shape = guide_legend(title = 'Cell genotype')) +
  viridis::scale_color_viridis(option = 'inferno', 'Maturity score') +
  annotate('text', x = 0, y = 0.12, label = paste('Pearson correlation:', round(cor(met_activity$metabolic_activity, met_activity$PC2), 3)))

p3 <- ggplot(met_activity, aes(x = PC3, y = metabolic_activity)) +
  geom_point(aes(shape = cell_genotype, color = maturity_score), size = 3) +
  geom_smooth(method = 'lm', color = 'black') +
  theme_classic() +
  labs(x = 'PC3 (4.52%)', y = 'Metabolic activity') +
  scale_fill_manual(values = c('tomato', 'slateblue1')) +
  guides(shape = guide_legend(title = 'Cell genotype')) +
  viridis::scale_color_viridis(option = 'inferno', 'Maturity score') +
  annotate('text', x = 0, y = 0.12, label = paste('Pearson correlation:', round(cor(met_activity$metabolic_activity, met_activity$PC3), 3)))

meta_rxn_consist.pca2_plot <- ggpubr::ggarrange(
  p1,
  p2,
  p3,
  common.legend = TRUE,
  legend = 'bottom',
  nrow = 1,
  ncol = 3
  )

rm(
  p1,
  p2,
  p3,
  met_activity,
  meta_rxn_consist.pca2_plot
  )

# Signed -log x spearman plots
meta_rxn_consist <- read.csv(here('compass', 'd14 null vs d14 cre', 'd14 null vs cre compass metareaction consistencies.csv'), header = TRUE)
meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'X')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))

met_activity <- read.csv(here('compass', 'd14 null vs d14 cre', 'metabolic activity calculations.csv'), header = TRUE) # PCA coords have already been added as columns

patho_cor <- as.data.frame(cor(meta_rxn_consist[-1], met_activity$patho_score, method = 'spearman'))
mature_cor <- as.data.frame(cor(meta_rxn_consist[-1], met_activity$maturity_score, method = 'spearman'))

patho_cor <- rownames_to_column(patho_cor, var = 'rowname')
mature_cor <- rownames_to_column(mature_cor, var = 'rowname')

colnames(patho_cor)[2] <- 'patho_spearman'
colnames(mature_cor)[2] <- 'mature_spearman'

combine.cor <- merge(patho_cor, mature_cor, by = 'rowname')

rm(
  patho_cor,
  mature_cor
)

compass.result <- read.csv(here('compass', 'd14 null vs d14 cre', 'd14 null vs cre compass analysis results.csv'), header = TRUE)
colnames(compass.result)[1] <- 'rowname'
compass.result <- select(compass.result, c(1, 4:8))
compass.result <- compass.result %>%
  mutate(neg_log_p = -log10(adjusted_pval))
compass.result <- compass.result %>%
  mutate(signed_neg_log_p = neg_log_p*sign(cohens_d))

compass.cor <- merge(compass.result, combine.cor, by = 'rowname')

rm(
  compass.result,
  met_activity,
  combine.cor
)

write.csv(compass.cor, here('compass', 'd14 null vs d14 cre', 'meta_rxn correlation with patho mature.csv'))

ggplot(compass.cor, aes(x = signed_neg_log_p, y = patho_spearman)) +
  geom_point(fill = 'gray', shape = 21, color = 'black', size = 3) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_classic() +
  labs(x = bquote(Signed~-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), y = 'Spearman correlation with pathogenicity score') +
  Seurat::NoLegend() +
  geom_point(data = subset(compass.cor, subset = subsystem == 'Glycolysis/gluconeogenesis'), fill = 'red', shape = 21, color = 'black', size = 3)

ggplot(compass.cor, aes(x = signed_neg_log_p, y = mature_spearman)) +
  geom_point(fill = 'gray', shape = 21, color = 'black', size = 3) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  geom_vline(xintercept = 0, linetype = 'dashed') +
  theme_classic() +
  labs(x = bquote(Signed~-log[10](BH-adjusted~Wilcoxon~rank~sum~p)), y = 'Spearman correlation with maturity score') +
  Seurat::NoLegend() +
  geom_point(data = subset(compass.cor, subset = subsystem == 'Glycolysis/gluconeogenesis'), fill = 'red', shape = 21, color = 'black', size = 3)

rm(
  compass.cor,
  meta_rxn_consist
)

# Generating marker gene x meta rxn heatmap
meta_rxn_consist <- read.csv(here('compass', 'd14 null vs d14 cre', 'd14 null vs cre compass metareaction consistencies.csv'), header = TRUE)
meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'X')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))
meta_rxn_consist <- rownames_to_column(meta_rxn_consist, var = 'rowname')

expression <- read.csv(here('compass', 'd14 null vs d14 cre', 'combined d14 micropooled_data.csv'), header = TRUE)
expression[892, 1] <- 'Marc2'
expression[13629, 1] <- 'March2'
expression <- column_to_rownames(expression, var = 'X')
expression <- t(expression)
expression <- as.data.frame(expression)
expression <- rownames_to_column(expression, var = 'rowname')
expression <- expression %>%
  select(c(
    'rowname',
    'Il1b',
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
    'Cd177'
  ))

combined.consist_express <- merge(expression, meta_rxn_consist, by = 'rowname')
combined.consist_express <- column_to_rownames(combined.consist_express, var = 'rowname')

meta_rxn_consist <- column_to_rownames(meta_rxn_consist, var = 'rowname')
meta_rxn_consist <- as.data.frame(t(meta_rxn_consist))

cor.table <- data.frame(reaction = rownames(meta_rxn_consist))
for (i in 1:22) {
  cor.result <- cor(combined.consist_express[23:ncol(combined.consist_express)], combined.consist_express[i], method = 'spearman')
  cor.table[ , ncol(cor.table) + 1] <- cor.result
  colnames(cor.table)[i+1] <- colnames(combined.consist_express)[i]
}

rm(
  cor.result,
  i,
  combined.consist_express,
  expression,
  meta_rxn_consist
)

write.csv(cor.table, here('compass', 'd14 null vs d14 cre', 'maturity genes correlation compass results heatmap.csv'))

compass.result <- read.csv(here('compass', 'd14 null vs d14 cre', 'd14 null vs cre compass analysis results.csv'), header = TRUE)
compass.result <- compass.result %>%
  filter(adjusted_pval < 0.1)

cor.table <- column_to_rownames(cor.table, var = 'reaction')
cor.table <- cor.table %>%
  filter(rownames(.) %in% compass.result$X)

heatmap.data <- as.matrix(cor.table)

rm(
  cor.table,
  compass.result
)

ha <- ComplexHeatmap::HeatmapAnnotation(
  'Gene markers' = c(rep('Pathogenic', 9), rep('Maturity', 13)),
  col = list('Gene markers' = c('Pathogenic' = 'tomato', 'Maturity' = 'slateblue1')),
  gp = grid::gpar(col = 'black'),
  show_annotation_name = FALSE
)

col_fun = circlize::colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))

ComplexHeatmap::Heatmap(
  heatmap.data,
  name = 'Spearman correlation',
  col = col_fun,
  column_names_side = 'top',
  cluster_columns = FALSE,
  cluster_column_slices = FALSE,
  cluster_rows = TRUE,
  cluster_row_slices = FALSE,
  show_parent_dend_line = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  show_row_names = FALSE,
  column_names_rot = 45,
  border = TRUE,
  top_annotation = ha
)

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
integrated$cell <- Idents(integrated)
integrated$cell_origin <- paste(Idents(integrated), integrated$sample_origin, sep = "_")
Idents(integrated) <- "cell_origin"
DefaultAssay(integrated) <- 'RNA'
integrated <- NormalizeData(integrated, verbose = TRUE)
integrated <- ScaleData(integrated, verbose = TRUE)

gran <- subset(integrated, idents = c(
  'Granulocytes 1',
  'Granulocytes 2',
  'Granulocytes 3',
  'Granulocytes 4',
  'Granulocytes 5',
  'Granulocytes 6',
  'Granulocytes 7',
  'Granulocytes 8',
  'Granulocytes 9'
))

gran <- subset(integrated, idents = c(
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
  "Granulocytes 8_d14_cre",
  
  "Granulocytes 9_d3_null",
  "Granulocytes 9_d3_cre",
  "Granulocytes 9_d14_null",
  "Granulocytes 9_d14_cre"
))

d3_null <- subset(gran, sample_origin =='d3_null')
d3_cre <- subset(gran, sample_origin =='d3_cre')
d14_null <- subset(gran, sample_origin =='d14_null')
d14_cre <- subset(gran, sample_origin =='d14_cre')

gran.avg <- AverageExpression(gran, return.seurat = TRUE)
d3_null.avg <- AverageExpression(d3_null, return.seurat = TRUE)
d3_cre.avg <- AverageExpression(d3_cre, return.seurat = TRUE)
d14_null.avg <- AverageExpression(d14_null, return.seurat = TRUE)
d14_cre.avg <- AverageExpression(d14_cre, return.seurat = TRUE)

avg.norm_counts <- gran.avg@assays$RNA@data
d3_null.avg.norm_counts <- d3_null.avg@assays$RNA@data
d3_cre.avg.norm_counts <- d3_cre.avg@assays$RNA@data
d14_null.avg.norm_counts <- d14_null.avg@assays$RNA@data
d14_cre.avg.norm_counts <- d14_cre.avg@assays$RNA@data

avg.norm_counts <- as.data.frame(avg.norm_counts)
d3_null.avg.norm_counts <- as.data.frame(d3_null.avg.norm_counts)
d3_cre.avg.norm_counts <- as.data.frame(d3_cre.avg.norm_counts)
d14_null.avg.norm_counts <- as.data.frame(d14_null.avg.norm_counts)
d14_cre.avg.norm_counts <- as.data.frame(d14_cre.avg.norm_counts)

avg.norm_counts <- rownames_to_column(avg.norm_counts, var = 'gene')
d3_null.avg.norm_counts <- rownames_to_column(d3_null.avg.norm_counts, var = 'gene')
d3_cre.avg.norm_counts <- rownames_to_column(d3_cre.avg.norm_counts, var = 'gene')
d14_null.avg.norm_counts <- rownames_to_column(d14_null.avg.norm_counts, var = 'gene')
d14_cre.avg.norm_counts <- rownames_to_column(d14_cre.avg.norm_counts, var = 'gene')

## Grabbing pathway sets
# GO.set <- msigdbr(species = 'Mus musculus', category = 'C5', subcategory = 'BP') # GO:BP
# CP.set <- msigdbr(species = 'Mus musculus', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
# HM.set <- msigdbr(species = 'Mus musculus', category = 'H') # Hallmark

## Grabbing gene lists for desired pathways
# pathways <- CP.set %>%
#   filter(gs_name == 'KEGG_GLYCOLYSIS_GLUCONEOGENESIS') # Set this to be whatever pathway you'd like

# gene_list <- pathways$human_gene_symbol
gene_list <- c(
  'Ilb1',
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
  'Cd177'
)

## Filtering average Seurat object for genes in pathway
heatmap.data <- avg.norm_counts %>%
  filter(gene %in% gene_list)  %>%
  filter_all(all_vars(. >= 0.1))
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]

colnames(heatmap.data) <- c(
  "G5_d3_null",  "G2_d3_null",  "G6_d3_null",  "G7_d3_null",  "G1_d3_null", 
  "G8_d3_null",  "G3_d3_null",  "G4_d3_null",  "G9_d3_null",  "G6_d3_cre",  
  "G1_d3_cre",   "G7_d3_cre",   "G2_d3_cre",   "G5_d3_cre",   "G3_d3_cre",  
  "G8_d3_cre",   "G4_d3_cre",   "G9_d3_cre",   "G1_d14_null", "G2_d14_null",
  "G8_d14_null", "G5_d14_null", "G3_d14_null", "G4_d14_null", "G7_d14_null",
  "G6_d14_null", "G9_d14_null", "G8_d14_cre",  "G2_d14_cre",  "G4_d14_cre", 
  "G3_d14_cre",  "G1_d14_cre",  "G5_d14_cre",  "G6_d14_cre",  "G7_d14_cre", 
  "G9_d14_cre"
)

col_order <- c(
  "G1_d3_null",
  "G2_d3_null",
  "G3_d3_null",
  "G4_d3_null",
  "G5_d3_null",
  "G6_d3_null",
  "G7_d3_null",
  "G8_d3_null",
  "G9_d3_null",
  
  "G1_d3_cre",
  "G2_d3_cre",
  "G3_d3_cre",
  "G4_d3_cre",
  "G5_d3_cre",
  "G6_d3_cre",
  "G7_d3_cre",
  "G8_d3_cre",
  "G9_d3_cre",
  
  "G1_d14_null",
  "G2_d14_null",
  "G3_d14_null",
  "G4_d14_null",
  "G5_d14_null",
  "G6_d14_null",
  "G7_d14_null",
  "G8_d14_null",
  "G9_d14_null",
  
  "G1_d14_cre",
  "G2_d14_cre",
  "G3_d14_cre",
  "G4_d14_cre",
  "G5_d14_cre",
  "G6_d14_cre",
  "G7_d14_cre",
  "G8_d14_cre",
  "G9_d14_cre"
)

heatmap.data <- heatmap.data[, col_order]

colnames(heatmap.data) <- c('G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'G7', 'G8', 'G9')

## Create vectors for column/row annotations
# mdsc_pmn <- c('')
# Check to make sure annotation order matches how the genes appear in the matrix
gene_anno <- c(
  'MDSC genes',
  'MDSC genes',
  'PMN genes',
  'PMN genes',
  'MDSC genes',
  'PMN genes',
  'PMN genes',
  'PMN genes',
  'MDSC genes',
  'MDSC genes',
  'PMN genes',
  'PMN genes',
  'PMN genes',
  'PMN genes',
  'MDSC genes',
  'MDSC genes',
  'PMN genes',
  'PMN genes',
  'MDSC genes',
  'PMN genes',
  'PMN genes'
)

col_anno <- c(
  "D3 Null",
  "D3 Null",
  "D3 Null",
  "D3 Null",
  "D3 Null",
  "D3 Null",
  "D3 Null",
  "D3 Null",
  "D3 Null",
  
  "D3 Cre",
  "D3 Cre",
  "D3 Cre",
  "D3 Cre",
  "D3 Cre",
  "D3 Cre",
  "D3 Cre",
  "D3 Cre",
  "D3 Cre",
  
  "D14 Null",
  "D14 Null",
  "D14 Null",
  "D14 Null",
  "D14 Null",
  "D14 Null",
  "D14 Null",
  "D14 Null",
  "D14 Null",
  
  "D14 Cre",
  "D14 Cre",
  "D14 Cre",
  "D14 Cre",
  "D14 Cre",
  "D14 Cre",
  "D14 Cre",
  "D14 Cre",
  "D14 Cre"
)

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = circlize::colorRamp2(c(0, htmp_range), c("white", "red"))

ComplexHeatmap::Heatmap(heatmap.data,
                        name = 'Average\nExpression',
                        col = col_fun,
                        rect_gp = grid::gpar(col = 'black', lwd = 2),
                        cluster_columns = TRUE,
                        cluster_column_slices = FALSE,
                        column_gap = grid::unit(5, 'mm'),
                        column_split = col_anno,
                        # row_split = gene_anno,
                        cluster_rows = TRUE,
                        cluster_row_slices = FALSE,
                        row_gap = grid::unit(5, 'mm'),
                        show_parent_dend_line = FALSE,
                        heatmap_legend_param = list(border = 'black'),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid::grid.text(sprintf("%.1f", heatmap.data[i, j]), x, y, gp = grid::gpar(fontsize = 10))
                        })

# Within-cluster pathways
pathway_list <- c(
  'REACTOME_RHO_GTPASES_ACTIVATE_NADPH_OXIDASES',
  'GOBP_HYDROGEN_PEROXIDE_METABOLIC_PROCESS',
  'REACTOME_RESPIRATORY_ELECTRON_TRANSPORT_ATP_SYNTHESIS_BY_CHEMIOSMOTIC_COUPLING_AND_HEAT_PRODUCTION_BY_UNCOUPLING_PROTEINS',
  'WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA',
  'REACTOME_DETOXIFICATION_OF_REACTIVE_OXYGEN_SPECIES',
  'KEGG_LEUKOCYTE_TRANSENDOTHELIAL_MIGRATION',
  'REACTOME_RESPIRATORY_ELECTRON_TRANSPORT',
  'KEGG_OXIDATIVE_PHOSPHORYLATION',
  'HALLMARK_OXIDATIVE_PHOSPHORYLATION',
  'REACTOME_RHO_GTPASE_EFFECTORS',
  'WP_OVERVIEW_OF_PROINFLAMMATORY_AND_PROFIBROTIC_MEDIATORS',
  'HALLMARK_INTERFERON_GAMMA_RESPONSE',
  'HALLMARK_TNFA_SIGNALING_VIA_NFKB',
  'REACTOME_INTERLEUKIN_10_SIGNALING',
  'NABA_SECRETED_FACTORS',
  'KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION',
  'HALLMARK_INFLAMMATORY_RESPONSE',
  'WP_PROSTAGLANDIN_SIGNALING',
  'PID_TNF_PATHWAY',
  'KEGG_APOPTOSIS',
  'BIOCARTA_NFKB_PATHWAY',
  'KEGG_RIG_I_LIKE_RECEPTOR_SIGNALING_PATHWAY',
  'WP_NOVEL_INTRACELLULAR_COMPONENTS_OF_RIGILIKE_RECEPTOR_RLR_PATHWAY',
  'BIOCARTA_IL1R_PATHWAY',
  'WP_SIGNAL_TRANSDUCTION_THROUGH_IL1R',
  'WP_PHOTODYNAMIC_THERAPYINDUCED_NFKB_SURVIVAL_SIGNALING',
  'KEGG_NOD_LIKE_RECEPTOR_SIGNALING_PATHWAY',
  'HALLMARK_KRAS_SIGNALING_UP',
  'PID_IL23_PATHWAY',
  'REACTOME_INTERFERON_GAMMA_SIGNALING',
  'REACTOME_INTERLEUKIN_1_FAMILY_SIGNALING',
  'REACTOME_CYTOKINE_SIGNALING_IN_IMMUNE_SYSTEM',
  'NABA_MATRISOME',
  'REACTOME_INTERLEUKIN_1_SIGNALING',
  'WP_INTERACTIONS_BETWEEN_IMMUNE_CELLS_AND_MICRORNAS_IN_TUMOR_MICROENVIRONMENT',
  'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION',
  'REACTOME_CLASS_A_1_RHODOPSIN_LIKE_RECEPTORS',
  'KEGG_HEMATOPOIETIC_CELL_LINEAGE',
  'REACTOME_PEPTIDE_LIGAND_BINDING_RECEPTORS',
  'REACTOME_DDX58_IFIH1_MEDIATED_INDUCTION_OF_INTERFERON_ALPHA_BETA',
  'REACTOME_SIGNALING_BY_INTERLEUKINS',
  'REACTOME_GPCR_LIGAND_BINDING',
  'REACTOME_INTERLEUKIN_4_AND_INTERLEUKIN_13_SIGNALING',
  'PID_AMB2_NEUTROPHILS_PATHWAY',
  'HALLMARK_COMPLEMENT',
  'REACTOME_DEATH_RECEPTOR_SIGNALLING',
  'WP_FOLATE_METABOLISM',
  'HALLMARK_IL6_JAK_STAT3_SIGNALING',
  'HALLMARK_INTERFERON_ALPHA_RESPONSE'
  )

## Filtering pathway data for pathways in pathway_list
heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- reshape2::dcast(heatmap.data, pathway ~ cluster_subset, value.var = 'NES', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'pathway')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- reshape2::dcast(heatmap2.data, pathway ~ cluster_subset, value.var = 'padj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'pathway')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

## Create vectors for column/row annotations for brain; Check to make sure annotation order matches how the genes appear in the matrix
celltype_anno <- c('N/A',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte',
                   'Granulocyte')

which(heatmap.data == max(heatmap.data), arr.ind = TRUE) # Use this to find the max value within the matrix so that you can set your upper bound
which(heatmap.data == min(heatmap.data), arr.ind = TRUE) # Use this to find the min value within the matrix so that you can set your lower bound

col_fun = circlize::colorRamp2(c(-4, 0, 4), c('blue', 'white', 'red')) # Check to make sure that the upper bound of this range isn't cutting off the expression of some genes in the matrix
ht <- ComplexHeatmap::Heatmap(heatmap.data,
                              name = 'NES',
                              column_title = 'Granulocytes 1',
                              col = col_fun,
                              rect_gp = grid::gpar(col = 'black', lwd = 2),
                              cluster_columns = TRUE,
                              cluster_column_slices = FALSE,
                              column_gap = grid::unit(5, 'mm'),
                              column_names_gp = grid::gpar(fontsize = 8),
                              cluster_rows = TRUE,
                              cluster_row_slices = FALSE,
                              row_gap = grid::unit(5, 'mm'),
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
                              })

ComplexHeatmap::draw(ht, padding = unit(c(70, 5, 2, 5), 'mm')) # bottom, left, top, right paddings

# Top 10 pathways by NES
pathways <- pathways %>%
  filter(str_detect(cluster_name, 'Granulocytes [:digit:]')) %>%
  filter(cluster_name != 'Granulocytes 9')

top_pathways_up <- pathways %>%
  arrange(desc(NES)) %>%
  group_by(cluster_num, cluster_subset) %>%
  do(head(., n = 5))

top_pathways_dn <- pathways %>%
  arrange(NES) %>%
  group_by(cluster_num, cluster_subset) %>%
  do(head(., n = 5))

top_pathways <- rbind(top_pathways_up, top_pathways_dn)

clust_num <- 7 # Choose your cluster of interest here

path_clust <- top_pathways %>% filter(cluster_num == clust_num)
clust_name <- unique(path_clust$cluster_name)
pathway_list <- unique(path_clust$pathway)

heatmap.data <- path_clust %>%
  filter(pathway %in% pathway_list)
heatmap.data <- reshape2::dcast(heatmap.data, pathway ~ cluster_subset, value.var = 'NES', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'pathway')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

heatmap2.data <- path_clust %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- reshape2::dcast(heatmap2.data, pathway ~ cluster_subset, value.var = 'padj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'pathway')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

heatmap.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap.data <- reshape2::dcast(heatmap.data, pathway ~ cluster_name + cluster_subset, value.var = 'NES', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'pathway')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]
heatmap.data <- t(heatmap.data)

heatmap2.data <- pathways %>%
  filter(pathway %in% pathway_list)
heatmap2.data <- reshape2::dcast(heatmap2.data, pathway ~ cluster_name + cluster_subset, value.var = 'padj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'pathway')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]
heatmap2.data <- t(heatmap2.data)

row_order <- c(
  'd3_null',
  'd3_cre',
  'd14_null',
  'd14_cre'
)

row_order <- c(
  'Granulocytes 1_d3_null',
  'Granulocytes 1_d3_cre',
  'Granulocytes 1_d14_null',
  'Granulocytes 1_d14_cre',
  
  'Granulocytes 2_d3_null',
  'Granulocytes 2_d3_cre',
  'Granulocytes 2_d14_null',
  'Granulocytes 2_d14_cre',
  
  'Granulocytes 3_d3_null',
  'Granulocytes 3_d3_cre',
  'Granulocytes 3_d14_null',
  'Granulocytes 3_d14_cre',
  
  'Granulocytes 4_d3_null',
  'Granulocytes 4_d3_cre',
  'Granulocytes 4_d14_null',
  'Granulocytes 4_d14_cre',
  
  'Granulocytes 5_d3_null',
  'Granulocytes 5_d3_cre',
  'Granulocytes 5_d14_null',
  'Granulocytes 5_d14_cre',
  
  'Granulocytes 6_d3_null',
  'Granulocytes 6_d3_cre',
  'Granulocytes 6_d14_null',
  'Granulocytes 6_d14_cre',
  
  'Granulocytes 7_d3_null',
  'Granulocytes 7_d3_cre',
  'Granulocytes 7_d14_null',
  'Granulocytes 7_d14_cre',
  
  'Granulocytes 8_d3_null',
  'Granulocytes 8_d3_cre',
  'Granulocytes 8_d14_null',
  'Granulocytes 8_d14_cre'
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
  'G1_D3 WT',
  'G1_D3 HIF-1a cKO',
  'G1_D14 WT',
  'G1_D14 HIF-1a cKO',
  
  'G2_D3 WT',
  'G2_D3 HIF-1a cKO',
  'G2_D14 WT',
  'G2_D14 HIF-1a cKO',
  
  'G3_D3 WT',
  'G3_D3 HIF-1a cKO',
  'G3_D14 WT',
  'G3_D14 HIF-1a cKO',
  
  'G4_D3 WT',
  'G4_D3 HIF-1a cKO',
  'G4_D14 WT',
  'G4_D14 HIF-1a cKO',
  
  'G5_D3 WT',
  'G5_D3 HIF-1a cKO',
  'G5_D14 WT',
  'G5_D14 HIF-1a cKO',
  
  'G6_D3 WT',
  'G6_D3 HIF-1a cKO',
  'G6_D14 WT',
  'G6_D14 HIF-1a cKO',
  
  'G7_D3 WT',
  'G7_D3 HIF-1a cKO',
  'G7_D14 WT',
  'G7_D14 HIF-1a cKO',
  
  'G8_D3 WT',
  'G8_D3 HIF-1a cKO',
  'G8_D14 WT',
  'G8_D14 HIF-1a cKO'
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
  "Reactome: Respiratory Electron Transport ATP Synthesis by Chemiosmoticc Coupling and Heat Production by Uncoupling Proteins",
  "Reactome: RHO GTPase Effectors",                                                                                            
  "Reactome: RHO GTPases Activate NADPH Oxidases",                                                                             
  "Reactome: RHOA GTPase Cycle",                                                                                               
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

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = circlize::colorRamp2(c(-htmp_range, 0, htmp_range), c('blue', 'white', 'red'))

ht <- ComplexHeatmap::Heatmap(heatmap.data,
                              name = 'NES',
                              # column_title = clust_name,
                              column_title = "Top enriched pathways per granulocyte cluster",
                              col = col_fun,
                              rect_gp = grid::gpar(col = 'black', lwd = 2),
                              cluster_columns = TRUE,
                              cluster_column_slices = FALSE,
                              column_gap = grid::unit(5, 'mm'),
                              column_names_gp = grid::gpar(fontsize = 8),
                              cluster_rows = FALSE,
                              cluster_row_slices = FALSE,
                              row_split = row_anno,
                              row_gap = grid::unit(5, 'mm'),
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
                              })

ComplexHeatmap::draw(ht, padding = unit(c(110, 5, 2, 5), 'mm')) # bottom, left, top, right paddings

# Top 10 expressed genes by cluster
# From DE data
top_genes <- de %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  do(head(., n = 10))

gene_list <- unique(top_genes$gene)

heatmap.data <- top_genes %>%
  filter(gene %in% gene_list)
heatmap.data <- reshape2::dcast(heatmap.data, gene ~ cluster, value.var = 'avg_log2FC', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]

heatmap2.data <- top_genes %>%
  filter(gene %in% gene_list)
heatmap2.data$p_val_adj <- as.numeric(heatmap2.data$p_val_adj)
heatmap2.data <- reshape2::dcast(heatmap2.data, gene ~ cluster, value.var = 'p_val_adj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'gene')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

col_fun = circlize::colorRamp2(c(-4, 0, 4), c('white', 'white', 'red')) # Check to make sure that the upper bound of this range isn't cutting off the expression of some genes in the matrix

ht <- ComplexHeatmap::Heatmap(heatmap.data,
                              name = 'Average\nNormalized\nCounts',
                              column_title = 'Top 10 Differentially Expressed Genes by Cluster',
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
                              )

ComplexHeatmap::draw(ht, padding = unit(c(2, 5, 2, 5), 'mm')) # bottom, left, top, right paddings



# Top 10 expressed genes by cluster (throw this into a loop later on for bonus points)
# From DE data
top_genes_up <- de %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  do(head(., n = 10))

top_genes_down <- de %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  do(tail(., n = 10))

top_genes <- rbind(top_genes_up, top_genes_down)

clust <- 10 # Choose your cluster of interest here

gene_clust <- top_genes %>% filter(cluster_num == clust)
clust_name <- unique(gene_clust$cluster_name)
gene_list <- unique(gene_clust$gene)

heatmap.data <- gene_clust %>%
  filter(gene %in% gene_list)
heatmap.data <- reshape2::dcast(heatmap.data, gene ~ cluster, value.var = 'avg_log2FC', fill = 0)
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]

heatmap2.data <- gene_clust %>%
  filter(gene %in% gene_list)
heatmap2.data$p_val_adj <- as.numeric(heatmap2.data$p_val_adj)
heatmap2.data <- reshape2::dcast(heatmap2.data, gene ~ cluster, value.var = 'p_val_adj', fill = 1)
heatmap2.data <- column_to_rownames(heatmap2.data, var = 'gene')
heatmap2.data <- as.matrix(heatmap2.data)
heatmap2.data <- heatmap2.data[rowSums(heatmap2.data) != 0, ]

col_order <- c(
  paste0(clust_name, '_d3_null'),
  paste0(clust_name, '_d3_cre'),
  paste0(clust_name, '_d14_null'),
  paste0(clust_name, '_d14_cre')
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

# Trajectory Analysis -----
# Finding trajectories =====
## Inferring trajectories
object_counts <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@counts), 'sparseMatrix'))
object_expression <- Matrix::t(as(as.matrix(neutro.integrated@assays$RNA@data), 'sparseMatrix'))
neutro.integrated_dyn <- wrap_expression(
  counts = object_counts, 
  expression = object_expression
)

rm(
  object_counts,
  object_expression
)

## Add a dimensionality reduction
neutro.integrated_dimred <- dyndimred::dimred_umap(neutro.integrated_dyn$expression)

## Infer the trajectory
neutro.integrated_model <- infer_trajectory(neutro.integrated_dyn, ti_slingshot(), verbose = TRUE)

## Plot trajectory & pseudotime
neutro.integrated_milestone.umap <- plot_dimred(neutro.integrated_model, label_milestones = TRUE, dimred = neutro.integrated_dimred, hex_cells = FALSE) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'UMAP 1', y = 'UMAP 2') ## Check this prior to rooting
neutro.integrated_milestone.pca <- plot_dimred(neutro.integrated_model, label_milestones = TRUE, dimred = 'pca', hex_cells = FALSE) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'PC 1', y = 'PC 2') ## Check this prior to rooting
neutro.integrated_traj.umap <- plot_dimred(neutro.integrated_model, dimred = neutro.integrated_dimred, grouping = neutro.integrated@active.ident, color_density = 'grouping', hex_cells = FALSE) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'UMAP 1', y = 'UMAP 2')
neutro.integrated_traj.pca <- plot_dimred(neutro.integrated_model, dimred = 'pca', grouping = neutro.integrated@active.ident, color_density = 'grouping', hex_cells = FALSE) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'PC 1', y = 'PC 2')
neutro.integrated_traj2.umap <- plot_dimred(neutro.integrated_model, dimred = neutro.integrated_dimred, grouping = neutro.integrated@meta.data$Sample.ID, color_density = 'grouping', hex_cells = FALSE) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'UMAP 1', y = 'UMAP 2')
neutro.integrated_traj2.pca <- plot_dimred(neutro.integrated_model, dimred = 'pca', grouping = neutro.integrated@meta.data$Sample.ID, color_density = 'grouping', hex_cells = FALSE) + theme_classic() + theme(legend.position = 'none', axis.text = element_blank(), axis.ticks = element_blank()) + labs(x = 'PC 1', y = 'PC 2')
neutro.integrated_pseudo.umap <- plot_dimred(neutro.integrated_model, "pseudotime", pseudotime = calculate_pseudotime(neutro.integrated_model), dimred = neutro.integrated_dimred, hex_cells = FALSE) + theme_classic() + labs(x = 'UMAP 1', y = 'UMAP 2') + theme(axis.text = element_blank(), axis.ticks = element_blank())
neutro.integrated_pseudo.pca <- plot_dimred(neutro.integrated_model, "pseudotime", pseudotime = calculate_pseudotime(neutro.integrated_model), dimred = 'pca', hex_cells = FALSE) + theme_classic() + labs(x = 'PC 1', y = 'PC 2') + theme(axis.text = element_blank(), axis.ticks = element_blank())

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

rm(
  ndim,
  optpoint,
  optpoint1,
  optpoint2,
  x,
  pca,
  proj,
  line,
  expression,
  clusterings,
  wh.cl,
  max_clusters
)

## Find trajectory
lineages <- getLineages(dimred, labels, start.clus = '4')
sling <- getCurves(lineages, shrink = 1L, reweight = TRUE, reassign = TRUE, thresh = 0.001, maxit = 10L, stretch = 2L, smoother = 'smooth.spline', shrink.method = 'cosine') ## Parameters taken from dynverse ti_slingshot code

rm(
  dimred,
  labels,
  lineages
)

# TradeSeq #####
## Find optimal number of knots to use for GAM fitting
counts <- neutro.integrated@assays$RNA@counts
filt <- rowSums(counts > 1) >= 120 # Filtering transcripts (genes w/count of at least two in at least 120 different cells)
counts.filt <- counts[filt, ]

rm(filt)

icMat <- evaluateK(counts = as.matrix(counts.filt), sds = sling, k = 3:15, nGenes = 300, verbose = TRUE, plot = TRUE)

## Fit GAM
sce <- fitGAM(as.matrix(counts.filt), sling, nknots = 10) # SCE method
converge <- table(rowData(sce)$tradeSeq$converged) # Check whether genes converged

# List method
# control <- gam.control()
# control$maxit <- 1000 # Set maximum number of iterations to 1,000
# gamList <- fitGAM(counts = as.matrix(counts.filt), pseudotime = slingPseudotime(sling, na = FALSE), cellWeights = slingCurveWeights(sling), control = control, sce = FALSE)
# pvalLineage <- getSmootherPvalues(gamList)
# statLineage <- getSmootherTestStats(gamList)

## Testing
assoRes <- associationTest(sce, lineages = TRUE) # Testing whether genes are significantly changed along pseudotime (independent lineages)
startRes <- startVsEndTest(sce, lineages = TRUE) # Testing whether genes are significantly changed between the start & end of pseudotime (independent lineages)
endRes <- diffEndTest(sce) # Testing whether genes are significantly changed at end of pseudotime
patternRes <- patternTest(sce) # Testing whether genes have significantly different expression patterns throughout pseudotime
earlyDE.plot <- plotGeneCount(curve = sling, counts = as.matrix(counts.filt), clusters = apply(slingClusterLabels(sling), 1, which.max), models = sce) # Visualize where knots are
earlyDERes <- earlyDETest(sce, knots = c(1, 3))

## Identifying significant changes in assoRes
assoRes.sig.lin1 <- rownames(assoRes)[which(p.adjust(assoRes$pvalue_1, "fdr") <= 0.05)]
assoRes.sig.lin2 <- rownames(assoRes)[which(p.adjust(assoRes$pvalue_2, "fdr") <= 0.05)]

assoRes.sig.upset <- upset(fromList(list('Lineage 1' = assoRes.sig.lin1 , 'Lineage 2' = assoRes.sig.lin2)))

yhatSmooth.1 <- predictSmooth(sce, gene = assoRes.sig.lin1, nPoints = 50)
heatSmooth.1 <- pheatmap(t(scale(t(yhatSmooth.1[, 1:50]))), cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)

yhatSmooth.2 <- predictSmooth(sce, gene = assoRes.sig.lin2, nPoints = 50)
heatSmooth.2 <- pheatmap(t(scale(t(yhatSmooth.2[, 1:50]))), cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE)

## Combining End & Pattern testing
patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes.comp <- patternRes[, c('Gene', 'pattern')]

endRes$Gene <- rownames(endRes)
endRes$end <- endRes$waldStat
endRes.comp <- endRes[, c('Gene', 'end')]

compare <- merge(patternRes.comp, endRes.comp, by = 'Gene', all = FALSE)
compare$transientScore <- rank(-compare$end, ties.method = 'min')^2 + rank(compare$pattern, ties.method = 'random')^2

rm(
  patternRes.comp,
  endRes.comp
)

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

rm(
  geneId,
  ii,
  xx,
  p,
  cId,
  nPointsClus,
  cUniq,
  clusterLabels
)

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

rm(
  geneSets,
  m_list,
  stats.lin1,
  stats.lin2
)

# Helper functions -----
files_to_read <- list.files(path = 'output/brain/GSEA/monomac', pattern = '\\.tsv$', full.names = TRUE)
all_files <- lapply(files_to_read, function(x) {
  read.table(file = x, 
             sep = '\t', 
             header = TRUE)
})

for(i in 1:length(all_files)) {
  write.csv(all_files[i], here('output', 'brain', 'GSEA', 'monomac', 'converted', paste0(list.files('output/brain/GSEA/monomac')[i], '.csv')))
}