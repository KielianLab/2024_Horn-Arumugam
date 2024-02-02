## Made by Christopher M. Horn, MS
## Kielian Lab data
## Created: 2022-11-30
## Updated: 2023-01-10

## Notes: Integrating patient data from three subjects: Subject 1, Subject 2, & Subject 5

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

rm(
  packages,
  new.packages
)

## Set options
options(future.globals.maxSize = 4000 * 1024^2)
set.seed(12345)

# Initializing objects -----

## Setting up Seurat objects
## Load in the data
# Subject 1
s1.blood.data <- Read10X(data.dir = here('raw_data', 'subject_1', 'blood'))
s1.tissue.data <- Read10X(data.dir = here('raw_data', 'subject_1', 'tissue'))

# Subject 2
s2.blood.data <- Read10X(data.dir = here('raw_data', 'subject_2', 'blood'))
s2.tissue.data <- Read10X(data.dir = here('raw_data', 'subject_2', 'tissue'))

# Subject 5
s5.blood.data <- Read10X(data.dir = here('raw_data', 'subject_5', 'blood'))
s5.blood.data <- s5.blood.data$`Gene Expression`
s5.tissue.data <- Read10X(data.dir = here('raw_data', 'subject_5', 'tissue'))
s5.tissue.data <- s5.tissue.data$`Gene Expression`

## Initialize Seurat object w/raw data
# Subject 1
s1.blood <- CreateSeuratObject(counts = s1.blood.data, project = 's1 blood', min.cells = 3, min.features = 200)
s1.tissue <- CreateSeuratObject(counts = s1.tissue.data, project = 's1 tissue', min.cells = 3, min.features = 200)

# Subject 2
s2.blood <- CreateSeuratObject(counts = s2.blood.data, project = 's2 blood', min.cells = 3, min.features = 200)
s2.tissue <- CreateSeuratObject(counts = s2.tissue.data, project = 's2 tissue', min.cells = 3, min.features = 200)

# Subject 5
s5.blood <- CreateSeuratObject(counts = s5.blood.data, project = 's5 blood', min.cells = 3, min.features = 200)
s5.tissue <- CreateSeuratObject(counts = s5.tissue.data, project = 's5 tissue', min.cells = 3, min.features = 200)

## Pre-processing and QC
# Subject 1
s1.blood[['percent.mt']] <- PercentageFeatureSet(s1.blood, pattern = '^MT-')
s1.blood <- subset(s1.blood, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

s1.tissue[['percent.mt']] <- PercentageFeatureSet(s1.tissue, pattern = '^MT-')
s1.tissue <- subset(s1.tissue, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Subject 2
s2.blood[['percent.mt']] <- PercentageFeatureSet(s2.blood, pattern = '^MT-')
s2.blood <- subset(s2.blood, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

s2.tissue[['percent.mt']] <- PercentageFeatureSet(s2.tissue, pattern = '^MT-')
s2.tissue <- subset(s2.tissue, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# Subject 5
s5.blood[['percent.mt']] <- PercentageFeatureSet(s5.blood, pattern = '^MT-')
s5.blood <- subset(s5.blood, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

s5.tissue[['percent.mt']] <- PercentageFeatureSet(s5.tissue, pattern = '^MT-')
s5.tissue <- subset(s5.tissue, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

## Create timing metadata
# Subject 1
sample.id_1 <- rep('s1 blood', length(s1.blood@meta.data$orig.ident))
sample.id_2 <- rep('s1 tissue', length(s1.tissue@meta.data$orig.ident))
names(sample.id_1) <- rownames(s1.blood@meta.data)
names(sample.id_2) <- rownames(s1.tissue@meta.data)
s1.blood <- AddMetaData(s1.blood, sample.id_1, col.name = 'Sample_origin')
s1.tissue <- AddMetaData(s1.tissue, sample.id_2, col.name = 'Sample_origin')

# Subject 1
sample.id_3 <- rep('s2 blood', length(s2.blood@meta.data$orig.ident))
sample.id_4 <- rep('s2 tissue', length(s2.tissue@meta.data$orig.ident))
names(sample.id_3) <- rownames(s2.blood@meta.data)
names(sample.id_4) <- rownames(s2.tissue@meta.data)
s2.blood <- AddMetaData(s2.blood, sample.id_3, col.name = 'Sample_origin')
s2.tissue <- AddMetaData(s2.tissue, sample.id_4, col.name = 'Sample_origin')

# Subject 1
sample.id_5 <- rep('s5 blood', length(s5.blood@meta.data$orig.ident))
sample.id_6 <- rep('s5 tissue', length(s5.tissue@meta.data$orig.ident))
names(sample.id_5) <- rownames(s5.blood@meta.data)
names(sample.id_6) <- rownames(s5.tissue@meta.data)
s5.blood <- AddMetaData(s5.blood, sample.id_5, col.name = 'Sample_origin')
s5.tissue <- AddMetaData(s5.tissue, sample.id_6, col.name = 'Sample_origin')

## Write individual object metadata to file
# Subject 1
write.csv(s1.blood@meta.data, here('output', 's1 blood_metadata.csv'))
write.csv(s1.tissue@meta.data, here('output', 's1 tissue_metadata.csv'))

# Subject 2
write.csv(s2.blood@meta.data, here('output', 's2 blood_metadata.csv'))
write.csv(s2.tissue@meta.data, here('output', 's2 tissue_metadata.csv'))

# Subject 5
write.csv(s5.blood@meta.data, here('output', 's5 blood_metadata.csv'))
write.csv(s5.tissue@meta.data, here('output', 's5 tissue_metadata.csv'))

## Write QC metrics
# Subject 1
s1.blood.QC_mets <- dim(s1.blood)
s1.tissue.QC_mets <- dim(s1.tissue)

write.csv(s1.blood.QC_mets, here('output', 'QC', 'subject_1', 's1 blood QC_mets.csv'))
write.csv(s1.tissue.QC_mets, here('output', 'QC', 'subject_1', 's1 tissue QC_mets.csv'))

s1.blood.QC_mets.plot  <- VlnPlot(s1.blood, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
s1.tissue.QC_mets.plot  <- VlnPlot(s1.tissue, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

ggsave('s1 blood QC_mets.plot.png', plot = s1.blood.QC_mets.plot, device = 'png', path = here('output', 'QC', 'subject_1'))
ggsave('s1 tissue QC_mets.plot.png', plot = s1.tissue.QC_mets.plot, device = 'png', path = here('output', 'QC', 'subject_1'))

# Subject 2
s2.blood.QC_mets <- dim(s2.blood)
s2.tissue.QC_mets <- dim(s2.tissue)

write.csv(s2.blood.QC_mets, here('output', 'QC', 'subject_2', 's2 blood QC_mets.csv'))
write.csv(s2.tissue.QC_mets, here('output', 'QC', 'subject_2', 's2 tissue QC_mets.csv'))

s2.blood.QC_mets.plot  <- VlnPlot(s2.blood, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
s2.tissue.QC_mets.plot  <- VlnPlot(s2.tissue, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

ggsave('s2 blood QC_mets.plot.png', plot = s2.blood.QC_mets.plot, device = 'png', path = here('output', 'QC', 'subject_2'))
ggsave('s2 tissue QC_mets.plot.png', plot = s2.tissue.QC_mets.plot, device = 'png', path = here('output', 'QC', 'subject_2'))

# Subject 5
s5.blood.QC_mets <- dim(s5.blood)
s5.tissue.QC_mets <- dim(s5.tissue)

write.csv(s5.blood.QC_mets, here('output', 'QC', 'subject_5', 's5 blood QC_mets.csv'))
write.csv(s5.tissue.QC_mets, here('output', 'QC', 'subject_5', 's5 tissue QC_mets.csv'))

s5.blood.QC_mets.plot  <- VlnPlot(s5.blood, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)
s5.tissue.QC_mets.plot  <- VlnPlot(s5.tissue, features = c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol = 3)

ggsave('s5 blood QC_mets.plot.png', plot = s5.blood.QC_mets.plot, device = 'png', path = here('output', 'QC', 'subject_5'))
ggsave('s5 tissue QC_mets.plot.png', plot = s5.tissue.QC_mets.plot, device = 'png', path = here('output', 'QC', 'subject_5'))

## Remove temp objects
rm(s1.blood.data,
   s1.tissue.data,
   s2.blood.data,
   s2.tissue.data,
   s5.blood.data,
   s5.tissue.data,
   sample.id_1,
   sample.id_2,
   sample.id_3,
   sample.id_4,
   sample.id_5,
   sample.id_6,
   s1.blood.QC_mets,
   s1.blood.QC_mets.plot,
   s1.tissue.QC_mets,
   s1.tissue.QC_mets.plot,
   s2.blood.QC_mets,
   s2.blood.QC_mets.plot,
   s2.tissue.QC_mets,
   s2.tissue.QC_mets.plot,
   s5.blood.QC_mets,
   s5.blood.QC_mets.plot,
   s5.tissue.QC_mets,
   s5.tissue.QC_mets.plot
)

gc()

# Annotating cell types -----
## Annotate the cells w/SingleR
# Subject 1
s1.blood <- NormalizeData(s1.blood)
s1.tissue <- NormalizeData(s1.tissue)

# Subject 2
s2.blood <- NormalizeData(s2.blood)
s2.tissue <- NormalizeData(s2.tissue)

# Subject 5
s5.blood <- NormalizeData(s5.blood)
s5.tissue <- NormalizeData(s5.tissue)

ref.se <- HumanPrimaryCellAtlasData() # snapshot date: 2020-10-27

# Subject 1
s1.blood_sce <- as.SingleCellExperiment(s1.blood)
s1.tissue_sce <- as.SingleCellExperiment(s1.tissue)

commonGenes.1 <- intersect(rownames(s1.blood_sce), rownames(ref.se))
commonGenes.2 <- intersect(rownames(s1.tissue_sce), rownames(ref.se))

ref.se_1 <- ref.se[commonGenes.1,]
ref.se_2 <- ref.se[commonGenes.2,]

s1.blood_sce <- s1.blood_sce[commonGenes.1,]
s1.tissue_sce <- s1.tissue_sce[commonGenes.2,]

pred.s1.blood <- SingleR(test = s1.blood_sce, ref = ref.se_1, labels = ref.se_1$label.main)
pred.s1.tissue <- SingleR(test = s1.tissue_sce, ref = ref.se_2, labels = ref.se_2$label.main)

s1.blood[['celltype']] <- pred.s1.blood$pruned.labels
s1.tissue[['celltype']] <- pred.s1.tissue$pruned.labels

write.csv(pred.s1.blood, here('output', 'singleR', 's1 blood singleR scores.csv'))
write.csv(pred.s1.tissue, here('output', 'singleR', 's1 tissue singleR scores.csv'))

# Subject 2
s2.blood_sce <- as.SingleCellExperiment(s2.blood)
s2.tissue_sce <- as.SingleCellExperiment(s2.tissue)

commonGenes.3 <- intersect(rownames(s2.blood_sce), rownames(ref.se))
commonGenes.4 <- intersect(rownames(s2.tissue_sce), rownames(ref.se))

ref.se_3 <- ref.se[commonGenes.3,]
ref.se_4 <- ref.se[commonGenes.4,]

s2.blood_sce <- s2.blood_sce[commonGenes.3,]
s2.tissue_sce <- s2.tissue_sce[commonGenes.4,]

pred.s2.blood <- SingleR(test = s2.blood_sce, ref = ref.se_3, labels = ref.se_3$label.main)
pred.s2.tissue <- SingleR(test = s2.tissue_sce, ref = ref.se_4, labels = ref.se_4$label.main)

s2.blood[['celltype']] <- pred.s2.blood$pruned.labels
s2.tissue[['celltype']] <- pred.s2.tissue$pruned.labels

write.csv(pred.s2.blood, here('output', 'singleR', 's2 blood singleR scores.csv'))
write.csv(pred.s2.tissue, here('output', 'singleR', 's2 tissue singleR scores.csv'))

# Subject 5
s5.blood_sce <- as.SingleCellExperiment(s5.blood)
s5.tissue_sce <- as.SingleCellExperiment(s5.tissue)

commonGenes.5 <- intersect(rownames(s5.blood_sce), rownames(ref.se))
commonGenes.6 <- intersect(rownames(s5.tissue_sce), rownames(ref.se))

ref.se_5 <- ref.se[commonGenes.5,]
ref.se_6 <- ref.se[commonGenes.6,]

s5.blood_sce <- s5.blood_sce[commonGenes.5,]
s5.tissue_sce <- s5.tissue_sce[commonGenes.6,]

pred.s5.blood <- SingleR(test = s5.blood_sce, ref = ref.se_5, labels = ref.se_5$label.main)
pred.s5.tissue <- SingleR(test = s5.tissue_sce, ref = ref.se_6, labels = ref.se_6$label.main)

s5.blood[['celltype']] <- pred.s5.blood$pruned.labels
s5.tissue[['celltype']] <- pred.s5.tissue$pruned.labels

write.csv(pred.s5.blood, here('output', 'singleR', 's5 blood singleR scores.csv'))
write.csv(pred.s5.tissue, here('output', 'singleR', 's5 tissue singleR scores.csv'))

## Remove temp objects
rm(ref.se,
   s1.blood_sce,
   s1.tissue_sce,
   s2.blood_sce,
   s2.tissue_sce,
   s5.blood_sce,
   s5.tissue_sce,
   commonGenes.1,
   commonGenes.2,
   commonGenes.3,
   commonGenes.4,
   commonGenes.5,
   commonGenes.6,
   ref.se_1,
   ref.se_2,
   ref.se_3,
   ref.se_4,
   ref.se_5,
   ref.se_6,
   pred.s1.blood,
   pred.s1.tissue,
   pred.s2.blood,
   pred.s2.tissue,
   pred.s5.blood,
   pred.s5.tissue
)

gc()

# Object integration -----
## Integrating all cells from all days
sample.list <- c(
  s1.blood,
  s1.tissue,
  s2.blood,
  s2.tissue,
  s5.blood,
  s5.tissue
)

names(sample.list) <- c(
  's1 blood',
  's1 tissue',
  's2 blood',
  's2 tissue',
  's5 blood',
  's5 tissue'
)

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
   s1.blood,
   s1.tissue,
   s2.blood,
   s2.tissue,
   s5.blood,
   s5.tissue,
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
new.cluster.ids <- c(
  'Granulocytes 1',
  'Granulocytes/Myelocytes 1',
  'Granulocytes 2',
  'Granulocytes 3',
  'Granulocytes 4',
  'Granulocytes 5',
  'Granulocytes 6',
  'T Cells',
  'Granulocytes 7',
  'Granulocytes/Myelocytes 2',
  'Granulocytes 8',
  'Granulocytes 9',
  'Granulocytes 10',
  'C13',
  'Granulocytes 11',
  'NK/T Cells',
  'Granulocytes 12'
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
    '#0073C2FF',
    '#EFC000FF',
    '#868686FF',
    '#CD534CFF',
    '#7AA6DCFF',
    '#003C67FF'
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
    '#0073C2FF',
    '#EFC000FF',
    '#868686FF',
    '#CD534CFF',
    '#7AA6DCFF',
    '#003C67FF'
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
c00 <- subset(ortho.integrated, idents = c('Granulocytes 1'))
c01 <- subset(ortho.integrated, idents = c('Granulocytes/Myelocytes 1'))
c02 <- subset(ortho.integrated, idents = c('Granulocytes 2'))
c03 <- subset(ortho.integrated, idents = c('Granulocytes 3'))
c04 <- subset(ortho.integrated, idents = c('Granulocytes 4'))
c05 <- subset(ortho.integrated, idents = c('Granulocytes 5'))
c06 <- subset(ortho.integrated, idents = c('Granulocytes 6'))
c07 <- subset(ortho.integrated, idents = c('T Cells'))
c08 <- subset(ortho.integrated, idents = c('Granulocytes 7'))
c09 <- subset(ortho.integrated, idents = c('Granulocytes/Myelocytes 2'))
c10 <- subset(ortho.integrated, idents = c('Granulocytes 8'))
c11 <- subset(ortho.integrated, idents = c('Granulocytes 9'))
c12 <- subset(ortho.integrated, idents = c('Granulocytes 10'))
c13 <- subset(ortho.integrated, idents = c('C13'))
c14 <- subset(ortho.integrated, idents = c('Granulocytes 11'))
c15 <- subset(ortho.integrated, idents = c('NK/T Cells'))
c16 <- subset(ortho.integrated, idents = c('Granulocytes 12'))

c.list <- list(
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

rm(
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

## Within cluster blood vs tissue DE Analysis
ortho.integrated$origin <- str_split(ortho.integrated$Sample_origin, ' ', simplify = TRUE)[,2]
ortho.integrated$cell <- Idents(ortho.integrated)

c00 <- subset(ortho.integrated, idents = c('Granulocytes 1'))
c01 <- subset(ortho.integrated, idents = c('Granulocytes/Myelocytes 1'))
c02 <- subset(ortho.integrated, idents = c('Granulocytes 2'))
c03 <- subset(ortho.integrated, idents = c('Granulocytes 3'))
c04 <- subset(ortho.integrated, idents = c('Granulocytes 4'))
c05 <- subset(ortho.integrated, idents = c('Granulocytes 5'))
c06 <- subset(ortho.integrated, idents = c('Granulocytes 6'))
c07 <- subset(ortho.integrated, idents = c('T Cells'))
c08 <- subset(ortho.integrated, idents = c('Granulocytes 7'))
c09 <- subset(ortho.integrated, idents = c('Granulocytes/Myelocytes 2'))
c10 <- subset(ortho.integrated, idents = c('Granulocytes 8'))
c11 <- subset(ortho.integrated, idents = c('Granulocytes 9'))
c12 <- subset(ortho.integrated, idents = c('Granulocytes 10'))
c13 <- subset(ortho.integrated, idents = c('C13'))
c14 <- subset(ortho.integrated, idents = c('Granulocytes 11'))
c15 <- subset(ortho.integrated, idents = c('NK/T Cells'))
c16 <- subset(ortho.integrated, idents = c('Granulocytes 12'))

c.list <- list(
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

rm(
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

gc()

for(i in 1:length(c.list)){
  DefaultAssay(c.list[[i]]) <- 'RNA'
  c.list[[i]] <- NormalizeData(c.list[[i]], verbose = TRUE)
  c.list[[i]] <- ScaleData(c.list[[i]], verbose = TRUE)
  c.list[[i]]$cell_origin <- paste(Idents(c.list[[i]]), c.list[[i]]$origin, sep = '_')
  Idents(c.list[[i]]) <- 'cell_origin'
}

for(i in 1:length(c.list)){
  DE <- FindAllMarkers(c.list[[i]], min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
  DE <- DE %>%
    arrange(desc(avg_log2FC))
  write.csv(DE, here('output', 'DE', 'within cluster', 'blood_v_tissue', paste0('DE_', i-1, '.csv')))
}

rm(
  i,
  c.list,
  DE
)

gc()

# Cluster-level GSEA -----
# Between cluster
GO.set <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Homo sapiens', category = 'C2') %>% filter(gs_subcat != 'CGP') # Canonical Pathways
HM.set <- msigdbr(species = 'Homo sapiens', category = 'H') # Hallmark

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
    eaRes <- fgsea(
      pathways = m_list,
      stats = stats,
      eps = 0.0,
      minSize = 10,
      maxSize = 500
    )
    
    eaRes <- arrange(eaRes, desc(NES))
    
    fwrite(
      eaRes,
      here(
        'output',
        'GSEA',
        unique(geneSet_list[[i]]$gs_cat),
        paste0(
          unique(geneSet_list[[i]]$gs_cat),
          '_eaRes_',
          ii-1,
          '.tsv'
        )
      ),
      sep = '\t',
      sep2 = c('', ' ', '')
    )
    
    temp <- read.table(
      here(
        'output',
        'GSEA',
        unique(geneSet_list[[i]]$gs_cat),
        paste0(
          unique(geneSet_list[[i]]$gs_cat),
          '_eaRes_',
          ii-1,
          '.tsv'
        )
      ),
      sep = '\t',
      header = TRUE
    )
    
    temp <- temp %>%
      mutate(
        cluster_num = rep(ii-1, length(nrow(temp))),
        cluster_name = rep(levels(ortho.integrated)[ii], length(nrow(temp))),
        pathway_db = rep(unique(geneSet_list[[i]]$gs_cat), length(nrow(temp)))
    )
    
    write.csv(
      temp,
      here(
        'output',
        'GSEA',
        'csv',
        unique(geneSet_list[[i]]$gs_cat),
        paste0(
          unique(geneSet_list[[i]]$gs_cat),
          '_eaRes_',
          ii-1,
          '.csv'
        )
      )
    )
  }
}

rm(
  geneSet_list,
  m_list,
  stats,
  glist,
  eaRes,
  DE,
  GO.set,
  HM.set,
  CP.set,
  i,
  ii,
  temp
)

gc()

c2_00 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_0.csv'))
c2_01 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_1.csv'))
c2_02 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_2.csv'))
c2_03 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_3.csv'))
c2_04 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_4.csv'))
c2_05 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_5.csv'))
c2_06 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_6.csv'))
c2_07 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_7.csv'))
c2_08 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_8.csv'))
c2_09 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_9.csv'))
c2_10 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_10.csv'))
c2_11 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_11.csv'))
c2_12 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_12.csv'))
c2_13 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_13.csv'))
c2_14 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_14.csv'))
c2_15 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_15.csv'))
c2_16 <- read.csv(here('output', 'GSEA', 'csv', 'C2', 'C2_eaRes_16.csv'))

h_00 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_0.csv'))
h_01 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_1.csv'))
h_02 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_2.csv'))
h_03 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_3.csv'))
h_04 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_4.csv'))
h_05 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_5.csv'))
h_06 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_6.csv'))
h_07 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_7.csv'))
h_08 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_8.csv'))
h_09 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_9.csv'))
h_10 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_10.csv'))
h_11 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_11.csv'))
h_12 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_12.csv'))
h_13 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_13.csv'))
h_14 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_14.csv'))
h_15 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_15.csv'))
h_16 <- read.csv(here('output', 'GSEA', 'csv', 'H', 'H_eaRes_16.csv'))

c5_00 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_0.csv'))
c5_01 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_1.csv'))
c5_02 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_2.csv'))
c5_03 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_3.csv'))
c5_04 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_4.csv'))
c5_05 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_5.csv'))
c5_06 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_6.csv'))
c5_07 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_7.csv'))
c5_08 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_8.csv'))
c5_09 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_9.csv'))
c5_10 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_10.csv'))
c5_11 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_11.csv'))
c5_12 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_12.csv'))
c5_13 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_13.csv'))
c5_14 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_14.csv'))
c5_15 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_15.csv'))
c5_16 <- read.csv(here('output', 'GSEA', 'csv', 'C5', 'C5_eaRes_16.csv'))

c2 <- rbind(c2_00, c2_01, c2_02, c2_03, c2_04, c2_05, c2_06, c2_07, c2_08, c2_09, c2_10, c2_11, c2_12, c2_13, c2_14, c2_15, c2_16)
h <- rbind(h_00, h_01, h_02, h_03, h_04, h_05, h_06, h_07, h_08, h_09, h_10, h_11, h_12, h_13, h_14, h_15, h_16)
c5 <- rbind(c5_00, c5_01, c5_02, c5_03, c5_04, c5_05, c5_06, c5_07, c5_08, c5_09, c5_10, c5_11, c5_12, c5_13, c5_14, c5_15, c5_16)

write.csv(c2, here('output', 'GSEA', 'full_C2.csv'))
write.csv(h, here('output', 'GSEA', 'full_H.csv'))
write.csv(c5, here('output', 'GSEA', 'full_C5.csv'))

rm(c2_00, c2_01, c2_02, c2_03, c2_04, c2_05, c2_06, c2_07, c2_08, c2_09, c2_10, c2_11, c2_12, c2_13, c2_14, c2_15, c2_16)
rm(h_00, h_01, h_02, h_03, h_04, h_05, h_06, h_07, h_08, h_09, h_10, h_11, h_12, h_13, h_14, h_15, h_16)
rm(c5_00, c5_01, c5_02, c5_03, c5_04, c5_05, c5_06, c5_07, c5_08, c5_09, c5_10, c5_11, c5_12, c5_13, c5_14, c5_15, c5_16)

# Within cluster
GO.set <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Homo sapiens', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
HM.set <- msigdbr(species = 'Homo sapiens', category = 'H') # Hallmark

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
      eaRes <- fgsea(
        pathways = m_list,
        stats = stats,
        eps = 0.0,
        minSize = 10,
        maxSize = 500
      )
      
      eaRes <- arrange(eaRes, desc(NES))
      
      fwrite(
        eaRes,
        here(
          'output',
          'GSEA',
          'within cluster',
          unique(geneSet_list[[i]]$gs_cat),
          paste0(
            unique(geneSet_list[[i]]$gs_cat),
            '_eaRes_',
            ii-1,
            '_',
            str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'),
            '.tsv'
          )
        ),
        sep = '\t',
        sep2 = c('', ' ', '')
      )
      
      temp <- read.table(
        here(
          'output',
          'GSEA',
          'within cluster',
          unique(geneSet_list[[i]]$gs_cat),
          paste0(
            unique(geneSet_list[[i]]$gs_cat),
            '_eaRes_',
            ii-1,
            '_',
            str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'),
            '.tsv')
          ),
        sep = '\t',
        header = TRUE
      )
      
      temp <- temp %>%
        mutate(
          cluster_num = rep(ii-1, length(nrow(temp))),
          cluster_name = rep(levels(ortho.integrated)[ii], length(nrow(temp))),
          pathway_db = rep(unique(geneSet_list[[i]]$gs_cat), length(nrow(temp))),
          sample_subset = rep(str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'), length(nrow(temp)))
      )
      
      write.csv(
        temp,
        here(
          'output',
          'GSEA',
          'within cluster',
          'csv',
          unique(geneSet_list[[i]]$gs_cat),
          paste0(
            unique(geneSet_list[[i]]$gs_cat),
            '_eaRes_',
            ii-1,
            '_',
            str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'),
            '.csv'
          )
        )
      )
    }
  }
}

rm(
  geneSet_list,
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
  iii,
  temp
)

gc()

c2_00_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_s1 blood.csv'))
c2_00_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_s2 blood.csv'))
c2_00_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_s5 blood.csv'))
c2_00_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_s1 tissue.csv'))
c2_00_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_s2 tissue.csv'))
c2_00_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_0_s5 tissue.csv'))

c2_01_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_s1 blood.csv'))
c2_01_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_s2 blood.csv'))
c2_01_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_s5 blood.csv'))
c2_01_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_s1 tissue.csv'))
c2_01_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_s2 tissue.csv'))
c2_01_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_1_s5 tissue.csv'))

c2_02_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_s1 blood.csv'))
c2_02_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_s2 blood.csv'))
c2_02_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_s5 blood.csv'))
c2_02_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_s1 tissue.csv'))
c2_02_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_s2 tissue.csv'))
c2_02_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_2_s5 tissue.csv'))

c2_03_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_s1 blood.csv'))
c2_03_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_s2 blood.csv'))
c2_03_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_s5 blood.csv'))
c2_03_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_s1 tissue.csv'))
c2_03_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_s2 tissue.csv'))
c2_03_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_3_s5 tissue.csv'))

c2_04_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_s1 blood.csv'))
c2_04_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_s2 blood.csv'))
c2_04_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_s5 blood.csv'))
c2_04_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_s1 tissue.csv'))
c2_04_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_s2 tissue.csv'))
c2_04_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_4_s5 tissue.csv'))

c2_05_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_s1 blood.csv'))
c2_05_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_s2 blood.csv'))
c2_05_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_s5 blood.csv'))
c2_05_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_s1 tissue.csv'))
c2_05_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_s2 tissue.csv'))
c2_05_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_5_s5 tissue.csv'))

c2_06_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_s1 blood.csv'))
c2_06_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_s2 blood.csv'))
c2_06_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_s5 blood.csv'))
c2_06_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_s1 tissue.csv'))
c2_06_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_s2 tissue.csv'))
c2_06_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_6_s5 tissue.csv'))

c2_07_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_s1 blood.csv'))
c2_07_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_s2 blood.csv'))
c2_07_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_s5 blood.csv'))
c2_07_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_s1 tissue.csv'))
c2_07_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_s2 tissue.csv'))
c2_07_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_7_s5 tissue.csv'))

c2_08_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_8_s1 blood.csv'))
c2_08_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_8_s2 blood.csv'))
c2_08_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_8_s5 blood.csv'))
c2_08_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_8_s1 tissue.csv'))
c2_08_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_8_s2 tissue.csv'))
c2_08_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_8_s5 tissue.csv'))

c2_09_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_9_s1 blood.csv'))
c2_09_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_9_s2 blood.csv'))
c2_09_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_9_s5 blood.csv'))
c2_09_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_9_s1 tissue.csv'))
c2_09_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_9_s2 tissue.csv'))
c2_09_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_9_s5 tissue.csv'))

c2_10_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_10_s1 blood.csv'))
c2_10_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_10_s2 blood.csv'))
c2_10_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_10_s5 blood.csv'))
c2_10_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_10_s1 tissue.csv'))
c2_10_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_10_s2 tissue.csv'))
c2_10_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_10_s5 tissue.csv'))

c2_11_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_11_s1 blood.csv'))
c2_11_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_11_s2 blood.csv'))
c2_11_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_11_s5 blood.csv'))
c2_11_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_11_s1 tissue.csv'))
c2_11_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_11_s2 tissue.csv'))
c2_11_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_11_s5 tissue.csv'))

c2_12_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_12_s1 blood.csv'))
c2_12_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_12_s2 blood.csv'))
c2_12_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_12_s5 blood.csv'))
c2_12_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_12_s1 tissue.csv'))
c2_12_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_12_s2 tissue.csv'))
c2_12_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_12_s5 tissue.csv'))

c2_13_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_13_s1 blood.csv'))
c2_13_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_13_s2 blood.csv'))
c2_13_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_13_s5 blood.csv'))
c2_13_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_13_s1 tissue.csv'))
c2_13_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_13_s2 tissue.csv'))
c2_13_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_13_s5 tissue.csv'))

c2_14_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_14_s1 blood.csv'))
c2_14_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_14_s2 blood.csv'))
c2_14_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_14_s5 blood.csv'))
c2_14_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_14_s1 tissue.csv'))
c2_14_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_14_s2 tissue.csv'))
c2_14_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_14_s5 tissue.csv'))

c2_15_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_15_s1 blood.csv'))
c2_15_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_15_s2 blood.csv'))
c2_15_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_15_s5 blood.csv'))
c2_15_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_15_s1 tissue.csv'))
c2_15_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_15_s2 tissue.csv'))
c2_15_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_15_s5 tissue.csv'))

c2_16_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_16_s1 blood.csv'))
c2_16_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_16_s2 blood.csv'))
c2_16_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_16_s5 blood.csv'))
c2_16_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_16_s1 tissue.csv'))
c2_16_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_16_s2 tissue.csv'))
c2_16_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C2', 'C2_eaRes_16_s5 tissue.csv'))

h_00_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_s1 blood.csv'))
h_00_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_s2 blood.csv'))
h_00_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_s5 blood.csv'))
h_00_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_s1 tissue.csv'))
h_00_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_s2 tissue.csv'))
h_00_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_0_s5 tissue.csv'))

h_01_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_s1 blood.csv'))
h_01_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_s2 blood.csv'))
h_01_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_s5 blood.csv'))
h_01_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_s1 tissue.csv'))
h_01_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_s2 tissue.csv'))
h_01_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_1_s5 tissue.csv'))

h_02_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_s1 blood.csv'))
h_02_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_s2 blood.csv'))
h_02_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_s5 blood.csv'))
h_02_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_s1 tissue.csv'))
h_02_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_s2 tissue.csv'))
h_02_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_2_s5 tissue.csv'))

h_03_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_s1 blood.csv'))
h_03_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_s2 blood.csv'))
h_03_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_s5 blood.csv'))
h_03_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_s1 tissue.csv'))
h_03_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_s2 tissue.csv'))
h_03_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_3_s5 tissue.csv'))

h_04_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_s1 blood.csv'))
h_04_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_s2 blood.csv'))
h_04_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_s5 blood.csv'))
h_04_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_s1 tissue.csv'))
h_04_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_s2 tissue.csv'))
h_04_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_4_s5 tissue.csv'))

h_05_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_s1 blood.csv'))
h_05_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_s2 blood.csv'))
h_05_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_s5 blood.csv'))
h_05_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_s1 tissue.csv'))
h_05_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_s2 tissue.csv'))
h_05_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_5_s5 tissue.csv'))

h_06_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_s1 blood.csv'))
h_06_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_s2 blood.csv'))
h_06_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_s5 blood.csv'))
h_06_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_s1 tissue.csv'))
h_06_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_s2 tissue.csv'))
h_06_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_6_s5 tissue.csv'))

h_07_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_s1 blood.csv'))
h_07_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_s2 blood.csv'))
h_07_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_s5 blood.csv'))
h_07_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_s1 tissue.csv'))
h_07_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_s2 tissue.csv'))
h_07_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_7_s5 tissue.csv'))

h_08_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_8_s1 blood.csv'))
h_08_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_8_s2 blood.csv'))
h_08_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_8_s5 blood.csv'))
h_08_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_8_s1 tissue.csv'))
h_08_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_8_s2 tissue.csv'))
h_08_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_8_s5 tissue.csv'))

h_09_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_9_s1 blood.csv'))
h_09_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_9_s2 blood.csv'))
h_09_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_9_s5 blood.csv'))
h_09_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_9_s1 tissue.csv'))
h_09_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_9_s2 tissue.csv'))
h_09_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_9_s5 tissue.csv'))

h_10_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_10_s1 blood.csv'))
h_10_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_10_s2 blood.csv'))
h_10_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_10_s5 blood.csv'))
h_10_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_10_s1 tissue.csv'))
h_10_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_10_s2 tissue.csv'))
h_10_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_10_s5 tissue.csv'))

h_11_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_11_s1 blood.csv'))
h_11_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_11_s2 blood.csv'))
h_11_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_11_s5 blood.csv'))
h_11_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_11_s1 tissue.csv'))
h_11_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_11_s2 tissue.csv'))
h_11_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_11_s5 tissue.csv'))

h_12_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_12_s1 blood.csv'))
h_12_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_12_s2 blood.csv'))
h_12_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_12_s5 blood.csv'))
h_12_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_12_s1 tissue.csv'))
h_12_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_12_s2 tissue.csv'))
h_12_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_12_s5 tissue.csv'))

h_13_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_13_s1 blood.csv'))
h_13_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_13_s2 blood.csv'))
h_13_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_13_s5 blood.csv'))
h_13_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_13_s1 tissue.csv'))
h_13_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_13_s2 tissue.csv'))
h_13_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_13_s5 tissue.csv'))

h_14_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_14_s1 blood.csv'))
h_14_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_14_s2 blood.csv'))
h_14_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_14_s5 blood.csv'))
h_14_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_14_s1 tissue.csv'))
h_14_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_14_s2 tissue.csv'))
h_14_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_14_s5 tissue.csv'))

h_15_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_15_s1 blood.csv'))
h_15_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_15_s2 blood.csv'))
h_15_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_15_s5 blood.csv'))
h_15_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_15_s1 tissue.csv'))
h_15_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_15_s2 tissue.csv'))
h_15_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_15_s5 tissue.csv'))

h_16_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_16_s1 blood.csv'))
h_16_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_16_s2 blood.csv'))
h_16_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_16_s5 blood.csv'))
h_16_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_16_s1 tissue.csv'))
h_16_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_16_s2 tissue.csv'))
h_16_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'H', 'H_eaRes_16_s5 tissue.csv'))



c5_00_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_s1 blood.csv'))
c5_00_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_s2 blood.csv'))
c5_00_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_s5 blood.csv'))
c5_00_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_s1 tissue.csv'))
c5_00_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_s2 tissue.csv'))
c5_00_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_0_s5 tissue.csv'))

c5_01_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_s1 blood.csv'))
c5_01_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_s2 blood.csv'))
c5_01_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_s5 blood.csv'))
c5_01_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_s1 tissue.csv'))
c5_01_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_s2 tissue.csv'))
c5_01_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_1_s5 tissue.csv'))

c5_02_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_s1 blood.csv'))
c5_02_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_s2 blood.csv'))
c5_02_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_s5 blood.csv'))
c5_02_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_s1 tissue.csv'))
c5_02_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_s2 tissue.csv'))
c5_02_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_2_s5 tissue.csv'))

c5_03_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_s1 blood.csv'))
c5_03_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_s2 blood.csv'))
c5_03_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_s5 blood.csv'))
c5_03_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_s1 tissue.csv'))
c5_03_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_s2 tissue.csv'))
c5_03_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_3_s5 tissue.csv'))

c5_04_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_s1 blood.csv'))
c5_04_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_s2 blood.csv'))
c5_04_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_s5 blood.csv'))
c5_04_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_s1 tissue.csv'))
c5_04_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_s2 tissue.csv'))
c5_04_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_4_s5 tissue.csv'))

c5_05_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_s1 blood.csv'))
c5_05_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_s2 blood.csv'))
c5_05_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_s5 blood.csv'))
c5_05_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_s1 tissue.csv'))
c5_05_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_s2 tissue.csv'))
c5_05_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_5_s5 tissue.csv'))

c5_06_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_s1 blood.csv'))
c5_06_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_s2 blood.csv'))
c5_06_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_s5 blood.csv'))
c5_06_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_s1 tissue.csv'))
c5_06_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_s2 tissue.csv'))
c5_06_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_6_s5 tissue.csv'))

c5_07_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_s1 blood.csv'))
c5_07_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_s2 blood.csv'))
c5_07_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_s5 blood.csv'))
c5_07_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_s1 tissue.csv'))
c5_07_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_s2 tissue.csv'))
c5_07_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_7_s5 tissue.csv'))

c5_08_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_8_s1 blood.csv'))
c5_08_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_8_s2 blood.csv'))
c5_08_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_8_s5 blood.csv'))
c5_08_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_8_s1 tissue.csv'))
c5_08_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_8_s2 tissue.csv'))
c5_08_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_8_s5 tissue.csv'))

c5_09_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_9_s1 blood.csv'))
c5_09_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_9_s2 blood.csv'))
c5_09_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_9_s5 blood.csv'))
c5_09_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_9_s1 tissue.csv'))
c5_09_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_9_s2 tissue.csv'))
c5_09_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_9_s5 tissue.csv'))

c5_10_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_10_s1 blood.csv'))
c5_10_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_10_s2 blood.csv'))
c5_10_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_10_s5 blood.csv'))
c5_10_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_10_s1 tissue.csv'))
c5_10_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_10_s2 tissue.csv'))
c5_10_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_10_s5 tissue.csv'))

c5_11_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_11_s1 blood.csv'))
c5_11_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_11_s2 blood.csv'))
c5_11_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_11_s5 blood.csv'))
c5_11_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_11_s1 tissue.csv'))
c5_11_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_11_s2 tissue.csv'))
c5_11_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_11_s5 tissue.csv'))

c5_12_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_12_s1 blood.csv'))
c5_12_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_12_s2 blood.csv'))
c5_12_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_12_s5 blood.csv'))
c5_12_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_12_s1 tissue.csv'))
c5_12_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_12_s2 tissue.csv'))
c5_12_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_12_s5 tissue.csv'))

c5_13_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_13_s1 blood.csv'))
c5_13_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_13_s2 blood.csv'))
c5_13_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_13_s5 blood.csv'))
c5_13_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_13_s1 tissue.csv'))
c5_13_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_13_s2 tissue.csv'))
c5_13_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_13_s5 tissue.csv'))

c5_14_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_14_s1 blood.csv'))
c5_14_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_14_s2 blood.csv'))
c5_14_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_14_s5 blood.csv'))
c5_14_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_14_s1 tissue.csv'))
c5_14_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_14_s2 tissue.csv'))
c5_14_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_14_s5 tissue.csv'))

c5_15_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_15_s1 blood.csv'))
c5_15_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_15_s2 blood.csv'))
c5_15_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_15_s5 blood.csv'))
c5_15_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_15_s1 tissue.csv'))
c5_15_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_15_s2 tissue.csv'))
c5_15_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_15_s5 tissue.csv'))

c5_16_s1.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_16_s1 blood.csv'))
c5_16_s2.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_16_s2 blood.csv'))
c5_16_s5.blood <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_16_s5 blood.csv'))
c5_16_s1.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_16_s1 tissue.csv'))
c5_16_s2.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_16_s2 tissue.csv'))
c5_16_s5.tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'csv', 'C5', 'C5_eaRes_16_s5 tissue.csv'))

s1.c2_blood <- rbind(
  c2_00_s1.blood,
  c2_01_s1.blood,
  c2_02_s1.blood,
  c2_03_s1.blood,
  c2_04_s1.blood,
  c2_05_s1.blood,
  c2_06_s1.blood,
  c2_07_s1.blood,
  c2_08_s1.blood,
  c2_09_s1.blood,
  c2_10_s1.blood,
  c2_11_s1.blood,
  c2_12_s1.blood,
  c2_13_s1.blood,
  c2_14_s1.blood,
  c2_15_s1.blood,
  c2_16_s1.blood
)

s1.c2_tissue <- rbind(
  c2_00_s1.tissue,
  c2_01_s1.tissue,
  c2_02_s1.tissue,
  c2_03_s1.tissue,
  c2_04_s1.tissue,
  c2_05_s1.tissue,
  c2_06_s1.tissue,
  c2_07_s1.tissue,
  c2_08_s1.tissue,
  c2_09_s1.tissue,
  c2_10_s1.tissue,
  c2_11_s1.tissue,
  c2_12_s1.tissue,
  c2_13_s1.tissue,
  c2_14_s1.tissue,
  c2_15_s1.tissue,
  c2_16_s1.tissue
)

s2.c2_blood <- rbind(
  c2_00_s2.blood,
  c2_01_s2.blood,
  c2_02_s2.blood,
  c2_03_s2.blood,
  c2_04_s2.blood,
  c2_05_s2.blood,
  c2_06_s2.blood,
  c2_07_s2.blood,
  c2_08_s2.blood,
  c2_09_s2.blood,
  c2_10_s2.blood,
  c2_11_s2.blood,
  c2_12_s2.blood,
  c2_13_s2.blood,
  c2_14_s2.blood,
  c2_15_s2.blood,
  c2_16_s2.blood
)

s2.c2_tissue <- rbind(
  c2_00_s2.tissue,
  c2_01_s2.tissue,
  c2_02_s2.tissue,
  c2_03_s2.tissue,
  c2_04_s2.tissue,
  c2_05_s2.tissue,
  c2_06_s2.tissue,
  c2_07_s2.tissue,
  c2_08_s2.tissue,
  c2_09_s2.tissue,
  c2_10_s2.tissue,
  c2_11_s2.tissue,
  c2_12_s2.tissue,
  c2_13_s2.tissue,
  c2_14_s2.tissue,
  c2_15_s2.tissue,
  c2_16_s2.tissue
)

s5.c2_blood <- rbind(
  c2_00_s5.blood,
  c2_01_s5.blood,
  c2_02_s5.blood,
  c2_03_s5.blood,
  c2_04_s5.blood,
  c2_05_s5.blood,
  c2_06_s5.blood,
  c2_07_s5.blood,
  c2_08_s5.blood,
  c2_09_s5.blood,
  c2_10_s5.blood,
  c2_11_s5.blood,
  c2_12_s5.blood,
  c2_13_s5.blood,
  c2_14_s5.blood,
  c2_15_s5.blood,
  c2_16_s5.blood
)

s5.c2_tissue <- rbind(
  c2_00_s5.tissue,
  c2_01_s5.tissue,
  c2_02_s5.tissue,
  c2_03_s5.tissue,
  c2_04_s5.tissue,
  c2_05_s5.tissue,
  c2_06_s5.tissue,
  c2_07_s5.tissue,
  c2_08_s5.tissue,
  c2_09_s5.tissue,
  c2_10_s5.tissue,
  c2_11_s5.tissue,
  c2_12_s5.tissue,
  c2_13_s5.tissue,
  c2_14_s5.tissue,
  c2_15_s5.tissue,
  c2_16_s5.tissue
)

s1.h_blood <- rbind(
  h_00_s1.blood,
  h_01_s1.blood,
  h_02_s1.blood,
  h_03_s1.blood,
  h_04_s1.blood,
  h_05_s1.blood,
  h_06_s1.blood,
  h_07_s1.blood,
  h_08_s1.blood,
  h_09_s1.blood,
  h_10_s1.blood,
  h_11_s1.blood,
  h_12_s1.blood,
  h_13_s1.blood,
  h_14_s1.blood,
  h_15_s1.blood,
  h_16_s1.blood
)

s1.h_tissue <- rbind(
  h_00_s1.tissue,
  h_01_s1.tissue,
  h_02_s1.tissue,
  h_03_s1.tissue,
  h_04_s1.tissue,
  h_05_s1.tissue,
  h_06_s1.tissue,
  h_07_s1.tissue,
  h_08_s1.tissue,
  h_09_s1.tissue,
  h_10_s1.tissue,
  h_11_s1.tissue,
  h_12_s1.tissue,
  h_13_s1.tissue,
  h_14_s1.tissue,
  h_15_s1.tissue,
  h_16_s1.tissue
)

s2.h_blood <- rbind(
  h_00_s2.blood,
  h_01_s2.blood,
  h_02_s2.blood,
  h_03_s2.blood,
  h_04_s2.blood,
  h_05_s2.blood,
  h_06_s2.blood,
  h_07_s2.blood,
  h_08_s2.blood,
  h_09_s2.blood,
  h_10_s2.blood,
  h_11_s2.blood,
  h_12_s2.blood,
  h_13_s2.blood,
  h_14_s2.blood,
  h_15_s2.blood,
  h_16_s2.blood
)

s2.h_tissue <- rbind(
  h_00_s2.tissue,
  h_01_s2.tissue,
  h_02_s2.tissue,
  h_03_s2.tissue,
  h_04_s2.tissue,
  h_05_s2.tissue,
  h_06_s2.tissue,
  h_07_s2.tissue,
  h_08_s2.tissue,
  h_09_s2.tissue,
  h_10_s2.tissue,
  h_11_s2.tissue,
  h_12_s2.tissue,
  h_13_s2.tissue,
  h_14_s2.tissue,
  h_15_s2.tissue,
  h_16_s2.tissue
)

s5.h_blood <- rbind(
  h_00_s5.blood,
  h_01_s5.blood,
  h_02_s5.blood,
  h_03_s5.blood,
  h_04_s5.blood,
  h_05_s5.blood,
  h_06_s5.blood,
  h_07_s5.blood,
  h_08_s5.blood,
  h_09_s5.blood,
  h_10_s5.blood,
  h_11_s5.blood,
  h_12_s5.blood,
  h_13_s5.blood,
  h_14_s5.blood,
  h_15_s5.blood,
  h_16_s5.blood
)

s5.h_tissue <- rbind(
  h_00_s5.tissue,
  h_01_s5.tissue,
  h_02_s5.tissue,
  h_03_s5.tissue,
  h_04_s5.tissue,
  h_05_s5.tissue,
  h_06_s5.tissue,
  h_07_s5.tissue,
  h_08_s5.tissue,
  h_09_s5.tissue,
  h_10_s5.tissue,
  h_11_s5.tissue,
  h_12_s5.tissue,
  h_13_s5.tissue,
  h_14_s5.tissue,
  h_15_s5.tissue,
  h_16_s5.tissue
)

s1.c5_blood <- rbind(
  c5_00_s1.blood,
  c5_01_s1.blood,
  c5_02_s1.blood,
  c5_03_s1.blood,
  c5_04_s1.blood,
  c5_05_s1.blood,
  c5_06_s1.blood,
  c5_07_s1.blood,
  c5_08_s1.blood,
  c5_09_s1.blood,
  c5_10_s1.blood,
  c5_11_s1.blood,
  c5_12_s1.blood,
  c5_13_s1.blood,
  c5_14_s1.blood,
  c5_15_s1.blood,
  c5_16_s1.blood
)

s1.c5_tissue <- rbind(
  c5_00_s1.tissue,
  c5_01_s1.tissue,
  c5_02_s1.tissue,
  c5_03_s1.tissue,
  c5_04_s1.tissue,
  c5_05_s1.tissue,
  c5_06_s1.tissue,
  c5_07_s1.tissue,
  c5_08_s1.tissue,
  c5_09_s1.tissue,
  c5_10_s1.tissue,
  c5_11_s1.tissue,
  c5_12_s1.tissue,
  c5_13_s1.tissue,
  c5_14_s1.tissue,
  c5_15_s1.tissue,
  c5_16_s1.tissue
)

s2.c5_blood <- rbind(
  c5_00_s2.blood,
  c5_01_s2.blood,
  c5_02_s2.blood,
  c5_03_s2.blood,
  c5_04_s2.blood,
  c5_05_s2.blood,
  c5_06_s2.blood,
  c5_07_s2.blood,
  c5_08_s2.blood,
  c5_09_s2.blood,
  c5_10_s2.blood,
  c5_11_s2.blood,
  c5_12_s2.blood,
  c5_13_s2.blood,
  c5_14_s2.blood,
  c5_15_s2.blood,
  c5_16_s2.blood
)

s2.c5_tissue <- rbind(
  c5_00_s2.tissue,
  c5_01_s2.tissue,
  c5_02_s2.tissue,
  c5_03_s2.tissue,
  c5_04_s2.tissue,
  c5_05_s2.tissue,
  c5_06_s2.tissue,
  c5_07_s2.tissue,
  c5_08_s2.tissue,
  c5_09_s2.tissue,
  c5_10_s2.tissue,
  c5_11_s2.tissue,
  c5_12_s2.tissue,
  c5_13_s2.tissue,
  c5_14_s2.tissue,
  c5_15_s2.tissue,
  c5_16_s2.tissue
)

s5.c5_blood <- rbind(
  c5_00_s5.blood,
  c5_01_s5.blood,
  c5_02_s5.blood,
  c5_03_s5.blood,
  c5_04_s5.blood,
  c5_05_s5.blood,
  c5_06_s5.blood,
  c5_07_s5.blood,
  c5_08_s5.blood,
  c5_09_s5.blood,
  c5_10_s5.blood,
  c5_11_s5.blood,
  c5_12_s5.blood,
  c5_13_s5.blood,
  c5_14_s5.blood,
  c5_15_s5.blood,
  c5_16_s5.blood
)

s5.c5_tissue <- rbind(
  c5_00_s5.tissue,
  c5_01_s5.tissue,
  c5_02_s5.tissue,
  c5_03_s5.tissue,
  c5_04_s5.tissue,
  c5_05_s5.tissue,
  c5_06_s5.tissue,
  c5_07_s5.tissue,
  c5_08_s5.tissue,
  c5_09_s5.tissue,
  c5_10_s5.tissue,
  c5_11_s5.tissue,
  c5_12_s5.tissue,
  c5_13_s5.tissue,
  c5_14_s5.tissue,
  c5_15_s5.tissue,
  c5_16_s5.tissue
)

write.csv(s1.c2_blood, here('output', 'GSEA', 'within cluster', 'full_C2_s1 blood.csv'))
write.csv(s1.c2_tissue, here('output', 'GSEA', 'within cluster', 'full_C2_s1 tissue.csv'))
write.csv(s2.c2_blood, here('output', 'GSEA', 'within cluster', 'full_C2_s2 blood.csv'))
write.csv(s2.c2_tissue, here('output', 'GSEA', 'within cluster', 'full_C2_s2 tissue.csv'))
write.csv(s5.c2_blood, here('output', 'GSEA', 'within cluster', 'full_C2_s5 blood.csv'))
write.csv(s5.c2_tissue, here('output', 'GSEA', 'within cluster', 'full_C2_s5 tissue.csv'))

write.csv(s1.h_blood, here('output', 'GSEA', 'within cluster', 'full_H_s1 blood.csv'))
write.csv(s1.h_tissue, here('output', 'GSEA', 'within cluster', 'full_H_s1 tissue.csv'))
write.csv(s2.h_blood, here('output', 'GSEA', 'within cluster', 'full_H_s2 blood.csv'))
write.csv(s2.h_tissue, here('output', 'GSEA', 'within cluster', 'full_H_s2 tissue.csv'))
write.csv(s5.h_blood, here('output', 'GSEA', 'within cluster', 'full_H_s5 blood.csv'))
write.csv(s5.h_tissue, here('output', 'GSEA', 'within cluster', 'full_H_s5 tissue.csv'))

write.csv(s1.c5_blood, here('output', 'GSEA', 'within cluster', 'full_C5_s1 blood.csv'))
write.csv(s1.c5_tissue, here('output', 'GSEA', 'within cluster', 'full_C5_s1 tissue.csv'))
write.csv(s2.c5_blood, here('output', 'GSEA', 'within cluster', 'full_C5_s2 blood.csv'))
write.csv(s2.c5_tissue, here('output', 'GSEA', 'within cluster', 'full_C5_s2 tissue.csv'))
write.csv(s5.c5_blood, here('output', 'GSEA', 'within cluster', 'full_C5_s5 blood.csv'))
write.csv(s5.c5_tissue, here('output', 'GSEA', 'within cluster', 'full_C5_s5 tissue.csv'))

# Within cluster blood vs tissue
GO.set <- msigdbr(species = 'Homo sapiens', category = 'C5', subcategory = 'BP') # GO:BP
CP.set <- msigdbr(species = 'Homo sapiens', category = 'C2') %>% filter(gs_subcat != 'CGP') # CP
HM.set <- msigdbr(species = 'Homo sapiens', category = 'H') # Hallmark

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
      eaRes <- fgsea(
        pathways = m_list,
        stats = stats,
        eps = 0.0,
        minSize = 10,
        maxSize = 500
      )
      
      eaRes <- arrange(eaRes, desc(NES))
      
      fwrite(
        eaRes,
        here(
          'output',
          'GSEA',
          'within cluster',
          'blood_v_tissue',
          unique(geneSet_list[[i]]$gs_cat),
          paste0(
            unique(geneSet_list[[i]]$gs_cat),
            '_eaRes_',
            ii-1,
            '_',
            str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'),
            '.tsv'
          )
        ),
        sep = '\t',
        sep2 = c('', ' ', '')
      )
      
      temp <- read.table(
        here(
          'output',
          'GSEA',
          'within cluster',
          'blood_v_tissue',
          unique(geneSet_list[[i]]$gs_cat),
          paste0(
            unique(geneSet_list[[i]]$gs_cat),
            '_eaRes_',
            ii-1,
            '_',
            str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'),
            '.tsv')
        ),
        sep = '\t',
        header = TRUE
      )
      
      temp <- temp %>%
        mutate(
          cluster_num = rep(ii-1, length(nrow(temp))),
          cluster_name = rep(levels(ortho.integrated)[ii], length(nrow(temp))),
          pathway_db = rep(unique(geneSet_list[[i]]$gs_cat), length(nrow(temp))),
          sample_subset = rep(str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'), length(nrow(temp)))
        )
      
      write.csv(
        temp,
        here(
          'output',
          'GSEA',
          'within cluster',
          'blood_v_tissue',
          'csv',
          unique(geneSet_list[[i]]$gs_cat),
          paste0(
            unique(geneSet_list[[i]]$gs_cat),
            '_eaRes_',
            ii-1,
            '_',
            str_extract(levels(as.factor(DE_list[[ii]]$cluster))[iii], '[^_]*$'),
            '.csv'
          )
        )
      )
    }
  }
}

rm(
  geneSet_list,
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
  iii,
  temp
)

gc()

c2_00_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_0_blood.csv'))
c2_00_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_0_tissue.csv'))

c2_01_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_1_blood.csv'))
c2_01_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_1_tissue.csv'))

c2_02_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_2_blood.csv'))
c2_02_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_2_tissue.csv'))

c2_03_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_3_blood.csv'))
c2_03_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_3_tissue.csv'))

c2_04_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_4_blood.csv'))
c2_04_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_4_tissue.csv'))

c2_05_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_5_blood.csv'))
c2_05_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_5_tissue.csv'))

c2_06_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_6_blood.csv'))
c2_06_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_6_tissue.csv'))

c2_07_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_7_blood.csv'))
c2_07_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_7_tissue.csv'))

c2_08_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_8_blood.csv'))
c2_08_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_8_tissue.csv'))

c2_09_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_9_blood.csv'))
c2_09_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_9_tissue.csv'))

c2_10_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_10_blood.csv'))
c2_10_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_10_tissue.csv'))

c2_11_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_11_blood.csv'))
c2_11_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_11_tissue.csv'))

c2_12_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_12_blood.csv'))
c2_12_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_12_tissue.csv'))

c2_13_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_13_blood.csv'))
c2_13_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_13_tissue.csv'))

c2_14_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_14_blood.csv'))
c2_14_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_14_tissue.csv'))

c2_15_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_15_blood.csv'))
c2_15_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_15_tissue.csv'))

c2_16_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_16_blood.csv'))
c2_16_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C2', 'C2_eaRes_16_tissue.csv'))



h_00_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_0_blood.csv'))
h_00_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_0_tissue.csv'))

h_01_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_1_blood.csv'))
h_01_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_1_tissue.csv'))

h_02_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_2_blood.csv'))
h_02_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_2_tissue.csv'))

h_03_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_3_blood.csv'))
h_03_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_3_tissue.csv'))

h_04_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_4_blood.csv'))
h_04_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_4_tissue.csv'))

h_05_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_5_blood.csv'))
h_05_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_5_tissue.csv'))

h_06_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_6_blood.csv'))
h_06_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_6_tissue.csv'))

h_07_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_7_blood.csv'))
h_07_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_7_tissue.csv'))

h_08_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_8_blood.csv'))
h_08_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_8_tissue.csv'))

h_09_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_9_blood.csv'))
h_09_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_9_tissue.csv'))

h_10_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_10_blood.csv'))
h_10_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_10_tissue.csv'))

h_11_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_11_blood.csv'))
h_11_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_11_tissue.csv'))

h_12_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_12_blood.csv'))
h_12_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_12_tissue.csv'))

h_13_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_13_blood.csv'))
h_13_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_13_tissue.csv'))

h_14_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_14_blood.csv'))
h_14_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_14_tissue.csv'))

h_15_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_15_blood.csv'))
h_15_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_15_tissue.csv'))

h_16_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_16_blood.csv'))
h_16_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'H', 'H_eaRes_16_tissue.csv'))

c5_00_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_0_blood.csv'))
c5_00_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_0_tissue.csv'))

c5_01_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_1_blood.csv'))
c5_01_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_1_tissue.csv'))

c5_02_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_2_blood.csv'))
c5_02_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_2_tissue.csv'))

c5_03_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_3_blood.csv'))
c5_03_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_3_tissue.csv'))

c5_04_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_4_blood.csv'))
c5_04_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_4_tissue.csv'))

c5_05_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_5_blood.csv'))
c5_05_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_5_tissue.csv'))

c5_06_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_6_blood.csv'))
c5_06_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_6_tissue.csv'))

c5_07_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_7_blood.csv'))
c5_07_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_7_tissue.csv'))

c5_08_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_8_blood.csv'))
c5_08_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_8_tissue.csv'))

c5_09_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_9_blood.csv'))
c5_09_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_9_tissue.csv'))

c5_10_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_10_blood.csv'))
c5_10_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_10_tissue.csv'))

c5_11_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_11_blood.csv'))
c5_11_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_11_tissue.csv'))

c5_12_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_12_blood.csv'))
c5_12_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_12_tissue.csv'))

c5_13_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_13_blood.csv'))
c5_13_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_13_tissue.csv'))

c5_14_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_14_blood.csv'))
c5_14_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_14_tissue.csv'))

c5_15_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_15_blood.csv'))
c5_15_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_15_tissue.csv'))

c5_16_blood <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_16_blood.csv'))
c5_16_tissue <- read.csv(here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'csv', 'C5', 'C5_eaRes_16_tissue.csv'))

c2_blood <- rbind(
  c2_00_blood,
  c2_01_blood,
  c2_02_blood,
  c2_03_blood,
  c2_04_blood,
  c2_05_blood,
  c2_06_blood,
  c2_07_blood,
  c2_08_blood,
  c2_09_blood,
  c2_10_blood,
  c2_11_blood,
  c2_12_blood,
  c2_13_blood,
  c2_14_blood,
  c2_15_blood,
  c2_16_blood
)

c2_tissue <- rbind(
  c2_00_tissue,
  c2_01_tissue,
  c2_02_tissue,
  c2_03_tissue,
  c2_04_tissue,
  c2_05_tissue,
  c2_06_tissue,
  c2_07_tissue,
  c2_08_tissue,
  c2_09_tissue,
  c2_10_tissue,
  c2_11_tissue,
  c2_12_tissue,
  c2_13_tissue,
  c2_14_tissue,
  c2_15_tissue,
  c2_16_tissue
)

h_blood <- rbind(
  h_00_blood,
  h_01_blood,
  h_02_blood,
  h_03_blood,
  h_04_blood,
  h_05_blood,
  h_06_blood,
  h_07_blood,
  h_08_blood,
  h_09_blood,
  h_10_blood,
  h_11_blood,
  h_12_blood,
  h_13_blood,
  h_14_blood,
  h_15_blood,
  h_16_blood
)

h_tissue <- rbind(
  h_00_tissue,
  h_01_tissue,
  h_02_tissue,
  h_03_tissue,
  h_04_tissue,
  h_05_tissue,
  h_06_tissue,
  h_07_tissue,
  h_08_tissue,
  h_09_tissue,
  h_10_tissue,
  h_11_tissue,
  h_12_tissue,
  h_13_tissue,
  h_14_tissue,
  h_15_tissue,
  h_16_tissue
)

c5_blood <- rbind(
  c5_00_blood,
  c5_01_blood,
  c5_02_blood,
  c5_03_blood,
  c5_04_blood,
  c5_05_blood,
  c5_06_blood,
  c5_07_blood,
  c5_08_blood,
  c5_09_blood,
  c5_10_blood,
  c5_11_blood,
  c5_12_blood,
  c5_13_blood,
  c5_14_blood,
  c5_15_blood,
  c5_16_blood
)

c5_tissue <- rbind(
  c5_00_tissue,
  c5_01_tissue,
  c5_02_tissue,
  c5_03_tissue,
  c5_04_tissue,
  c5_05_tissue,
  c5_06_tissue,
  c5_07_tissue,
  c5_08_tissue,
  c5_09_tissue,
  c5_10_tissue,
  c5_11_tissue,
  c5_12_tissue,
  c5_13_tissue,
  c5_14_tissue,
  c5_15_tissue,
  c5_16_tissue
)

write.csv(c2_blood, here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'full_C2_blood.csv'))
write.csv(c2_tissue, here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'full_C2_tissue.csv'))

write.csv(h_blood, here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'full_H_blood.csv'))
write.csv(h_tissue, here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'full_H_tissue.csv'))

write.csv(c5_blood, here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'full_C5_blood.csv'))
write.csv(c5_tissue, here('output', 'GSEA', 'within cluster', 'blood_v_tissue', 'full_C5_tissue.csv'))

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
# Cluster-averaged normalized counts
## Average single cell data & pull out normalized counts
## Counts are normalized by dividing the counts for a given feature by the total counts per cell, multiplying by a scale factor (default == 10,000), and then taking the natural log using log1p()
s1.blood <- subset(ortho.integrated, Sample_origin == 's1 blood')
s1.tissue <- subset(ortho.integrated, Sample_origin == 's1 tissue')
s2.blood <- subset(ortho.integrated, Sample_origin == 's2 blood')
s2.tissue <- subset(ortho.integrated, Sample_origin == 's2 tissue')
s5.blood <- subset(ortho.integrated, Sample_origin == 's5 blood')
s5.tissue <- subset(ortho.integrated, Sample_origin == 's5 tissue')

ortho.avg <- AverageExpression(ortho.integrated, return.seurat = TRUE)
s1.blood.avg <- AverageExpression(s1.blood, return.seurat = TRUE)
s1.tissue.avg <- AverageExpression(s1.tissue, return.seurat = TRUE)
s2.blood.avg <- AverageExpression(s2.blood, return.seurat = TRUE)
s2.tissue.avg <- AverageExpression(s2.tissue, return.seurat = TRUE)
s5.blood.avg <- AverageExpression(s5.blood, return.seurat = TRUE)
s5.tissue.avg <- AverageExpression(s5.tissue, return.seurat = TRUE)

ortho_anc <- ortho.avg@assays$RNA@data
s1.blood_anc <- s1.blood.avg@assays$RNA@data
s1.tissue_anc <- s1.tissue.avg@assays$RNA@data
s2.blood_anc <- s2.blood.avg@assays$RNA@data
s2.tissue_anc <- s2.tissue.avg@assays$RNA@data
s5.blood_anc <- s5.blood.avg@assays$RNA@data
s5.tissue_anc <- s5.tissue.avg@assays$RNA@data

ortho_anc <- as.data.frame(ortho_anc)
s1.blood_anc <- as.data.frame(s1.blood_anc)
s1.tissue_anc <- as.data.frame(s1.tissue_anc)
s2.blood_anc <- as.data.frame(s2.blood_anc)
s2.tissue_anc <- as.data.frame(s2.tissue_anc)
s5.blood_anc <- as.data.frame(s5.blood_anc)
s5.tissue_anc <- as.data.frame(s5.tissue_anc)

ortho_anc <- rownames_to_column(ortho_anc, var = 'gene')
s1.blood_anc <- rownames_to_column(s1.blood_anc, var = 'gene')
s1.tissue_anc <- rownames_to_column(s1.tissue_anc, var = 'gene')
s2.blood_anc <- rownames_to_column(s2.blood_anc, var = 'gene')
s2.tissue_anc <- rownames_to_column(s2.tissue_anc, var = 'gene')
s5.blood_anc <- rownames_to_column(s5.blood_anc, var = 'gene')
s5.tissue_anc <- rownames_to_column(s5.tissue_anc, var = 'gene')

rm(
  s1.blood,
  s1.tissue,
  s2.blood,
  s2.tissue,
  s5.blood,
  s5.tissue,
  ortho.avg,
  s1.blood.avg,
  s1.tissue.avg,
  s2.blood.avg,
  s2.tissue.avg,
  s5.blood.avg,
  s5.tissue.avg
)

gc()

# MDSC v PMN genes based on Alshetaiwi et al (Sci Immunol 2020)
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

# PMN genes based on Lung et al (Cell Metab 2022)
gene_list <- c(
  'S100a8',
  'Cxcr2',
  'Mmp8',
  'Ltf'
)

# G-MDSC Type I genes based on Lung et al (Cell Metab 2022)
g_mdsc.1 <- c(
  'NGP',
  'LTF',
  'CD177',
  'ANXA1',
  'MMP8',
  'S100A8',
  'S100A9',
  'CEBPE',
  'LTB4R1',
  'CYBB'
)

# G-MDSC Type II genes based on Lung et al (Cell Metab 2022)
g_mdsc.2 <- c(
  'CCL4',
  'CCL3',
  'CXCL2',
  'CXCL3',
  'SPP1',
  'IL1B',
  'NFKBIA',
  'SOCS3',
  'MIF',
  'KLF6',
  'ATF3',
  'PTGS2',
  'XBP1'
)

# M-MDSC genes based on Lung et al (Cell Metab 2022)
m_mdsc <- c(
  'ITGAM',
  'LY6G',
  'ARG1',
  'NOS2',
  'LY6C1'
)

gene_list <- c(g_mdsc.1, g_mdsc.2, m_mdsc)

rm(
  g_mdsc.1,
  g_mdsc.2,
  m_mdsc
)

gc()

# PMN-MDSC based on Veglia et al (Nat Rev Immunol 2021)
pmn_mdsc <- c(
  'STAT1',
  'STAT3',
  'STAT6',
  'IRF1',
  'S100A8',
  'S100A9',
  'ANXA1',
  'LYZ2',
  'CXCL1',
  'CXCL2',
  'CXCR1',
  'CXCR2',
  'IL8',
  'LILRA3',
  'TREM1',
  'PTGS2',
  'ARG1',
  'ARG2',
  'TGFB1',
  'VEGF',
  'IL6',
  'CSF1',
  'IL1B',
  'WFDC17',
  'IL4R',
  'OLR1',
  'CD84'
)

# M-MDSC based on Veglia et al (Nat Rev Immunol 2021)
m_mdsc <- c(
  'S100A9',
  'S100A8',
  'ARG1',
  'ARG2',
  'NOS2',
  'IL10',
  'VEGFA',
  'WFDC17',
  'CD14',
  'TGFB1',
  'TNF',
  'STAT3',
  'IL6'
)

gene_list <- m_mdsc

rm(
  pmn_mdsc,
  m_mdsc
)

gc()

## Filtering average Seurat object for genes in pathway
anc <- ortho_anc # Set which set of counts you'd like to build a heatmap for
htmp.title <- 'Integrated Data'

heatmap.data <- anc %>%
  filter(gene %in% toupper(gene_list)) %>%
  filter_all(all_vars(. >= 0.1))
heatmap.data <- column_to_rownames(heatmap.data, var = 'gene')
heatmap.data <- as.matrix(heatmap.data)
heatmap.data <- heatmap.data[rowSums(heatmap.data) != 0, ]

# Row annotations
# Alshetaiwi et al (Sci Immunol 2020)
gene_anno <- c(
  'MDSC',
  'PMN',
  'PMN',
  'PMN',
  'MDSC',
  'PMN',
  'PMN',
  'MDSC',
  'PMN',
  'MDSC',
  'PMN',
  'PMN',
  'PMN',
  'MDSC',
  'PMN',
  'MDSC',
  'PMN'
)

# Lung et al (Cell Metab 2022) MDSC
gene_anno <- c(
  'G-MDSC (Type I)',
  'G-MDSC (Type I)',
  'G-MDSC (Type II)',
  'G-MDSC (Type II)',
  'G-MDSC (Type II)',
  'G-MDSC (Type I)',
  'G-MDSC (Type II)',
  'G-MDSC (Type II)',
  'G-MDSC (Type II)',
  'M-MDSC',
  'G-MDSC (Type I)',
  'G-MDSC (Type II)',
  'G-MDSC (Type I)',
  'G-MDSC (Type I)',
  'G-MDSC (Type II)',
  'M-MDSC',
  'G-MDSC (Type II)',
  'G-MDSC (Type II)',
  'G-MDSC (Type II)',
  'G-MDSC (Type I)',
  'G-MDSC (Type II)',
  'G-MDSC (Type II)',
  'G-MDSC (Type I)'
)

# Blood v Tissue
col_anno <- c(
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Blood',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue',
  'Tissue'
)

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = circlize::colorRamp2(c(0, htmp_range), c("white", "red"))

ht <- ComplexHeatmap::Heatmap(
  heatmap.data,
  column_title = htmp.title,
  name = 'Average\nNormalized\nCounts',
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
    }
)

ComplexHeatmap::draw(ht, padding = unit(c(10, 5, 2, 5), 'mm')) # bottom, left, top, right paddings

rm(
  gene_list,
  gene_anno,
  htmp_range,
  htmp.title,
  col_fun,
  anc,
  heatmap.data
)

gc()

# Top pathways
top_pathways <- pathways %>%
  arrange(desc(NES)) %>%
  group_by(cluster_num, sample_subset) %>%
  do(head(., n = 10))

clust_num <- 3 # Choose your cluster of interest here

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
  's1 blood',
  's1 tissue',
  's2 blood',
  's2 tissue',
  's5 blood',
  's5 tissue'
)

heatmap.data <- heatmap.data[row_order, ]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = circlize::colorRamp2(
  c(-htmp_range, 0, htmp_range),
  c('blue', 'white', 'red')
)

ht <- ComplexHeatmap::Heatmap(
  heatmap.data,
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

ComplexHeatmap::draw(ht, padding = unit(c(70, 5, 2, 5), 'mm')) # bottom, left, top, right paddings

# Top DE genes
top_genes_up <- de %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  do(head(., n = 10))

top_genes_down <- de %>%
  arrange(desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  do(tail(., n = 10))

top_genes <- rbind(top_genes_up, top_genes_down)

clust <- 3 # Choose your cluster of interest here

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
  paste0(clust_name, '_s1 blood'),
  paste0(clust_name, '_s1 tissue'),
  paste0(clust_name, '_s2 blood'),
  paste0(clust_name, '_s2 tissue'),
  paste0(clust_name, '_s5 blood'),
  paste0(clust_name, '_s5 tissue')
)

heatmap.data <- heatmap.data[, col_order]

col_order <- colnames(heatmap.data)
row_order <- rownames(heatmap.data)

heatmap2.data <- heatmap2.data[row_order, col_order]

htmp_range <- c(max(heatmap.data), -min(heatmap.data))
htmp_range <- max(htmp_range)
htmp_range <- ceiling(htmp_range)

col_fun = circlize::colorRamp2(c(-htmp_range, 0, htmp_range), c('blue', 'white', 'red'))

ht <- ComplexHeatmap::Heatmap(
  heatmap.data,
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
sling <- getCurves(lineages, shrink = 1L, reweight = T, reassign = T, thresh = 0.001, maxit = 10L, stretch = 2L, smoother = 'smooth.spline', shrink.method = 'cosine') ## Parameters taken from dynverse ti_slingshot code

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

rm(
  clust,
  i
)

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