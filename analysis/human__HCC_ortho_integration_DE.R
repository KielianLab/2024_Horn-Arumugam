## Load in packages
packages <- c(
  'tidyverse',
  'Seurat',
  'SingleCellExperiment',
  'MAST',
  'future'
)

install.packages(packages, repos = 'http://cran.us.r-project.org')
invisible(lapply(packages, library, character.only = TRUE))

rm(packages)

## Set options
options(future.globals.maxSize = 4000 * 1024^5)
set.seed(12345)

plan('multicore')

## Bring in objects
brain.integrated <- readRDS('/work/kielianlab/cmhorn/ortho.integrated.rds')

## Differential expression analysis
DefaultAssay(brain.integrated) <- 'RNA'
brain.integrated <- NormalizeData(brain.integrated, verbose = TRUE)
brain.integrated <- ScaleData(brain.integrated, verbose = TRUE)

DE <- FindAllMarkers(brain.integrated, min.pct = 0, logfc.threshold = 0, test.use = 'MAST')

write.csv(DE, '/work/kielianlab/cmhorn/output/DE/full DE.csv')

for (i in seq_along(levels(DE$cluster))) {
  if(i != length(DE$cluster)) {
    single.DE <- DE %>% filter(DE$cluster %in% c(levels(DE$cluster)[i]))
    write.csv(single.DE, paste0('/work/kielianlab/cmhorn/output/DE/DE_', i-1, '.csv'))
  }
}

rm(
  i,
  single.DE,
  DE
)

gc()

## Within cluster DE Analysis
c00 <- subset(cellplex.integrated, idents = c('Granulocytes 1'))
c01 <- subset(cellplex.integrated, idents = c('Granulocytes/Myelocytes 1'))
c02 <- subset(cellplex.integrated, idents = c('Granulocytes 2'))
c03 <- subset(cellplex.integrated, idents = c('Granulocytes 3'))
c04 <- subset(cellplex.integrated, idents = c('Granulocytes 4'))
c05 <- subset(cellplex.integrated, idents = c('Granulocytes 5'))
c06 <- subset(cellplex.integrated, idents = c('Granulocytes 6'))
c07 <- subset(cellplex.integrated, idents = c('T Cells'))
c08 <- subset(cellplex.integrated, idents = c('Granulocytes 7'))
c09 <- subset(cellplex.integrated, idents = c('Granulocytes/Myelocytes 2'))
c10 <- subset(cellplex.integrated, idents = c('Granulocytes 8'))
c11 <- subset(cellplex.integrated, idents = c('Granulocytes 9'))
c12 <- subset(cellplex.integrated, idents = c('Granulocytes 10'))
c13 <- subset(cellplex.integrated, idents = c('C13'))
c14 <- subset(cellplex.integrated, idents = c('Granulocytes 11'))
c15 <- subset(cellplex.integrated, idents = c('NK/T Cells'))
c16 <- subset(cellplex.integrated, idents = c('Granulocytes 12'))

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
  c.list[[i]]$cell_origin <- paste(Idents(c.list[[i]]), c.list[[i]]$Sample_ID, sep = '_')
  c.list[[i]]$cell <- Idents(c.list[[i]])
  Idents(c.list[[i]]) <- 'cell_origin'
}

for(i in 1:length(c.list)){
  DE <- FindAllMarkers(c.list[[i]], min.pct = 0, logfc.threshold = 0, test.use = 'MAST')
  DE <- DE %>%
    arrange(desc(avg_log2FC))
  write.csv(DE, paste0('/work/kielianlab/cmhorn/output/DE/within_cluster/DE_', i-1, '.csv'))
}

rm(
  i,
  c.list,
  DE
)

gc()