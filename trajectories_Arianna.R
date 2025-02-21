sce <- readRDS('~/Downloads/peds_labelled_data3_Cells3.0.rds')
require(scuttle)
sce <- logNormCounts(sce)

# (1) just relevant markers
sub <- sce[c('CD28','Lag3', 'Tim3', '41BB', 'CD45', 'CD8a', 'CD279', 'FoxP3', 'Ki67', 'CD3', 'HLADR', 'GranzymeB', 'CD274'), 
           which(sce$Clusters == 'CD8+ T-cells')]
# OR
# (2) all markers
sub <- sce[, which(sce$Clusters == 'CD8+ T-cells')]

require(BiocSingular)
pca <- runPCA(t(assay(sub,'logcounts')), rank=nrow(sub))
plot(pca$sdev^2)
#pairs(pca$x[,1:4], asp=1, col=colorby(sub$Subcluster))

require(uwot)
umap <- umap(pca$x)
plot(umap, asp=1, col=colorby(sub$Subcluster), main='All 29 Markers')
legendby(sub$Subcluster, pos='left')




