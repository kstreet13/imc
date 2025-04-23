# CD4 T cells
# might be two starting points:
# - *resident memory T cell
# - generic CD3/CD4

# 2 or 3 endpoints, similar to CD8 (+/- Proliferating)

sce <- readRDS('~/Downloads/peds_labelled_data3_Cells3.0.rds')
require(scuttle)
require(slingshot)
sce <- logNormCounts(sce)

# (1) just relevant markers (don't know these for CD4)
# sub <- sce[c('CD28','Lag3', 'Tim3', '41BB', 'CD45', 'CD8a', 'CD279', 'FoxP3', 'Ki67', 'CD3', 'HLADR', 'GranzymeB', 'CD274'),
#            which(sce$Clusters == 'CD8+ T-cells')]
# OR
# (2) all markers
sub <- sce[, which(sce$Clusters == 'CD4+ T-cells')]

require(BiocSingular)
pca <- runPCA(t(assay(sub,'logcounts')), rank=nrow(sub))
plot(pca$sdev^2)
#pairs(pca$x[,1:4], asp=1, col=colorby(sub$Subcluster))

require(uwot)
set.seed(1)
umap <- umap(pca$x)
plot(umap, asp=1, cex=.5, col=colorby(sub$Subcluster), main='All 29 Markers')
legendby(sub$Subcluster, pos='topright', cex=.75)

plot(pca$sdev^2)
K <- 5

require(mclust)
mc <- Mclust(pca$x[,1:K], G=2:8)
#
set.seed(1)
km <- kmeans(pca$x[,1:K], centers = 7, nstart = 100)
set.seed(1)
km <- kmeans(umap, centers = 6, nstart = 100)



# pick starting cluster 
# clus <- mc$classification
clus <- km$cluster
pct <- sapply(sort(unique(clus)), function(clID){
    sum(clus==clID & sub$Subcluster=="Resident memory CD4+ T-cells") / sum(clus==clID)
})
st <- which.max(pct)

print(slingLineages(getLineages(pca$x[,1:K], clus, start.clus = st)))
pto <- slingshot(umap, clus, start.clus = st)
emb <- embedCurves(pto, umap)

# PLOTS #
#########
#png(filename = '~/Desktop/CD8_celltype.png', width = 6, height = 6, units = 'in', res = 300)
plot(umap, asp=1, col=colorby(sub$Subcluster), main='CD4+ T Cells', pch = 16, cex = .5,
     xlab = 'UMAP-1', ylab = 'UMAP-2')
lines(SlingshotDataSet(pto))
#dev.off()


