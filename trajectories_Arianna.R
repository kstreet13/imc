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


require(mclust)
mc <- Mclust(pca$x)
# pick starting cluster (which cluster has highest % of CD3+ CD8+ T-cells?)
pct <- sapply(sort(unique(mc$classification)), function(clID){
    sum(mc$classification==clID & sub$Subcluster=="CD3+ CD8+ T-cells") / sum(mc$classification==clID)
})
st <- which.max(pct)


require(slingshot)

for(k in 5:29){
    cat('k =',k,'\n')
    #print(slingLineages(getLineages(pca$x[,1:k], sub$Subcluster, start.clus = "CD3+ CD8+ T-cells")))
    print(slingLineages(getLineages(pca$x[,1:k], mc$classification, start.clus = st)))
    cat('\n')
}

pto <- slingshot(pca$x[,1:20], mc$classification, start.clus = st)

emb <- embedCurves(pto, umap)

plot(umap, asp=1, col=colorby(sub$Subcluster), main='All 29 Markers')
lines(SlingshotDataSet(emb))

plot(umap, asp=1, col=colorby(assay(pto)[,1], alpha = slingCurveWeights(pto)[,1]), main='All 29 Markers')

plot(umap, asp=1, col= rgb(slingCurveWeights(pto,as.probs=TRUE)[,1], slingCurveWeights(pto,as.probs=TRUE)[,2], slingCurveWeights(pto,as.probs=TRUE)[,3]), main='All 29 Markers')



