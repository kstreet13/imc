sce <- readRDS('~/Downloads/peds_labelled_data3_Cells3.0.rds')
require(scuttle)
require(slingshot)
sce <- logNormCounts(sce)

# (1) just relevant markers
# sub <- sce[c('CD28','Lag3', 'Tim3', '41BB', 'CD45', 'CD8a', 'CD279', 'FoxP3', 'Ki67', 'CD3', 'HLADR', 'GranzymeB', 'CD274'),
#            which(sce$Clusters == 'CD8+ T-cells')]
# OR
# (2) all markers
sub <- sce[, which(sce$Clusters == 'CD8+ T-cells')]

require(BiocSingular)
pca <- runPCA(t(assay(sub,'logcounts')), rank=nrow(sub))
plot(pca$sdev^2)
#pairs(pca$x[,1:4], asp=1, col=colorby(sub$Subcluster))

require(uwot)
set.seed(1)
umap <- umap(pca$x)
plot(umap, asp=1, col=colorby(sub$Subcluster), main='All 29 Markers')
legendby(sub$Subcluster, pos='left')

plot(pca$sdev^2)
K <- 8

# remove disconnected group
drop <- which(umap[,1] < -2 &
                  umap[,2] > 2 &
                  umap[,2] < 4)
sub2 <- sub[,-drop]
umap2 <- umap[-drop,]
pca$x <- pca$x[-drop, ]

#require(mclust)
# mc <- Mclust(pca$x[,1:K], G=2:6)
# pick starting cluster (which cluster has highest % of CD3+ CD8+ T-cells?)
set.seed(1)
km <- kmeans(pca$x[,1:K], centers = 7, nstart = 100)

# clus <- mc$classification
clus <- km$cluster
pct <- sapply(sort(unique(clus)), function(clID){
    sum(clus==clID & sub2$Subcluster=="CD3+ CD8+ T-cells") / sum(clus==clID)
})
st <- which.max(pct)

print(slingLineages(getLineages(pca$x[,1:K], clus, start.clus = st)))
pto <- slingshot(pca$x[,1:K], clus, start.clus = st)
emb <- embedCurves(pto, umap2)

# PLOTS #
#########
png(filename = '~/Desktop/CD8_celltype.png', width = 6, height = 6, units = 'in', res = 300)
plot(umap, asp=1, col=colorby(sub$Subcluster), main='CD8+ T Cells', pch = 16, cex = .5,
     xlab = 'UMAP-1', ylab = 'UMAP-2')
lines(SlingshotDataSet(emb))
dev.off()

png(filename = '~/Desktop/CD8_celltype_legend.png', width = 6, height = 6, units = 'in', res = 300)
plot.new()
legendby(sub$Subcluster, pos='left')
dev.off()


plot(umap2, asp=1, col=colorby(assay(pto)[,1], alpha = slingCurveWeights(pto)[,1]), main='All 29 Markers')
plot(umap2, asp=1, col=colorby(assay(pto)[,2], alpha = slingCurveWeights(pto)[,2]), main='All 29 Markers')
plot(umap2, asp=1, col=colorby(assay(pto)[,3], alpha = slingCurveWeights(pto)[,3]), main='All 29 Markers')

png(filename = '~/Desktop/CD8_lineages.png', width = 6, height = 6, units = 'in', res = 300)
plot(umap, asp=1, col = 'grey75', main='CD8+ T Cells', xlab = 'UMAP-1', ylab = 'UMAP-2', pch = 16, cex = .5)
points(umap2, asp=1, col = rgb(slingCurveWeights(pto,as.probs=TRUE)[,1], slingCurveWeights(pto,as.probs=TRUE)[,2], slingCurveWeights(pto,as.probs=TRUE)[,3]), pch = 16, cex = .5)
lines(SlingshotDataSet(emb))
dev.off()

png(filename = '~/Desktop/CD8_lineages_legend.png', width = 6, height = 6, units = 'in', res = 300)
plot.new()
legend('left', legend = c("Uncommitted","Lineage 1","Lineage 2","Lineage 3","Excluded"),
       pch = 16, col = c(rgb(1/3,1/3,1/3),rgb(1,0,0),rgb(0,1,0),rgb(0,0,1),'grey75'), bty = 'n')
dev.off()

# umap colored by pseudotime
png(filename = '~/Desktop/CD8_pseudotime.png', width = 6, height = 6, units = 'in', res = 300)
plot(umap, asp=1, col = 'grey75', main='CD8+ T Cells', xlab = 'UMAP-1', ylab = 'UMAP-2', pch=16, cex=.5)
points(umap2, asp=1, col = colorby(slingAvgPseudotime(pto)), pch=16, cex=.5)
lines(SlingshotDataSet(emb))
dev.off()

png(filename = '~/Desktop/CD8_pseudotime_legend.png', width = 6, height = 6, units = 'in', res = 300)
plot(0:1,0:1, col = 'white', axes=FALSE, xlab='', ylab='')
li <- as.raster(matrix(colorby(seq(1,0,length.out=50)), ncol=1))
rasterImage(li, 0, .25, .1, .75)
text(c(.1,.1),c(.25,.75), c('Initial','Terminal'), pos = 4)
dev.off()

# umap colored by group
png(filename = '~/Desktop/CD8_group.png', width = 6, height = 6, units = 'in', res = 300)
plot(umap, asp=1, col = groupcol[sub$group], main='CD8+ T Cells', xlab = 'UMAP-1', ylab = 'UMAP-2', pch=16, cex=.5)
lines(SlingshotDataSet(emb))
dev.off()

png(filename = '~/Desktop/CD8_group_legend.png', width = 6, height = 6, units = 'in', res = 300)
plot.new()
legend('left', legend = names(groupcol)[1:2], col = groupcol[1:2], pch = 16, bty='n')
dev.off()


# density plots
lineage_density <- function(sub2, pto, lineage){
    groupcol <- c(NR = '#E0B0FF', TCMR = '#DA70D6', CR = '#800080')
    dl <- lapply(unique(sub2$group), function(g){
        d <- density(slingPseudotime(pto, na = FALSE)[which(sub2$group == g), lineage],
                     weights = slingCurveWeights(pto)[which(sub2$group == g), lineage])
        d$x <- c(min(d$x), d$x, max(d$x))
        d$y <- c(0, d$y, 0) * mean(sub$group == g)
        return(d)
    })
    names(dl) <- unique(sub2$group)
    
    plot(range(c(dl[[1]]$x, dl[[2]]$x)), range(c(dl[[1]]$y, dl[[2]]$y)), col='white',
         main = paste('Lineage',lineage), xlab = 'Pseudotime', ylab='density')
    polygon(dl[[1]], col = groupcol[names(dl)[1]])
    polygon(dl[[2]], col = groupcol[names(dl)[2]])
    polygon(dl[[1]], col = 'transparent', lty = 2)
    polygon(dl[[2]], col = 'transparent', lty = 2)
    legend('topleft', legend=names(groupcol), fill=groupcol, bty='n')
}

lineage_density(sub2, pto, 1)
lineage_density(sub2, pto, 2)
lineage_density(sub2, pto, 3)
