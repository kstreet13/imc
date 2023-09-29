
require(SingleCellExperiment)

sce <- readRDS('data/adult_IMC_data.rds')
assay(sce,'unscaled') <- asinh(assay(sce,'counts')/5)

{
    require(uwot)
    sixPlotSummary <- function(sub, file='~/Desktop/summary.png'){
        subumap <- umap(t(assay(sub,'exprs')))
        subpca <- prcomp(t(assay(sub,'exprs')))
        
        png(filename = file, width = 4000, height = 2000, res = 150)
        layout(matrix(1:6, nrow=2))
        barplot(table(sub$group), col = brewer.pal(3,'Set1'))
        plot(subumap, pch=16, col = colorby(sub$Subcluster, alpha=.5), asp=1, main='UMAP - subcluster')
        legendby(sub$Subcluster, pos='bottomleft')
        plot(subpca$x[,2:1], col = colorby(sub$group, alpha=.5), asp=1, main='PCA (1 vs. 2)', xlab='PC-2', ylab='PC-1')
        plot(subumap, pch=16, col = colorby(sub$group, alpha=.5), asp=1, main='UMAP - group')
        plot(subpca$x[,c(3,1)], col = colorby(sub$group, alpha=.5), asp=1, main='PCA (1 vs. 3)', xlab='PC-3', ylab='PC-1')
        plot(subpca$x[,3:2], col = colorby(sub$group, alpha=.5), asp=1, main='PCA (2 vs. 3)', xlab='PC-3', ylab='PC-2')
        dev.off()
    }
} # sixPlotSummary

barplot(table(sce$Clusters), las=2)
barplot(table(sce$group), las=2)

require(uwot)
umap <- umap(t(assay(sce,'exprs')))
umap2 <- umap(t(assay(sce,'unscaled')))
pca <- prcomp(t(assay(sce,'exprs')))

png(filename = '~/Desktop/umap2.png', width = 1500, height = 1500, res = 150)
plot(umap, asp=1, col=colorby(sce$Clusters))
legendby(sce$Clusters, pos='topleft')
dev.off()

png(filename = '~/Desktop/pca2.png', width = 3000, height = 3000, res = 150)
pairs(pca[,1:4], asp=1, col=colorby(sce$Clusters))
dev.off()


##########################
# cluster-level analysis #
##########################
# Maybe do each metacluster individually and analyze their respective
# subclusters? And maybe possibly also look at 2 metaclusters together such as
# CD4+ and CD8+ T cells, or T cells and macrophages/monocytes? I think the
# metaclusters we are more interested in are those that have subclusters (CD4,
# CD8, macrophages, monocytes, B cells), since it would be interesting to see if
# there’s some sort of progression these cells have when going from no rejection
# to rejection.

# summaries 

# CD8+ T-cells (scaled counts look OK in this subgroup)
sub <- sce[, which(sce$Clusters == 'CD8+ T-cells')]
sixPlotSummary(sub, file='~/Desktop/cd8.png')
# CD4+ T-cells (scaled counts have big outliers in this subgroup)
sub <- sce[, which(sce$Clusters == 'CD4+ T-cells')]
sixPlotSummary(sub, file='~/Desktop/cd4.png')
# B cells (???)
sub <- sce[, which(sce$Clusters == 'B cells')]
sixPlotSummary(sub, file='~/Desktop/Bcell.png')
# Monocytes (???)
sub <- sce[, which(sce$Clusters == 'Monocytes')]
sixPlotSummary(sub, file='~/Desktop/monocyte.png')
# Macrophages (???)
sub <- sce[, which(sce$Clusters == 'Macrophages')]
sixPlotSummary(sub, file='~/Desktop/macrophage.png')

# combos
sub <- sce[, which(sce$Clusters %in% c('CD8+ T-cells','CD4+ T-cells'))]
sixPlotSummary(sub, file='~/Desktop/cd4_cd8.png')
sub <- sce[, which(sce$Clusters %in% c('CD8+ T-cells','CD4+ T-cells','Monocytes'))]
sixPlotSummary(sub, file='~/Desktop/tcell_monocyte.png')
sub <- sce[, which(sce$Clusters %in% c('CD8+ T-cells','CD4+ T-cells','Macrophages'))]
sixPlotSummary(sub, file='~/Desktop/tcell_macrophage.png')




# subclusters with specific markers
# CD8+ T-cells (scaled counts look OK in this subgroup)
sub <- sce[c('CD28', 'CD16', 'CD11b', 'CD45', 'CD8', 'CD279', 'Foxp3', 'Ki67', 'CD3', 'HLADR', 'GranzymeB'), 
           which(sce$Clusters == 'CD8+ T-cells')]
sixPlotSummary(sub, file='~/Desktop/cd8.png')
# CD4+ T-cells (scaled counts have big outliers in this subgroup)
sub <- sce[c('CD28', 'CD16', 'CD11b', 'CD45', 'CD4', 'CD279', 'Foxp3', 'Ki67', 'CD3', 'HLADR'), 
           which(sce$Clusters == 'CD4+ T-cells')]
sixPlotSummary(sub, file='~/Desktop/cd4.png')
# B cells (???)
sub <- sce[c('CD45', 'CD20', 'CD279', 'Foxp3', 'Ki67', 'HLADR'), which(sce$Clusters == 'B cells')]
sixPlotSummary(sub, file='~/Desktop/Bcell.png')
# Monocytes (???)
sub <- sce[c('CD16', 'CD11b', 'CD45', 'Foxp3', 'CD163', 'CD68', 'Ki67', 'HLADR'), 
           which(sce$Clusters == 'Monocytes')]
sixPlotSummary(sub, file='~/Desktop/monocyte.png')
# Macrophages (???)
sub <- sce[c('CD16', 'CD11b', 'CD45', 'Foxp3', 'CD163', 'CD68', 'Ki67', 'HLADR'), 
           which(sce$Clusters == 'Macrophages')]
sixPlotSummary(sub, file='~/Desktop/macrophage.png')








sub <- sce[, which(sce$Clusters == 'CD8+ T-cells')]

subumap <- umap(t(assay(sub,'exprs')))
subumap3 <- umap(t(assay(sub,'exprs')), n_components = 3)
subpca <- prcomp(t(assay(sub,'exprs')))
# markers of interest:
# CD28, CD16, CD11b, CD45, CD8, CD279, Foxp3, Ki67, CD3, HLADR, GranzymeB
intumap <- umap(t(assay(sub,'exprs'))[,c('CD28', 'CD16', 'CD11b', 'CD45', 'CD8', 'CD279', 'Foxp3', 'Ki67', 'CD3', 'HLADR', 'GranzymeB')])

plot(subumap, col = colorby(sub$Subcluster, alpha=.5), asp=1)
plot(subumap, col = colorby(sub$group, alpha=.5), asp=1)
pairs(subpca$x[,1:4], asp=1, col=colorby(sub$group))

require(rgl)
plot3d(subpca$x[,2:4], aspect = 'iso', col = colorby(sub$group), size=5)
plot3d(subumap3, aspect = 'iso', col = colorby(sub$Subcluster), size=5)





# CD4+ T-cells (scaled counts have big outliers in this subgroup)
sub <- sce[, which(sce$Clusters == 'CD4+ T-cells')]
sixPlotSummary(sub, file='~/Desktop/cd4.png')

subumap1 <- umap(t(assay(sub,'exprs')))
subumap2 <- umap(t(assay(sub,'unscaled')))
subpca <- prcomp(t(assay(sub,'exprs')))

barplot(table(sub$group), col = brewer.pal(3,'Set1'))
plot(subumap, col = colorby(sub$Subcluster, alpha=.5), asp=1)
plot(subumap, col = colorby(sub$group, alpha=.5), asp=1)
pairs(subpca$x[,1:4], asp=1, col=colorby(sub$group))


png(filename = '~/Desktop/cd4.png', width = 5000, height = 2500, res = 150)
layout(matrix(1:6, nrow=2))
barplot(table(sub$group), col = brewer.pal(3,'Set1'))
plot(subumap2, col = colorby(sub$group, alpha=.5), asp=1, main='UMAP - *Unscaled*')
plot(subpca$x[,2:1], col = colorby(sub$group, alpha=.5), asp=1, xlab='PC-2', ylab='PC-1')
plot(subumap1, col = colorby(sub$group, alpha=.5), asp=1, main='UMAP - Scaled')
plot(subpca$x[,c(3,1)], col = colorby(sub$group, alpha=.5), asp=1, xlab='PC-3', ylab='PC-1')
plot(subpca$x[,3:2], col = colorby(sub$group, alpha=.5), asp=1, xlab='PC-3', ylab='PC-2')
dev.off()









require(cluster)
p2 <- pam(subpca$x[,1:10], k = 2)
p3 <- pam(subpca$x[,1:10], k = 3)
p4 <- pam(subpca$x[,1:10], k = 4)
p5 <- pam(subpca$x[,1:10], k = 5)
p6 <- pam(subpca$x[,1:10], k = 6)
p7 <- pam(subpca$x[,1:10], k = 7)
p8 <- pam(subpca$x[,1:10], k = 8)
p9 <- pam(subpca$x[,1:10], k = 9)
mean(p2$silinfo$widths[,'sil_width'])
mean(p3$silinfo$widths[,'sil_width'])
mean(p4$silinfo$widths[,'sil_width'])
mean(p5$silinfo$widths[,'sil_width'])
mean(p6$silinfo$widths[,'sil_width'])
mean(p7$silinfo$widths[,'sil_width'])
mean(p8$silinfo$widths[,'sil_width'])
mean(p9$silinfo$widths[,'sil_width'])

#







# trying to get clever by taking neighbors of a cluster

require(BiocNeighbors)
mnn <- findMutualNN(t(assay(sce,'exprs')), t(assay(sce,'exprs')), k1 = 21)

# snn <- bluster::makeSNNGraph(t(assay(sce,'exprs')))

knn <- findKNN(t(assay(sce,'exprs')), 10)$index
mnn <- apply(knn,2,function(k){
    keep <- rowMins(abs(knn[k,] - 1:nrow(knn))) == 0
    out <- k
    out[!keep] <- NA
    return(out)
})


unique(sce$Clusters)
group <- "CD8+ T-cells"

ind <- which(sce$Clusters == group)
neighbors <- unique(as.numeric(mnn[ind, ]))
ind <- unique(c(ind, neighbors))
ind <- ind[!is.na(ind)]

# X <- t(assay(sce,'exprs'))[ind,]
# X <- cbind(X, log1p(colData(sce)$area)[ind], colData(sce)$eccentricity[ind])

g.umap <- umap(t(assay(sce,'exprs'))[ind,])


ind0 <- ind %in% which(sce$Clusters == group)
ind2 <- !ind0

png(filename = '~/Desktop/gumap.png', width = 1500, height = 1500, res = 150)
plot(range(g.umap[,1]), range(g.umap[,2]), asp=1, col='white', main = group)
points(g.umap[which(ind0),], pch=16, col=1)
points(g.umap[which(ind2),], pch=1, col=colorby(sce$Clusters[which(ind2)], alpha=.5))
legendby(as.character(sce$Clusters[which(ind2)]), pos='bottomright')
dev.off()




# NR: #E0B0FF
# TCMR: #DA70D6
# CR: #800080

