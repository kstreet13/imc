require(SingleCellExperiment); require(slingshot)

sce <- readRDS('~/Downloads/adult_IMC_data.rds')


set.seed(1)
sub <- sce[c('CD16', 'CD11b', 'CD45', 'Foxp3', 'CD163', 'CD68', 'Ki67', 'HLADR'), 
           which(sce$Clusters == 'Monocytes')]
# Monocyte subclusters: CD16, CD11b, CD45, Foxp3, CD163, CD68, Ki67, HLADR
reducedDim(sub,'umap') <- uwot::umap(t(assay(sub,'exprs')))
reducedDim(sub,'umap') <- cbind(-reducedDim(sub,'umap')[,1], reducedDim(sub,'umap')[,2])
# remove outliers
drop <- which(reducedDim(sub,'umap')[,1] > -2 &
                  reducedDim(sub,'umap')[,1] < 2 &
                  reducedDim(sub,'umap')[,2] > 4)
sub <- sub[,-drop]

km <- kmeans(reducedDim(sub,'umap'), centers = 3)
st <- which.max(km$centers[,1])
sub <- slingshot(sub, clusterLabels = km$cluster, start.clus = st)



png(filename = '~/Desktop/monocyte.png', width = 3000, height = 3000, res = 150)

layout(matrix(1:4, 2,2, byrow = TRUE))
par(mar = c(5,4,4,2))


groupcol <- c(NR = '#E0B0FF', TCMR = '#DA70D6', CR = '#800080')
plot(reducedDim(sub,'umap'), asp=1, pch=16, col=groupcol[sub$group], xlab = 'UMAP-1', ylab='UMAP-2')
legend('bottomleft', legend=names(groupcol), pch=16, col=groupcol, bty='n')

subcluscol <- c('#FFEA00', '#bf812d', '#ae017e', '#0096FF')
names(subcluscol) <- c('Non-classical monocytes', 'Classical monocytes', 'Intermediate monocytes', 'Activated monocytes')
plot(reducedDim(sub,'umap'), asp=1, pch=16, col=subcluscol[sub$Subcluster], xlab = 'UMAP-1', ylab='UMAP-2')
legend('bottomleft', legend=names(subcluscol), pch=16, col=subcluscol, bty='n')

plot(reducedDim(sub,'umap'), asp=1, pch=16, col=colorby(slingAvgPseudotime(sub)), xlab = 'UMAP-1', ylab='UMAP-2')
lines(SlingshotDataSet(sub), lwd=3, col='white')
lines(SlingshotDataSet(sub), lwd=2, col='black')


dl <- lapply(c('TCMR','CR','NR'), function(g){
    d <- density(slingAvgPseudotime(sub)[which(sub$group == g)])
    d$x <- c(min(d$x), d$x, max(d$x))
    d$y <- c(0, d$y, 0) * mean(sub$group == g)
    return(d)
})

plot(range(c(dl[[1]]$x, dl[[2]]$x, dl[[3]]$x)), range(c(dl[[1]]$y, dl[[2]]$y, dl[[3]]$y)), col='white',
     , xlab = 'Pseudotime', ylab='density')
polygon(dl[[1]], col = groupcol['TCMR'])
polygon(dl[[2]], col = groupcol['CR'])
polygon(dl[[3]], col = groupcol['NR'])
polygon(dl[[1]], col = 'transparent', lty = 2)
polygon(dl[[2]], col = 'transparent', lty = 2)
polygon(dl[[3]], col = 'transparent', lty = 2)
legend('topleft', legend=names(groupcol), fill=groupcol, bty='n')

dev.off()






par(mar = c(5,4,4,2)+.1)

