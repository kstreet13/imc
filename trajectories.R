
# My best analysis of each metacluster

require(SingleCellExperiment)
sce <- readRDS('data/adult_IMC_data.rds')
group_col <- c(NR = '#E0B0FF', TCMR = '#DA70D6', CR = '#800080')
sce$group_col <- '#E0B0FF'
sce$group_col[sce$group == 'CR'] <- '#800080'
sce$group_col[sce$group == 'TCMR'] <- '#DA70D6'
sce$group <- as.character(sce$group)
require(uwot)
require(slingshot)



################
# CD8+ T-cells #
################
set.seed(1)
sub <- sce[c('CD28', 'CD16', 'CD11b', 'CD45', 'CD8', 'CD279', 'Foxp3', 'Ki67', 'CD3', 'HLADR', 'GranzymeB'), 
           which(sce$Clusters == 'CD8+ T-cells')]
reducedDim(sub,'umap') <- umap(t(assay(sub,'exprs')))
# remove outliers
keep <- which(reducedDim(sub,'umap')[,2] > -6 | (reducedDim(sub,'umap')[,1] > -4.5 & reducedDim(sub,'umap')[,1] < 2))
sub <- sub[, keep]

km <- kmeans(reducedDim(sub,'umap'), centers = 5)
start.clus <- which.max(km$centers[,1])
end.clus <- c(which.min(km$centers[,1]), which.min(km$centers[,2]))
sub <- slingshot(sub, start.clus = start.clus, end.clus = end.clus, clusterLabels = km$cluster)

{
    png(filename = '~/Desktop/trajectory_cd8.png', width = 2000, height = 2000, res = 100)
    layout(matrix(1:4, 2,2, byrow = TRUE))
    plot(reducedDim(sub,'umap'), col = alpha(sub$group_col, alpha=.5), asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', pch=16, legend = names(group_col), col = group_col, bty='n')
    
    subcluscols <- c('#FFEA00', '#ae017e', '#bf812d', '#0096FF', '#ff77ff')
    subclusnames <- c('CD3+ CD8+ T-cells', 'Cytotoxic T-cells', 'PD1+ CD8+ T-cells', 'PD1+CD28+ CD8+ T-cells', 'Proliferating CD8+ T-cells')
    names(subcluscols) <- subclusnames
    plot(reducedDim(sub,'umap'), col = subcluscols[sub$Subcluster], asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', legend = sort(subclusnames), pch=16, col = subcluscols[order(subclusnames)], bty='n')
    
    plot(reducedDim(sub,'umap'), col = colorby(slingAvgPseudotime(sub)),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    lines(SlingshotDataSet(sub), lwd=3, col='white')
    lines(SlingshotDataSet(sub), lwd=2)
    # group-level densities
    d <- lapply(names(group_col), function(g){
        out <- density(slingAvgPseudotime(sub)[which(sub$group==g)], na.rm = TRUE)
        out$y <- c(0,out$y,0) * mean(sub$group == g)
        out$x <- c(min(out$x), out$x, max(out$x))
        return(out)
    })
    names(d) <- names(group_col)
    plot(c(min(c(d[[1]]$x,d[[2]]$x,d[[3]]$x)), max(c(d[[1]]$x,d[[2]]$x,d[[3]]$x))), 
         c(0, max(c(d[[1]]$y,d[[2]]$y,d[[3]]$y))), col='white',
         xlab='Pseudotime', ylab='density')
    ord <- names(sort(table(sub$group), decreasing = TRUE))
    for(g in ord){
        polygon(d[[g]], col=group_col[g])
    }
    for(g in ord){
        polygon(d[[g]], col='transparent', lty=2)
    }
    legend('topleft', legend = names(group_col), fill = group_col, bty='n')
    dev.off()
    layout(1)
}

{
    png(filename = '~/Desktop/features_cd8.png', width = 2000, height = 2500, res = 100)
    layout(matrix(1:24, nrow=6,ncol=4, byrow = TRUE))
    sub2 <- sce[, colnames(sub)]
    for(j in 1:nrow(sub2)){
        plot(slingAvgPseudotime(sub), assay(sub2,'exprs')[j,], main = rownames(sub2)[j],
             xlab='Pseudotime', col = alpha(sub$group_col, alpha=.5))
    }
    dev.off()
    layout(1)
}


################
# CD4+ T-cells #
################
set.seed(1)
sub <- sce[, which(sce$Clusters == 'CD4+ T-cells')]
reducedDim(sub,'umap') <- umap(t(assay(sub,'exprs')))
# remove outliers
keep <- which(reducedDim(sub,'umap')[,1] > -6)
sub <- sub[, keep]

km <- kmeans(reducedDim(sub,'umap'), centers = 2)
start.clus <- which.min(km$centers[,2])
sub <- slingshot(sub, start.clus = start.clus, clusterLabels = km$cluster)



{
    png(filename = '~/Desktop/trajectory_cd4.png', width = 2000, height = 2000, res = 100)
    layout(matrix(1:4, 2,2, byrow = TRUE))
    plot(reducedDim(sub,'umap'), col = alpha(sub$group_col, alpha=.5), asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', pch=16, legend = names(group_col), col = group_col, bty='n')
    
    subcluscols <- c('#FFEA00', 'orange', '#ff77ff', '#6a87a3', '#0096FF', '#bf812d', '#00ffff', '#ae017e', '#7F00FF')
    subclusnames <- c('CD3+ CD4+ T-cells', 'Resident memory CD4+ T-cells', 'Proliferating CD4+ T-cells', 'Naive CD4+ T-cells', 'PD1+ CD4+ T-cells', 'Activated CD4+ T-cells', 'CD16+ CD4+ T-cells', 'HLADR+ CD4+ Tregs', 'HLADR- CD4+ Tregs')
    names(subcluscols) <- subclusnames
    plot(reducedDim(sub,'umap'), col = subcluscols[sub$Subcluster], asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', legend = sort(subclusnames), pch=16, col = subcluscols[order(subclusnames)], bty='n')
    
    plot(reducedDim(sub,'umap'), col = colorby(slingAvgPseudotime(sub)),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    lines(SlingshotDataSet(sub), lwd=3, col='white')
    lines(SlingshotDataSet(sub), lwd=2)
    # group-level densities
    d <- lapply(names(group_col), function(g){
        out <- density(slingAvgPseudotime(sub)[which(sub$group==g)], na.rm = TRUE)
        out$y <- c(0,out$y,0) * mean(sub$group == g)
        out$x <- c(min(out$x), out$x, max(out$x))
        return(out)
    })
    names(d) <- names(group_col)
    plot(c(min(c(d[[1]]$x,d[[2]]$x,d[[3]]$x)), max(c(d[[1]]$x,d[[2]]$x,d[[3]]$x))), 
         c(0, max(c(d[[1]]$y,d[[2]]$y,d[[3]]$y))), col='white',
         xlab='Pseudotime', ylab='density')
    ord <- names(sort(table(sub$group), decreasing = TRUE))
    for(g in ord){
        polygon(d[[g]], col=group_col[g])
    }
    for(g in ord){
        polygon(d[[g]], col='transparent', lty=2)
    }
    legend('topleft', legend = names(group_col), fill = group_col, bty='n')
    dev.off()
    layout(1)
}

{
    png(filename = '~/Desktop/features_cd4.png', width = 2000, height = 2500, res = 100)
    layout(matrix(1:24, nrow=6,ncol=4, byrow = TRUE))
    sub2 <- sce[, colnames(sub)]
    for(j in 1:nrow(sub2)){
        plot(slingAvgPseudotime(sub), assay(sub2,'exprs')[j,], main = rownames(sub2)[j],
             xlab='Pseudotime', col = alpha(sub$group_col, alpha=.5))
    }
    dev.off()
    layout(1)
}



###########
# B cells #
###########
set.seed(1)
sub <- sce[c('CD45', 'CD20', 'CD279', 'Foxp3', 'Ki67', 'HLADR'), which(sce$Clusters == 'B cells')]
reducedDim(sub,'umap') <- umap(t(assay(sub,'exprs')))
# remove outliers
keep <- which(reducedDim(sub,'umap')[,2] > -5)
sub <- sub[, keep]

km <- kmeans(reducedDim(sub,'umap'), centers = 2)
start.clus <- which.min(km$centers[,2])
sub <- slingshot(sub, start.clus = start.clus, clusterLabels = km$cluster)


{
    png(filename = '~/Desktop/trajectory_Bcell.png', width = 2000, height = 2000, res = 100)
    layout(matrix(1:4, 2,2, byrow = TRUE))
    plot(reducedDim(sub,'umap'), col = alpha(sub$group_col, alpha=.5), asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', pch=16, legend = names(group_col), col = group_col, bty='n')
    plot(reducedDim(sub,'umap'), col = colorby(sub$Subcluster),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legendby(sub$Subcluster, pos='bottomleft')
    plot(reducedDim(sub,'umap'), col = colorby(slingAvgPseudotime(sub)),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    lines(SlingshotDataSet(sub), lwd=3, col='white')
    lines(SlingshotDataSet(sub), lwd=2)
    # group-level densities
    d <- lapply(names(group_col), function(g){
        out <- density(slingAvgPseudotime(sub)[which(sub$group==g)], na.rm = TRUE)
        out$y <- c(0,out$y,0) * mean(sub$group == g)
        out$x <- c(min(out$x), out$x, max(out$x))
        return(out)
    })
    names(d) <- names(group_col)
    plot(c(min(c(d[[1]]$x,d[[2]]$x,d[[3]]$x)), max(c(d[[1]]$x,d[[2]]$x,d[[3]]$x))), 
         c(0, max(c(d[[1]]$y,d[[2]]$y,d[[3]]$y))), col='white',
         xlab='Pseudotime', ylab='density')
    ord <- names(sort(table(sub$group), decreasing = TRUE))
    for(g in ord){
        polygon(d[[g]], col=group_col[g])
    }
    for(g in ord){
        polygon(d[[g]], col='transparent', lty=2)
    }
    legend('topleft', legend = names(group_col), fill = group_col, bty='n')
    dev.off()
    layout(1)
}

{
    png(filename = '~/Desktop/features_Bcell.png', width = 2000, height = 2500, res = 100)
    layout(matrix(1:24, nrow=6,ncol=4, byrow = TRUE))
    sub2 <- sce[, colnames(sub)]
    for(j in 1:nrow(sub2)){
        plot(slingAvgPseudotime(sub), assay(sub2,'exprs')[j,], main = rownames(sub2)[j],
             xlab='Pseudotime', col = alpha(sub$group_col, alpha=.5))
    }
    dev.off()
    layout(1)
}



#############
# Monocytes #
#############
set.seed(1)
sub <- sce[c('CD16', 'CD11b', 'CD45', 'Foxp3', 'CD163', 'CD68', 'Ki67', 'HLADR'), 
           which(sce$Clusters == 'Monocytes')]
reducedDim(sub,'umap') <- umap(t(assay(sub,'exprs')))
# remove outliers
keep <- which(reducedDim(sub,'umap')[,2] < 5 | reducedDim(sub,'umap')[,1] < -2)
sub <- sub[, keep]

km <- kmeans(reducedDim(sub,'umap'), centers = 2)
start.clus <- which.max(km$centers[,1])
sub <- slingshot(sub, start.clus = start.clus, clusterLabels = km$cluster)

{
    png(filename = '~/Desktop/trajectory_monocyte.png', width = 2000, height = 2000, res = 100)
    layout(matrix(1:4, 2,2, byrow = TRUE))
    plot(reducedDim(sub,'umap'), col = alpha(sub$group_col, alpha=.5), asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', pch=16, legend = names(group_col), col = group_col, bty='n')
    plot(reducedDim(sub,'umap'), col = colorby(sub$Subcluster),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legendby(sub$Subcluster, pos='bottomleft')
    plot(reducedDim(sub,'umap'), col = colorby(slingAvgPseudotime(sub)),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    lines(SlingshotDataSet(sub), lwd=3, col='white')
    lines(SlingshotDataSet(sub), lwd=2)
    # group-level densities
    d <- lapply(names(group_col), function(g){
        out <- density(slingAvgPseudotime(sub)[which(sub$group==g)], na.rm = TRUE)
        out$y <- c(0,out$y,0) * mean(sub$group == g)
        out$x <- c(min(out$x), out$x, max(out$x))
        return(out)
    })
    names(d) <- names(group_col)
    plot(c(min(c(d[[1]]$x,d[[2]]$x,d[[3]]$x)), max(c(d[[1]]$x,d[[2]]$x,d[[3]]$x))), 
         c(0, max(c(d[[1]]$y,d[[2]]$y,d[[3]]$y))), col='white',
         xlab='Pseudotime', ylab='density')
    ord <- names(sort(table(sub$group), decreasing = TRUE))
    for(g in ord){
        polygon(d[[g]], col=group_col[g])
    }
    for(g in ord){
        polygon(d[[g]], col='transparent', lty=2)
    }
    legend('topleft', legend = names(group_col), fill = group_col, bty='n')
    dev.off()
    layout(1)
}

{
    png(filename = '~/Desktop/features_monocyte.png', width = 2000, height = 2500, res = 100)
    layout(matrix(1:24, nrow=6,ncol=4, byrow = TRUE))
    sub2 <- sce[, colnames(sub)]
    for(j in 1:nrow(sub2)){
        plot(slingAvgPseudotime(sub), assay(sub2,'exprs')[j,], main = rownames(sub2)[j],
             xlab='Pseudotime', col = alpha(sub$group_col, alpha=.5))
    }
    dev.off()
    layout(1)
}



##################
# M1 Macrophages #
##################
set.seed(1)
sub <- sce[, grep('M1', sce$Subcluster)]
reducedDim(sub,'umap') <- umap(t(assay(sub,'exprs')))
# remove outliers
keep <- which(reducedDim(sub,'umap')[,1] < 5 & reducedDim(sub,'umap')[,2] > -5)
sub <- sub[, keep]

km <- kmeans(reducedDim(sub,'umap'), centers = 2)
start.clus <- which.max(km$centers[,2])
sub <- slingshot(sub, start.clus = start.clus, clusterLabels = km$cluster)

{
    png(filename = '~/Desktop/trajectory_M1macrophage.png', width = 2000, height = 2000, res = 100)
    layout(matrix(1:4, 2,2, byrow = TRUE))
    plot(reducedDim(sub,'umap'), col = alpha(sub$group_col, alpha=.5), asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', pch=16, legend = names(group_col), col = group_col, bty='n')
    plot(reducedDim(sub,'umap'), col = colorby(sub$Subcluster),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legendby(sub$Subcluster, pos='bottomleft')
    plot(reducedDim(sub,'umap'), col = colorby(slingAvgPseudotime(sub)),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    lines(SlingshotDataSet(sub), lwd=3, col='white')
    lines(SlingshotDataSet(sub), lwd=2)
    # group-level densities
    d <- lapply(names(group_col), function(g){
        out <- density(slingAvgPseudotime(sub)[which(sub$group==g)], na.rm = TRUE)
        out$y <- c(0,out$y,0) * mean(sub$group == g)
        out$x <- c(min(out$x), out$x, max(out$x))
        return(out)
    })
    names(d) <- names(group_col)
    plot(c(min(c(d[[1]]$x,d[[2]]$x,d[[3]]$x)), max(c(d[[1]]$x,d[[2]]$x,d[[3]]$x))), 
         c(0, max(c(d[[1]]$y,d[[2]]$y,d[[3]]$y))), col='white',
         xlab='Pseudotime', ylab='density')
    ord <- names(sort(table(sub$group), decreasing = TRUE))
    for(g in ord){
        polygon(d[[g]], col=group_col[g])
    }
    for(g in ord){
        polygon(d[[g]], col='transparent', lty=2)
    }
    legend('topleft', legend = names(group_col), fill = group_col, bty='n')
    dev.off()
    layout(1)
}

{
    png(filename = '~/Desktop/features_M1macrophage.png', width = 2000, height = 2500, res = 100)
    layout(matrix(1:24, nrow=6,ncol=4, byrow = TRUE))
    sub2 <- sce[, colnames(sub)]
    for(j in 1:nrow(sub2)){
        plot(slingAvgPseudotime(sub), assay(sub2,'exprs')[j,], main = rownames(sub2)[j],
             xlab='Pseudotime', col = alpha(sub$group_col, alpha=.5))
    }
    dev.off()
    layout(1)
}

# highlight CRs
{
    png(filename = '~/Desktop/trajectory_M1macrophage_CR.png', width = 2000, height = 2000, res = 100)
    layout(matrix(1:4, 2,2, byrow = TRUE))
    plot(reducedDim(sub,'umap'), col = 'grey70', asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    points(reducedDim(sub,'umap')[which(sub$group=='CR'),], col = group_col['CR'], pch=16)
    legend('bottomleft', pch=16, legend = c('Non-CR','CR'), col = c('grey70',group_col['CR']), bty='n')
    plot(reducedDim(sub,'umap'), col = 'grey70', asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    points(reducedDim(sub,'umap')[which(sub$group=='CR'),], col = colorby(sub$Subcluster)[which(sub$group=='CR')], pch=16)
    legendby(sub$Subcluster, pos='bottomleft')
    plot(reducedDim(sub,'umap'), col = 'grey70', asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    points(reducedDim(sub,'umap')[which(sub$group=='CR'),], 
           col = colorby(slingAvgPseudotime(sub))[which(sub$group=='CR')],pch=16)
    lines(SlingshotDataSet(sub), lwd=3, col='white')
    lines(SlingshotDataSet(sub), lwd=2)
    # group-level densities
    d <- list()
    out <- density(slingAvgPseudotime(sub)[which(sub$group=='CR')], na.rm = TRUE)
    out$y <- c(0,out$y,0) * mean(sub$group == g)
    out$x <- c(min(out$x), out$x, max(out$x))
    d[[1]] <- out
    out <- density(slingAvgPseudotime(sub)[which(sub$group!='CR')], na.rm = TRUE)
    out$y <- c(0,out$y,0) * mean(sub$group == g)
    out$x <- c(min(out$x), out$x, max(out$x))
    d[[2]] <- out
    names(d) <- c('CR','Non-CR')
    plot(c(min(c(d[[1]]$x,d[[2]]$x)), max(c(d[[1]]$x,d[[2]]$x))), 
         c(0, max(c(d[[1]]$y,d[[2]]$y))), col='white',
         xlab='Pseudotime', ylab='density')
    polygon(d[[2]], col='grey70')
    polygon(d[[1]], col=group_col['CR'])
    polygon(d[[2]], col='transparent', lty=2)
    legend('topleft', legend = c('Non-CR','CR'), fill = c('grey70',group_col['CR']), bty='n')
    dev.off()
    layout(1)
}

{
    png(filename = '~/Desktop/features_M1macrophage_CR.png', width = 2000, height = 2500, res = 100)
    layout(matrix(1:24, nrow=6,ncol=4, byrow = TRUE))
    sub2 <- sce[, colnames(sub)]
    for(j in 1:nrow(sub2)){
        plot(slingAvgPseudotime(sub), assay(sub2,'exprs')[j,], main = rownames(sub2)[j],
             xlab='Pseudotime', col = 'grey70')
        points(slingAvgPseudotime(sub)[which(sub$group=='CR')], assay(sub2,'exprs')[j,which(sub$group=='CR')],
               col = group_col['CR'])
    }
    dev.off()
    layout(1)
}

##################
# M2 Macrophages #
##################
set.seed(1)
sub <- sce[, grep('M2', sce$Subcluster)]
reducedDim(sub,'umap') <- umap(t(assay(sub,'exprs')))
# remove outliers
keep <- which(reducedDim(sub,'umap')[,1] < 4.75)
sub <- sub[, keep]

km <- kmeans(reducedDim(sub,'umap'), centers = 2)
start.clus <- which.max(km$centers[,2])
sub <- slingshot(sub, start.clus = start.clus, clusterLabels = km$cluster)

{
    png(filename = '~/Desktop/trajectory_M2macrophage.png', width = 2000, height = 2000, res = 100)
    layout(matrix(1:4, 2,2, byrow = TRUE))
    plot(reducedDim(sub,'umap'), col = alpha(sub$group_col, alpha=.5), asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legend('bottomleft', pch=16, legend = names(group_col), col = group_col, bty='n')
    plot(reducedDim(sub,'umap'), col = colorby(sub$Subcluster),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    legendby(sub$Subcluster, pos='bottomleft')
    plot(reducedDim(sub,'umap'), col = colorby(slingAvgPseudotime(sub)),asp=1, xlab='UMAP-1', ylab='UMAP-2', pch=16)
    lines(SlingshotDataSet(sub), lwd=3, col='white')
    lines(SlingshotDataSet(sub), lwd=2)
    # group-level densities
    d <- lapply(names(group_col), function(g){
        out <- density(slingAvgPseudotime(sub)[which(sub$group==g)], na.rm = TRUE)
        out$y <- c(0,out$y,0) * mean(sub$group == g)
        out$x <- c(min(out$x), out$x, max(out$x))
        return(out)
    })
    names(d) <- names(group_col)
    plot(c(min(c(d[[1]]$x,d[[2]]$x,d[[3]]$x)), max(c(d[[1]]$x,d[[2]]$x,d[[3]]$x))), 
         c(0, max(c(d[[1]]$y,d[[2]]$y,d[[3]]$y))), col='white',
         xlab='Pseudotime', ylab='density')
    ord <- names(sort(table(sub$group), decreasing = TRUE))
    for(g in ord){
        polygon(d[[g]], col=group_col[g])
    }
    for(g in ord){
        polygon(d[[g]], col='transparent', lty=2)
    }
    legend('topleft', legend = names(group_col), fill = group_col, bty='n')
    dev.off()
    layout(1)
}

{
    png(filename = '~/Desktop/features_M2macrophage.png', width = 2000, height = 2500, res = 100)
    layout(matrix(1:24, nrow=6,ncol=4, byrow = TRUE))
    sub2 <- sce[, colnames(sub)]
    for(j in 1:nrow(sub2)){
        plot(slingAvgPseudotime(sub), assay(sub2,'exprs')[j,], main = rownames(sub2)[j],
             xlab='Pseudotime', col = alpha(sub$group_col, alpha=.5))
    }
    dev.off()
    layout(1)
}






