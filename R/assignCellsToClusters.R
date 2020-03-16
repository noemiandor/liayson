assignCellsToClusters <- function(outc, xps, similarity = T) {
    fr = plyr::count(outc$sps)
    fr = fr[sort(fr$freq, decreasing = T, index.return = T)$ix, ]
    fr$x = as.character(fr$x)
    cnvs = t(.grpstats(t(outc$cnps), outc$sps, "mean")$mean)
    #devFromInt = apply(abs(round(cnvs) - cnvs), 1, median)
    ### Don't use segments with unstable CN states within given clone
    seg = rownames(cnvs); #[devFromInt < maxDevFromInt]
    
    if (similarity) {
        method = "correlation"
        d = as.matrix(proxy::simil(t(cbind(cnvs[seg, ], xps[seg, ])), method = method))
        d = d[colnames(cnvs), colnames(xps)]
        sps = as.numeric(rownames(d)[apply(d, 2, which.max)])
        # diff=apply(d,2,max) - apply(d,2, function(x) sort(x,decreasing = T)[2]) sps[diff<1E-4]=NA
    } else {
        method = "Euclidean"
        d = as.matrix(dist(t(cbind(cnvs[seg, ], xps[seg, ])), method = method))
        d = d[colnames(cnvs), colnames(xps)]
        sps = as.numeric(rownames(d)[apply(d, 2, which.min)])
    }
    names(sps) = colnames(d)
    
    outc$cnps = cbind(outc$cnps, xps[rownames(cnvs), !is.na(sps)])
    outc$sps = c(outc$sps, sps[!is.na(sps)])
    
    x = sort(outc$sps, index.return = T, decreasing = T)
    outc$cnps = outc$cnps[, x$ix]
    outc$sps = x$x
    
    
    ## Visualize Cap
    minV = quantile(as.numeric(d), 0.025, na.rm = T)
    d[d < minV] = minV
    # Get distance matrix in shape
    d = sweep(d, MARGIN = 2, apply(d, 2, max), FUN = "/")
    
    HMCOLS = fliplr(brewer.pal(11, "RdBu"))
    col = rainbow(max(outc$sps))
    tmp = outc$sps[names(outc$sps) %in% colnames(d)]
    hm = try(heatmap.2(t(d[, names(tmp)]), margins = c(13, 6), cexCol = 0.85, col = HMCOLS, colRow = col[tmp], Rowv = NULL, Colv = NULL, trace = "n", symm = F, hclustfun = function(x) hclust(x, 
        method = "ward.D")))
    
    ## Visualize CN differences between clones
    ix = sort(cnvs[seg, fr$x[1]], index.return = T)$ix
    sz = 1 + 0.25 * nrow(fr)
    par(mfrow = c(2, 1))
    plot(0.05 + cnvs[seg[ix], fr$x[1]], pch = 20, cex = sz, col = 1, xlab = "locus", ylab = "copy number", log = "y")
    for (i in 2:nrow(fr)) {
        points(0.05 + cnvs[seg[ix], fr$x[i]], pch = 20, cex = sz - (i * 0.25), col = i, log = "y")
    }
    legend("topleft", paste("Clone", as.character(fr$x)), fill = 1:nrow(fr))
    ## Compare G0G1 freq to cycling freq
    fr2 = plyr::count(outc$sps[names(outc$sps) %in% colnames(xps)])
    rownames(fr2) = as.character(fr2$x)
    y = fr2[fr$x, ]$freq
    y[is.na(y)] = 0
    plot(fr$freq, y, xlab = "G0G1", ylab = "cycling", pch = 20, cex = 2)
    
    return(outc)
}
