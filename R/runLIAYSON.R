runLIAYSON <- function(X, S, sName, mingps = 20, GRCh = 37, h = 0.2, minSegLength = 1e+06, outD = NULL) {
    message(paste("Running LIAYSON on", sName, "..."))
    ## Genes expressed per cell:
    gpc = apply(X > 0, 2, sum)
    names(gpc) = colnames(X)
    
    ## Expression per segment per cell:
    cps <- matrix(NA, nrow(S), ncol(X))
    rownames(cps) = rownames(S)
    colnames(cps) = colnames(X)
    eps <- gps <- cps
    tmp = aggregateSegmentExpression(X, S, GRCh = GRCh, mingps = 1)
    if(is.null(tmp)){
        print("BiomaRt Web service for annotation of gene locations is not available. Aborting.")
        return()
    }
    eps[rownames(tmp$eps), ] = tmp$eps
    gps[rownames(tmp$gps), ] = tmp$gps
    
    ## Copy number per segment per cell:
    iK = rownames(tmp$gps)[apply(tmp$gps, 1, median) >= mingps]
    cps[iK, ] = segmentExpression2CopyNumber(tmp$eps[iK, ], gpc, cn = S[iK, "CN_Estimate"], seed = 0, nCores = ceil(length(iK)/2), stdOUT = paste0("log.", sName))
    
    ## Cluster cells
    cps = cps[apply(!is.na(cps), 1, any), ]
    outc = clusterCells(cps, h = h, minSegLength = minSegLength)
    
    ## Save output
    if (!is.null(outD)) {
        saveClusteredCells(outc, outD, sName)
    }
    return(outc)
}
