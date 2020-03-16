## Aggregates segment specific gene expression
aggregateSegmentExpression <- function(epg, segments, mingps = 20, GRCh = 37) {
    LOCCOLS = c("chr", "startpos", "endpos")
    cells = colnames(epg)
    
    ## Map gene IDs to genomic coordinates
    ensembl = NULL
    if (GRCh == 36) {
        ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "may2009.archive.ensembl.org"))
    } else if (GRCh == 37) {
        ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "feb2014.archive.ensembl.org"), silent = T)
        if (class(ensembl) == "try-error") {
            ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch37.ensembl.org"))
        }
    } else if (GRCh == 38) {
        # ensembl = useEnsembl(biomart='ensembl',GRCh=NULL)
        ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "dec2016.archive.ensembl.org"), silent = T)
        if (class(ensembl) == "try-error") {
            ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "grch38.ensembl.org"))
        }
    } 
    if (is.null(ensembl) || class(ensembl) == "try-error") {
        warning(paste("GRCh", GRCh, "not supported. Only supporting GRCh=36, GRCh=37 or GRCh=38. Abborting."))
        return()
    }
    mart = biomaRt::useDataset("hsapiens_gene_ensembl", mart = ensembl)
    genes <- biomaRt::getBM(c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"), mart = mart)
    colnames(genes) = gsub("end_position", "endpos", gsub("start_position", "startpos", gsub("chromosome_name", "chr", colnames(genes))))
    
    x = .intersect_MatlabV(rownames(epg), genes[, "hgnc_symbol"])
    epg = epg[x$ia, , drop = F]
    genes = genes[x$ib, , drop = F]
    epg = cbind(genes[, LOCCOLS], epg)
    epg$chr = as.numeric(gsub("Y", 24, gsub("X", 23, epg$chr)))
    epg = as.matrix(epg[!is.na(epg$chr), ])
    rownames(epg) = paste0(epg[, "chr"], ":", epg[, "startpos"], "-", epg[, "endpos"])
    
    ## Aggregate gene expression per each segment
    map = .assignCBSToMutation(epg, segments)
    ## Genes that could be mapped to a segment
    iC = which(!is.na(map[, "quantityID"]))
    ## Segments
    s = rownames(segments)[map[iC, "quantityID"]]
    ## Expression per copy number segment
    eps = .grpstats(epg[iC, cells], s, "mean")$mean
    
    ## Keep only 'valid' segments, containing on average >=mingps expressed genes Count genes expressed per copy number segment
    gps = .grpstats(epg[iC, cells] > 0, s, "sum")$sum
    iK = which(apply(gps, 1, median) >= mingps)
    eps = eps[iK, ]
    gps = gps[iK, ]
    message(paste(">=", mingps, "genes expressed in", length(iK), "segments for", length(cells), "cells"))
    
    return(list(eps = eps, gps = gps))
}



.intersect_MatlabV <- function(a, b) {
    x = intersect(a, b)
    ia = match(x, a)
    ib = match(x, b)
    return(list(a = x, ia = ia, ib = ib))
}
