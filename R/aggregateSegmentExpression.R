## Aggregates segment specific gene expression
aggregateSegmentExpression <- function(epg, segments, dataset="hsapiens_gene_ensembl", mingps = 20, GRCh = 37, host=NULL) {
    join_id = list(mmusculus_gene_ensembl="mgi_symbol", hsapiens_gene_ensembl="hgnc_symbol")
    LOCCOLS = c("chr", "startpos", "endpos")
    cells = colnames(epg)
    
    ## Map gene IDs to genomic coordinates
    ensembl = NULL
    if(!is.null(host)){
        ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host = host), silent = T)
    }else{
        if (GRCh == 36) {
            ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host = "may2009.archive.ensembl.org"), silent = T)
        } else if (GRCh == 37) {
            ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host = "feb2014.archive.ensembl.org"), silent = T)
            if (class(ensembl) == "try-error") {
                ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host = "grch37.ensembl.org"), silent = T)
            }
        } else if (GRCh == 38) {
            # ensembl = useEnsembl(biomart='ensembl',GRCh=NULL)
            ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host = "may2021.archive.ensembl.org"), silent = T)
            if (class(ensembl) == "try-error") {
                ensembl = try(useMart("ENSEMBL_MART_ENSEMBL", dataset = dataset, host = "grch38.ensembl.org"), silent = T)
            }
        } 
    }
    
    if(class(ensembl) == "try-error"){
        print(paste("BiomaRt service currently not available. Try again later."))
        return()
    }else if (is.null(ensembl)) {
        warning(paste("GRCh", GRCh, "not supported. Only supporting GRCh=36, GRCh=37 or GRCh=38. Abborting."))
        return()
    }
    mart = biomaRt::useDataset(dataset, mart = ensembl)
    genes <- biomaRt::getBM(c("ensembl_gene_id",join_id[[dataset]], "chromosome_name", "start_position", "end_position"), mart = mart)
    colnames(genes) = gsub("end_position", "endpos", gsub("start_position", "startpos", gsub("chromosome_name", "chr", colnames(genes))))
    
    x = .intersect_MatlabV(toupper(rownames(epg)), toupper(genes[,join_id[[dataset]]]))
    
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
