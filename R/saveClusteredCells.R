saveClusteredCells <- function(outc, outD, sName) {
    write.table(outc$cnps, file = paste0(outD, filesep, sName, ".sc.cbs"), sep = "\t", quote = F)
    ## Save results: spstats
    spstats = plyr::count(outc$sps)
    spstats$freq = round(spstats$freq/sum(spstats$freq), getNumRes())
    colnames(spstats) = c("ID", "Mean Weighted")
    write.table(spstats, file = paste0(outD, filesep, sName, ".spstats"), sep = "\t", quote = F, row.names = F)
    ## Save results: sps profiles
    cnvs = t(.grpstats(t(outc$cnps), outc$sps, "mean")$mean)
    cnvs = cnvs[, rownames(spstats), drop = F]
    colnames(cnvs) = paste0("Clone", spstats$ID, "_", round(spstats$`Mean Weighted`, getNumRes()))
    cnvs = as.data.frame(cbind(.parseLOCUS(rownames(cnvs)), cnvs))
    cnvs$CN_Estimate = apply(outc$cnps[rownames(cnvs), ], 1, mean, na.rm = T)
    write.table(cbind(rownames(cnvs), cnvs), file = paste0(outD, filesep, sName, ".sps.cbs"), sep = "\t", quote = F, row.names = F, col.names = c("LOCUS", colnames(cnvs)))
    write.tree(as.phylo(outc$tree), file = paste0(outD, filesep, sName, ".tree"))
}
