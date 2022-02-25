saveClusteredCells <- function(outc, expression, ccMembership, sName, outD="~/Downloads", NUMRES=getNumRes()) {
    if(!"G1Malignant" %in% names(ccMembership)){
        warning("List 'ccMembership' must include entry 'G1Malignant', containing the IDs of G0/G1 tumor cells. No data was saved.", immediate. = T)
        return()
    }
    
    ## G1 clones
    idx=match(ccMembership$G1Malignant,names(outc$sps))
    outc_G1=list(cnps=outc$cnps[,idx], sps=outc$sps[idx], tree=outc$tree)
    
    ## Save G1 clone frequencies only
    write.table(outc_G1$cnps, file = paste0(outD, filesep, sName, ".sc.cbs"), sep = "\t", quote = F)
    ## Save results: spstats
    spstats = plyr::count(outc_G1$sps)
    spstats$freq = round(spstats$freq/sum(spstats$freq), getNumRes())
    colnames(spstats) = c("ID", "Mean Weighted")
    write.table(spstats, file = paste0(outD, filesep, sName, ".spstats"), sep = "\t", quote = F, row.names = F)
    ## Save results: sps profiles
    cnvs = t(.grpstats(t(outc_G1$cnps), outc_G1$sps, "mean")$mean)
    cnvs = cnvs[, rownames(spstats), drop = F]
    colnames(cnvs) = paste0("Clone", spstats$ID, "_", round(spstats$`Mean Weighted`, getNumRes()))
    cnvs = as.data.frame(cbind(.parseLOCUS(rownames(cnvs)), cnvs))
    cnvs$CN_Estimate = apply(outc_G1$cnps[rownames(cnvs), ], 1, mean, na.rm = T)
    write.table(cbind(rownames(cnvs), cnvs), file = paste0(outD, filesep, sName, ".sps.cbs"), sep = "\t", quote = F, row.names = F, col.names = c("LOCUS", colnames(cnvs)))
    write.tree(as.phylo(outc_G1$tree), file = paste0(outD, filesep, sName, ".tree"))
    
    ##Save clone membership for each Clone, including S and G2M cells
    for(id in spstats$ID){
        sz=spstats$`Mean Weighted`[spstats$ID==id]
        iC=names(outc$sps)[outc$sps==id];         
        dm=rbind(outc$cnps[,outc$sps==id], as.matrix(expression[,iC])); ####copy number followed by Gene expresssion
        colnames(dm)=iC
        ##@TODO: identical profiles will not be saved multiple times in DB CLONEID -- fix by counting!
        cols=  rep(1/ncol(dm),ncol(dm)) + .addNoise(ncol(dm))
        alias= gsub("_","",colnames(dm));
        state= rep(""        ,ncol(dm));    
        for(thisstate in names(ccMembership)){
            state[ colnames(dm) %in% ccMembership[[thisstate]]]=thisstate
        }
        write.table(cbind(rownames(dm),dm),row.names = F, col.names = c("LOCUS", paste0("Clone_",cols,"_",state,"_",alias)),
                    sep="\t",quote=F, file=paste0(outD,filesep,sName,".",round(sz,NUMRES),".sps.cbs")); 
    }
}

.addNoise <- function(n) {
    sample(1E5,n)/1E8; ##@TODO: should be a function of min value entry rather than array length
}