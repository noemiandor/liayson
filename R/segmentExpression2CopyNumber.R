segmentExpression2CopyNumber <- function(eps, gpc, cn, seed = 0, outF = NULL, maxPloidy = 8, nCores = 2, stdOUT = "log.applyAR2seg") {
    if (nCores > round(detectCores() * 0.5)) {
        warning(paste("Attempted to use more than 50% of the available cores."))
    }
    visfit <- function(i) {
        plot(gpc, as.numeric(eps[i, ]), pch = 20, xlab = "Genes per cell", ylab = "Mean expression of segment-genes", cex = 0.15, cex.lab = 1.5, cex.axis = 1.4)
        abline(m[[i]], col = "red")
    }
    HMCOLS = fliplr(brewer.pal(11, "RdBu"))
    MINCELLS = 5
    rcells = colnames(eps)
    # ploidy=round(median(cn))
    
    # ############################################## ## Test gpc for bimodal ploidy distribution ## r=quantile(gpc,c(0,0.95)); p=density(gpc[gpc>=r[1] & gpc<=r[2]], bw = 200, adjust = 1.2, kernel =
    # c('gaussian')); #bw = 'SJ'; ##bw ~ number of genes per cell par(mfrow=c(3,1)); plot(p$x,p$y); second_deriv <- diff(sign(diff(p$y))) if(sum(second_deriv == -2)==2){ ##Currently ignores >2
    # modes ## Shift 2nd peak onto first peak library(mclust) em=Mclust(gpc,G=2); iMax=which.max(em$parameters$mean) shift=max(em$parameters$mean)-min(em$parameters$mean)
    # hist(gpc,200,main=paste('before shift peaks:', paste(round(em$parameters$mean),collapse = ', '))) gpc[em$classification==iMax]=gpc[em$classification==iMax]-shift hist(gpc,200,main='after
    # shift') }
    
    
    #################################################################### Numerical optimization: integer copy number per segment per cell## Approximates # of genes per segment
    gps = apply(eps, 1, function(x) nrow(plyr::count(x)))
    ## Median expression per segment across cells
    x = apply(eps, 1, median)
    ## Step 1: How much of inter-cell differences in expression per segment can be explained by differential gene coverage per cell? Linear regression fit
    m = apply(eps, 1, function(x) lm(as.numeric(x) ~ gpc))
    # par(mfrow=c(2,2)); tmp=lapply(c(1:4), visfit ); ##Visualize first 4 segments' fits Rest of fit: expression/segment ~ gene coverage
    eps = t(sapply(m, function(x) x$residuals))
    ## Shift by median expression per segment
    eps = sweep(eps, 1, x, FUN = "+")
    colnames(eps) = rcells
    ## Step 2: Assumption --> CN segments have variable average copy numbers across cells mlr=mlv(x,method='mfv')$M/ploidy; ##Ratio: maximum likelihood UMI count/ maximum likelihood ploidy
    ## eps=eps/mlr ##Shift to known maximum likelihood ploidy align to population copy number per segment
    eps = sweep(eps, 1, cn/apply(eps[names(cn), ], 1, mean), FUN = "*")
    
    
    ## Define CN states
    bg = as.numeric(eps)
    cnStates = plyr::count(round(bg))
    cnStates = cnStates[cnStates$freq > MINCELLS, ]
    ## Step 3: Cap outliers @TODO: min ploidy should be function parameter
    cnStates = cnStates[cnStates$x >= 0 & cnStates$x <= maxPloidy, ]
    
    
    
    out = eps
    out[T] = 0
    ######################################################### Step 3: Dependency on cell i's read count at locus j### a-priori likelihood: initiate
    PRIR = list()
    PRIR[as.character(cnStates$x)] = rep(list(out), nrow(cnStates))
    for (seg in rownames(eps)) {
        rcl = .cnp(cn[seg], eps[seg, ], 1)
        for (cnState in intersect(names(rcl), as.character(cnStates$x))) {
            PRIR[[cnState]][seg, ] = rcl[[cnState]]
        }
    }
    ## Avoid 0 prior
    for (cn_ in names(PRIR)) {
        PRIR[[cn_]][PRIR[[cn_]] < 9e-224] = 9e-224
    }
    ## Set seed:
    out[T] = NA
    minP = max(sapply(PRIR, max, na.rm = T), na.rm = T)
    pc = 0
    while (pc < seed & minP > -0.01) {
        minP = minP - 0.005
        out = .pickML(out, PRIR, cnStates, minP = minP)
        pc = .progressReport(out, paste("a-priori probabilities @ minRP >=", round(minP, 2) ), verbose = F)
    }
    pc = .progressReport(out, paste("a-priori probabilities @ minRP >=", round(minP, 2) ))
    pc1 = pc
    ## Visualize
    if (!is.null(outF)) {
        pdf(paste0(outF, ".pdf"))
        ## A-priori distribution of CNs for first 4 segments across cells
        par(mfrow = c(4, 2))
        tmp = lapply(c(1:8), function(i) hist(eps[i, ], 100, main = "", col = "cyan", xlim = quantile(cnStates$x, c(0, 1)), xlab = "Copy number", ylab = "# Cells", cex.lab = 1.5, cex.axis = 1.4))
        ## Distributions for each copy number state
        par(mfrow = c(3, 2))
        for (cnState in cnStates$x) {
            try(hist(PRIR[[as.character(cnState)]], 40, xlim = c(0, 1), col = "cyan", main = paste("CN state:", cnState), xlab = "Probability of CN state", ylab = "# of cell-segment pairs"), silent = T)
        }
        ## Seed cells for each segment
        for (seg in rownames(out)) {
            ii = which(!is.na(out[seg, ]))
            hi = hist(eps[seg, ], 100, col = "gray", xlim = c(0, maxPloidy), main = paste0(seg, ": ", cn[seg]), freq = F, ylim = c(0, 1))
            points(eps[seg, ii], rep(0, length(ii)), col = "red", pch = 8)
        }
        ii = sample(ncol(out), min(500, ncol(out)))
        ## Heatmap
        tmp = try(heatmap.2(eps[, ii], cexCol = 0.75, trace = "none", na.color = "black", symm = F, col = HMCOLS, margins = c(13, 10)), silent = T)
        try(heatmap.2(out[, ii], cexCol = 0.75, trace = "none", margins = c(13, 10), na.color = "black", Rowv = tmp$rowDendrogram, Colv = tmp$colDendrogram, symm = F, col = HMCOLS), silent = T)
    }
    
    
    ############################################################ Step 4: Dependency on cell i's read count at other loci:##
    PSTR = PRIR
    ## % cells with >=2 non-NaN
    s_cells = sum(apply(!is.na(out), 2, sum) >= 2)/ncol(out)
    #### Association rule mining### cbind makes sure that output of apply is list
    L = apply(cbind(out, rep(NA, nrow(out))), 2, function(x) paste(names(x[!is.na(x)]), x[!is.na(x)]))
    L = L[sapply(L, length) > 1]
    if (length(L) > 0) {
        L = as(L, "transactions")
        RULES = apriori(L, parameter = list(support = s_cells * 0.1, confidence = 0.75, maxlen = 4), control = list(verbose = F))
        RULES = RULES[apply(RULES@lhs@data, 2, sum) > 0]
        # RULES=sort(RULES,by='lift'); # RULES=RULES[quality(RULES)$lift>1]
        message(paste("Mined", length(RULES), "rules. Now applying rules..."))
        
        # ####Apply association rules###
        message(paste("Using ", nCores, " cores to apply ARs to segments."))
        if (!is.null(stdOUT)) {
            message(paste("Stdout and stderr connections will be redirected to ", stdOUT))
        }
        cl <- makeCluster(nCores, outfile = stdOUT)
        try(on.exit(stopCluster(cl)))
        ## Split expression profile by sample ID
        input = list()
        for (seg in rownames(out)) {
            ## Find association rules where right hand side contains this segment
            sI = grep(seg, RULES@rhs@itemInfo[, 1])
            rules = RULES[apply(RULES@rhs@data[sI, , drop = F], 2, any)]
            if (length(rules) == 0) {
                next
            }
            pstr = t(sapply(PRIR, function(x) x[seg, , drop = F]))
            rownames(pstr) = names(PRIR)
            colnames(pstr) = colnames(PRIR[[1]])
            input[[seg]] = list(seg = seg, out = out, rules = rules, PSTR = pstr)
        }
        ## Distribute jobs
        results = clusterApply(cl, input, .applyAR2seg)
        # Gather results
        for (i in 1:length(results)) {
            pstr = results[[i]]
            pstr = log(pstr + 1)
            pstr = pstr/max(pstr, na.rm=T)
            # pstr = sweep(pstr, 2, apply(pstr, 2, max, na.rm=T), FUN = '/')
            seg = rownames(out)[i]
            for (cn_ in rownames(pstr)) {
                PSTR[[cn_]][seg, colnames(pstr)] = pstr[cn_, ] + PRIR[[cn_]][seg, colnames(pstr)]
            }
        }
    } else {
        warning("No association rule mining was performed.")
        PSTR = PRIR
    }
    # ###############
    
    
    ## Pick maximum likelihood among PSTR-entries:
    out[T] = NA
    out = .pickML(out, PSTR, cnStates)
    pc = .progressReport(out, "posteriori probabilities")
    
    
    ## Visualize phylogeny ###
    if (!is.null(outF)) {
        ## cells
        cI = which(apply(!is.na(out), 2, sum) > 0.5 * nrow(out))
        ## segments
        gI = which(apply(!is.na(out), 1, sum) > 0.5 * ncol(out))
        if (length(cI) > 1000) {
            cI = sample(cI, 1000, replace = F)
        }
        tmp = out[gI, cI]
        hm = try(heatmap.2(tmp, trace = "none", symm = F, hclustfun = function(x) hclust(x, method = "ward.D2"), distfun = function(x) dist(x, method = "euclidean"), na.color = "black", main = paste0(nrow(tmp), 
                                                                                                                                                                                                        " segments x ", length(cI), " cells"), col = HMCOLS))
        dev.off()
    }
    
    return(out)
}

######################################################################################################### 

## Pick maximum likelihood among Copy number states:
.pickML <- function(out, prob_CNstate, cnStates, minP = 0) {
    # prob_CNstate - list with fields corresponding to copy number state. Each field contains a cell-by-segment matrix of likelihoods minP - will only assign copy number to cell-segment pairs of
    # >=minP certainty
    for (x in cnStates$x) {
        p = prob_CNstate[[as.character(x)]]
        ## Index of cell segment pairs to which we will assign copy number state x
        cnI = which(is.na(out) & p >= minP)
        for (y in setdiff(cnStates$x, x)) {
            ## Clear-cut winner: use >=
            cnI = setdiff(cnI, which(prob_CNstate[[as.character(y)]] >= p))
        }
        out[cnI] = x
    }
    return(out)
}


# ##Bayes: a-priori probability of copy number state x given measurement x_ .cnp<-function(x_, x, l, ploidy){ px = sum( round(x_)==x ) / length(x_) SD= px * 1; #*l ##Penalize deprature from
# local CN state ##and reward segment precision l # SD= SD / (1+abs(x-ploidy)); ##Penalize deprature from global ploidy state px_Gx = dnorm(x_,mean = x,sd = SD) px_Gx = px_Gx/ dnorm(x,mean =
# x,sd = SD) return(px_Gx) }

## Assume only 2 CN states per segment
.cnp <- function(cns, x, pctCells = 1) {
    ## cns - average CN for this segment x - expression of each cell in this segment pctCells - % cells to be assigned a CN state
    
    # normal state := most common state
    ns = round(cns)
    ## mutated state := adjacent state
    ms = ns + sign(cns - ns)
    
    ## Kernel fit in case CN states are not adjacent:
    p = density(x, bw = "SJ", adjust = 0.5, kernel = c("gaussian"))
    second_deriv <- diff(sign(diff(p$y)))
    ## Maxima
    ii = which(second_deriv == -2) + 1
    ii = ii[p$y[ii] >= 0.1]
    ii = ii[sort(p$y[ii], index.return = T, decreasing = T)$ix]
    ## Reassign normal and mutated state only if kernel fit was successfull
    if (length(ii) > 1 && round(p$x[ii[1]]) != round(p$x[ii[2]])) {
        ns = round(p$x[ii[1]])
        ms = round(p$x[ii[2]])
    }
    ## clust = mixtools::normalmixEM(x, k = 2) ns = round(clust$mu[which.max(clust$lambda)] ); ms = round(clust$mu[which.min(clust$lambda)] ); plot(p$x,p$y, main=paste(names(cns),ns,',',ms))
    
    ## fraction mutated cells
    sp = (cns - ns)/(ms - ns)
    ## prevent NA if only one CN state exists for this segment
    sp = max(1/length(x), sp, na.rm = T)
    ms_cells = max(1, round(sp * length(x) * pctCells))
    ns_cells = max(1, round((1 - sp) * length(x) * pctCells))
    
    ## Assign ms and ns states to ms_cells and ns_cells respectively:
    ms_idx = sort(abs(x - ms), index.return = T)$ix[1:ms_cells]
    ns_idx = sort(abs(x - ns), index.return = T)$ix[1:ns_cells]
    
    ## Calculate probabilities
    PRIR = list()
    px_Gx = dnorm(x, mean = ms, sd = 2 * sd(x[ms_idx]))
    PRIR[[as.character(ms)]] = px_Gx/dnorm(ms, mean = ms, sd = 2 * sd(x[ms_idx]))
    
    px_Gx = dnorm(x, mean = ns, sd = 2 * sd(x[ns_idx]))
    PRIR[[as.character(ns)]] = px_Gx/dnorm(ns, mean = ns, sd = 2 * sd(x[ns_idx]))
    
    return(PRIR)
}

.progressReport <- function(out, step, verbose = T) {
    v = sum(!is.na(as.numeric(out)))/length(out)
    if (verbose) {
        message(paste0("After ", step, ": ", round(100 * v, 2), "% copy numbers inferred"))
    }
    return(v)
}




# .getLocusEntropy<-function(out){ tmp=cbind(out,repmat(min(out,na.rm=T):max(out,na.rm=T),nrow(out),1));##Ensure every CN state is represented for every segment l=apply(tmp,1,function(x)
# plyr::count(round(x[!is.na(x)]))$freq) ent=apply(l,2,function(x) vegan::diversity(x/sum(x), index = 'simpson')) return(ent) }
