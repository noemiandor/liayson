clusterCells <- function(cnps, k = NA, h = NA, weights = NULL, minSegLength = 1e+06, chrOrder = NULL, HFUN = "ward.D2", ...) {
    ## CN state colorcode for heatmaps
    HMCOLS = fliplr(brewer.pal(11, "RdBu"))
    MAXK = 30
    cnps = as.matrix(cnps)
    
    ######################### Cells & loci to keep ##
    cI = which(apply(!is.na(cnps), 2, sum) > 0.5 * nrow(cnps))
    ## cells:
    cnps = cnps[, cI]
    gI = which(apply(!is.na(cnps), 1, sum) > 0.75 * ncol(cnps))
    ## segments:
    cnps = cnps[gI, ]
    ## Impute missing values
    cnps = t(e1071::impute(t(cnps), "mean"))
    
    colCol = repmat("black", ncol(cnps), 1)
    names(colCol) = colnames(cnps)
    rowCol = repmat("black", nrow(cnps), 1)
    names(rowCol) = rownames(cnps)
    
    ############################################### Segment weights --> pairwise cell distances##
    loci = as.data.frame(.parseLOCUS(rownames(cnps)))
    rownames(loci) = rownames(cnps)
    loci = loci[loci$seglength >= minSegLength, ]
    cnps = cnps[rownames(loci), ]
    w = loci[, "seglength"]/1e+06
    if (!is.null(weights)) {
        w = as.numeric(weights[rownames(cnps)])
    }
    # d=distances::distances(t( cnps ), weights =w )
    hd = .weightedHamming(t(cnps), w)
    d = hd$dist
    # # ##Exclude singletons and recalculate distance # Z=hclust(d, method = HFUN); # # Z$height=round(Z$height); # TC=cutree(Z,h=0.05); # singletons=plyr::count(TC); # singletons=names(TC)[TC %in%
    # singletons$x[singletons$freq<3]] # message(paste('Excluding',length(singletons),'singletons out of ',ncol(cnps),'cells')) # cnps=cnps[,!colnames(cnps) %in% singletons] # hd=.weightedHamming(t(
    # cnps ), w); d=hd$dist;
    
    ################################ Build maximum parsimony tree## trinary=hd$trinary tree <- bionjs(d) o=apply(d,1,median) tree=root(tree, outgroup = names(o[o==max(o)]),resolve.root =T) ## root
    ################################ tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE); ## make ultrametric # treeRatchet <- pratchet(trinary, start=tree, maxit=100, k=5, trace=0) # treeRatchet <- acctran(treeRatchet, trinary)
    ################################ # d=as.dist(cophenetic(tree))
    
    ############################ Find number of clones:k ##
    if (is.na(k) && is.na(h)) {
        message("Neither k nor h is set.")
        # message('Using indices from NbClust-package to decide number of clusters') # idxs=c('kl', 'ch', 'hartigan','cindex', 'db', 'silhouette', 'ratkowsky', 'ball', 'ptbiserial', 'gap', 'frey',
        # 'mcclain', 'dunn','sdindex', 'sdbw', 'gamma', 'gplus', 'tau') idxs=c('frey', 'mcclain', 'cindex', 'sihouette','dunn') ks=matrix(NA,length(idxs),1); colnames(ks)='K'; rownames(ks)=idxs for(idx
        # in setdiff(idxs,c('gamma', 'gplus', 'tau'))) { # nbc=try(NbClust::NbClust(data = t(cnps), distance = 'euclidean', min.nc = 1, max.nc = min(MAXK+1,nrow(cnps)-1), method = HFUN, index = idx),silent=T)
        # nbc=try(NbClust::NbClust( diss=as.dist(d), distance = NULL, min.nc = 1, max.nc = min(MAXK+1,nrow(cnps)-1), method = HFUN, index = idx),silent=T) if(class(nbc)!='try-error'){
        # ks[idx,'K']=nbc$Best.nc['Number_clusters'] } } message(ks) ks=ks[ks<MAXK+1,]; ##Don't trust results at max of range k=round(mean(ks[is.finite(ks)]));
        # #round(modeest::mlv(ks[is.finite(ks)],method='mfv')$M)
        
        # message('Using Calinski-Harabasz criterion to decide number of clusters') k=fpc::pamk(as.dist(d), krange=1:MAXK, criterion = 'ch')$nc # message('Using average silhouette criterion to decide
        # number of clusters') # k=fpc::pamk(as.dist(d), krange=1:MAXK, criterion='asw')$nc
        
        message("Using Akaike information criterion to decide number of clusters...")
        ks = c()
        ## Repeat for robustness
        for (i in 1:min(25, ncol(cnps)-1) ) {
            fit <- sapply(1:MAXK, function(x) .kmeansAIC(kmeans(t(cnps), centers = x))$AIK)
            ks[i] = which.min(fit)
        }
        ks = plyr::count(ks)
        k = ks$x[which.max(ks$freq)]
        message(paste("Looking for", k, "clusters..."))
    }
    
    ############################# Cluster into k clones ###
    d = as.matrix(d)
    rownames(d) <- colnames(d) <- colnames(cnps)
    d = as.dist(d)
    Z = hclust(d, method = HFUN)
    # Z=reorder(tree, order = 'cladewise')
    if (!is.na(k)) {
        TC = cutree(Z, k = k)
    } else if (!is.na(h)) {
        TC = .cutreeAtCnvAbudnance(Z, t(cnps), h, loci)
    }
    Z = as.phylo(Z)
    
    # ############################################################# ##For each segment, assign ML CN-state to each clone member## if(consolidate<1){ out_ml=cnps; ##will contain more homogeneous
    # clones (consolidate=0 <=> strictly homogeneous, consolidate=1 <=> not enforcing homogeneity) for(cl in unique(TC)){ ii=find(TC==cl) for(s in rownames(cnps)){ fr=plyr::count(cnps[s,ii]);
    # fr=fr[!is.na(fr$x),,drop=F]; fr$freq=fr$freq/sum(fr$freq) if(max(fr$freq)>consolidate){ out_ml[s,ii]=fr$x[which.max(fr$freq)]; ##A cell's ML CN-state assigned based on clone-membership } } }
    # if(sum(sum(out_ml!=cnps,na.rm=T))>0*sum(sum(!is.na(cnps),na.rm=T))){ ##Repeat clustering if ML assignment effectively changed at least one cell's CN state
    # outL=clusterCells(out_ml,k=k,h=h,weights=weights,chrOrder=chrOrder,consolidate = 1) return(outL) } }
    
    ############## Phylogeny ###
    ucol = rainbow(length(unique(TC)))
    colI = ucol[TC]
    names(colI) = names(TC)
    tmp = t(cnps)
    Colv = T
    if (!is.null(chrOrder)) {
        Colv = NULL
        chrOrder = chrOrder[chrOrder %in% colnames(tmp)]
        ii = sort(match(colnames(tmp), chrOrder), index.return = T)$ix
        tmp = tmp[, ii]
    }
    ## Visualize tree
    ii = fliplr(Z$tip.label[Z$edge[Z$edge[, 2] <= length(Z$tip.label), 2]])
    plot(Z, show.tip.label = T, tip.color = colI, cex = 0.2)
    hm = try(heatmap.2(tmp[ii, ], dendrogram = "row", margins = c(13, 6), cexCol = 0.85, Rowv = NULL, Colv = Colv, col = HMCOLS, colRow = colCol[ii], colCol = rowCol[colnames(tmp)], RowSideColors = colI[ii], 
        symm = F, trace = "none"))
    # , hclustfun = function(x) hclust(x,method =HFUN), distfun = function(x) dist(x,method ='euclidean')))
    
    colnames(cnps) = paste0("SP", TC, "_", colnames(cnps))
    return(list(cnps = cnps, sps = TC, tree = Z))
}


.weightedHamming = function(cellXseg, w) {
    # o=cellXseg
    w = w - min(w)
    w = 2000 * w/max(w)
    w = round(w + 1)
    ## Rep according to weight
    o = lapply(1:length(w), function(i) repmat(cellXseg[, i, drop = F], 1, w[i]))
    o = do.call(cbind, o)
    
    ## Calc hamming distance
    u <- o
    # ploidy=median(round(o)) u[round(o)>ploidy]='A' u[round(o)<ploidy]='D' u[round(o)==ploidy]='N'
    trinary = phyDat(u, type = "USER", levels = unique(as.character(u)))
    d <- dist.hamming(trinary)
    return(list(dist = d, trinary = trinary))
}


## A subtree was defined as a clone if the maximum distance between its cell members was less than XX% of the genome Dendextend version
.cutreeAtCnvAbudnance <- function(Z, out, h, loci) {
    ## Build tree
    Z = as.phylo(Z)
    genome = sum(loci$seglength)
    subtrees <- subtrees(Z)
    subtrees = sapply(subtrees, function(x) x$tip.label)
    clones = list()
    while (any(sapply(subtrees, length) >= 1)) {
        ## Sort in ascending order of length for efficacy
        subtrees = subtrees[sort(sapply(subtrees, length), index.return = T)$ix]
        subtrees = subtrees[sapply(subtrees, length) >= 1]
        clone = subtrees[[1]]
        if (any(sapply(subtrees, length) > 1)) {
            tree0 = subtrees[[1]]
            ## All supertrees of current tree
            for (tree in subtrees[sapply(subtrees, function(x) length(intersect(tree0, x)) == length(tree0))]) {
                # la=proxy::dist(out[tree,],method = function(x,y) sum(x!=y)/length(x)) ii=which( as.matrix(la) == max(la), arr.ind = TRUE)[1,]
                la = proxy::dist(out[tree, ], method = function(x, y) sum(loci$seglength[which(x != y)])/genome)
                if (max(la, na.rm = T) > h) {
                  break
                }
                clone = tree
            }
        }
        clones[[length(clones) + 1]] = clone
        ## exclude this clone's ancestors from further processing to prevent scambled clone membership
        subtrees = subtrees[sapply(subtrees, function(x) length(intersect(clone, x)) == 0)]
        ## exclude this clone's members from further processing subtrees=sapply(subtrees, function(x) setdiff(x,clone));
    }
    ## Add singletons
    clones = c(clones, subtrees[sapply(subtrees, length) == 1])
    ## Clean up
    clones = clones[sapply(clones, length) > 0]
    ## Add singletons
    clones = c(clones, as.list(setdiff(Z$tip.label, unlist(clones))))
    if (sum(sapply(clones, length)) != length(Z$tip.label)) {
        warning("Cells lost while defining clones!", immediate. = T)
    }
    TC = rep(0, length(Z$tip.label))
    names(TC) = Z$tip.label
    for (i in 1:length(clones)) {
        TC[clones[[i]]] = i
    }
    return(TC)
}

.kmeansAIC = function(fit) {
    m = ncol(fit$centers)
    n = length(fit$cluster)
    k = nrow(fit$centers)
    D = fit$tot.withinss
    return(list(AIK = D + 2 * m * k, BIC = D + log(n) * m * k))
}


## A subtree was defined as a clone if the maximum distance between its cell members was less than XX% of the genome Phylobase version .cutreeAtCnvAbudnance<-function(Z, out, h, loci){
## genome=sum(loci$seglength) clones=list() phy=phylo4(as.phylo(Z)) tips=Z$labels while(!isempty(tips)){ d=0 tree0_1=list(tips[1], tips[1]) while(d<h){ A=ancestor(phy, tree0_1[[2]])
## tree=names(descendants(phy,A,type='tips')) la=proxy::dist(out[tree,],method = function(x,y) sum(x!=y)) ii=which( as.matrix(la) == max(la), arr.ind = TRUE)[1,] la=proxy::dist(out[tree[ii],],
## method=function(x,y) sum(loci$seglength[which(x!=y)])/genome ) d=max(la) tree0_1=list(tree0_1[[2]],A) if(is.na(ancestor(phy, A))){ tree0_1[[1]]=tree0_1[[2]] break } }
## clone=names(descendants(phy,tree0_1[[1]],type='tips')) clones[[length(clones)+1]]=clone tips=setdiff(tips,clone) } if(sum(sapply(clones,length))!=length(Z$labels)){ warning('Cells lost while
## defining clones!', immediate. = T) } TC=rep(0, length(Z$labels)); names(TC)=Z$labels; for(i in 1:length(clones)) { TC[clones[[i]]]=i; } return(TC) }
