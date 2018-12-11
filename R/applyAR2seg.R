.applyAR2seg <- function(varargin) {
    seg = varargin$seg
    out = varargin$out
    rules = varargin$rules
    PSTR = varargin$PSTR
    
    sI = grep(seg, rules@rhs@itemInfo[, 1])
    print(paste("Applying association rules on segment", seg, "(", which(seg == rownames(out)), "out of", nrow(out), ")"))
    
    cn_status = rep(NA, length(rules))
    for (cn_ in unique(out[!is.na(out)])) {
        eqI = intersect(grep(paste0(" ", cn_), rules@rhs@itemInfo[, 1]), sI)
        if (!isempty(eqI)) {
            cn_status[rules@rhs@data[eqI, ]] = cn_
        }
    }
    
    
    ## Match each left hand rule to set of cells to which they apply
    cmap = matrix(F, length(rules), ncol(out))
    colnames(cmap) = colnames(out)
    for (i in 1:length(rules)) {
        seg1_EQ = rules@lhs@itemInfo[rules@lhs@data[, i], ]
        cI = sapply(seg1_EQ, .testARstatement, out)
        cI = apply(cI, 1, all)
        cI[is.na(cI)] = F
        ## Cells for which one can conclude copy number state of @seg
        cmap[i, names(cI)] = cI
    }
    print(paste("Found", length(rules), "rules applying to", sum(apply(cmap, 2, any)), "cells."))
    
    ## Calculate cummulative support for CNstate == ploidy vs. CNstate != ploidy
    cellsWithRules = c()
    for (cell in colnames(out)) {
        ## Calculate posteriori:= a-priori + support
        for (cn_ in unique(cn_status)) {
            # All rules come to conclusion that seg has CNstate == ploidy
            eq_ii = which(cmap[, cell] & cn_status == cn_)
            if (!isempty(eq_ii)) {
                cellsWithRules = c(cellsWithRules, cell)
            }
            PSTR[as.character(cn_), cell] = sum(quality(rules[eq_ii])$confidence)
        }
    }
    cellsWithRules = unique(cellsWithRules)
    print(paste("Applied rules for", length(cellsWithRules), "out of", ncol(out), "cells for segment", seg))
    
    return(PSTR)
}
