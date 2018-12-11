.testARstatement <- function(seg_eqp, out) {
    seg_eqp = strsplit(seg_eqp, " ")[[1]]
    seg = seg_eqp[1]
    EQp = as.numeric(seg_eqp[2])
    return(out[seg, ] == EQp)
}
