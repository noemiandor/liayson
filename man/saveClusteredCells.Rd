\name{saveClusteredCells}
\alias{saveClusteredCells}
\title{
Saving clones to user-defined output.
}
\description{
Writes clone statistics, dendrogram and clone-specific mutation profiles.
}
\usage{
saveClusteredCells(outc, outD, sName)
}
\arguments{
  \item{outc}{
    Output of \code{clusterCells} or \code{assignCellsToClusters}: list containing segment-by-cell matrix, clone membership of each cell and the underlying dendrogram.
}
  \item{outD}{
    The output directory.
}
  \item{sName}{
    Prefix for the output files (typically the sample name).
}
}
\details{
Writes each of the following aspects of a sample's clonal composition into an output file:\cr
1. The clone membership of each cell (*.spstats)\cr
2. The segment-by-cell matrix of copy number states (*.sc.cbs)\cr
3. The consensus copy number profile of each detected clone, calculated as the average profile of cells that are members of the respective clone (*.sps.cbs)\cr
4. The cell dendrogram (*.tree).
}

\author{
Noemi Andor
}


