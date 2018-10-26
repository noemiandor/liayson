\name{saveClusteredCells}
\alias{saveClusteredCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Saving clones to user-defined output.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
saveClusteredCells(outc, outD, sName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{outc}{
    Output of \code{clusterCells} or \code{assignCellsToClusters}: list containing segment-by-cell matrix, clone membership of each cell and phylogenetic tree.
}
}
\details{
Saves:\cr
1. The clone membership of each cell (*.spstats)\cr
2. The segment-by-cell matrix of copy number states (*.sc.cbs)\cr
3. The consensus copy number profile of each detected clone, calculated as the average profile of cells that are members of the respective clone (*.sps.cbs)\cr
4. The phylogeny (*.tree).
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

