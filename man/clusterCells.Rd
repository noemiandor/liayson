\name{clusterCells}
\alias{clusterCells}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Grouping cells into clones.
}
\description{
Clusters cells according to their copy number profile.
}
\usage{
clusterCells(cnps, k=NA, h=NA, consolidate=1, weights=NULL, minSegLength=1E6, 
   colCol=NULL, rowCol=NULL, main="", chrOrder=NULL, cexRow = 0.07, trace ="none")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cnps}{
    Segment-by-cell matrix of copy number states (output of \code{segmentExpression2CopyNumber}).
  }
  \item{consolidate}{
    Whether or not to homogenize the profiles of cells with same clone membership (consolidate=0 <=> all cells within a clone will be identical; consolidate=1 <=> intra-clone homogeneity not enforced).
  }
  \item{minSegLength}{
    Minimum number of base pairs below which a segment is to be excluded when defining clones.
  }
}
\details{
Let CNF be the matrix of copy number states per non-private segment per cell, with entries \eqn{(i, j)} pointing to the copy number state of cell  \eqn{i} at locus  \eqn{j}. 
Pairwise distances between cells are calculated in Euclidean space of their segmental copy number profiles (rows in CNF), weighted by segment length. Hierarchical clustering with the ward.D2 agglomeration method is used to build a phylogenetic tree of the cells. A subtree is defined as a clone if the maximum distance between its cell members is less than a user-defined fraction of the genome (\eqn{h}.\cr
Alternatively, if \eqn{k} is set, the tree is cut to obtain \eqn{k} clones.\cr
If neither \eqn{h} nor \eqn{k} are set, Akaike information criterion is used to decide on anywhere between 1 and 30 clones.
}
\value{List with three fields:
	  \item{cnps}{
	        Segment-by-cell matrix of copy number states.
    }         
	  \item{sps}{
          The clone membership of each cell (that is, columns in cnps).
    }
    \item{tree}{
          An object of class hclust which describes the tree produced by the clustering process.
    }
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
\examples{
data(cnps)
outc = clusterCells(cnps, h=0.15)
}
