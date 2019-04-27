\name{clusterCells}
\alias{clusterCells}
\title{
Grouping cells into clones.
}
\description{
Clusters cells according to their copy number profile.
}
\usage{
clusterCells(cnps, k=NA, h=NA, weights=NULL, minSegLength=1E6, 
                       chrOrder=NULL,  HFUN="ward.D2",...)
}
\arguments{
  \item{cnps}{
    Segment-by-cell matrix of copy number states (output of \code{segmentExpression2CopyNumber}).
  }
  \item{k}{
    Desired number of clusters (see also \code{\link{cutree}}).
  }
  \item{h}{
    Threshold used to define clones from hierarchical clustering result. A subtree is defined as a clone if the maximum distance between its cell members is less than 100*\eqn{h}\% of the genome.
  }
  \item{weights}{
    Vector of weights assigning differential importance to segments (typically calculated based on segment lengths).
  }
  \item{minSegLength}{
    Minimum number of base pairs below which a segment is to be excluded when defining clones.
  }
  \item{chrOrder}{
    Specifies order in which chromosomes should be plotted.
  }
  \item{HFUN}{
    Agglomeration method used to compute the hierarchical clustering (see also \code{\link{hclust}}).
  }
  \item{...}{additional arguments passed on to \code{\link{heatmap.2}} }
}
\details{
Let CNF be the matrix of copy number states per non-private segment per cell, with entries \eqn{(i, j)} pointing to the copy number state of cell  \eqn{j} at locus  \eqn{i}. 
Pairwise distances between cells are calculated in Hamming space of their segmental copy number profiles (rows in CNF), weighted by segment length. Hierarchical clustering is used to build a tree of the cells from the distance matrix. A subtree is defined as a clone if the maximum distance between its cell members is less than a user-defined fraction of the genome (\eqn{h}).\cr
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
          An object of class hclust.
    }
}
\references{
Andor, N.*, Lau, B.*, Catalanotti, C., Kumar, V., Sathe, A., Belhocine, K., Wheeler, T., et al. (2018) Joint single cell DNA-Seq and RNA-Seq of gastric cancer reveals subclonal signatures of genomic instability and gene expression. doi: https://doi.org/10.1101/445932
}
\author{
Noemi Andor
}
\examples{
data(cnps)
set.seed(3)
rcells = sample(colnames(cnps), 120)
outc = clusterCells(cnps[apply(cnps, 1, var)>0, rcells])
}
