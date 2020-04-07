\name{assignCellsToClusters}
\alias{assignCellsToClusters}
\title{Assigns cells to previously defined clones.}
\description{Cells that have not been used to define clones (such as cycling or apoptotic cells) can retrospectively be assigned a clone membership.}
\usage{assignCellsToClusters(outc, xps, similarity=T)}
\arguments{
  \item{outc}{
	  List containing segment-by-cell matrix and clone membership of each cell. See \code{clusterCells}).
}
  \item{xps}{
	  Segment-by-cell matrix of expression- or copy number states. Columns represent new cells to be assigned to existing clones.
}
  \item{similarity}{
	  Whether to use similarity (similarity=T) or distance (similarity=F), when comparing cells to existing clones. Default similarity metric is "correlation. Default distance metric is "Euclidean". 
}
}
\details{
Let S := \{ \eqn{S_1, S_2, ... S_n} \} be the set of \eqn{n} genomic segments obtained from bulk DNA-sequencing. Further, let \eqn{S_I \in S} be the subset of segments for which cells within a clone have a well defined copy number state. 
Pearson Correlation Coefficients are calculated as similarity metric between each new cell and the consensus profile of each clone, based on segments \eqn{s \in S_I}. Each cell is assigned to the clone to which it is most similar.\cr
Alternatively, if similarity is set to false, the Euclidean distance metric is used instead of the Pearson Correlation.
}
\value{List with same components as input:
	  \item{cnps}{
	     Segment-by-cell matrix of copy number states, with new cells added as columns.
    }         
	  \item{sps}{
       The clone membership of each cell (that is columns in cnps).
    }
}
\author{Noemi Andor}
\examples{
  data(cnps)
  data(eps)
  set.seed(3)
  rcells1 = sample(colnames(cnps), 120)
  rcells2 = setdiff(colnames(eps), rcells1)
  outc = clusterCells(cnps[apply(cnps, 1, var)>0, rcells1])
  outc = assignCellsToClusters(outc, eps[,rcells2])
}
