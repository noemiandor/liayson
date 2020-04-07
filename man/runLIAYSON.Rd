\name{runLIAYSON}
\alias{runLIAYSON}
\title{
Main Function.
}
\description{
Given an RNA-seq derived cell-by-gene matrix and an DNA-seq derived copy number segmentation, LIAYSON predicts the number of clones present in a tumor, their size, the copy number profile of each clone and the clone membership of each single cell.
}
\usage{
runLIAYSON(X, S, sName, mingps = 20, GRCh = 37, h = 0.2, minSegLength=1E6, outD = NULL)
}
\arguments{
 \item{X}{
    Gene-by-cell matrix of expression. Recommendation is to cap extreme UMI counts (e.g. at the 99\% quantile) and to include only cells expressing at least 1,000 genes.
  }
  \item{S}{
    Matrix in which each row corresponds to a copy number segment as calculated by a circular binary segmentation algorithm. Has to contain at least the following column names:\cr
    	\bold{chr} - chromosome;\cr
    	\bold{startpos} - the first genomic position of a copy number segment; \cr
    	\bold{endpos} - the last genomic position of a copy number segment; \cr
	    \bold{CN_Estimate} - the copy number estimated for each segment.
  }
  \item{sName}{
    Sample name.
  }
  \item{mingps}{
    Minimum number of expressed genes a segment needs to contain in order to be included in output.
  }
  \item{GRCh}{
    Human reference genome version to be used for annotating gene coordinates.
  }
  \item{h}{
    Height at which the tree should be cut (see also \code{\link{cutree}}).
  }
  \item{minSegLength}{
    Minimum number of base pairs below which a segment is to be excluded when defining clones.
  }
  \item{outD}{
    The output directory.
  }
}
\author{
Noemi Andor
}

\seealso{
\code{\link{clusterCells}}
\code{\link{segmentExpression2CopyNumber}}
}
\examples{
data(epg)
data(segments)
\donttest{
	out = runLIAYSON(epg, segments, sName="NCI-N87", GRCh = 38, h=0.05)
}
}
