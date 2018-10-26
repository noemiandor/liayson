\name{aggregateSegmentExpression}
\alias{aggregateSegmentExpression}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Aggregating genes across copy number segments.
}
\description{
Calculates average expression of genes grouped by common segment membership.
}
\usage{
aggregateSegmentExpression(epg, segments, mingps, GRCh=37)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{epg}{
    Gene-by-cell matrix of expression. Recommendation is to cap UMI counts at the 99% quantile and to include only cells expressing at least 1,000 genes.
  }
  \item{segments}{
    Matrix in which each row corresponds to a copy number segment as calculated by a circular binary segmentation algorithm. Has to contain at least the following column names:\cr
    	\bold{chr} - chromosome;\cr
    	\bold{startpos} - the first genomic position of a copy number segment; \cr
    	\bold{endpos} - the last genomic position of a copy number segment; \cr
	    \bold{CN_Estimate} - the copy number estimated for each segment. 
	  Recommendation is to exclude genomic segments in which the average number of expressed genes per cell is below 20.
  }
  \item{mingps}{
    Minimum number of expressed genes a segment needs to contain in order to be included in output.
  }
}
\details{
Let S := \{ \eqn{S_1, S_2, ... S_n} \} be the set of \eqn{n} genomic segments that have been obtained from DNA-sequencing a given sample (e.g. from bulk exome-sequencing, scDNA-sequencing, etc.). Genes are mapped to their genomic coordinates using the biomaRt package and assigned to a segment based on their coordinates. Genes are grouped by their segment membership, to obtain the average number of UMIs and the number of expressed genes per segment \eqn{S_j} per cell i. 
}
\value{List with fields:
	  \item{eps}{
	        Segment-by-cell matrix of expression values.
    }         
	  \item{gps}{
          Segment-by-cell matrix of the number of expressed genes.
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
data(epg)
data(segments)
X=aggregateSegmentExpression(epg, segments, mingps=20, GRCh=38)

##View output: 
library(futile.matrix)
lapply(X, peek)
}

