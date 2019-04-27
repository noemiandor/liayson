\name{aggregateSegmentExpression}
\alias{aggregateSegmentExpression}
\title{
Aggregating genes across copy number segments.
}
\description{
Calculates average expression of genes grouped by common segment membership.
}
\usage{
aggregateSegmentExpression(epg, segments, mingps = 20, GRCh=37)
}
\arguments{
  \item{epg}{
    Gene-by-cell matrix of expression. Recommendation is to cap extreme UMI counts (e.g. at the 99\% quantile) and to include only cells expressing at least 1,000 genes.
  }
  \item{segments}{
    Matrix in which each row corresponds to a copy number segment as calculated by a circular binary segmentation algorithm. Has to contain at least the following column names:\cr
    	\bold{chr} - chromosome;\cr
    	\bold{startpos} - the first genomic position of a copy number segment; \cr
    	\bold{endpos} - the last genomic position of a copy number segment; \cr
	    \bold{CN_Estimate} - the copy number estimated for each segment.
  }
  \item{mingps}{
    Minimum number of expressed genes a segment needs to contain in order to be included in output.
  }
  \item{GRCh}{
    Human reference genome version to be used for annotating gene coordinates.
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

\author{
Noemi Andor
}

\examples{
data(epg)
data(segments)
#X=aggregateSegmentExpression(epg, segments, mingps=20, GRCh=38)

}

