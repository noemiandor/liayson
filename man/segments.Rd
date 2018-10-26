\name{segments}
\alias{segments}
\docType{data}
\title{Bulk copy number profile of NCI-N87 cell line.}
\description{Copy number segmentation matrix obtained as average among G0G1 cells.}
\usage{data(segments)}
\format{
	Matrix in which each row corresponds to a copy number segment as calculated by a circular binary segmentation algorithm. Has to contain at least the following column names:\cr
  \bold{chr} - chromosome;\cr
  \bold{startpos} - the first genomic position of a copy number segment; \cr
  \bold{endpos} - the last genomic position of a copy number segment; \cr
  \bold{CN_Estimate} - the copy number estimated for each segment. 
}

\source{Data obtained from Ji lab at Stanford.}

\keyword{datasets}
