\name{segmentExpression2CopyNumber}
\alias{segmentExpression2CopyNumber}
\title{
Calling CNVs.
}
\description{
Maps single cell expression profiles to copy number profiles.
}
\usage{
segmentExpression2CopyNumber(eps, gpc, cn, seed=0, outF=NULL, maxPloidy=8, 
                                       nCores=2, stdOUT="log.applyAR2seg")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eps}{
    Segment-by-cell matrix of expression.
  }
  \item{gpc}{
    Number of genes expressed per cell.
  }
  \item{cn}{
    Average copy number across cells for each segment (i.e. row in eps).
  }
  \item{seed}{
    The fraction of entries in a-priori segment-by-cell copy number matrix to be used as seed for association rule mining.
  }
  \item{outF}{
    Output file prefix in which to print intermediary heatmaps and histograms, or NULL (default) if no print.
  }
  \item{maxPloidy}{
    The maximum ploidy to accept as solution.
  }
  \item{nCores}{
    The numbers of threads used.
  }
  \item{stdOUT}{
    Log-file to which standard output is redirected during parallel processing.
  }
}
\details{
Let S := \{ \eqn{S_1, S_2, ... S_n} \} be the set of \eqn{n} genomic segments obtained from bulk DNA-sequencing. Let \eqn{E_{ij}} and \eqn{G_{ij}} be the average number of UMIs and the number of expressed genes per segment \eqn{i} per cell \eqn{j}. The segment-by-cell expression matrix is first normalized by gene coverage. For each \eqn{x \in S}, the linear regression model: \cr
              \deqn{E_{x*} \sim \sum_{i \in S}G_{i*} } \cr
, fits the average segment expression per cell onto the cell's overall gene coverage. The model's residuals \eqn{R_{ij}} reflect inter-cell differences in expression per segment that cannot be explained by differential gene coverage per cell. A first approximation of the segment-by-cell copy number matrix CN is given by:\cr
              \deqn{CN_{ij} := R_{ij} * (cn_i /  \bar{R_{i*}}  )}\cr
, where \eqn{cn_i} is the population-average copy number of segment \eqn{i} derived from DNA-seq. Above transformation of \eqn{E_{ij}} into \eqn{CN_{ij}} is in essence a numerical optimization, shifting the distribution of each segment to the average value expected from bulk DNA-seq.
\cr
\cr
Let  \eqn{x' \in CN} be the measured copy number of a given segment-cell pair, and \eqn{x} its corresponding true copy number state. The probability of assigning copy number  \eqn{x} to a cell  \eqn{j} at locus  \eqn{i} depends on:
\cr
\bold{A. Cell  \eqn{j}'s read count at locus \eqn{i}}, calculated conditional on the measurement  \eqn{x'}. 
Using a Gaussian smoothing kernel, we compute the kernel density estimate of the read counts at locus \eqn{i} across cells to identify the major (\eqn{M}) and the minor (\eqn{m}) copy number states of \eqn{i} as the highest and second highest peak of the fit respectively. Then we calculate the proportion of cells expected at state  \eqn{m} as  \eqn{f = \frac{cn_i - M}{m - M} }. The probability of assigning copy number \eqn{x} to a cell \eqn{j} at locus \eqn{i} is calculated as:\cr
\eqn{  P_A(x|x') \sim }\cr
\bold{:}  \eqn{  0, if x \notin {m,M}}\cr
\bold{:}  \eqn{  P_{ij}(x'|N(m, sd = f)), if x == m}\cr
\bold{:}  \eqn{  P_{ij}(x'|N(M, sd = 1-f)), if x == M}\cr 
\cr
\bold{B. Cell  \eqn{j}'s read count at other loci}, i.e. how similar the cell is to other cells that have copy number \eqn{x} at locus \eqn{i}. We use Apriori - an algorithm for association rule mining - to find groups of loci that tend to have correlated copy number states across cells. Let \eqn{V_{i,K \to x}} be the set of rules concluding copy number \eqn{x} for locus \eqn{i}, where \eqn{k \in K} are copy number profiles of up to \eqn{n=4} loci in the form \{ \eqn{S_1=x_1, S_2=x_2, ... S_n=x_n} \}. Further let \eqn{C_r} be the confidence of a rule \eqn{r \in V_{i,K \to x}}. For each cell \eqn{j \in J} matching any of the copy number profiles in \eqn{K}, we calculate:\cr
\eqn{  P_B(x) \sim \sum_{r \in V_{i,K \to x}}C_r } \cr
, the cumulative confidence of the rules in support of \eqn{x} at \eqn{i}.
\cr
\cr
We first obtain a seed of cell-segment pairs by assigning a-priori copy number states only when \eqn{argmax_{x \in [1,8]}  P_A (x|x') > t}. We use this seed as input to \bold{B}. Finally, a-posteriori copy number for segment \eqn{i} in cell \eqn{j} is calculated as:\cr
\deqn{argmax_{x \in [1,8]} P_A(x|x') + P_B(x) }
}
\value{
  Segment-by-cell matrix of copy number states.
}
\references{
Andor, N.*, Lau, B.*, Catalanotti, C., Kumar, V., Sathe, A., Belhocine, K., Wheeler, T., et al. (2018) Joint single cell DNA-Seq and RNA-Seq of gastric cancer reveals subclonal signatures of genomic instability and gene expression. doi: https://doi.org/10.1101/445932
\cr
\cr
Borgelt C & Kruse R. (2002) Induction of Association Rules: Apriori Implementation.
}
\author{
Noemi Andor
}
\seealso{
\code{\link{apriori}}
}
\examples{
##Calculate number of genes expressed per each cell:
data(epg)
gpc = apply(epg>0, 2, sum)

##Call function:
data(eps)
data(segments)
cn=segments[rownames(eps),"CN_Estimate"]
\donttest{
	cnps = segmentExpression2CopyNumber(eps, gpc, cn, seed=0.5, nCores=2, stdOUT="log")
	head(eps[,1:3]); ##Expression of first three cells
	head(cnps[,1:3]); ##Copy number of first three cells
}
}
