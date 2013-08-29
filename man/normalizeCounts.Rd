\name{normalizeCounts}
\alias{normalizeCounts}

\title{
Count data normalization
}
\description{
Normalize count data to remove systematic technical effects.
}

\usage{
normalizeCounts(counts, group=rep.int(1,ncol(counts)), method=c("TMM", "cqn"),
                common.disp = FALSE, prior.df=8, annot=NULL, lib.sizes=NULL, verbose=TRUE)
}

\arguments{
  \item{counts}{
    numeric data.frame or matrix containing the count data.
  }
  \item{group}{
    vector giving the experimental group/condition for each
    sample/library. This argument is only relevant
    when \code{method="cqn"}.
  }
  \item{method}{
    specific method to use in order to normalize the input matrix of counts. By default
    this is set to TMM (Robinson and Oshlack, 2010) using the implementation available
    in the edgeR package. The other option is cqn (Hansen, Irizarry and Wu, 2012).
  }
  \item{common.disp}{
    logical indicating whether a common or tagwise (default) dispersions
    should be estimated and employed when adjusting counts. This argument is only relevant
    when \code{method="TMM"}.
  }
  \item{prior.df}{
    argument provided to the call of \code{\link[edgeR]{estimateTagwiseDisp}}
    which defines the prior degrees of freedom. It is used in
    calculating 'prior.n' which, in turn, defines the amount of shrinkage of the estimated tagwise
    dispersions to the common one. By default \code{prior.df=8} thus
    assumming no shrinkage toward that common dispersion. This
    argument is not used if \code{common.disp=TRUE}. This argument is only relevant
    when \code{method="TMM"}.
  }
  \item{annot}{
    matrix or data frame with row names matching at least part of the row names in the
    \code{counts} input matrix, containing feature/tag/gene lengths in bp on its first column,
    and a second covariate, such as G+C content, on its second column. These two pieces
    of information are provided to arguments \code{lengths} and \code{x} when calling
    \code{\link[cqn]{cqn}}. This argument is only relevant when \code{method="TMM"}.
  }
  \item{lib.sizes}{
    vector of the total number of reads to be considered per sample/library. If
    \code{lib.sizes=NULL} (default) then these quantities are estimated as the column
    sums in the input matrix of counts.
  }
  \item{verbose}{
    logical indicating whether progress should be reported.
  }
}

\details{
  This function encapsulates calls to RNA-seq normalization procedures
  available in
  the \code{\link[edgeR]{edgeR}} and \code{\link[cqn]{cqn}} packages
  in order to try
  to remove systematic technical effects from raw counts.By default,
  the
  TMM method described in Robinson and Oshlack (2010) is employed
  to calculate normalization factors which are applied to
  estimate effective library sizes, then common and tagwise
  (only when the
  argument common.disp=TRUE) dispersions are calculated
  (Robinson and Smyth,
  Bioinformatics 2007) and finally counts are adjusted so
  that library sizes
  are approximately equal for the given dispersion values
  (Robinson and
  Smyth, Biostatistics 2008).Setting the argument
  \code{method="cqn"}, conditional
  quantile normalization (Hansen, Irizarry and Wu,
  2012) is applied which aims at
  adjusting for tag/feature/gene length and other
  covariate such as G+C content. This
  information should be provided through the
  \code{annot} argument. This procedure
  calculates, for every gene and every sample,
  an offset to apply to the log2 reads per
  million (RPM) and the function
  \code{normalizeCounts()} adds this offset to the
  the log2 RPM values calculated from the
  input count data matrix, unlogs them and rolls
  back these normalized RPM values into
  integer counts. Details on these two normalization
  procedures are given in the
  documentation of the \code{edgeR} and \code{cqn} Bioconductor
  packages.
  
}

\value{
  A matrix of normalized counts.
}

\references{
  K.D. Hansen, R.A. Irizarry and Z. Wu. Removing technical variability in RNA-seq
  data using conditional quantile normalization. \emph{Biostatistics}, 2012.

  M.D. Robinson and A. Oshlack. A scaling normalization method for
  differential expression analysis of RNA-seq data. \emph{Genome Biol},
  11:R25, 2010.

  Robinson MD and Smyth GK (2007). Moderated statistical tests for
  assessing differences in tag abundance. _Bioinformatics_ 23,
  2881-2887

  Robinson MD and Smyth GK (2008). Small-sample estimation of
  negative binomial dispersion, with applications to SAGE data.
  _Biostatistics_, 9, 321-332
}

\author{
  J.R. Gonzalez and R. Castelo 
}

\seealso{
  \code{\link{filterCounts}}
}

\examples{
# Generate a random matrix of counts
counts <- matrix(rPT(n=1000, a=0.5, mu=10, D=5), ncol = 40)

colSums(counts)
counts[1:5, 1:5]

# Normalize counts
normCounts <- normalizeCounts(counts, rep(c(1,2), 20))

colSums(normCounts)
normCounts[1:5, 1:5]
}

\keyword{misc}
 