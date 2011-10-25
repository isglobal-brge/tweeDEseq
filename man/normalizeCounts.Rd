\name{normalizeCounts}
\alias{normalizeCounts}

\title{
Count data normalization
}
\description{
Normalize count data to remove systematic technical effects.
}

\usage{
normalizeCounts(counts, group, common.disp = FALSE, prior.n=1)
}

\arguments{
  \item{counts}{
    numeric data.frame or matrix containing the count data.
  }
  \item{group}{
    vector giving the experimental group/condition for each
    sample/library.
  }
  \item{common.disp}{
    logical indicating whether a common or tagwise (default) dispersions
    should be estimated and employed when adjusting counts.
  }
  \item{prior.n}{
    argument provided to the call of \code{\link[edgeR]{estimateTagwiseDisp}}
    which defines the amount of shrinkage of the estimated tagwise
    dispersions to the common one. By default \code{prior.n=1} thus
    assumming no shrinkage toward that common dispersion. This
    argument is not used if \code{common.disp=TRUE}.
  }

}

\details{
  This function aims at removing systematic technical effects from raw counts
  by using normalization procedures available in the \code{\link[edgeR]{edgeR}}
  package. In particular, the TMM method described in Robinson and Oshlack
  (2010) is employed to calculate normalization factors which are applied to
  estimate effective library sizes, then common and tagwise (only when the
  argument \code{common.disp=TRUE}) dispersions are calculated and finally
  counts are adjusted so that library sizes are approximately equal for the
  given dispersion values.
}

\value{
  A matrix of normalized counts.
}

\references{
  M.D. Robinson and A. Oshlack. A scaling normalization method for
  differential expression analysis of RNA-seq data. \emph{Genome Biol},
  11:R25, 2010.
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
