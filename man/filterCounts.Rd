\name{filterCounts}
\alias{filterCounts}

\title{
Count data filtering
}
\description{
Filter count data to remove lowly expressed genes.
}

\usage{
filterCounts(counts, cpm.cutoff=0.5, n.samples.cutoff=2)
}

\arguments{
  \item{counts}{
    numeric data.frame or matrix containing the count data.
  }
  \item{cpm.cutoff}{
    expression level cutoff defined as the minimum number of counts per million.    By default this is set to 0.5 counts per million.
  }
  \item{n.samples.cutoff}{
    minimum number of samples where a gene should meet the counts per million cutoff (\code{cpm.cutoff}) in order to be kept as part of the count data matrix.
  }
}

\details{
  This function removes genes with very low expression level defined in terms of a minimum
  number of coutns per million occurring in a minimum number of samples. Such a policy was
  described by Davis McCarthy in a message at the bioc-sig-sequencing mailing list.
  By default, this function keeps genes that are expressed at a level of 0.5 counts per
  million or greater in at least two samples.
}

\value{
  A matrix of filtered genes.
}

\references{
  Davis McCarthy, \url{https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2011-June/002072.html}.
}

\author{
  J.R. Gonzalez and R. Castelo
}

\seealso{
  \code{\link{normalizeCounts}}
  }

\examples{
# Generate a random matrix of counts
counts <- matrix(rPT(n=1000, a=0.5, mu=10, D=5), ncol = 40)

dim(counts)

# Filter genes with requiring the minimum expression level on every sample
filteredCounts <- filterCounts(counts, n.samples.cutoff=dim(counts)[2])

dim(filteredCounts)
}

\keyword{misc}
