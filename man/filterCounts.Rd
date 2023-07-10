\name{filterCounts}
\alias{filterCounts}

\title{
Count data filtering
}
\description{
Filter count data to remove lowly expressed genes.
}

\usage{
filterCounts(counts, 
             cpm.cutoff=0.5, 
             n.samples.cutoff=2, 
             mean.cpm.cutoff=0, 
             lib.sizes=NULL)
}

\arguments{
  \item{counts}{
    numeric data.frame or matrix containing the count data.
  }
  \item{cpm.cutoff}{
    expression level cutoff defined as the minimum number of counts per million. By 
    default this is set to 0.5 counts per million.
  }
  \item{n.samples.cutoff}{
    minimum number of samples where a gene should meet the counts per million cutoff 
    (\code{cpm.cutoff}) in order to be kept as part of the count data matrix. When 
    \code{n.samples.cutoff} is a number between 0 and 1, then it is interpreted as the 
    fraction of samples that should meet the counts per million cutoff (\code{cpm.cutoff}).
  }
  \item{mean.cpm.cutoff}{
    minimum mean of counts per million cutoff that a gene should meet in order to be kept.
    When the value of this argument is larger than 0 then it overrules the other arguments
    \code{cpm.cutoff} and \code{n.samples.cutoff}.
  }
  \item{lib.sizes}{
    vector of the total number of reads to be considered per sample/library. If
    \code{lib.sizes=NULL} (default) then these quantities are estimated as the column
    sums in the input matrix of counts.
  }
}

\details{
  This function removes genes with very low expression level defined in terms of a minimum
  number of counts per million occurring in a minimum number of samples. Such a policy was
  described by Davis McCarthy in a message at the bioc-sig-sequencing mailing list.
  By default, this function keeps genes that are expressed at a level of 0.5 counts per
  million or greater in at least two samples. Alternatively, one can use the 
  \code{mean.cpm.cutoff} to set a minimum mean expression level through all the samples.
}

\value{
  A matrix of filtered genes.
}

\references{
  Davis McCarthy, 
  \url{https://stat.ethz.ch/pipermail/bioc-sig-sequencing/2011-June/002072.html}.
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
