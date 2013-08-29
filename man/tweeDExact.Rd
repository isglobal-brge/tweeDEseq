\name{tweeDExact}
\alias{tweeDExact}
\alias{exactTestPT}

\title{
  Exact test for differences between two Poisson-Tweedie groups
}
\description{
  Carry out an exact test for differences between two Poisson-Tweedie populations.
}
\usage{
tweeDExact(counts, group, tol = 1e-15, mc.cores = 1)
exactTestPT(counts, group, tol = 1e-15, threshold = 150e3)
}

\arguments{
  \item{counts}{
    The RNA-seq counts. An object of type 'matrix' or 'data.frame' for 'tweeDExact', or an
    object of type 'vector' for 'exactTest'. 
  }
  \item{group}{
    vector giving the experimental group/condition for each
    sample/library.
  }
  \item{tol}{
    Tolerance for the Poisson-Tweedie probability computations. The
    probabilities under the 'tol' value will automatically considered as 0.
  }
  \item{threshold}{
    an integer (default is 50e3). If the sum of all counts in a
  certain gene excedes this value 'testPoissonTweedie' will be called
  instead of 'exactTest'. Larger values will result in a longer
  computing time.}
  \item{mc.cores}{
    number of cpu cores to be used. This option is only available
    when the 'multicore' package is installed and loaded
    first. In such a case, if the default value of
    'mc.cores=1' is not changed, all available cores will be
    used.
  }
}
\details{
  'exactTest' performs the exact test for a vector of counts.

  'tweeDExact' performs the test for a whole 'data.frame'. The P-values
  are then corrected using the Benjamini and Hochberg method.
}
\value{
  'exactTest' returns the p-value resulting from the exact test between
  two different Poisson-Tweedie populations, as well as the method that
  was used to compute it. 

  'tweeDExact' returns a 'data.frame'. Each row corresponds to a gene
  and it contains the following information:

  - In the first columns the mean of counts in each of the subgroups.

  - In the third column the p-value of the test for differential
  expression between the two subgroups.

  - In the fourth column the p-value corrected for multiple comparisons
  using the Benjamini-Hochberg FDR procedure.

  - In the last (fifth) column the method that was used to compute the p-value.
}
\references{
  P. Hougaard, M.L. Ting Lee, and G.A. Whitmore (1997). Analysis of
  overdispersed count data by mixtures of poisson variables and
  poisson processes. Biometrics 53, pages 1225-1238.
}
\author{
  Mikel Esnaola
}

\seealso{
  \code{\link{testPoissonTweedie}}
  \code{\link{tweeDExact}}
}
\examples{
counts <- matrix(rPT(n = 1000, a = 0.5, mu = 10, D = 5), ncol = 40)

tweeDExact(counts, group = rep(c(1,2),20))
}
\keyword{htest}
