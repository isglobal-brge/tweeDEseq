\name{compareCountDist}
\alias{compareCountDist}

\title{
Compare count data distributions
}
\description{
  Compares the empirical and estimated distributions for different count data models
}
\usage{
  compareCountDist(x, plot=TRUE, ...)
}

\arguments{
  \item{x}{numeric vector containing the read counts.
  }
  \item{plot}{If \code{TRUE} (the default) then the plot with the ECDF function for the counts and the three different Poisson-Tweedie distributions is produced, otherwise no graphical output is given and this only makes sense if one is interested in the returned value (see value section below).
  }
  \item{...}{Further arguments to be passed to the plot function.
}
}
\details{
  This function serves the purpose of comparing a empirical distribution of counts with three Poisson-Tweedie distributions arising from estimating mean, dispersion and setting \eqn{a=1} for comparing against a Poisson, \eqn{a=0} for comparing against a negative binomial and estimating the shape parameter a from data too. The legend shows the values of the \eqn{a} parameter and the P-value of the likelihood ratio test on whether the expression profile follows a negative binomial distribution (\eqn{H_0:a=0}).
}
\value{
  List with the following components:
  \item{a}{shape parameter estimated from the input data \code{x}.}
  \item{p.value}{P-value for the test that the data follows a negative binomial distribution, i.e., \eqn{H_0:a=0}.}
}

\references{
  Esnaola M, Puig P, Gonzalez D, Castelo R and Gonzalez JR (2013). A
  flexible count data model to fit the wide diversity of expression
  profiles arising from extensively replicated RNA-seq experiments. BMC Bioinformatics 14: 254
    
}

\seealso{
  \code{\link{qqchisq}}
  \code{\link{testShapePT}}
  }

\examples{
# Generate 500 random counts following a Poisson Inverse Gaussian
# distribution with mean = 20 and dispersion = 5
randomCounts <- rPT(n = 500, mu = 20, D = 5, a = 0.5)

xx <- compareCountDist(randomCounts, plot=FALSE)
xx
}

\keyword{models}
