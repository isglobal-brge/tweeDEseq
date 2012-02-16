\name{testShapePT}
\alias{testShapePT}

\title{
Test shape parameter of PT
}
\description{
Function to test whether the shape parameter is equal to a given value. 
}
\usage{
testShapePT(x, a = 0)
}

\arguments{
  \item{x}{
    object of class 'mlePT'.
}
  \item{a}{
    numeric scalar smaller than 1. The function will test whether the
    shape parameter is equal to the introduced 'a' (default is 0).
}
}
\details{
  By default a = 0, and therefore the function tests whether the count
  data follows a Negative-Binomial distribution or not. In this case, a
  Likelihood Ratio Test is performed. When the given value for 'a' is
  different from 0, a Wald test is performed.

  If a = 1, the function tests whether the count data follows a
  Poisson distribution or not.

  If a = 0.5, the function tests whether the count data follows a
  Poisson-inverse Gaussian distribution or not.

  If a = -1, the function tests whether the count data follows a
  Polya-Aeppli distribution or not.
}
\value{
  numeric p-value of the test.
}
\references{
  M. Esnaola, P. Puig, D. Gonzalez, R. Castelo, J.R. Gonzalez.
  A flexible count data model to fit the wide diversity of expression
  profiles arising from extensively replicated RNA-seq
  experiments. Submitted.

  A.H. El-Shaarawi, R. Zhu, H. Joe (2010). Modelling species abundance
  using the Poisson-Tweedie family. Environmetrics 22, pages 152-164.
  
  P. Hougaard, M.L. Ting Lee, and G.A. Whitmore (1997). Analysis of
  overdispersed count data by mixtures of poisson variables and poisson
  processes. Biometrics 53, pages 1225-1238.
}

\seealso{
  \code{\link{gofTest}}
  \code{\link{mlePoissonTweedie}}
  \code{\link{compareCountDist}}
  }

\examples{
# Generate a random matrix of counts
counts <- rPT(n=1000, a=0.5, mu=10, D=5)

# Maximum likelihood estimation of the Poisson-Tweedie parameters
mleEstimate <- mlePoissonTweedie(x = counts, a.ini = 0, D.ini
= 10)

# Test whether data comes from Negative-Binomial distribution
testShapePT(mleEstimate)

# Test whether data comes from Poisson-inverse Gaussian
testShapePT(mleEstimate, a = 0.5)
}

\keyword{htest}