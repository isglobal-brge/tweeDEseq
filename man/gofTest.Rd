\name{gofTest}
\alias{gofTest}

\title{
Test the goodness of fit of every row in a matrix of counts
}
\description{
Function to test the goodness of fit of every row in a matrix of counts
}
\usage{
gofTest(counts, a = 0, mc.cores = 1)
}

\arguments{
  \item{counts}{
    matrix of counts
}
  \item{a}{
    numeric scalar smaller than 1. The function will test whether the
    shape parameter is equal to the introduced 'a' (default is 0).
  }
  \item{mc.cores}{
    number of cpu cores to be used. This option is only
    available when the 'multicore' package is installed and loaded first.
    In such a case, if the default value of \code{mc.cores=1} is not changed,
    all available cores will be used.
  }
}
\details{
  By default a = 0, and therefore the function tests for every row of the
  input matrix of counts whether the count data follows a Negative-Binomial
  distribution. In this case, a Likelihood Ratio Test is performed. When
  the given value for 'a' is different from 0, a Wald test is performed. This
  function calls \code{\link{testShapePT}}.
}
\value{
  a vector of statistics that follows a \eqn{\chi^2} distribution with one
  degree of freedom under the null hypothesis.
}
\references{
  Esnaola M, Puig P, Gonzalez D, Castelo R and Gonzalez JR (2013). A
  flexible count data model to fit the wide diversity of expression
  profiles arising from extensively replicated RNA-seq experiments. BMC Bioinformatics 14: 254

  A.H. El-Shaarawi, R. Zhu, H. Joe (2010). Modelling species abundance
  using the Poisson-Tweedie family. Environmetrics 22, pages 152-164.
  
  P. Hougaard, M.L. Ting Lee, and G.A. Whitmore (1997). Analysis of
  overdispersed count data by mixtures of poisson variables and poisson
  processes. Biometrics 53, pages 1225-1238.
}

\seealso{
  \code{\link{testShapePT}}
  }

\examples{
## Generate a random matrix of counts
counts <- matrix(rPT(n=2000, a=0.5, mu=10, D=5), nrow=20)

## Perform the goodness-of-fit tests for every row in the matrix
chi2gof <- gofTest(counts)

## Calculate and sort the corresponding P-values for the
## null hypothesis that counts follow a negative binomial distribution
sort(pchisq(chi2gof, df=1, lower.tail=FALSE))
}

\keyword{htest}
