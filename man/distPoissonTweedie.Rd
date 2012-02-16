\name{dPT}
\alias{dPT}
\alias{rPT}

\title{
The Poisson-Tweedie family of distributions
}
\description{
  Density function and random generation for the Poisson-Tweedie family of distributions. 
}
\usage{
dPT(x, mu, D, a, tol = 1e-15)
rPT(n, mu, D, a, max = 10*sqrt(mu*D), tol = 1e-4)
}

\arguments{
  \item{x}{
    an object of class 'mlePT' or a non-negative vector 
    containing the integers in which the distribution
    should be evaluated.
}
  \item{mu}{
    numeric positive scalar giving the mean of the distribution.
}
  \item{D}{
    numeric positive scalar giving the dispersion of the distribution.
}
  \item{a}{
    numeric scalar smaller than 1 giving the shape parameter of the distribution. 
}
  \item{tol}{
    numeric scalar giving the tolerance.
  }
  \item{n}{
    integer scalar giving number of random values to return.
  }
  \item{max}{
    numeric scalar containing the maximum number of counts to be used in
    the sampling process.
  }
}

\value{
  If 'x' is of class 'mlePT', 'dPT' will return the Poisson-Tweedie
  distribution with parameters equal to the ones estimated by
  'mlePoissonTweedie' evaluated on the data that was used to estimate
  the parameters. If 'x' is a numeric vector, 'dPT' will return the density of the
  specified Poisson-Tweedie distribution evaluated on 'x'.
  
  'rPT' generates random deviates.
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
  \code{\link{compareCountDist}}
  \code{\link{testShapePT}}
  }

\examples{
# To compute the density function in 1:100 of the Polya-Aeppli
# distribution with mean = 20 and dispersion = 5
dPT(x = 1:100, mu = 20, D = 5, a = -1)

# To generate 100 random counts of the same distribution with same
# parameters
rPT(n = 100, mu = 20, D = 5, a = -1) 
}

\keyword{distribution}