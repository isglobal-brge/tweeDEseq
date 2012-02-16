\name{mlePoissonTweedie}
\alias{mlePoissonTweedie}
\alias{getParam}
\alias{mlePT}

\title{
Maximum likelihood estimation of the Poisson-Tweedie parameters
}
\description{
Maximum likelihood estimation of the Poisson-Tweedie parameters using
L-BFGS-B quasi-Newton method.
}
\usage{
mlePoissonTweedie(x, a, D.ini, a.ini, maxit = 100, loglik=TRUE, maxCount=20000, ...)
getParam(object)
}

\arguments{
  \item{x}{numeric vector containing the read counts.
  }
  \item{a}{numeric scalar smaller than 1, if specified the PT shape parameter will be fixed.
  }
  \item{D.ini}{numeric positive scalar giving the initial value for the dispersion.
  }
  \item{a.ini}{numeric scalar smaller than 1 giving the initial value for the shape
    parameter (ignored if 'a' is specified).
  }
  \item{maxit}{numeric scalar providing the maximum number of 'L-BFGS-B'
    iterations to be performed (default is '100').
  }
  \item{loglik}{
    is log-likelihood computed? The default is TRUE
  }
  \item{object}{
    an object of class 'mlePT'.
  }
  \item{maxCount}{
    if max(x) > maxCount, then moment method is used to estimate model parameters to reduce computation time. The default is 20000.
  }
  \item{\dots}{additional arguments to be passed to the 'optim'
    'control' options.
  }
}
\details{
The L-BFGS-B quasi-Newton method is used to calculate iteratively the
maximum likelihood estimates of the three Poisson-Tweedie parameters. If
'a' argument is specified, this parameter will be fixed and the method
will only estimate the other two. 
}
\value{
  An object of class 'mlePT' containing the following information:

  par: numeric vector giving the estimated mean ('mu'), dispersion ('D') and shape parameter 'a'.

  se: numeric vector containing the standard errors of the estimated parameters 'mu', 'D' and 'a'.

  loglik: numeric scalar providing the value of the loglikelihod for the
  estimated parameters.

  iter: numeric scalar giving the number of performed iterations.

  paramZhu: numeric vector giving the values of the estimated parameters
  in the Zhu parameterization 'a', 'b' and 'c'.

  paramHou: numeric vector giving the values of the estimated parameters
  in the Hougaard parameterization 'alpha', 'delta' and 'theta'.

  skewness: numeric scalar providing the estimate of the skewness given
  the estimated parameters.

  x: numeric vector containing the count data introduced as the 'x'
  argument by the user.

  convergence:  A character string giving any additional information returned
  by the optimizer, or 'NULL'.
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
  \code{\link{testShapePT}}
  \code{\link{print.mlePT}}
  }
  
\examples{
# Generate 500 random counts following a Poisson Inverse Gaussian
# distribution with mean = 20 and dispersion = 5
randomCounts <- rPT(n = 500, mu = 20, D = 5, a = 0.5)

# Estimate all three parameters
res1 <- mlePoissonTweedie(x = randomCounts, a.ini = 0, D.ini
= 10)
res1
getParam(res1)

#Fix 'a = 0.5' and estimate the other two parameters
res2 <- mlePoissonTweedie(x = randomCounts, a = 0.5, D.ini
= 10)
res2
getParam(res2)
}

\keyword{models}
