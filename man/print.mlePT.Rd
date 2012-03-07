\name{Methods for objects of class 'mlePT'}
\alias{print.mlePT}
\alias{logLik.mlePT}
\alias{confint.mlePT}

\title{
Methods for objects of class 'mlePT'
}
\description{
print, extract loglikelihood or compute confidence interval for an
object of class 'mlePT'.
}
\usage{
\method{print}{mlePT}(x, digits = 3, ...)
\method{logLik}{mlePT}(object, ...)
\method{confint}{mlePT}(object, parm, level = 0.95, ...)
}

\arguments{
  \item{x}{
    object of class 'mlePT'.
  }
  \item{object}{
    object of class 'mlePT'.
    }
  \item{digits}{
    integer scalar giving the number of digits to be rounded the solution.
  }
  \item{parm}{
    a specification of which parameters are to be given
    confidence intervals, either a vector of numbers or a vector
    of names.  If missing, all parameters are considered.
  }
  \item{level}{
    the confidence level required (default is 0.95).
    }
  \item{\dots}{
    additional arguments.
  }
}

\value{
  'logLik' returns the loglikelihood of the selected model.
  
  'confint' returns a matrix (or vector) with columns giving lower and upper
  confidence limits for each parameter. 
}

\seealso{
  \code{\link{mlePoissonTweedie}}
  }

\examples{
# Load and aggregate the 'seizure' database
data(seizure)
aggCounts <- aggregate(x = cbind(seizure$count, seizure$trx), by =
list(seizure$id), FUN = sum)

# Estimate the parameters
mleSeizure <- mlePoissonTweedie(x = aggCounts[,2], a.ini = 0, D.ini =
10)

# Print
mleSeizure

# Extract loglikelihood
logLik(mleSeizure)

# Compute confidence inerval
confint(mleSeizure)
}

\keyword{methods}