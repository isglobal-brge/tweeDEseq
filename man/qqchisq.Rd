\name{qqchisq}
\alias{qqchisq}

\title{
Chi-square quantile-quantile plot
}
\description{
Make a chi-square quantile-quantile plot.
}
\usage{
qqchisq(stat, df=1, normal=FALSE, rangeExpected=FALSE, obsQuantiles=c(0.50, 0.75, 0.95), ...)
}

\arguments{
  \item{stat}{
    vector of \eqn{\chi^2} statistics.
}
  \item{df}{
    degrees of freedom of \code{stat}.
}
  \item{normal}{
    logical; set to \code{TRUE} if the \eqn{\chi^2} statistics in \code{stat}
    should be transform into normal z-scores in order to improve the display
    of lower quantiles. For this purpose, this function uses the
    \code{\link{zscoreGamma}} function from the limma package. Default is set
    to \code{FALSE}.
}
  \item{rangeExpected}{
    logical; set to \code{TRUE} if the displayed range of the observed \eqn{\chi^2}
    statistics is restricted to the range of their expected values. Default is
    set to \code{FALSE}.
}
  \item{obsQuantiles}{
    observed quantiles to indicate by horizontal dash lines. By default, these
    are set to 50\%, 75\% and 95\%.
}
  \item{...}{
    further arguments to pass to the \code{\link{plot}} function.
}
}
\details{
  The main purpose of this function in the \code{tweeDEseq} package is
  to provide means to assess the goodness of fit of count data to the
  negative binomial distribution. The main input argument \code{stats}
  should be the output of \code{\link{gofTest}}.
}
\value{
  it returns invisibly a list with two components \code{x} and \code{y}
  corresponding to the coordinates of the plotted statistics.
}
\references{
  M. Esnaola, P. Puig, D. Gonzalez, R. Castelo, J.R. Gonzalez.
  A flexible count data model to fit the wide diversity of expression
  profiles arising from extensively replicated RNA-seq
  experiments. Submitted.
}

\seealso{
  \code{\link{compareCountDist}}
  \code{\link{testShapePT}}
  }


\examples{
## Generate a random matrix of counts
counts <- matrix(rPT(n=2000, a=0.5, mu=10, D=5), nrow=20)

## Perform the goodness-of-fit tests for every row in the matrix
chi2gof <- gofTest(counts)

\dontrun{
qqchisq(chi2gof)
}
}

\keyword{htest}
