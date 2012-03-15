\name{tweeDEseq-internal}
\alias{tweeDEseq-internal}

\alias{loglikPoissonTweedie}
\alias{loglikPoissonTweedie2}
\alias{momentEstimates}
\alias{shapeTrend}
\alias{loglikGlmPT}
\alias{glmPT.fit}

\title{Internal 'tweeDEseq' functions}
\description{Internal tweeDEseq functions}

\usage{
loglikPoissonTweedie(p, x, mu, verbose = FALSE, tol = 1e-15, probs =
FALSE, w = NULL)
loglikPoissonTweedie2(p, a, x, mu, verbose = FALSE, tol = 1e-15, probs =
FALSE, w = NULL)
momentEstimates(x, w)
shapeTrend(x)
loglikGlmPT(par, X, Y, offset = NULL, allFactors=FALSE, a=NULL,
tol=1e-300, maxCount, verbose = FALSE)
glmPT.fit(X, Y, offset=NULL, allFactors=FALSE, a = NULL, maxCount =
2000, ...)
}

\details{These are not to be called by the user}
\keyword{internal}
