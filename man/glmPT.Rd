\name{glmPT}
\alias{glmPT}
\alias{tweeDEglm}
\title{
Fit Poisson-Tweedie generalized linear model.
}
\description{
'glmPT' is used to fit generalized linear models for the Poisson-Tweedie
family of distributions.
}
\usage{
tweeDEglm(formula, counts, data, mc.cores = 1, a = NULL, offset = NULL,  ...)
glmPT(formula, data, offset = NULL, a = NULL, ...)
}
\arguments{
  \item{formula}{
    an object of class 'formula': a symbolic description of the model to
    be fitted.
  }
  \item{counts}{
    Matrix or data.frame of counts for the 'tweeDEglm'.
  }
  \item{data}{
    an optional data frame, list or environment containing the
    variables in the model. If not found in 'data', the
    variables are taken from 'environment(formula)', typically
    the environment from which 'glm' is called.
  }
  \item{mc.cores}{
    number of cpu cores to be used. This option is only available
    when the 'multicore' package is installed and loaded first.
    In such a case, if the default value of 'mc.cores=1' is not
    changed, all available cores will be used.
  }
  \item{offset}{
    this can be used to specify an _a priori_ known component to
    be included in the linear predictor during fitting.
  }
  \item{a}{
    numeric vector (for 'tweeDEglm') or numeric scalar (for 'glmPT') smaller than 1. If specified the PT shape
    parameter will be fixed. For 'tweeDEglm', if the provided 'a' is a
    scalar this value will be used for all rows of 'counts' (genes).
  }
  \item{\dots}{
    additional arguments to be passed to the 'optim' 'control'
    options.
  }
}

\value{
  An object of class 'glmPT' containing the following information:
  \item{call}{
    the matched call.
  }
  \item{contrasts}{
    (where relevant) the contrasts used.
  }
  \item{convergence}{
    A character string giving any additional information
    returned by the optimizer, or 'NULL'.
  }
  \item{counts}{
    A two-element integer vector giving the number of calls to
    'fn' and 'gr' respectively. This excludes those calls needed
    to compute the Hessian, if requested, and any calls to 'fn'
    to compute a finite-difference approximation to the gradient.
  }
  \item{df}{
    Number of estimated parameters.
  }
  \item{fitted.values}{
    The fitted mean values, obtained by transforming the
    linear predictors by the inverse of the link function.
  }
  \item{hessian}{
    A symmetric matrix giving
    an estimate of the Hessian at the solution found.
  }
  \item{message}{
    A character string giving any additional information returned
    by the optimizer, or 'NULL'.
  }
  \item{ncov}{
    An integer giving the number of variables in the model.
  }
  \item{par}{
    A vector giving the estimated parameters.
  }
  \item{residuals}{
    The residuals in the final iteration of the IWLS fit.
  }
  \item{se}{
    A vector giving the standard error of the estimated parameters.
  }
  \item{value}{
    Value of the log-likelihood of the model in the last iteration.
  }
  \item{...}{
    Further arguments to be passed to the 'glm.fit' function.
  }
}
\author{
  Mikel Esnaola <mesnaola@creal.cat>
}

\examples{
counts <- matrix(rPT(n = 1000, a = 0.5, mu = 10, D = 5), ncol = 40)
g <- factor(rep(c(0,1), 20))
mod1 <- glmPT(counts[1,]~g)
mod1
summary(mod1)
anova(mod1)

mod2 <- tweeDEglm(~ g, counts) 
mod2
}
