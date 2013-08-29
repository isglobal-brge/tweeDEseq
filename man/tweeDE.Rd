\name{tweeDE}
\alias{tweeDE}
\alias{testPoissonTweedie}
\alias{print.tweeDE}
\alias{MAplot}
\alias{MAplot.tweeDE}
\alias{Vplot}
\alias{Vplot.tweeDE}
\title{
Score test for differences between two Poisson-Tweedie groups
}
\description{
Carry out a score test for differences between two Poisson-Tweedie groups.
}
\usage{
tweeDE(object, group, mc.cores = 1, pair = NULL, a = NULL, \dots)
testPoissonTweedie(x, group, saveModel = FALSE, a = NULL, log = FALSE,  \dots)
MAplot(x, \dots)
Vplot(x, \dots)
\method{print}{tweeDE}(x, n=6L, sort.by="pval", log2fc.cutoff=0, pval.adjust.cutoff=1, print=TRUE, \dots)
\method{MAplot}{tweeDE}(x, log2fc.cutoff=0, highlight=NULL, \dots)
\method{Vplot}{tweeDE}(x, log2fc.cutoff=0, pval.adjust.cutoff=1, highlight=NULL, ylab=expression(paste(-log[10], " Raw P-value")), \dots)
}

\arguments{
  \item{object}{
    a \code{data.frame} or a \code{matrix} of RNA-seq counts.
  }
  \item{group}{
    vector giving the experimental group/condition for each
    sample/library.
  }
  \item{mc.cores}{
    number of cpu cores to be used. This option is only
    available when the 'multicore' package is installed and loaded first.
    In such a case, if the default value of \code{mc.cores=1} is not changed,
    all available cores will be used.
  }
  \item{pair}{
    vector of two elements containing the representants of each of the
    two groups (default is 'NULL').
  }
  \item{a}{
    for 'tweeDE' function, numerical vector with values strictly less than 1. It allows to fix the shape
    parameter 'a' during differential expression tests for each of the
    genes (rows of 'object').
    for 'testPoissonTweedie' function, numeral value strictly less than
    1 which fixes the shape parameter 'a' for the test.
  }
  \item{n}{
    maximum number of genes printed.
  }
  \item{sort.by}{
    character string, indicating whether genes should be ranked by
    their P-value (\code{pval}), which is the default setting, or
    by absolute log2 fold-change (\code{log2fc}).
  }
  \item{log2fc.cutoff}{
    cutoff on the minimum value of the log2 fold change.
  }
  \item{pval.adjust.cutoff}{
    cutoff on the maximum adjusted P-value (FDR).
  }
  \item{print}{
    logical; it indicates whether the output should be printed on the terminal.
  }
  \item{highlight}{
    list of arguments to the \code{points()} plotting function in order to highlight
    genes in the MA or volcano plots. A component called \code{genes} is expected to
    have the identifiers of the genes to be higlighted.
  }
  \item{ylab}{
    label on the y-axis of the volcano plot set by default to -log10 of the raw
     P-value which is what this plot displays on that axis.
  }
  \item{x}{
    object returned by the function \code{tweeDE} in the case of \code{print}
    and vector of count data in the case of \code{testPoissonTweedie}.
  }
  \item{saveModel}{
    logical indicating whether the results of fitting the model should
    be saved or not (default is 'FALSE').
  }
  \item{log}{
    logical (default is FALSE). If FALSE, the tested Null Hypothesis states that
    the difference between the means is 0 while, if TRUE, it states that
    the quotient between the logarithm of means is equal to 1. For
    this last case, the standard error is computed using the Delta
    Method.
  }
  \item{\dots}{
    additional arguments.
  }
}
\details{
  'testPoissonTweedie' performs the test for a vector of counts.
  
  'tweeDE' performs the test for a whole 'data.frame'. The P-values are
  then corrected using the Benjamini and Hochberg method.
}
\value{
  'testPoissonTweedie' returns a list with:
  
     'mean': means for each group
     'pvalue':  p-value for the test

  'tweeDE' returns a 'data.frame' with columns

  'overallMean': overall mean counts 
  'meanA': mean counts of the first group
  'meanB': mean counts of the second group
  'log2fc': logarigthm (base 2) of the fold-change (second group vs. first group)
  'pval': p-value for the test
  'pval.adjust': adjusted p-value using Benjamini-Hochberg method
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
  \code{\link{normalizeCounts}}
  \code{\link{mlePoissonTweedie}}
  }

\examples{
# Generate a random matrix of counts
counts <- matrix(rPT(n = 1000, a = 0.5, mu = 10, D = 5), ncol = 40)

# Test for differences between the two groups
tweeDE(counts, group = rep(c(1,2),20))
}

\keyword{htest}
