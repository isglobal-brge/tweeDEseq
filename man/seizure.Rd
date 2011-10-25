\name{seizure}
\alias{seizure}
\docType{data}
\title{
Epileptic seizure counts
}
\description{
Data on seizure counts for 59 epileptics.
}
\usage{data(seizure)}
\format{
  A data frame with 236 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{a numeric vector, identification number for each patient}
    \item{\code{count}}{a numeric vector, seizure counts}
    \item{\code{visit}}{a numeric vector, visit number}
    \item{\code{trx}}{a numeric vector, treatment: progabide (1) or placebo (0)}
    \item{\code{baseline}}{a numeric vector, baseline 8 week seizure count}
    \item{\code{age}}{a numeric vector, age of patient}
  }
}
\details{
  The data are from a placebo-controlled clinical trial of 59
  epileptics.
  Patients with partial seizures were enrolled in a randomized
  clinical trial of the anti-epileptic drug, progabide. Participants
  in the study were randomized to either progabide or a placebo, as
  an adjuvant to the standard anti-epileptic chemotherapy. Progabide
  is an anti-epileptic drug whose primary mechanism of action is to
  enhance gamma-aminobutyric acid (GABA) content; GABA is the
  primary inhibitory neurotransmitter in the brain.
  Prior to receiving treatment, baseline data on the number of
  epileptic seizures during the preceding 8-week interval were
  recorded. Counts of epileptic seizures during  2-week intervals
  before each of four successive post-randomization clinic visits
  were recorded.
}
\source{
  P.F Thall, and S.C. Vail (1990). Some covariance models for
  longitudinal count data with overdispersion. Biometrics, 46, 657-671,
}
\references{
P. Hougaard, M.L. Ting Lee, and G.A. Whitmore (1997): Analysis of
overdispersed count data by mixtures of poisson variables and poisson
processes. Biometrics 53, pages 1225-1238.
}
\seealso{
  \code{\link{testPoissonTweedie}}
  \code{\link{mlePoissonTweedie}}
  }
\examples{
# Although this is not a differential expression dataset, it is appropriate
# to illustrate the application of the Poisson-Tweedie in
# epidemiological studies

data(seizure)
summary(seizure)

# Aggregate
aggCounts <- aggregate(x = cbind(seizure$count, seizure$trx), by =
list(seizure$id), FUN = sum)

# Estimation of the three parameters for all individuals
mleSeizure <- mlePoissonTweedie(x = aggCounts[,2], a.ini = 0, D.ini = 10)
mleSeizure

#Poisson-Tweedie test
testPoissonTweedie(x = aggCounts[,2], group = aggCounts[,3])
}
\keyword{datasets}
