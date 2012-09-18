exactTestPT <- function(counts, group, tol = 1e-15, threshold = 150e3){
  if (!is.factor(group))
    group <- as.factor(group)
  groups <- levels(group)
  ngroups <- length(groups)
  if (!ngroups >= 2)
    stop("two groups are required")
  if (ngroups > 2)
    stop("more than two groups is not yet supported")
  if (length(counts) != length(group))
    stop("'group' and 'counts' must have the same length")

  countsA <- counts[group == groups[1]]
  countsB <- counts[group == groups[2]]
  kA <- sum(countsA)
  kB <- sum(countsB)
  k <- kA + kB
  nA <- length(countsA)
  nB <- length(countsB)

  ## Estimate the parameters for each of the groups
  mle <- mlePoissonTweedie(counts)
  par <- mle$paramZhu

  info <- NA
  
  if(is.na(par[1])){
    p <- testPoissonTweedie(counts, group)$pvalue
    info <- "testPoissonTweedie - NA parameters"
  }
  else{
    if(par[1] == 0){
      p <- testPoissonTweedie(counts, group)$pvalue
      info <- "testPoissonTweedie - a=0"
    }
    else if(par[1] >= 0.95){
      ## Poisson
      lambda <- mle$par[1]
      lambdaA <- nA * lambda
      lambdaB <- nB * lambda
      pA <- dpois(0:k, lambdaA)
      pB <- dpois(0:k, lambdaB)
      pObsA <- pA[kA+1]
      pObsB <- pB[kB+1]
      probs <- pA*pB[length(pB):1]
      ix <- which(probs <= pObsA*pObsB)
      num <- sum(probs[ix])
      den <- sum(probs)
      p <- num/den
      info <- "Poisson"
    }
    else{
      if (k>threshold){
        p <- testPoissonTweedie(counts, group)$pvalue
        info <- "testPoissonTweedie - k>threshold"
      }
      else{
        ## Poisson - Tweedie
        aA <- par[1]
        bA <- par[2]*nA
        cA <- par[3]
        aB <- par[1]
        bB <- par[2]*nB
        cB <- par[3]
        ## Reparametrization (from Zhu to Hougaard)
        alphaA <- aA
        deltaA <- bA*cA^aA
        thetaA <- (1-cA)/cA
        alphaB <- aB
        deltaB <- bB*cB^aB
        thetaB <- (1-cB)/cB
        ##    muA <- (bA*cA)/(1-cA)^(1-aA)
        ##    muB <- (bB*cB)/(1-cB)^(1-aB)
        ##    cat("aA =", aA, ", bA =", bA, ", cA =", cA, "\n")
        ##    cat("aB =", aB, ", bB =", bB, ", cB =", cB, "\n")
        ##    cat("kA =", kA, ", muA =", muA, ", kB =", kB, ", muB =", muB, "\n")
        ## Probabilites for the estimated distribution of the sum
        pAl <- .Call("logprobs", as.integer(k), alphaA, deltaA, thetaA, tol)
        if(nA == nB)
          pBl <- pAl
        else
          pBl <- .Call("logprobs", as.integer(k), alphaB, deltaB, thetaB, tol)
        ## Transform (probabilities are in the logarithm scale)
        pA <- ifelse(pAl == 0, 0, exp(pAl))
        pB <- ifelse(pBl == 0, 0, exp(pBl))
        probs <- pA*pB[length(pB):1]
        ## Probability of observed values
        pObsAl <- ifelse(pA[kA+1] == 0, -Inf, pAl[kA+1])
        pObsBl <- ifelse(pB[kB+1] == 0, -Inf, pBl[kB+1])
        pObsl <- pObsAl + pObsBl
        ## Test statistic
        probsl <- pAl + pBl[length(pBl):1]
        ix <- which(probsl <= pObsl)
        num <- sum(probs[ix])
        den <- sum(probs)
        p <- num/den
        info <- "exactTest"
      }
    }
  }
  c(p,info)
}


