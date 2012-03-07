loglikGlmPT <- function(par, X, Y, offset=NULL, allFactors=FALSE, a=NULL, tol=1e-300, verbose = FALSE){
  ## WHILE IS NOT IMPLEMENTED IN C FOR ALL FACTORS:
  allFactors <- FALSE

  l <- 0
  ncov <- ncol(X)
  nobs <- length(Y)
  c <- par[ncov+1]
  if (is.null(a))
    a <- par[ncov+2]
  if(allFactors && is.null(offset)){
    offset <- rep(0, nobs)
    Xcol <- apply(X,1,paste,collapse="")
    occ <- table(Xcol)
    occI <- as.list(occ)
    for(i in 1:length(occ))
      occI[[i]] <- which(Xcol==names(occ[i]))
    valsChar <- strsplit(names(occ), "")
    vals <- lapply(valsChar, as.integer)
    for(i in 1:length(vals)){
      occItab <- table(Y[occI[[i]]])
      uniq <- as.integer(names(occItab))
      m <- exp(vals[[i]]%*%par[1:ncov] + offset[i])
      b <- (m*(1-c)^(1-a))/c
      if(abs(a)<0.001)
        pUniq <- dnbinom(uniq, mu = m, size = b)
      else if(a<0.999){
        probs <- .Call("zhuprobs", as.integer(max(uniq)), a, b, c, tol)
        pUniq <- probs[uniq+1]
      }
      else
        pUniq <- dpois(uniq,m)
      if(any(pUniq<=0) || any(is.nan(pUniq)))
        l <- l - occ[i]*1e7
      else
        l <- l + log(pUniq)%*%occItab
    }
  }
  else{
    if(is.null(offset))
      offset <- rep(0, nobs)

    X2 <- as.list(as.data.frame(t(X)))
    X2 <- lapply(X2, as.numeric)
    OFFSET <- offset
    A <- a
    C <- c
    NOBS <- nobs
    NCOV <- ncov
    PAR <- par

    res <- .Call("loglikGlm", as.integer(NOBS), as.integer(NCOV), A,  C, PAR,
                 X2, as.integer(Y), 1e-300, OFFSET)
    l <- res
  }
  if(verbose)
    cat("loglik=",l,", a=",a,", c=",c,"\n", sep="")

  l
}




