loglikGlmPT <- function(par, X, Y, offset=NULL, allFactors=FALSE, a=NULL, tol=1e-300){
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
      if(a!=1){
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
    for(i in 1:nobs){
      m <- exp(X[i,]%*%par[1:ncov] + offset[i])
      b <- (m*(1-c)^(1-a))/c
      if(a!=1){
        probs <- .Call("zhuprobs", as.integer(Y[i]), a, b, c, tol)
        probs <- probs[length(probs)]
      }
      else
        probs <- dpois(Y[i],m)
      if(probs<=0 || is.nan(probs))
        l <- l - 1e7
      else
        l <- l + log(probs)
    }
  }
  l
}
