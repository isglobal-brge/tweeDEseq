
#
# Zhu reparameterization mu, D, a (mu estimated using mean sample)
#

loglikPoissonTweedie <- function(p, x, mu, verbose=FALSE, tol=1e-15, probs = FALSE, w = NULL){
  x.t <- table(x)
  x.unique <- as.numeric(names(x.t))  
  mm <- max(x.unique)
  nn <- length(x.unique)

  if(!is.null(w))
    w.t <- aggregate(w~x, FUN=sum)$w
  
  # p = c(D, a) -> reparam pag. XXX Zhou

  D <- p[1]  
  a <- p[2]
  b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))   
  c <- (D-1)/(D-a)

  if(abs(a)<0.001)
    prx <- dnbinom(0:mm, mu=mu, size=b)
  else  if(a<=1-1e-3)
    prx <- .Call("zhuprobs", as.integer(mm), a, b, c, tol)
  else
    prx <- dpois(0:mm, b)
  q <- prx[x.unique+1]
  if(any(q==0) || any(is.na(q)))
    loglik <- -1e20
  else{
    if(is.null(w))
      loglik <- sum(log(q)*x.t)
    else
      loglik <- sum(log(q)*w.t)#/sum(w.t)
  }    
  if (verbose)
   # cat("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")
      message("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")

  if(probs)
    attr(loglik, "probs") <- prx[x+1]
  
  loglik
}


#
# Zhu 2-parameters    
#


loglikPoissonTweedie2 <- function(p, a, x, mu, verbose=FALSE, tol=1e-15, probs = FALSE, w = NULL){

  x.t <- table(x)
  x.unique <- as.numeric(names(x.t))  
  mm <- max(x.unique)
  nn <- length(x.unique)

  if(!is.null(w))
    w.t <- aggregate(w~x, FUN=sum)$w
  
  # p = c(D) -> reparam pag. XXX Zhou

  D <- p[1]  
  b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))  
  c <- (D-1)/(D-a)


  if(abs(a)<=0.001)
    prx <- dnbinom(0:mm, mu=mu, size=b)
  else  if(a<=1-1e-3)
    prx <- .Call("zhuprobs", as.integer(mm), a, b, c, tol)
  else
    prx <- dpois(0:mm, b)
  
  q <- prx[x.unique+1]
  if(any(q==0) || any(is.na(q)))
    loglik <- -1e20
  else{
    if(is.null(w))
      loglik <- sum(log(q)*x.t)
    else
      loglik <- sum(log(q)*w.t)#/sum(w.t)
  }    
  if (verbose)
   # cat("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")
      message("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")

  if(probs)
    attr(loglik, "probs") <- prx[x+1]
  
  loglik
}


loglikPoissonTweedie3 <- function(p, D, x, mu, verbose=FALSE, tol=1e-15, probs = FALSE, w = NULL){

  x.t <- table(x)
  x.unique <- as.numeric(names(x.t))  
  mm <- max(x.unique)
  nn <- length(x.unique)

  if(!is.null(w))
    w.t <- aggregate(w~x, FUN=sum)$w
  
  # p = c(D) -> reparam pag. XXX Zhou

  a <- p[1]  
  b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))  
  c <- (D-1)/(D-a)

  if(abs(a)<=0.001)
    prx <- dnbinom(0:mm, mu=mu, size=b)
  else  if(a<=1-1e-3)
    prx <- .Call("zhuprobs", as.integer(mm), a, b, c, tol)
  else
    prx <- dpois(0:mm, b)

  q <- prx[x.unique+1]
  if(any(q==0) || any(is.na(q)))
    loglik <- -1e20
  else{
    if(is.null(w))
      loglik <- sum(log(q)*x.t)
    else
      loglik <- sum(log(q)*w.t)#/sum(w.t)
  }    
  if (verbose)
   # cat("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")
      message("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")

  if(probs)
    attr(loglik, "probs") <- prx[x+1]
  
  loglik
}
