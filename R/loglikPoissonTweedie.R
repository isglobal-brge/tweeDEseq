
#
# Zhu reparameterization mu, D, a (mu estimated using mean sample)
#


loglikPoissonTweedie <- function(p, x, mu, verbose=FALSE, tol=1e-15){
  x.t <- table(x)
  x.unique <- as.numeric(names(x.t))  
  mm <- max(x.unique)
  nn <- length(x.unique)

  # p = c(D, a) -> reparam pag. XXX Zhou

  D <- p[1]  
  a <- p[2]
  b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))   
  c <- (D-1)/(D-a)


  if(a!=1)
    prx <- .Call("zhuprobs", as.integer(mm), a, b, c, tol)
  else
    prx <- dpois(0:mm, b)
  q <- prx[x.unique+1]
  if(any(q==0) || any(is.na(q)))
   loglik <- -1e20
  else
  loglik <- sum(log(q)*x.t)

  if (verbose)
   cat("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")

  loglik
}


#
# Zhu 2-parameters    
#


loglikPoissonTweedie2 <- function(p, a, x, mu, verbose=FALSE, tol=1e-15){

  x.t <- table(x)
  x.unique <- as.numeric(names(x.t))  
  mm <- max(x.unique)
  nn <- length(x.unique)

  # p = c(D) -> reparam pag. XXX Zhou

  D <- p[1]  
  b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))  
  c <- (D-1)/(D-a)


  if(a==0)
    prx <- dnbinom(0:mm, mu=mu, size=b)
  else if(a==1)
    prx <- dpois(0:mm, b)
  else
    prx <- .Call("zhuprobs", as.integer(mm), a, b, c, tol)

  q <- prx[x.unique+1]
  if(any(q==0) || any(is.na(q)))
   loglik <- -1e20
  else
  loglik <- sum(log(q)*x.t)

  if (verbose)
   cat("loglik=",loglik,",D=", D,",a=", a,",b=",b,",c=",c,"\n")

  loglik
}
