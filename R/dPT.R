dPT <- function (x, mu, D, a, tol=1e-15)
{
  if (inherits(x, "mlePT"))
    {
      a <- x$paramZhu[1]
      b <- x$paramZhu[2]
      c <- x$paramZhu[3]
      mu <- x$par[1]
      obs <- x$x
    }
  else
    {
      if(any(x<0))
        stop("'x' must be non-negative")
      if (missing(mu) || missing(D) || missing(a))
        stop("parameters mu, D and a are required")
      a <- a
      b <- (mu*(1-a)^(1-a))/((D-1)*(D-a)^(-a))
      c <- (D-1)/(D-a)
      obs <- x
    }
  
  x.t <- table(obs)
  x.unique <- as.numeric(names(x.t))
  mm <- max(x.unique)
  
  if (a==0)
    prx <- dnbinom(0:mm, mu=mu, size=b)
  else if (a==1)
    prx <- dpois(0:mm, b)
  else
    prx <- .Call("zhuprobs", as.integer(mm), a, b, c, tol)
  
  res <- prx[obs+1]
  res
}
