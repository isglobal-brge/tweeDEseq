getZhuParam <- function(x)
       {
         a <- x[1]
         b <- x[2]
         c <- x[3]
         mu <- (b*c)/((1-c)^(1-a))
         D <- (1-a*c)/(1-c)
         a <- a 
         ans <- c(mu,D,a)
         names(ans) <- c("mu", "D", "a")
         ans
       }


puig2hou <- function(x)
  {
   mu <- x[1]
   d <- x[2]
   k <- x[3]
   alpha <- (k-2)/(k-1)
   delta <- mu*((k-1)*(d-1))^(1/(k-1))
   theta <- (k-1)*(d-1)
   ans <- c(delta, theta, alpha) 
   names(ans) <- c("delta", "theta", "alpha")
   ans
  }


hou2puig <- function(x)
 {
   delta <- x[1]
   theta <- x[2]
   alpha <- x[3]
   mu <- delta*theta^(alpha-1)
   d <- 1+(1-alpha)*theta
   k <- (alpha-2)/(alpha-1)
   ans <- c(mu, d, k)
   names(ans) <- c("mu", "d", "k")
   ans
 }


zhu2hou <-function (x) 
{
    if (attr(x, "method") == "MLE") 
        x <- x$par
    else x <- c(x$par[2], x$b, x$par[3])
    alpha <- x[1]
    delta <- x[2] * (x[3])^x[1]
    theta <- (1 - x[3])/x[3]
    ans <- c(alpha, delta, theta)
    names(ans) <- c("alpha", "delta", "theta")
    ans
}
