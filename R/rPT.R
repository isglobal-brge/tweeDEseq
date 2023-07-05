rPT <- function(n, mu, D, a, max=10*sqrt(mu*D), tol=1e-4)
 {
   dd <- cumsum(dPT(0:max, mu, D, a))
   if (1-max(dd) > tol)
    warning("Increase 'max' argument to obtain better simulated data")
   prob <- runif(n, 0, max(dd))

   ans <- vapply(prob, function(x, dd, max) c(0:max)[dd>=x][1], dd=dd, max=max)
   ans
 } 
