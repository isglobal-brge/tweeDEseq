print.mlePT <- function(x, digits=3, ...)
 {

  if(!inherits(x, "mlePT"))
    stop("x should be an object of class 'mlePT'")

  out <- cbind (x$par, x$se)
  colnames(out) <- c("estimate", "s.e.")
  rownames(out) <- names(x$par)


  cat("\n")
  
  cat("Poisson-Tweedie parameter estimates (", x$method, ")", "\n", sep="")

  cat("\n")
  print(out, digits=2)
  cat("Skewness:", round(x$skewness,2))

  cat("\n")
  cat("\n Zhu parameterization \n")
  print(round(x$paramZhu, digits=digits))

  cat("\n Hougaard parameterization \n")
  print(round(x$paramHou, digits=digits))

  cat("\n log-likelihood:", x$loglik)
  cat("\n number of iterations:", x$iter)  
  cat("\n")
 }
