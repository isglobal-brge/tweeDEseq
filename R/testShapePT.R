testShapePT <- function (x, a=0){
  stat <- pval <- NA
  if(a==0){
    xx <- x$x
    paramNB <- fitdistr(xx, densfun="negative binomial")
    verNB <- sum(log(dnbinom(xx, size=paramNB$estimate[1], mu=paramNB$estimate[2])))
    stat <- 2*(x$log - verNB)
    if (stat >= 0)
      pval <- 1-pchisq(stat, 1)
    else
      stat <- NA
  }
  else{
    stat <- abs(x$par[3]-a)/x$se[3]
    pval <- unname(2*pnorm(stat, lower.tail=FALSE))
  }
  ans <- list(statistic = stat, pvalue = pval)
  ans
 }
