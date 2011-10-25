compareCountDist <- function(x, plot=TRUE, ...) {

  ## Empirical distribution
  ans <- ecdf(x)

  ##
  ## Poisson-Tweedie
  ##

  ans1 <- mlePoissonTweedie(x)
  paramPT <- round(ans1$par[3],2)
  pp1 <- dPT(0:max(x), mu=ans1$par[1], D=ans1$par[2], a=ans1$par[3])

  ##
  ## Negative Binomial
  ##

  paramNB <- MASS::fitdistr(x,densfun="negative binomial")
  pp2 <- dnbinom(0:max(x), size=paramNB$estimate[1], mu=paramNB$estimate[2])
  verNB <- sum(log(dnbinom(x, size=paramNB$estimate[1], mu=paramNB$estimate[2])))

  ##
  ## Poisson
  ##
  paramPoi <- MASS::fitdistr(x,densfun="Poisson")
  pp3 <- dpois(0:max(x), lambda=paramPoi$estimate[1])
  p <- 1-pchisq(2*(ans1$log - verNB), 1)
  pNB <- formatC(p, dig=2)

  if (plot) {
    ss <- 0:max(x)
    col <- c("red", "blue", "darkgreen")

    plot(ans, cex=1.5, cex.lab=1.5, cex.axis=1.5, xlab="Counts", ...)
    if (min(x)!=0)  points(0,0,pch=19)
    lines(ss, cumsum(pp1), lwd=4, col=col[1])
    lines(ss, cumsum(pp2), lwd=4, col=col[2])
    lines(ss, cumsum(pp3), lwd=4, col=col[3])
    legend("bottomright", inset=0.01, bg="white",
           c(paste("P=", pNB, "(H0: a=0)"), 
           bquote(paste("PT(D, ", mu, ", ", .(paramPT), ")", sep="")), 
           expression(paste("PT(D, ",mu, ", 0) - NB",  sep="")),  
           expression(paste("PT(D, ",mu,", 1) - Poisson", sep=""))), 
           col=c("white", col), lwd=4, cex=1.4)
  }

  invisible(list(a=as.numeric(paramPT), p.value=p))
}
