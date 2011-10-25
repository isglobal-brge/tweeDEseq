qqchisq <- function(stat, df=1, normal=FALSE, rangeExpected=FALSE, obsQuantiles=c(0.50, 0.75, 0.95), ...) {
  stat <- stat[!is.na(stat)]
  expctd <- qchisq(ppoints(length(stat)), df=df)[order(order(stat))]
  obsvd <- stat

  if (normal) {
    stat[stat==0] <- stat[stat==0] + .Machine$double.eps
    z <- zscoreGamma(stat, shape=df/2, scale=2)
    rng <- range(z)
    if (rangeExpected)
      rng <- range(qnorm(ppoints(length(z))))
    qq <- qqnorm(z, pch=".", cex=4, xlab="Expected Z",
                 ylab="Observed Z", ylim=rng, ...)
    abline(0, 1, lwd=2)
    expctd <- qq$x
    obsvd <- qq$y
  } else {
    rng <- range(obsvd)
    if (rangeExpected)
      rng <- range(expctd)
    plot(expctd, obsvd, xlab=expression(paste("Expected ", chi^2)),
         ylab=expression(paste("Observed ", chi^2)), pch=".", cex=4, ylim=rng, ...)
    abline(0, 1, lwd=2)
  }

  if (length(obsQuantiles) > 0) {
    abline(h=quantile(obsvd, prob=obsQuantiles), lty=3, lwd=3, col=grey(0.75))
    axis(4, at=quantile(obsvd, prob=obsQuantiles),
         labels=paste(ceiling(obsQuantiles*100), "%", sep=""), las=1)
  }

  invisible(list(x=expctd, y=obsvd))
}
