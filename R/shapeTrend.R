shapeTrend <- function(x){
  ll <- apply(x, 1, momentEstimates)
  sx <- log2(ll[1,])
  sy <- ll["k",]
  l <- lowess(sx, sy, f = 0.2)
  f <- approxfun(l, rule = 2)
  fx <- f(sx)
  fitted.a <- ((fx-2)/(fx-1))*(fx>1) + (-4)*(fx<=1)
  param <- list(mu=ll["mu",], D=ll["D",], a=fitted.a)
  param
}
