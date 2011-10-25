gofTest <- function(counts, a=0) {
  gof <- apply(counts, 1,
               function(x) {
                 ans <- NA
                 suppressWarnings(ans <- try(testShapePT(mlePoissonTweedie(x), a=a), TRUE))
                 if (!inherits(ans, "try-error")) {
                   ans <- ans$statistic
                 } else {
                   ans <- NA
                 }
                 ans
               })
  gof
}
