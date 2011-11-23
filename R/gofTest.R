gofTest <- function(counts, a=0, mc.cores=1) {
  ff <- function(x) {
    ans <- NA
    suppressWarnings(ans <- try(testShapePT(mlePoissonTweedie(x), a=a), TRUE))
    if (!inherits(ans, "try-error")) {
      ans <- ans$statistic
    } else {
      ans <- NA
    }
    ans
  }
  haveMulticore <- tweeDEseq:::.isPackageLoaded("multicore")
  if (mc.cores > 1){
    if(!haveMulticore)
      stop("In order to run calculations in parallel with multiple cores the 'multicore' library should be loaded first")
    else{
      mclapp <- get('mclapply', envir=getNamespace('multicore'))
      countsList <- as.list(as.data.frame(t(counts)))
      gof <- unlist(mclapp(countsList, ff), use.names=FALSE)
    }
  }
  else
    gof <- apply(counts, 1, ff)
  gof
}
