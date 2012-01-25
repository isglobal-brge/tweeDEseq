tweeDEglm <- function(formula, counts, data, mc.cores = 1,...){
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf, contrasts)
  nas <- attr(mf, "na.action")
  countData <- data.frame(t(counts))
  test.i <- function(x, mat, nc, ...){
    if (nc > 1) {
      masterDesc <- get('masterDescriptor', envir=getNamespace('parallel'))
      if (masterDesc() == 4) {
        aux <<- aux + 1
        setTxtProgressBar(pb, nc * aux)
      }
    }
    else {
      aux <<- aux + 1
      setTxtProgressBar(pb, nc * aux)
    }
    modFull <- glmPT.fit(X=mat, Y=x, ...)
    modNull <- glmPT.fit(X=matrix(rep(1,length(x)),ncol=1), Y=x, ...)
    pval <- anova.glmPT(modelNull=modNull, object=modFull)
    pval
  }
  ngenes <- nrow(counts)
  aux <- 0
  pb <- txtProgressBar(min = 0, max = ngenes, initial = 0,
                       style = 3)
  if (mc.cores > 1)
    res <- t(data.frame(mclapply(countData, test.i, mat=X, nc = mc.cores, mc.cores=mc.cores, ...)))
  else
    res <- lapply(countData, test.i, mat=X, nc = 1, ...)
  close(pb)
  res
}
