tweeDExact <- function(counts, group, tol = 1e-15, mc.cores = 1){
  if (!is.factor(group))
    group <- as.factor(group)
  groups <- levels(group)
  ngroups <- length(groups)
  if (!ngroups >= 2)
    stop("two groups are required")
  if (ngroups > 2)
    stop("more than two groups is not yet supported")
  if (ncol(counts) != length(group))
    stop("'group' and 'counts' must have the same length")
  
  test.i <- function(countsi, group, tol){
    aux <<- aux + 1
    setTxtProgressBar(pb, aux)
    suppressWarnings(p <- try(exactTestPT(counts = countsi, group = group, tol = tol)))
    if (!inherits(p, "try-error"))
      ans <- p
    else
      ans <- NA
    ans
  }

  data <- as.data.frame(t(counts))

  test.i.mc <- function(countsi, group, col, nc, coreID, tol){
    masterDesc <- get("masterDescriptor", envir = getNamespace("parallel"))
    if (masterDesc() == coreID) {
      aux <<- aux + 1
      setTxtProgressBar(pb, nc * aux)
    }
    suppressWarnings(p <- try(exactTestPT(counts = countsi, group = group, tol = tol)))
    if (!inherits(p, "try-error"))
      ans <- p
    else
      ans <- NA
    ans
  }
  
  ngenes <- nrow(counts)
  aux <- 0
  pb <- txtProgressBar(min = 0, max = ngenes, initial = 0, 
                       style = 3)
  masterDesc <- try(get("masterDescriptor", envir = getNamespace("parallel")), 
                    TRUE)
  if (mc.cores > 1) {
    if (class(masterDesc) == "try-error") 
      stop("It appears you are trying to use multiple cores from Windows, this is not possible")
    nAvailableCores <- detectCores()
    coreID <- mclapply(as.list(1:mc.cores), function(x) masterDesc(), 
                       mc.cores = mc.cores)[[1]]
    res <- data.frame(mclapply(data, test.i.mc, mc.cores = mc.cores, 
                               group = group, nc = mc.cores, coreID = coreID, tol = tol))
    setTxtProgressBar(pb, ngenes)
  }
  else res <- data.frame(lapply(data, test.i, group = group, 
                                  tol = tol))
  close(pb)

  means <- apply(counts, 1, function(x, gg){
    levs <- unique(gg)
    res <- c(mean(x[gg==levs[1]]), mean(x[gg==levs[2]]))
    res
  }
                 , gg = group)
  res2 <- data.frame(means[1,], means[2,], as.numeric(t(res[1,])),p.adjust(t(res[1,]), "BH"),t(res[2,]))
  rownames(res2) <- rownames(counts)
  colnames(res2) <- c(paste("mean.",unique(group),sep=""), "pval", "pval.adjust", "info")
  res2
}
