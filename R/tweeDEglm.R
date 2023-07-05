tweeDEglm <- function(formula, counts, data, mc.cores = 1, a = NULL, offset = NULL, ...){
  if( !is.null(offset) && !all(dim(counts) == dim(offset)) )
    stop("If provided, 'offset' must have the same dimensions as 'counts'")
  if(is.null(a))
    suppressWarnings(countData <- data.frame(t(counts)))
  else{
    if(length(a)==1){
      a <- rep(a, nrow(counts))
      warning("Only one 'a' value provided. Using this value for all genes")
    }
    if(length(a)!=nrow(counts))
      stop("If provided, 'a' must be a numeric vector of length equal to the number of rows of 'counts' (number of genes)")
    suppressWarnings(countData <- data.frame(t(cbind(counts,a))))
  }

  # countData <- lapply(as.list(1:ncol(countData)), function(i, counts, offset)
    countData <- lapply(as.list(seq_len(ncol(countData))), function(i, counts, offset)
                      list(countsi = counts[,i], offseti = offset[i,]),
                      counts = countData, offset = offset)
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
                                        #  shapes <- shapeTrend(counts)$a
                                        #  temp <- cbind(counts, shapes)
                                        #  if(is.null(names(counts)))
                                        #    countData <- data.frame(t(temp), row.names = paste("X",1:ncol(temp), sep=""))
                                        #  else
                                        #    countData <- data.frame(t(temp))
  test.i <- function(x, mat, nc, ...){
    if (nc > 1) {
      masterDesc <- try(get('masterDescriptor', envir=getNamespace('parallel')), TRUE)
      if (masterDesc() == 4) {
        aux <<- aux + 1
        setTxtProgressBar(pb, nc * aux)
      }
    }
    else {
      aux <<- aux + 1
      setTxtProgressBar(pb, nc * aux)
    }
    offseti <- x$offseti
    x <- x$countsi
    if(length(x)==nrow(mat))
      a <- NULL
    else if(length(x)==nrow(mat)+1){
      a <- x[length(x)]
      x <- x[-length(x)]
    }
    else
      stop("Dimension error. Please be sure the provided parameters are correctly specified")
#    modFull <- tryCatch(glmPT.fit(X=mat, Y=x, a=a, ...), warning = function(w) w)
#    modNull <- tryCatch(glmPT.fit(X=matrix(rep(1,length(x)),ncol=1), Y=x, a=a, ...), warning = function(w) w)
    modFull <- glmPT.fit(X=mat, Y=x, a=a, offset = offseti, ...)
    modNull <- glmPT.fit(X=matrix(rep(1,length(x)),ncol=1), Y=x, a=a, offset = offseti, ...)
    if( "warning"%in%class(modFull) && "warning"%in%class(modNull) )
      ans <- c("AICfull" = NA, "AICnull" = NA, "pval" = NA)
    else if ("warning"%in%class(modFull)){
      class(modNull) <- "glmPT"
      ans <- c("AICfull" = NA, "AICnull" = AIC(modNull), "pval" = NA)
    }
    else if ("warning"%in%class(modNull)){
      class(modFull) <- "glmPT"
      ans <- c("AICfull" = AIC(modFull), "AICnull" = NA, "pval" = NA)
    }
    else{
      class(modNull) <- "glmPT"
      class(modFull) <- "glmPT"
      pval <- anova.glmPT(modelNull=modNull, object=modFull)
      ans <- c("AICfull" = AIC(modFull), "AICnull" = AIC(modNull), "pval" = pval)
    }

    ans
  }
  ngenes <- nrow(counts)
  aux <- 0
  pb <- txtProgressBar(min = 0, max = ngenes, initial = 0,
                       style = 3)
  if (mc.cores > 1)
    temp <- t(data.frame(mclapply(countData, test.i, mat=X, nc = mc.cores, mc.cores=mc.cores, ...)))
  else
    temp <- t(data.frame(lapply(countData, test.i, mat=X, nc = 1, ...))) 
  setTxtProgressBar(pb, ngenes)
  pval.adjust <- p.adjust(temp[,"pval"], "BH")
  res <- data.frame(temp, pval.adjust)
  rownames(res) <- rownames(counts)
  close(pb)
  res
}
