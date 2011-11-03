tweeDE <- function(object, group, mc.cores=1, pair=NULL, ...)
  {
    x <- object
    totalTime <- proc.time()

    if (!is.matrix(object) && !is.data.frame(object))
        stop("'object' must be a matrix or a data.frame")
      
    if (!is.factor(group))
      group <- as.factor(group)

    groups <- levels(group)

   
    if (is.null(pair))
      pair <- groups

    if (length(pair) != 2) 
        stop("'pair' (number of groups) must be of length 2.")


    if (is.numeric(pair)) 
      pair <- levels(group)[pair]
    else
      pair <- as.character(pair)

    this.pair <- (group %in% pair)

    if (!all(pair %in% groups)) 
        stop("At least one element of given pair is not a group.\n Groups are: ", 
            paste(groups, collapse = " "), "\n")

    x <- x[ , this.pair] 
    group <- factor(as.vector(group[this.pair]))

    cat("Comparing groups:", as.vector(pair[2]), "-", as.vector(pair[1]), "\n")
    
    test.i <- function(x, g, cont, nc, ...)
      {
        aux <<- aux + 1
        setTxtProgressBar(pb, nc*aux)

        suppressWarnings( mod <- try(testPoissonTweedie(x, g, loglik=FALSE, maxCount=6000, ...), TRUE) )
        if (!inherits(mod, "try-error")) {
          ans <- c(mod$mean, mod$stat, mod$pval)
        }
        else {
          mm <- aggregate(x, by=list(g), FUN=mean) 
          ans <- c(unlist(mm[,2]), NA, NA)
        }
        ans
      }
    
    test.i.mc <- function(x, g, cont, nc, coreID, ...)
      {
        masterDesc <- get('masterDescriptor', envir=getNamespace('multicore'))
        if(masterDesc() == coreID){
          aux <<- aux + 1
          setTxtProgressBar(pb, nc*aux)
        }
        
        suppressWarnings( mod <- try(testPoissonTweedie(x, g, loglik=FALSE, maxCount=6000, ...), TRUE) )
        
        if (!inherits(mod, "try-error")) {
          ans <- c(mod$mean, mod$stat, mod$pval)
        }
        else {
          mm <- aggregate(x, by=list(g), FUN=mean) 
          ans <- c(unlist(mm[,2]), NA, NA)
        }
        ans
      }
    
    data <- data.frame(t(x))
    ngenes <- nrow(x)
    
    aux <- 0
    pb <- txtProgressBar(min=0,max=ngenes,initial=0,style=3)

    haveMulticore <- tweeDEseq:::.isPackageLoaded("multicore")
    if (mc.cores > 1 && !haveMulticore)
      stop("In order to run calculations in parallel with multiple cores the 'multicore' library should be loaded first")
    
    if (haveMulticore) {

      mclapp <- get('mclapply', envir=getNamespace('multicore'))
      detCor <- get('detectCores', envir=getNamespace('multicore'))
      masterDesc <- get('masterDescriptor', envir=getNamespace('multicore'))

      nAvailableCores <- detCor()
      if (mc.cores == 1)
        mc.cores <- nAvailableCores

      coreID <- mclapp(as.list(1:mc.cores), function(x) masterDesc(), mc.cores=mc.cores)[[1]]
      
      res <- t(data.frame(mclapp(data, test.i.mc, mc.cores = mc.cores, 
                                 g = group, nc = mc.cores, coreID = coreID)))
      setTxtProgressBar(pb, ngenes)
    }
    else
      res <-t(data.frame(lapply(data, test.i, g = group, nc = 1)))
    
    close(pb)
    colnames(res)[1:2] <- groups
    
    pval.bh <- p.adjust(res[,4], "BH")
    if(identical(pair, groups))
       pair.ind <- 2:1
    else
       pair.ind <- 1:2
    
    log2fc <- log2(res[,pair.ind[1]]/res[,pair.ind[2]])

    ans <- data.frame(overallMean=round(rowMeans(x),3),
                      round(res[,1:2],3), log2fc=log2fc, stat=res[,3], pval=res[,4],
                      pval.adjust=pval.bh)
 
    class(ans) <- c("tweeDE", "data.frame")
    attr(ans, "comparison") <- paste(as.vector(pair[2]), "-", as.vector(pair[1]))

    ans
  }

## from https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## function: isPackageLoaded
## purpose: to check whether the package specified by the name given in
##          the input argument is loaded. this function is borrowed from
##          the discussion on the R-help list found in this url:
##          https://stat.ethz.ch/pipermail/r-help/2005-September/078974.html
## parameters: name - package name
## return: TRUE if the package is loaded, FALSE otherwise

.isPackageLoaded <- function(name) {
  ## Purpose: is package 'name' loaded?
  ## --------------------------------------------------
  (paste("package:", name, sep="") %in% search()) ||
  (name %in% loadedNamespaces())
}

