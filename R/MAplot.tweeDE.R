MAplot <- function(x, ...) UseMethod("MAplot")

MAplot.tweeDE <- function(x, log2fc.cutoff=0, highlight=NULL, ...)
 {
  if(!inherits(x, "tweeDE"))
   stop("x must be an object of class 'tweeDE'")
  
  A <- log2(x$overallMean)
  M <- x$log2fc
  names(A) <- names(M) <- rownames(x)

  plot(A, M, col="grey", pch=21, bg="grey", ...)
  abline(h=0, lwd=3, col="black")
  if (log2fc.cutoff > 0)
    abline(h=c(-1, 1)*log2fc.cutoff, lwd=3, lty=2, col=grey(0.85))
  if (!is.null(highlight)) {
    if (!is.null(names(highlight))) highlight <- list(highlight)
    for (i in 1:length(highlight)) {
      if (any(is.na(match(highlight[[i]]$genes, rownames(x)))))
        warning(paste("Some genes in the", i, "'highlight' argument are not in 'x'"))
      args.points <- list(x=A[highlight[[i]]$genes],
                          y=M[highlight[[i]]$genes])
      args.points[names(highlight[[i]])[-match("genes", names(highlight[[i]]))]] <- highlight[[i]][-match("genes", names(highlight[[i]]))]
      do.call("points", args.points)
    }
  }
  invisible(ifelse(is.null(highlight),list(A=A, M=M),list(A=A, M=M, args.points=args.points)))
 }
