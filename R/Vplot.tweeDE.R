Vplot <- function(x, ...) UseMethod("Vplot")

Vplot.tweeDE <- function(x, log2fc.cutoff=0, pval.adjust.cutoff=1,
                         highlight=NULL, ylab=expression(paste(-log[10], " Raw P-value")), ...)
 {
  if(!inherits(x, "tweeDE"))
   stop("x must be an object of class 'tweeDE'")
  
  M <- x$log2fc
  P <- -log10(x$pval)
  names(M) <- names(P) <- rownames(x)

  plot(M, P, col="grey", pch=21, bg="grey", ylab=ylab, ...)
  if (log2fc.cutoff > 0)
    abline(v=c(-1, 1)*log2fc.cutoff, col=grey(0.85), lwd=3, lty=2)
  if (pval.adjust.cutoff < 1) {
    abline(h=-log10(max(x$pval[x$pval.adjust <= pval.adjust.cutoff])),
           col=grey(0.85), lwd=3, lty=2)
    text(max(x$log2fc)*0.9, 0.9*-log10(max(x$pval[x$pval.adjust <= pval.adjust.cutoff])), labels=sprintf("%.0f%% FDR", pval.adjust.cutoff*100), pos=3)
  }
  if (!is.null(highlight)) {
    if (!is.null(names(highlight))) highlight <- list(highlight)
    # for (i in 1:length(highlight)) {
    for (i in seq_len(length(highlight))) {
      if (any(is.na(match(highlight[[i]]$genes, rownames(x))))){
            warning("Some genes in the ", i, "'highlight' argument are not in 'x'")
       }
      args.points <- list(x=M[highlight[[i]]$genes],
                          y=P[highlight[[i]]$genes])
      args.points[names(highlight[[i]])[-match("genes", names(highlight[[i]]))]] <- highlight[[i]][-match("genes", names(highlight[[i]]))]
      do.call("points", args.points)
    }
  }
  invisible(list(M=M, P=P, args.points=args.points))
 }
