print.tweeDE <- function(x, n=6L, sort.by="pval", log2fc.cutoff=0,
                         pval.adjust.cutoff=1, print=TRUE, ...)
 {
  if(!inherits(x, "tweeDE"))
   stop("x must be an object of class 'tweeDE'")

  sort.by <- match.arg(sort.by, c("pval", "log2fc"))
  
  if (print)
    cat(paste("Comparison of groups:", attr(x, "comp"), "\n"))
  orderedIdx <- NA
  sortstr <- ""
  switch(sort.by,
         "log2fc"={orderedIdx <- order(abs(x[, "log2fc"])); sortstr <- "log2 fold-change"},
         "pval"={orderedIdx <- order(x[, "pval"]); sortstr <- "P-value"}
        )
  x <- x[orderedIdx, ]
  x <- x[abs(x$log2fc) > log2fc.cutoff & x$pval.adjust <= pval.adjust.cutoff, ]

  if (print) {
    topstr <- "Showing"
    if (nrow(x) > n)
      topstr <- paste("Showing top", n)
    cat(paste(topstr, "genes ranked by", sortstr, "\n"))
    cat(paste("Minimum absolute log2 fold-change of", round(log2fc.cutoff, digits=2), "\n"))
    cat(paste("Maximum adjusted P-value of", pval.adjust.cutoff, "\n"))
    print.data.frame(as.data.frame(x)[1:min(nrow(x), n), ], ...)
  }
  invisible(as.data.frame(x)[1:min(nrow(x), n), ])
 }
