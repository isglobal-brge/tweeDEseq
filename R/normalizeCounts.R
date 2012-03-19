## the code for this function has been derived from the steps
## taken by the function exactTest() from the edgeR package
## right before differential expression is tested within that function

normalizeCounts <- function(counts, group=rep.int(1,ncol(counts)), method=c("TMM", "cqn"),
                            common.disp=FALSE, prior.n=1, annot=NULL, lib.sizes=NULL, verbose=TRUE) {
  method <- match.arg(method)

  if (is.null(lib.sizes))
    lib.sizes <- colSums(counts)

  if (method == "TMM") {
    if (verbose)
      message("Using edgeR-TMM normalization.")
    d <- edgeR::DGEList(counts=counts, group=group, lib.size=lib.sizes)
    if (verbose)
      message("Calculating normalization factors with the TMM method.")
    d <- edgeR::calcNormFactors(d)
    if (verbose)
      message("Estimating common dispersion.")
    d <- edgeR::estimateCommonDisp(d)
    if (!common.disp) {
      if (verbose)
        message("Estimating tagwise dispersions.")
      d <- edgeR::estimateTagwiseDisp(d, prior.n=prior.n)
    }
    if (verbose)
      message("Calculating effective library sizes.")
    lib.size <- d$samples$lib.size * d$samples$norm.factors
    d2 <- edgeR::DGEList(counts=counts, group=group, lib.size=lib.size)
    if (common.disp) {
      dispersion <- rep(d$common.dispersion, nrow(counts))
      if (verbose)
        message("Adjusting counts to effective library sizes using a common dispersion.")
    } else {
      dispersion <- d$tagwise.dispersion
      if (verbose)
        message("Adjusting counts to effective library sizes using tagwise dispersions.")
    }
    dispersion <- pmax(dispersion, 1e-06)
    counts <- ceiling(edgeR::equalizeLibSizes(d2, disp=dispersion, null.hypothesis=TRUE)$pseudo-0.5)
  }

  if (method == "cqn") {
    if (verbose)
      message("Using cqn normalization.")
    if (!is.matrix(annot) && !is.data.frame(annot) || is.null(dim(annot)) || dim(annot)[2] != 2)
      stop("cqn requires a two-column matrix or data frame in the argument 'annot'. Do help(normalizeCounts) for further information.")
    common <- intersect(rownames(counts), rownames(annot[!is.na(annot[, 1]) & !is.na(annot[, 2]), ]))
    if (length(common) == 0)
      stop("No row names in 'counts' match any row name in 'annot'.")

    cqnNorm <- cqn::cqn(counts[common, ],
                        lengths=as.numeric(annot[common, 1]),
                        x=as.numeric(annot[common, 2]),
                        sizeFactors=lib.sizes, verbose=verbose)
    log2rpm <- cqnNorm$y + cqnNorm$offset
    log2r <- sweep(cqnNorm$y + cqnNorm$offset, 2, log2(lib.sizes/1e6), FUN="+")
    counts <- ceiling(2^log2r-0.5)
  }

  counts
}

filterCounts <- function(counts, cpm.cutoff=0.5, n.samples.cutoff=2, mean.cpm.cutoff=0, lib.sizes=NULL) {
  if (is.null(lib.sizes))
    lib.sizes <- colSums(counts)

  sel <- rep(TRUE, nrow(counts))
  if (mean.cpm.cutoff <= 0) {
    cpm <- 1e06 * t( t(counts) / lib.sizes)
    f <- 1
    if (n.samples.cutoff > 0 && n.samples.cutoff < 1)
      f <- ncol(counts)
    sel <- rowSums(cpm > cpm.cutoff) >= floor(n.samples.cutoff * f)
  } else {
    cpm <- sweep(log2(counts+1), 2, log2(lib.sizes/10^6))
    mean.cpm <- rowMeans(cpm)
    sel <- mean.cpm > log2(mean.cpm.cutoff)
  }

  counts[sel, ]
}
