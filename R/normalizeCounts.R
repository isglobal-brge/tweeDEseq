## the code for this function has been derived from the steps
## taken by the function exactTest() from the edgeR package
## right before differential expression is tested within that function

normalizeCounts <- function(counts, group, common.disp=FALSE, prior.n=1) {
  message("Using edgeR normalization methods.")
  d <- DGEList(counts=counts, group=group)
  message("Calculating normalization factors with the TMM method.")
  d <- calcNormFactors(d)
  message("Estimating common dispersion.")
  d <- estimateCommonDisp(d)
  if (!common.disp) {
    message("Estimating tagwise dispersions.")
    d <- estimateTagwiseDisp(d, prior.n=prior.n)
  }
  message("Calculating effective library sizes.")
  lib.size <- d$samples$lib.size * d$samples$norm.factors
  d2 <- DGEList(counts=counts, group=group, lib.size=lib.size)
  if (common.disp) {
    dispersion <- rep(d$common.dispersion, nrow(counts))
    message("Adjusting counts to effective library sizes using a common dispersion.")
  } else {
    dispersion <- d$tagwise.dispersion
    message("Adjusting counts to effective library sizes using tagwise dispersions.")
  }
  dispersion <- pmax(dispersion, 1e-06)
  counts <- ceiling(equalizeLibSizes(d2, disp=dispersion, null.hypothesis=TRUE)$pseudo-0.5)
  counts
}

filterCounts <- function(counts, cpm.cutoff=0.5, n.samples.cutoff=2) {
  cpm <- 1e06 * t( t(counts) / colSums(counts))
  sel <- rowSums(cpm > cpm.cutoff) >= n.samples.cutoff
  counts[sel, ]
}
