#' Function to format percentages (from the stats package)
#'
#' @param probs the percentages
#' @param digits the number of digits
format.perc <- function (probs, digits) {
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%");
}