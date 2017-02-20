#' Make a manhattan plot from GWAS results.
#'
#' @param res A dataframe with at least three named columns: a numeric column of
#'   chromosomes (\code{Chromosome}), a numeric column of base-pair positions
#'   (\code{Position}), and a numeric column of marker p-values
#'   (\code{p.value}).
#' @param colors A character vector of length two where elements are R color
#'   names.
#' @param cutoff A significance cutoff (default = \code{NULL}).
#' @param ylab A label for the y-axis (default =
#'   \code{expression(paste(-log[10], "p-value"))}).
#' @param ... Other graphics arguments.
#'
#' @author Aaron Kusmec
#'
#' @export

manhattan_plot <- function(res, colors = c("black", "grey80"), cutoff = NULL,
                           ylab = expression(paste(-log[10], "p-value")), ...) {
  # Check for appropriate column names
  if (!("Chromosome" %in% names(res))) stop("No column 'Chromosome.'")
  if (!("Position" %in% names(res))) stop("No column 'Position'.")
  if (!("p.value" %in% names(res))) stop("No column 'p.value'.")

  # Check the colors vector
  if (length(colors) < 2) stop("'colors' must have 2 elements.")
  if (length(colors) > 2)
    warning("length(colors) > 2; only the first 2 elements will be used.")

  res <- res[, c("Chromosome", "Position", "p.value")]
  res <- subset(res, !is.na(p.value))
  res <- res[order(res$Chromosome, res$Position, decreasing = FALSE), ]
  res$Color <- factor(res$Chromosome %% 2 + 1)
  res$Order <- 1:nrow(res)

  breakpoints <- split(res, res$Chromosome)
  breakpoints <- sapply(breakpoints, function(x) {
    (max(x$Order) + min(x$Order))/2
  })
  breakpoints <- data.frame(Chromosome = 1:length(breakpoints),
                            Position = breakpoints)

  ylim <- c(0, max(-log10(res$p.value*0.1)))
  plot(1, type = "n", xlim = c(1, max(res$Order)), ylim = ylim,
       xlab = "", ylab = ylab, xaxt = "n", las = 2, ...)
  axis(1, at = breakpoints$Position, labels = paste("chr", breakpoints$Chromosome),
       las = 2)
  points(res$Order, -log10(res$p.value), col = colors[res$Color], ...)
  if (!is.null(cutoff)) {
    abline(h = -log10(cutoff), lwd = 2, lty = 2, col = "red")
  }
}
