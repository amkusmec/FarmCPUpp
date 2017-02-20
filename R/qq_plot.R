#' Make a QQ plot from GWAS results.
#'
#' @param res A dataframe with at least three named columns: a numeric column of
#'   chromosomes (\code{Chromosome}), a numeric column of base-pair positions
#'   (\code{Position}), and a numeric column of marker p-values
#'   (\code{p.value}).
#' @param ... Other graphics arguments.
#'
#' @author Aaron Kusmec
#'
#' @export

qq_plot <- function(res, ...) {
  p_values <- res$p.value
  p_values <- p_values[!is.na(p_values) & p_values > 0 & p_values <= 1]

  if (length(p_values) < 1) return(NULL)

  p_values <- p_values[order(p_values)]
  quantiles <- (1:length(p_values))/(length(p_values) + 1)
  log_p_values <- -log10(p_values)
  log_quantiles <- -log10(quantiles)

  ci95 <- rep(NA, length(log_quantiles))
  ci05 <- rep(NA, length(log_quantiles))
  for (i in 1:length(log_quantiles)) {
    j <- ceiling((10^-log_quantiles[i])*length(log_quantiles))
    if (j == 0) j <- 1
    ci95[i] <- qbeta(0.95, j, length(log_quantiles) - j + 1)
    ci05[i] <- qbeta(0.05, j, length(log_quantiles) - j + 1)
  }

  plot(1, type = "n", xlim = c(0, max(log_quantiles)), ylim = c(0, max(log_p_values)),
       xlab = expression(Expected~~-log[10](p-value)),
       ylab = expression(Observed~~-log[10](p-value)), ...)
  polygon(c(rev(log_quantiles), log_quantiles),
          c(rev(-log10(ci05)), -log10(ci95)), col = "gray", border = NA)
  abline(a = 0, b = 1, col = "red", lwd = 2)
  points(log_quantiles, log_p_values, ...)
}
