#' Substitute p-values for pseudo-QTNs.
#'
#' Replace the p-values calculated for pseudo-QTNs according to one of several
#' rules.
#'
#' @param GM A dataframe with three columns containing the SNP ID, chromosome,
#'   and base-pair position of each marker.
#' @param GLM A \code{FarmCPUpp} object returned by \code{\link{quick_lm}}.
#' @param QTN A three column matrix containing the SNP ID, chromosome, and
#'   base-pair position of each pseudo-QTN.
#' @param method A character scalar, one of "penalty", "reward", "mean",
#'   "median", or "onsite".
#' @param npc An integer scalar. The number of user-defined covariates.
#'
#' @return The input \code{GLM} object with the pseudo-QTN p-values in
#'   \code{GLM[["P"]]} substituted.
#'
#' @details Fill in details.
#'
#' @author Xiaolei Liu
#' @author Zhiwu Zhang
#' @author Aaron Kusmec

substitute_pvalue <- function(GM, GLM, QTN = NULL, method = "mean", npc = 0) {

  # pseudo-QTNs are required
  if (nrow(QTN) < 1) return(GLM)

  # `GLM$P` will have more than 4 columns if p-values exist for pseudo-QTNs.
  # If execution makes it here (`QTN` is non-NULL) and there are 4 or fewer
  # columns, an error has occurred elsewhere, so exit.
  stopifnot(ncol(GLM$P) > 4)

  # Match the pseudo-QTNs with their position in the p-value matrix
  position <- match(QTN[, 1], GM[, 1], nomatch = 0)
  position <- position[position != 0]
  nqtn <- length(position)

  if (nqtn > 0) {
    # Marker p-values are always in column 4; following columns contain
    # pseudo-QTN p-values.
    index <- 5:ncol(GLM$P)

    if (length(index) > 1) {
      if (method == "penalty") P.QTN <- apply(GLM$P[, index], 2, max, na.rm = TRUE)
      if (method == "reward") P.QTN <- apply(GLM$P[, index], 2, min, na.rm = TRUE)
      if (method == "mean") P.QTN <- apply(GLM$P[, index], 2, mean, na.rm = TRUE)
      if (method == "median") P.QTN <- apply(GLM$P[, index], 2, median, na.rm = TRUE)
      if (method == "onsite") P.QTN <- GLM$P0[(npc + 2):length(GLM$P0)]
    } else {
      if (method == "penalty") P.QTN <- max(GLM$P[, index], na.rm = TRUE)
      if (method == "reward") P.QTN <- min(GLM$P[, index], na.rm = TRUE)
      if (method == "mean") P.QTN <- mean(GLM$P[, index], na.rm = TRUE)
      if (method == "median") P.QTN <- median(GLM$P[, index], median, na.rm = TRUE)
      if (method == "onsite") P.QTN <- GLM$P0[npc + 2]
    }

    GLM$P[position, "pvalue"] <- P.QTN
  }

  return(GLM)
}
