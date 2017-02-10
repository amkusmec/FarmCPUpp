#' Set marker p-values based on prior expectations.
#'
#' Modifies marker p-values calculated by the GLM based on prior expectations.
#'
#' @param GM A dataframe with three columns containing the SNP ID, chromosome,
#'   and base-pair position of each marker.
#' @param P A numeric vector containing marker p-values.
#' @param Prior A dataframe with four columns. The first three correspond to the
#'   columns of \code{GM} and the fourth contains prior probabilities.
#'
#' @return A numeric vector of updated marker p-values.
#'
#' @details If \code{Prior} is \code{NULL}, p-values are not modified. Marker
#'   p-values are multiplied by the user-provided prior p-values. A value of 1
#'   will not modify the p-values calculated by the GLM.
#'
#' @author Zhiwu Zhang
#' @author Aaron Kusmec

prior <- function(GM , P = NULL, Prior = NULL) {

  if (is.null(Prior)) return(P)

  index <- match(Prior[, 1], GM[, 1], nomatch = 0)
  P[index] <- P[index]*Prior[, 4]

  return(P)
}
