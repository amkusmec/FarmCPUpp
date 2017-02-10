#' Remove highly correlated pseudo-QTNs.
#'
#' Removes highly correlated pseudo-QTNs based on the Pearson correlation
#' coefficients between marker scores. In each pair of highly correlated
#' markers, the marker with the lower p-value is retained.
#'
#' @param GDP A pointer to a big.matrix object containing the marker scores.
#' @param GM A data frame with three columns: marker ID, chromosome, and
#'   base-pair position.
#' @param seqQTN An integer vector of columns in GDP that have been selected as
#'   pseudo-QTNs.
#' @param seqQTN.p A numeric vector of p-values for the pseudo-QTNs.
#' @param threshold A numeric scalar. Pairs of pseudo-QTNs with correlation
#'   coefficients above this value will be pruned.
#'
#' @return An integer vector specifying the columns in \code{GDP} that are
#'   selected as pseudo-QTNs.
#'
#' @author Zhiwu Zhang
#' @author Aaron Kusmec

remove_qtns <- function(GDP = NULL, GM = NULL, seqQTN = NULL, seqQTN.p = NULL,
                        threshold = 0.99) {
  # Return NULL if we don't have any pseudo-QTNs
  if (is.null(seqQTN)) return(NULL)

  # Sort pseudo-QTN by p-value and get their map information
  seqQTN <- seqQTN[order(seqQTN.p, decreasing = FALSE)]
  binmap <- GM[seqQTN, ]
  hugeNum <- 1e10
  cb <- binmap[, 2]*hugeNum + binmap[, 3]
  cb.unique <- unique(cb)
  index <- match(cb.unique, cb, nomatch = 0)
  seqQTN <- seqQTN[index]

  # Number of pseudo-QTNs after ID construction
  n <- length(seqQTN)

  # Get marker scores for the pseudo-QTNs.
  X <- bigmemory::as.matrix(bigmemory::deepcopy(GDP, cols = seqQTN))

  # Remove pseudo-QTNs that are too highly correlated with other pseudo-QTNs
  R <- cor(X, method = "pearson")
  index <- abs(R) > threshold
  b <- R*0
  b[index] <- 1
  C <- 1 - b
  C[lower.tri(C, diag = TRUE)] <- 1
  bd <- apply(C, 2, prod)
  position <- bd == 1
  seqQTN <- seqQTN[position]
  seqQTN <- seqQTN[!is.na(seqQTN)]

  return(seqQTN)
}
