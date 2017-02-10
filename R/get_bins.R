#' Bin markers and select representative markers based on p-values.
#'
#' Create bins of an arbitrary size and select an arbitrary number of
#' representative markers across bins with the smallest p-values.
#'
#' @param GI A data frame with three columns: marker name, chromosome, and
#'   base-pair position.
#' @param GP A data frame with four columns: marker name, chromosome, base-pair
#'   position, and p-value.
#' @param bin.size An integer scalar. The size (in base-pairs) of each bin.
#' @param nqtn An integer scalar. The number of pseudo-QTNs to select.
#' @param MaxBP An integer scalar. A scaling factor used to generate unique
#'   marker IDs.
#'
#' @return A logical vector indicating whether or not a marker is a pseudo-QTN.
#'
#' @details Markers are binned into bins of arbitrary size and represented by
#'   the marker with the smallest p-value within each bin.
#'
#' @author Zhiwu Zhang
get_bins <- function(GI = NULL, GP = NULL, bin.size = 1e7, nqtn = NULL,
                     MaxBP = 1e10) {

  if (is.null(GP)) return(NULL)


  # Bin SNPs --------------------------------------------------------------
  # Create unique SNP and bin IDs
  ID.GP <- as.numeric(GP[, 3]) + as.numeric(GP[, 2])*MaxBP
  bin.GP <- floor(ID.GP/bin.size)

  # Create a table with bin ID, SNP ID, and p-value. Columns 2 and 3 are used as
  # indicator variables
  binP <- data.frame(bin.GP = bin.GP, X1 = NA, X2 = NA, ID.GP = ID.GP,
                     pvalue = as.numeric(GP[, 4]))
  n <- nrow(binP)

  # Sort by p-value within bin ID
  binP <- binP[order(binP$bin.GP, binP$pvalue, decreasing = FALSE), ]

  # Set the indicator columns
  binP[, "X1"] <- c(0, binP[1:(n - 1), "bin.GP"])
  binP[, "X2"] <- binP[, "bin.GP"] - binP[, "X1"]

  # ???
  ID.GP <- binP[binP[, "X2"] > 0, ]


  # Choose the most influential bins as pseudo-QTNs -----------------------
  ID.GP <- ID.GP[order(ID.GP[, "pvalue"], decreasing = FALSE), ]

  index <- !is.na(ID.GP[, "ID.GP"])
  ID.GP <- ID.GP[index, "ID.GP"]
  num_bins <- length(ID.GP)
  if (!is.null(nqtn)) {
    if (!is.na(nqtn)) {
      available <- min(nqtn, num_bins)
      if (available == 0) {
        ID.GP <- -1
      } else {
        ID.GP <- ID.GP[1:available] #keep the top ones selected
      }
    }
  }

  # Create an index in GI
  theIndex <- NULL
  if (!is.null(GI)) {
    ID.GI <- GI[, 3] + GI[, 2]*MaxBP
    theIndex <- ID.GI %in% ID.GP
  }

  return(theIndex)
}
