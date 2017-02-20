#' Select pseudo-QTNs.
#'
#' Use p-values estimated by \code{\link{quick_lm}} to select pseudo-QTNs for
#' use as covariates in subsequent single-marker regression tests. Pseudo-QTNs
#' may be selected using a set number of bins and a set number of pseudo-QTNs or
#' may be optimized based on a set of optional sizes.
#'
#' @param Y A data frame with two columns. The first column contains taxa names.
#'   The second column contains numeric phenotypic values.
#' @param GM A data frame with three columns: marker IDs, chromosome, and
#'   base-pair positions.
#' @param P A numeric vector of marker p-values.
#' @param method A character string, one of \code{c("static", "optimum")}.
#' @param bin.sizes A vector of possible bin sizes in base-pairs.
#' @param nqtn A vector of possible numbers of pseudo-QTNs to select.
#' @param CV A numeric matrix of user-specified covariates.
#' @param GDP A pointer to a big.matrix object containing marker scores.
#' @param theLoop The current iteration of the overal FarmCPU process. Controls
#'   bin size selection for \code{method = "static"}.
#'
#' @return A logical vector specifying which markers have been chosen as
#'   pseudo-QTNs.
#'
#' @author Zhiwu Zhang
#' @author Aaron Kusmec

get_qtns <- function(Y, GM, P, method, bin.sizes, nqtn, CV = NULL, GDP = NULL,
                     theLoop = NULL) {
  # NULL p-values means we cannot select pseudo-QTNs
  if (is.null(P)) return(NULL)

  # Remove missing values from Y and set the upper bound for number of
  # pseudo-QTNs
  missing <- is.na(unlist(Y[, 2], use.names = FALSE))
  Y <- Y[!missing, ]
  n <- nrow(Y)
  bound <- round(sqrt(n)/sqrt(log10(n)))

  # Keep only numbers of pseudo-QTNs within the upper bound
  nqtn[nqtn > bound] <- bound
  nqtn <- unique(nqtn[nqtn <= bound])

  # If there is only one combination of bin.sizes and nqtn, set the method to
  # static bin selection.
  if (length(bin.sizes)*length(nqtn) == 1 & method == "optimum") {
    cat("Warning: Bin selection method changed to 'static'.\n")
    method <- "static"
  }

  # Static bin selection ----------------------------------------------------
  if (method == "static") {
    if (theLoop >= length(bin.sizes)) {
      bin.sizes <- bin.sizes[1]
    } else {
      bin.sizes <- bin.sizes[length(bin.sizes) - theLoop + 2]
    }

    nqtn <- bound

    cat("Static selection of", nqtn, "QTNs in", bin.sizes, "length bins.\n")

    mySpecify <- get_bins(GI = GM, GP = cbind(GM, P), bin.size = bin.sizes,
                          nqtn = nqtn)
    seqQTN <- which(mySpecify == TRUE)
  }

  # Optimum bin selection ---------------------------------------------------
  if (method == "optimum") {
    # Remove missing values from CV; Y has already had missing data removed
    CV <- CV[!missing, ]

    # Compute -2*LL for each combination of bin size and pseudo-QTN number in
    # parallel
    cat("Optimizing possible QTNs...\n")
    GD_desc <- bigmemory::describe(GDP)
    params <- mapply(function(x, y) c(x, y),
                     rep(bin.sizes, each = length(nqtn)),
                     rep(nqtn, times = length(bin.sizes)),
                     SIMPLIFY = FALSE)
    reml <- parallel::parLapply(X = params, fun = function(x) {
      local_GDP <- bigmemory::attach.big.matrix(GD_desc)
      mySpecify <- get_bins(GI = GM, GP = cbind(GM, P), bin.size = x[1],
                            nqtn = x[2])
      seqQTN <- which(mySpecify == TRUE)
      GK <- bigmemory::as.matrix(bigmemory::deepcopy(local_GDP, rows = !missing, cols = seqQTN))
      myREML <- bin_reml(Y = matrix(Y[, 2], ncol = 1), CV = CV, GK = GK)
      list(res = c(x[1], x[2], myREML$LL, myREML$vg, myREML$ve), seqQTN = seqQTN)
    })

    # Reduce the REML results
    reml.table <- unlist(lapply(reml, function(l) { return(l$res) }))
    reml.table <- matrix(reml.table, ncol = 5, byrow = TRUE)

    # Nicely print the results
    cat(sprintf("%9s %13s %10s %10s %10s \n", "bin.size", "bin.selection",
                "-2LL", "VG", "VE"))
    for (i in 1:nrow(reml.table)) {
      cat(sprintf("%9s %13s %10s %10s %10s \n", reml.table[i, 1],
                  reml.table[i, 2], reml.table[i, 3], reml.table[i, 4],
                  reml.table[i, 5]))
    }

    # Identify the optimal bin size and number of pseudo-QTNs
    opt <- which.min(reml.table[, 3])
    cat("Optimal bin size =", reml.table[opt, 1], "with", reml.table[opt, 2],
        "pseudo-QTNs.\n")
    seqQTN <- reml[[opt]]$seqQTN
  }

  return(seqQTN)
}
