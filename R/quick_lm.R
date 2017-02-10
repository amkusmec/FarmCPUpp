#' Perform single marker regression tests.
#'
#' Returns single marker regression test results for all markers in \code{GDP},
#' effect estimates for covariates (user-specified and pseudo-QTNs), and
#' p-values for pseudo-QTNs used by \code{\link{get_qtns}}.
#'
#' @param Y A numeric vector of phenotype values. May include missing values.
#' @param GDP A pointer to a big.matrix object.
#' @param CV A numeric matrix where each column contains values for a covariate.
#' @param npc An integer scalar. Number of columns of \code{CV} that contain
#'   user-defined covariates.
#' @param seqQTN An integer vector denoting the column in \code{GDP} of each
#'   pseudo-QTN.
#'
#' @return A list containing the following components:
#'   \describe{
#'     \item{\code{P}}{a numeric matrix whose first four columns contain
#'   the effect estimate, standard error, t statistic, and p-value for each
#'   single marker regression test. Subsequent columns contain p-values for
#'   pseudo-QTNs when used as a covariate when testing each
#'   marker.}
#'     \item{\code{P0}}{a numeric vector containing p-values for the mean,
#'   user-specified covariates, and pseudo-QTNs in the null
#'   model.}
#'     \item{\code{betapc}}{a numeric vector of effect estimates for
#'   user-specified covariates.}}
#'
#' @author Aaron Kusmec

quick_lm <- function(Y = NULL, GDP = NULL, CV = NULL, npc = 0, seqQTN = NULL) {

  # Y must be provided
  if (is.null(Y)) return(NULL)

  # GDP or CV must be provided
  if (is.null(GDP) & is.null(CV)) return(NULL)

  # Add a column for the mean to an existing CV or create a design matrix with a
  # mean only if CV does not exist
  if (!is.null(CV)) {
    ccv <- cbind(rep(1, length(Y)), as.matrix(CV))
    nf <- ncol(ccv)
  } else {
    ccv <- matrix(rep(1, length(Y)), ncol = 1)
    nf <- 1
  }

  cat("Number of covariates in current loop:", nf - 1, "\n")
  cat("Scanning...\n")

  # Remove missing data points
  missing <- which(!is.na(Y))
  Y <- Y[missing]
  ccv <- as.matrix(ccv[missing, ])

  # Adjust indices of missing data points for use in C++ code
  missing <- missing - 1

  res <- QuickLM(Y, ccv, GDP@address, missing, npc, length(seqQTN))
  tvalue <- res$coefficients/res$stderr
  pvalue <- 2*pt(abs(tvalue), res$df, lower.tail = FALSE)

  # Reformat the results
  if (length(seqQTN) == 0) {
    P <- cbind(res$coefficients, res$stderr, tvalue, pvalue)
  } else {
    res$seqQTN <- apply(abs(res$seqQTN), 2, pt, df = res$df, lower.tail = FALSE)
    P <- cbind(res$coefficients, res$stderr, tvalue, pvalue, 2*res$seqQTN)
  }

  # Effect estimates, standard errors, t statistics, and p-values for
  # pseudo-QTNs in P will be incorrect due to multicollinearity with marker data
  P[seqQTN, 1:4] <- NA

  # Add dimension names
  if (ncol(P) > 4) {
    colnames(P) <- c("estimate", "stderr", "tvalue", "pvalue",
                    paste0("QTN", 1:length(seqQTN)))
  } else {
    colnames(P) <- c("estimate", "stderr", "tvalue", "pvalue")
  }
  rownames(P) <- colnames(GDP)

  # Calculate effects of covariates and seqQTN in the null model
  fast.lm <- RcppEigen::fastLmPure(y = Y, X = ccv)
  tvalue <- fast.lm$coefficients/fast.lm$se
  P0 <- 2*pt(abs(tvalue), fast.lm$df.residual, lower.tail = FALSE)
  if (npc != 0) {
    betapc <- fast.lm$coefficients[2:(npc + 1)]
    P[seqQTN, "estimate"] <- fast.lm$coefficients[-c(1:(npc + 1))]
    P[seqQTN, "stderr"] <- fast.lm$se[-c(1:(npc + 1))]
    P[seqQTN, "tvalue"] <- tvalue[-c(1:(npc + 1))]
  } else {
    betapc <- NULL
    P[seqQTN, "estimate"] <- fast.lm$coefficients[-1]
    P[seqQTN, "stderr"] <- fast.lm$se[-1]
    P[seqQTN, "tvalue"] <- tvalue[-1]
  }

  myLM <- list(P = P, P0 = P0, betapc = betapc)
  gc()

  return(myLM)
}
