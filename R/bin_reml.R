#Evaluate the maximum log-likelihood.
#
#Estimate genetic and error variance components for a phenotype using a set of
#markers as explanatory variables.
#
#@param Y A numeric matrix with one column containing phenotypic values.
#@param CV A numeric matrix containing any user-specified covariates.
#@param GK A numeric matrix containing marker scores.
#
#@return A list containing the following components: \describe{
#  \item{\code{LL}}{A numeric scalar. The -2*LL estimated by REML.}
#  \item{\code{vg}}{A numeric scalar. The estimated genetic variance.}
#  \item{\code{ve}}{A numeric scalar. The estimated error variance.}}
#
#@author Xiaolei Liu, Jiabo Wang, Qishan Wang, Feng Tian, Zhiwu Zhang, Aaron Kusmec

bin_reml <- compiler::cmpfun(function(Y = NULL, CV = NULL, GK = NULL) {
  # Add a column for the mean to CV
  if (!is.null(CV)) {
    CV <- as.matrix(CV)
    X <- as.matrix(cbind(rep(1, nrow(CV)), CV))
  } else {
    X <- matrix(rep(1, nrow(Y)), nrow = nrow(Y), ncol = 1)
  }


  # Calculate the log-likelihood ------------------------------------------
  if (!is.null(GK) && any(var(GK) == 0)) {
    deltaExpStart <- 100
    deltaExpEnd <- deltaExpStart
  } else {
    deltaExpStart <- -5
    deltaExpEnd <- 5
  }

  # Singular value decomposition on genotypes
  K.GK.svd <- svd(GK)
  d <- K.GK.svd$d[K.GK.svd$d > 1e-08]
  d <- d^2
  U1 <- K.GK.svd$u[, 1:length(d)]

  # Handle a GK that contains one SNP
  if (is.null(dim(U1))) U1 <- matrix(U1, ncol = 1)

  n <- nrow(U1)
  U1TX <- crossprod(U1, X)
  U1TY <- crossprod(U1, Y)
  yU1TY <- Y - U1 %*% U1TY
  XU1TX <- X - U1 %*% U1TX
  IU <- -tcrossprod(U1, U1)
  diag(IU) <- rep(1, n) + diag(IU)
  IUX <- crossprod(IU, X)
  IUY <- crossprod(IU, Y)

  # Iteration on delta
  count <- 0
  for (m in seq(deltaExpStart, deltaExpEnd, by = 0.1)) {
    count <- count + 1
    delta <- exp(m)

    # Calculate beta -------------------------------------------------------
    beta1 <- 0
    for (i in 1:length(d)) {
      one <- matrix(U1TX[i, ], nrow = 1)
      beta1 <- beta1 + crossprod(one, (one/(d[i] + delta)))
    }

    beta2 <- 0
    for (i in 1:nrow(U1)) {
      one <- matrix(IUX[i, ], nrow = 1)
      beta2 <- beta2 + crossprod(one, one)
    }
    beta2 <- beta2/delta

    beta3 <- 0
    for (i in 1:length(d)) {
      one1 <- matrix(U1TX[i, ], nrow = 1)
      one2 <- matrix(U1TY[i, ], nrow = 1)
      beta3 <- beta3 + crossprod(one1, (one2/(d[i] + delta)))
    }

    beta4 <- 0
    for (i in 1:nrow(U1)) {
      one1 <- matrix(IUX[i, ], nrow = 1)
      one2 <- matrix(IUY[i, ], nrow = 1)
      beta4 <- beta4 + crossprod(one1, one2)
    }
    beta4 <- beta4/delta

    zw1 <- try(solve(beta1 + beta2), silent = TRUE)
    if (inherits(zw1, "try-error")) {
      zw1 <- MASS::ginv(beta1 + beta2)
    }

    zw2 <- beta3 + beta4
    beta <- crossprod(zw1, zw2)

    # Calculate LL ----------------------------------------------------------
    part11 <- n*log(2*3.14)
    part12 <- 0
    for (i in 1:length(d)) {
      part12 <- part12 + log(d[i] + delta)
    }
    part13 <- (nrow(U1) - length(d))*log(delta)
    part1 <- -(1/2)*(part11 + part12 + part13)

    part21 <- nrow(U1)
    part221 <- 0
    for (i in 1:length(d)) {
      one1 <- matrix(U1TX[i, ], nrow = 1)
      one2 <- matrix(U1TY[i, ], nrow = 1)
      part221 <- part221 + ((one2 - one1 %*% beta)^2)/(d[i] + delta)
    }

    part222 <- 0
    for (i in 1:n) {
      one1 <- matrix(XU1TX[i, ], nrow = 1)
      one2 <- matrix(yU1TY[i, ], nrow = 1)
      part222 <- part222 + ((one2 - one1 %*% beta)^2)/delta
    }
    part22 <- n*log((1/n)*(part221 + part222))
    part2 <- -(1/2)*(part21 + part22)

    LL <- part1 + part2
    part1 <- 0
    part2 <- 0

    # Save the optimum ----------------------------------------------------
    if (count == 1) {
      beta.save <- beta
      delta.save <- delta
      LL.save <- LL
    } else {
      if (LL > LL.save) {
        beta.save <- beta
        delta.save <- delta
        LL.save <- LL
      }
    }
  }

  names(delta.save) <- NULL
  names(LL.save) <- NULL

  # Calculate Vg and Ve ---------------------------------------------------
  sigma_a1 <- 0
  for (i in 1:length(d)) {
    one1 <- matrix(U1TX[i, ], nrow = 1)
    one2 <- matrix(U1TY[i, ], nrow = 1)
    sigma_a1 <- sigma_a1 + ((one2 - one1 %*% beta.save)^2)/(d[i] + delta.save)
  }

  sigma_a2 <- 0
  for (i in 1:nrow(U1)) {
    one1 <- matrix(IUX[i, ], nrow = 1)
    one2 <- matrix(IUY[i, ], nrow = 1)
    sigma_a2 <- sigma_a2 + (one2 - one1 %*% beta.save)^2
  }
  sigma_a2 <- sigma_a2/delta.save
  sigma_a <- (1/n)*(sigma_a1 + sigma_a2)

  sigma_e <- delta.save*sigma_a


  return(list(LL = -2*LL.save, vg = sigma_a, ve = sigma_e))
})
