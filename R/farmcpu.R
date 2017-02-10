#' Perform GWAS using the FarmCPU model.
#'
#' @param Y A dataframe with two columns: a character column of sample names and
#'   a numeric column of phenotypic values. This dataframe may contain missing
#'   values.
#' @param GD A \code{big.matrix} object.
#' @param GM A dataframe with three columns: a character column of sample names,
#'   a numeric column of chromosomes, and a numeric column of base-pair
#'   positions.
#' @param CV An optional numeric matrix of covariates.
#' @param GP A numeric vector of p-values to select pseudo-QTNs before
#'   single-marker regression in the first iteration.
#' @param method.sub Character string indicating the method used to substitute
#'   p-values for pseudo-QTNs. One of \code{c('reward', 'mean', 'median',
#'   'penalty', 'onsite')}.
#' @param method.sub.final Character string indicating the method used to
#'   substitute p-values for pseudo-QTNs on the final iteration.
#' @param method.bin Character string indicating the method used to select
#'   pseudo-QTNs. One of \code{c('static', 'optimum')}.
#' @param bin.size A numeric vector specifying bin sizes in base-pairs for
#'   pseudo-QTN selection.
#' @param bin.selection Numeric vector of numbers of pseudo-QTNs to select.
#' @param memo Character string added to file names.
#' @param Prior A dataframe with the same format as \code{GM} with a fourth
#'   numeric column containing prior marker probabilities.
#' @param ncores.glm Integer scalar number of cores to use for single-marker
#'   regression.
#' @param maxLoop Integer scalar maximum number of iterations to perform.
#' @param converge Numeric scalar such that \code{0 < converge < 1}. Controls
#'   the percentage of overlapping pseudo-QTNs between two iterations required
#'   to terminate the GWAS.
#' @param iteration.output Logical scalar. Whether to include the output of each
#'   iteration or not.
#' @param p.threshold Numeric scalar. P-value threshold for inclusion of
#'   pseudo-QTNs in the model on the first iteration. Defaults to a 0.01
#'   Bonferroni correction.
#' @param ncores.reml Integer scalar number of cores to use for pseudo-QTN selection.
#' @param threshold Numeric scalar maximum correlation allowed between
#'   pseudo-QTNs.
#'
#' @return A list. If \code{iteration.output = FALSE}, the list has one element.
#'   Otherwise, the list additionally has one element for each iteration. Each
#'   element has one required element and one optional element. The required
#'   element is a dataframe called \code{GWAS} with seven columns: the marker
#'   ID, the chromosome, the base-pair position, the marker p-value, the effect
#'   estimate, the standard error, and the t-statistic. If there are
#'   user-specified covariates, the results list contains a matrix called
#'   \code{CV} that contains the effect estimates for the covariates.
#'
#' @author Xiaolei Liu
#' @author Zhiwu Zhang
#' @author Aaron Kusmec
#'
#' @export

farmcpu <- function(Y, GD, GM, CV = NULL, GP = NULL, method.sub = "reward",
                    method.sub.final = "reward", method.bin = "static",
                    bin.size = c(5e5, 5e6, 5e7), bin.selection = seq(10, 100, 10),
                    memo = NULL, Prior = NULL, ncores.glm = 1, maxLoop = 10,
                    converge = 1, iteration.output = FALSE, p.threshold = NA,
                    ncores.reml = 1, threshold = 0.7) {

  cat("------------------------------ Welcome to FarmCPU ------------------------------\n")
  cat("Version 1.0.0 02.10.2017\n\n")

  # Argument checking -----------------------------------------------------
  sub.methods <- c("reward", "mean", "median", "penalty", "onsite")
  if (!(method.sub %in% sub.methods)) stop("Invalid value for method.sub.")
  if (!(method.sub.final %in% sub.methods)) stop("Invalid value for method.sub.final.")
  if (!(method.bin %in% c("static", "optimum", "none"))) stop("Invalid value for method.bin.")
  if (nrow(GM) == nrow(GD)) stop("GD is not in column orientation.")
  if (sum(is.na(Y[, 2])) > nrow(Y)) stop("No non-missing values in Y.")

  # Handle optional input -------------------------------------------------
  # Get the number of user-specified covariates
  if (!is.null(CV)) {
    CV <- as.matrix(CV)
    npc <- ncol(CV)
  } else {
    npc <- 0
  }

  # Handle prior probabilities for selecting pseudo-QTNs in the first iteration
  P <- NULL
  if (!is.null(GP)) P <- GP[, 4]

  # Control variables for the main loop -----------------------------------
  name.of.trait <- colnames(Y)[2]
  if (!is.null(memo)) name.of.trait <- paste(memo, name.of.trait, sep = ".")
  if (method.bin == "none") maxLoop <- 1 # Force exit after one loop (GLM)
  seqQTN.save <- 0
  seqQTN.pre <- -1
  theLoop <- 0
  theConverge <- 0
  isDone <- FALSE
  myResults <- vector(mode = "list")

  # Configure parallel processing for single-marker regression tests ------
  RcppParallel::setThreadOptions(numThreads = ncores.glm)

  # Configure parallel processing for bin optimization --------------------
  snow::setDefaultClusterOptions(outfile = "/dev/null")
  cl <- snow::makeMPIcluster(count = ncores.reml)
  invisible(snow::clusterEvalQ(cl, library(bigmemory)))
  invisible(snow::clusterEvalQ(cl, library(MASS)))
  invisible(snow::clusterEvalQ(cl, library(FarmCPUpp)))
  doSNOW::registerDoSNOW(cl)

  # Main loop -------------------------------------------------------------
  while (!isDone) {
    theLoop <- theLoop + 1
    cat(paste("\nCurrent loop:", theLoop, "/", maxLoop, "\n"))

    # Modify output file name for iteration output
    if (iteration.output) {
      spacer <- "0"
      if (theLoop > 9) spacer <- ""
      name.of.trait.2 <- paste0("Iteration_", spacer, theLoop, ".", name.of.trait)
    }

    # Step 1: Set prior marker p-values -----------------------------------
    myPrior <- prior(GM = GM, P = P, Prior = Prior)

    # Step 2: Set bins ----------------------------------------------------
    if (theLoop <= 2) { # Without pseudo-QTNs
      myBin <- get_qtns(Y = Y, GM = GM, P = myPrior,
                                   method = method.bin, bin.sizes = bin.size,
                                   nqtn = bin.selection, CV = CV,
                                   GDP = GD, theLoop = theLoop)
    } else { # Will have pseudo-QTNs after iteration 2
      myBin <- get_qtns(Y = Y, GM = GM, P = myPrior,
                                   method = method.bin, bin.sizes = bin.size,
                                   nqtn = bin.selection, CV = theCV,
                                   GDP = GD, theLoop = theLoop)
    }

    # Step 3: Check for significant pseudo-QTNs --------------------------
    seqQTN <- myBin

    # If no pseudo-QTN is significant enough, set all pseudo-QTNs to NULL
    if (theLoop == 2) {
      if (!is.na(p.threshold)) {
        if (min(myPrior, na.rm = TRUE) > p.threshold) {
          seqQTN <- NULL
          cat("Top SNPs have too small effects. Setting seqQTN to NULL.\n")
        }
      } else {
        # Use a Bonferroni correction at the 0.01 level is no threshold is
        # provided
        if (min(myPrior, na.rm = TRUE) > 0.01/nrow(GM)) {
          seqQTN <- NULL
          cat("Top SNPs have too small effects. Setting seqQTN to NULL.\n")
        }
      }
    }

    # If FarmCPU cannot find any significant pseudo-QTNs, exit.
    if (theLoop == 2 && is.null(seqQTN)) {
      GWAS <- cbind(GM, P, myGLM$P[, 1:3])
      colnames(GWAS) <- c(colnames(GM), "p.value", "estimate", "stderr", "tvalue")
      myIteration <- list(GWAS = GWAS)

      if (npc != 0) {
        betapc <- cbind(c(1:npc), myGLM$betapc)
        colnames(betapc) <- c("CV", "Effect")
        myIteration$CV <- betapc
      }

      myResults[[name.of.trait]] <- myIteration
      break
    }

    # Combine old and new pseudo-QTNs
    if (!is.null(seqQTN.save) && theLoop > 1) {
      if (seqQTN.save[1] != 0 & seqQTN.save[1] != -1 & !is.null(seqQTN)) {
        seqQTN <- union(seqQTN, seqQTN.save)
      }
    }

    # Screen combined pseudo-QTNs for significance
    if (theLoop != 1 & !is.null(seqQTN)) {
      seqQTN.p <- myPrior[seqQTN]

      # Screen using a user-specified threshold, or a default of 0.01
      if (!is.na(p.threshold)) {
        index.p <- seqQTN.p < p.threshold
      } else {
        index.p <- seqQTN.p < 0.01
      }

      # Regardless of significance, keep pseudo-QTNs that have been picked in
      # multiple iterations
      if (theLoop > 2) index.p[seqQTN %in% seqQTN.save] <- TRUE

      seqQTN.p <- seqQTN.p[index.p]
      seqQTN <- seqQTN[index.p]
      seqQTN.p <- seqQTN.p[!is.na(seqQTN)]
      seqQTN <- seqQTN[!is.na(seqQTN)]
    }

    # Step 4: Remove highly correlated pseudo-QTNs --------------------------
    cat("Removing highly correlated pseudo-QTNs.\n")
    seqQTN <- remove_qtns(GDP = GD, GM = GM, seqQTN = seqQTN,
                                     seqQTN.p = seqQTN.p, threshold = threshold)

    # Step 5: Check exit conditions -----------------------------------------
    theConverge <- length(intersect(seqQTN, seqQTN.save))/length(union(seqQTN, seqQTN.save))
    circle <- length(union(seqQTN, seqQTN.pre)) == length(intersect(seqQTN, seqQTN.pre))

    # Handle the first iteration
    if (is.null(seqQTN.pre)) {
      circle <- FALSE
    } else {
      if (seqQTN.pre[1] == 0) circle <- FALSE
      if (seqQTN.pre[1] == -1) circle <- FALSE
    }

    isDone <- (theLoop >= maxLoop) | (theConverge >= converge) | circle

    seqQTN.pre <- seqQTN.save
    seqQTN.save <- seqQTN

    # Step 6: Screen markers using the GLM ----------------------------------
    rm(myBin)
    gc()

    cat("pseudo-QTNs:\n", seqQTN, "\n")
    if (theLoop == maxLoop) {
      cat(paste("Total number of possible QTNs in the model is:", length(seqQTN), "\n"))
    }

    # Add pseudo-QTNs to the covariate matrix
    theCV <- CV
    if (!is.null(seqQTN)) {
      theCV <- cbind(CV, bigmemory::as.matrix(bigmemory::deepcopy(GD, cols = seqQTN)))
    }

    # Screen using the general linear model
    start.time <- proc.time()[3]
    myGLM <- quick_lm(Y = unlist(Y[, 2], use.names = FALSE), GDP = GD,
                                 CV = theCV, npc = npc, seqQTN = seqQTN)
    end.time <- proc.time()[3]
    cat("GLM time:", end.time - start.time, "(s)\n")

    # Step 7: Substitute p-values for pseudo-QTNs --------------------------
    if (!isDone) {
      cat(paste0("Substituting pseudo-QTN p-values with method ", method.sub, ".\n"))
      myGLM <- substitute_pvalue(GM = GM, GLM = myGLM,
                                            QTN = GM[seqQTN, ],
                                            method = method.sub, npc = npc)
    } else {
      cat(paste0("Substituting pseudo-QTN p-values with method ", method.sub.final, ".\n"))
      myGLM <- substitute_pvalue(GM = GM, GLM = myGLM,
                                            QTN = GM[seqQTN, ],
                                            method = method.sub.final, npc = npc)
    }
    P <- myGLM$P[, "pvalue"]
    P[P == 0] <- min(P[P != 0], na.rm = TRUE)*0.01

    # Exit and/or report
    if (isDone | iteration.output) {
      GWAS <- cbind(GM, P, myGLM$P[, 1:3])
      colnames(GWAS) = c(colnames(GM), "p.value", "estimate", "stderr", "tvalue")
      myIteration <- list(GWAS = GWAS)

      if (npc != 0) {
        betapc <- cbind(c(1:npc), myGLM$betapc)
        colnames(betapc) <- c("CV", "Effect")
        myIteration$CV <- betapc
      }

      if (isDone) {
        myResults[[name.of.trait]] <- myIteration
      } else {
        myResults[[name.of.trait.2]] <- myIteration
      }
    }
  }

  # Clean-up and exit
  snow::stopCluster(cl)
  gc()
  cat("*********************** FarmCPU Accomplished Successfully **********************\n")
  return(myResults)
}
