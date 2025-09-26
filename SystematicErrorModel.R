# Use multivariate normal distribution to generalize old systematic error model to multiple
# databases.

library(mvtnorm)
library(dplyr)
library(tidyr)

fitSystematicErrorModel <- function(data) {

  constructCovarianceMatrix <- function(params, nDatabases) {
    choleskyParams <- params[(nDatabases + 1):length(params)]
    choleskyMatrix <- matrix(0, nrow = nDatabases, ncol = nDatabases)
    choleskyMatrix[lower.tri(choleskyMatrix, diag = TRUE)] <- choleskyParams
    covarianceMatrix <- choleskyMatrix %*% t(choleskyMatrix)
    return(covarianceMatrix)
  }

  computeNegativeLogLikelihood <- function(params, logRrMatrix, seLogRrMatrix, nDatabases) {
    meanVector <- params[1:nDatabases]
    covarianceMatrix <- constructCovarianceMatrix(params, nDatabases)
    totalLogLikelihood <- 0
    for (i in 1:nrow(logRrMatrix)) {
      logRr <- logRrMatrix[i, ]
      seLogRr <- seLogRrMatrix[i, ]

      # Check for missing data for this outcome
      validIndices <- !is.na(seLogRr)

      if (sum(validIndices) < 2) {
        next # Cannot compute likelihood for a single point
      }
      # Subset matrices to non-missing data
      logRrSubset <- logRr[validIndices]
      seLogRrSubset <- seLogRr[validIndices]
      meanVectorSubset <- meanVector[validIndices]
      covarianceMatrixSubset <- covarianceMatrix[validIndices, validIndices]

      # Within-database variance matrix (measurement error)
      varianceMatrixSubset <- diag(seLogRrSubset^2)

      # Total covariance for this outcome
      totalCovarianceSubset <- covarianceMatrixSubset + varianceMatrixSubset

      logLikelihood <- tryCatch({
        mvtnorm::dmvnorm(logRrSubset, mean = meanVectorSubset, sigma = totalCovarianceSubset, log = TRUE)
      }, error = function(e) {
        # Return a very small number if the covariance matrix is not invertible
        return(-1e10)
      })
      totalLogLikelihood <- totalLogLikelihood + logLikelihood
    }
    return(-totalLogLikelihood)
  }

  # Get database IDs and count
  databaseIds <- unique(data$databaseId)
  nDatabases <- length(databaseIds)

  dataWide <- data %>%
    pivot_wider(
      id_cols = outcomeId,
      names_from = databaseId,
      values_from = c(logRr, seLogRr),
      names_sort = TRUE
    )

  databaseIds <- sort(databaseIds)
  logRrCols <- paste0("logRr_", databaseIds)
  seLogRrCols <- paste0("seLogRr_", databaseIds)
  logRrMatrix <- as.matrix(dataWide[, logRrCols])
  seLogRrMatrix <- as.matrix(dataWide[, seLogRrCols])

  initialMean <- rep(0, nDatabases)
  initialCholesky <- t(chol(diag(nDatabases) * 0.1))
  initialCholeskyParams <- initialCholesky[lower.tri(initialCholesky, diag = TRUE)]
  initialParams <- c(initialMean, initialCholeskyParams)
  optimResult <- optim(
    par = initialParams,
    fn = computeNegativeLogLikelihood,
    logRrMatrix = logRrMatrix,
    seLogRrMatrix = seLogRrMatrix,
    nDatabases = nDatabases,
    method = "L-BFGS-B",
    control = list(maxit = 2000)
  )
  if (optimResult$convergence != 0) {
    warning(sprintf("Optim failed to converge with covergence = %s", optimResult$convergence))
  }
  finalParams <- optimResult$par
  finalMean <- finalParams[1:nDatabases]
  names(finalMean) <- databaseIds

  finalCovarianceMatrix <- constructCovarianceMatrix(finalParams, nDatabases)
  rownames(finalCovarianceMatrix) <- databaseIds
  colnames(finalCovarianceMatrix) <- databaseIds
  return(list(
    mean = finalMean,
    covarianceMatrix = finalCovarianceMatrix
  ))
}

calibrateCiRandomEffects <- function(fitResult, newData) {
  modelDbIds <- names(fitResult$mean)
  availableDbIds <- intersect(modelDbIds, newData$databaseId)

  if (length(availableDbIds) == 0) {
    stop("None of the databases in newData are present in the systematic error model.")
  }

  # Subset the model parameters and new data to the common set of databases
  meanSys <- fitResult$mean[availableDbIds]
  covSys <- fitResult$covarianceMatrix[availableDbIds, availableDbIds, drop = FALSE]
  newData <- newData[newData$databaseId %in% availableDbIds, ]
  newData <- newData[match(availableDbIds, newData$databaseId), ] # Ensure order

  nDatabases <- length(availableDbIds)

  if (nDatabases == 1) {
    # --- Case 1: Only one database estimate is available ---
    # Heterogeneity (tau^2) is not applicable.
    # The model simplifies to a fixed-effect calibration.

    # De-bias the estimate by subtracting the database-specific mean systematic error.
    muHatFinal <- newData$logRr - meanSys

    # The total variance is the sum of the systematic error variance (for this DB)
    # and the random error variance (from the new estimate).
    muVarFinal <- covSys[1, 1] + newData$seLogRr^2

    # Compute 95% CI and exponentiate to HR scale
    ciLowerLog <- muHatFinal - 1.96 * sqrt(muVarFinal)
    ciUpperLog <- muHatFinal + 1.96 * sqrt(muVarFinal)

    return(data.frame(
      estimate = exp(as.numeric(muHatFinal)),
      ciLower = exp(as.numeric(ciLowerLog)),
      ciUpper = exp(as.numeric(ciUpperLog)),
      tau2 = NA,
      nDatabases = nDatabases,
      logRr = as.numeric(muHatFinal),
      seLogRr = sqrt(muVarFinal)
    ))

  } else {
    # --- Case 2: More than one database estimate is available ---
    # Proceed with the random-effects meta-analysis.

    # De-bias the vector of estimates by subtracting the mean systematic error vector.
    # This gives us y', the vector of estimates cleared of their average bias.
    yPrime <- newData$logRr - meanSys

    # Calculate the covariance matrix of the OBSERVED, DE-BIASED data.
    # This combines the correlated systematic error with the independent random error.
    sigmaTotal <- covSys + diag(newData$seLogRr^2)
    ones <- rep(1, nDatabases)

    # Profile the likelihood to find the MLE of tau^2 (between-database heterogeneity).
    # This function calculates the negative log-likelihood for a given tau^2,
    # having already integrated out the mean effect (mu).
    negLogLike <- function(logTau2) {
      tau2 <- exp(logTau2)
      if (is.infinite(tau2)) return(1e10)

      # *** THIS IS THE KEY STEP ***
      # V is the final covariance matrix of the de-biased data (yPrime).
      # It combines three components:
      # 1. covSys: Correlated systematic error (off-diagonal).
      # 2. diag(seLogRr^2): Independent random error (diagonal).
      # 3. diag(tau2): Independent between-database heterogeneity (diagonal).
      # By adding tau^2 as a diagonal matrix, we correctly model the true effects
      # as being independent, while allowing the systematic errors to be correlated.
      V <- sigmaTotal + diag(tau2, nDatabases)
      VInv <- try(solve(V), silent = TRUE)
      if (inherits(VInv, "try-error")) return(1e10)

      # Calculate the weighted average for this value of tau^2
      muHat <- (t(ones) %*% VInv %*% yPrime) / (t(ones) %*% VInv %*% ones)

      # Calculate the likelihood of observing yPrime given this mu and V
      ll <- mvtnorm::dmvnorm(yPrime, mean = rep(as.numeric(muHat), nDatabases), sigma = V, log = TRUE)
      return(-ll)
    }

    # Find the value of tau^2 that maximizes the likelihood.
    opt <- optim(par = log(0.01), fn = negLogLike, method = "BFGS", hessian = TRUE)
    estimatedTau2 <- exp(opt$par)

    # Use the Hessian matrix to approximate the variance of log(tau^2).
    # The variance is the inverse of the Fisher information, approximated by the Hessian.
    fisherInfo <- opt$hessian[1, 1]
    seLogTau2 <- NA
    tau2CiLower <- NA
    tau2CiUpper <- NA
    if (fisherInfo > 0) { # Check if Hessian is positive definite
      seLogTau2 <- sqrt(1 / fisherInfo)
      ciLowerLogTau2 <- opt$par - 1.96 * seLogTau2
      ciUpperLogTau2 <- opt$par + 1.96 * seLogTau2
      tau2CiLower <- exp(ciLowerLogTau2)
      tau2CiUpper <- exp(ciUpperLogTau2)
    }

    # Now, compute the final pooled estimate and its variance using the estimated tau^2.
    VFinal <- sigmaTotal + diag(estimatedTau2, nDatabases)
    VFinalInv <- solve(VFinal)

    # The final estimate is a weighted average, where the weights are derived
    # from the inverse of the full final covariance matrix (VFinalInv).
    # This correctly accounts for variance AND covariance.
    muHatFinal <- (t(ones) %*% VFinalInv %*% yPrime) / (t(ones) %*% VFinalInv %*% ones)

    # The variance of the final estimate is also derived from the inverse covariance matrix.
    muVarFinal <- 1 / (t(ones) %*% VFinalInv %*% ones)

    # Compute 95% CI and exponentiate to HR scale
    ciLowerLog <- muHatFinal - 1.96 * sqrt(muVarFinal)
    ciUpperLog <- muHatFinal + 1.96 * sqrt(muVarFinal)

    return(data.frame(
      estimate = exp(as.numeric(muHatFinal)),
      ciLower = exp(as.numeric(ciLowerLog)),
      ciUpper = exp(as.numeric(ciUpperLog)),
      tau2 = estimatedTau2,
      seLogTau2 = seLogTau2,
      tau2CiLower = tau2CiLower,
      tau2CiUpper = tau2CiUpper,
      nDatabases = nDatabases,
      logRr = as.numeric(muHatFinal),
      seLogRr = sqrt(muVarFinal)
    ))
  }
}
