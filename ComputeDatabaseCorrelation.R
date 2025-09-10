# Code for computing correlation between databases based on their estimates, using a frequentist
# approach. This doesn't really seem to work.

library(mvtnorm)
library(dplyr)
library(tidyr)

reconstructMatrices <- function(params, nDatabases) {
  nCorrelationParams <- nDatabases * (nDatabases - 1) / 2

  # Reconstruct the correlation matrix
  offDiagonals <- tanh(params[1:nCorrelationParams])
  correlationMatrix <- diag(nDatabases)
  correlationMatrix[lower.tri(correlationMatrix)] <- offDiagonals
  correlationMatrix <- correlationMatrix + t(correlationMatrix) - diag(nDatabases)

  # Reconstruct the standard deviations (taus)
  logTaus <- params[(nCorrelationParams + 1):(nCorrelationParams + nDatabases)]
  taus <- exp(logTaus)

  return(list(
    correlationMatrix = correlationMatrix,
    taus = taus
  ))
}

computeNegativeLogLikelihood <- function(params, logRrMatrix, seLogRrMatrix, nDatabases) {
  # p <<- params
  # params <- p
  reconstructed <- reconstructMatrices(params, nDatabases)
  correlationMatrix <- reconstructed$correlationMatrix
  taus <- reconstructed$taus

  # A correlation matrix must be positive definite. The Cholesky decomposition
  # will fail if it is not. We use try() to catch this failure.
  cholTest <- try(chol(correlationMatrix), silent = TRUE)
  if (inherits(cholTest, "try-error")) {
    return(1e9)
  }

  # Between-database covariance matrix (of true effects)
  covarianceMatrix <- diag(taus) %*% correlationMatrix %*% diag(taus)

  totalLogLikelihood <- 0

  # Iterate over each outcome (row)
  for (i in 1:nrow(logRrMatrix)) {
    y <- logRrMatrix[i, ]
    se <- seLogRrMatrix[i, ]

    # Check for missing data for this outcome
    validIndices <- !is.na(y)
    if (sum(validIndices) < 2) {
      next # Cannot compute likelihood for a single point
    }

    # Subset matrices to non-missing data
    ySub <- y[validIndices]
    seSub <- se[validIndices]
    covSub <- covarianceMatrix[validIndices, validIndices]

    # Within-database variance matrix (measurement error)
    varianceMatrix <- diag(seSub^2)

    # Total covariance for this outcome
    totalCovariance <- covSub + varianceMatrix

    logLikelihood <- tryCatch({
      invTotalCov <- solve(totalCovariance)
      muHat <- sum(invTotalCov %*% ySub) / sum(invTotalCov)
      dmvnorm(ySub, mean = rep(muHat, length(ySub)), sigma = totalCovariance, log = TRUE)

    }, error = function(e) {
      print("error", e$message)
      -1e9
    })
    # print(paste(i, logLikelihood))
    totalLogLikelihood <- totalLogLikelihood + logLikelihood
  }

  # print(paste(totalLogLikelihood, paste(params, collapse = ",")))
  print(totalLogLikelihood)
  # Return the negative log-likelihood for minimization
  return(-totalLogLikelihood)
}

estimateCorrelationMatrix <- function(data) {
  # Get database IDs and count
  databaseIds <- unique(data$databaseId)
  nDatabases <- length(databaseIds)

  dataWide <- data %>%
    pivot_wider(
      id_cols = outcomeId,
      names_from = databaseId,
      values_from = c(logRr, seLogRr),
      names_sort = TRUE # Ensure consistent column order
    )

  databaseIds <- sort(databaseIds)
  logRrCols <- paste0("logRr_", databaseIds)
  seLogRrCols <- paste0("seLogRr_", databaseIds)

  logRrMatrix <- as.matrix(dataWide[, logRrCols])
  seLogRrMatrix <- as.matrix(dataWide[, seLogRrCols])

  # Set initial parameter values for optimization
  nCorrelationParams <- nDatabases * (nDatabases - 1) / 2
  initialCorrelationParams <- rep(0, nCorrelationParams)
  initialLogTaus <- rep(log(0.1), nDatabases)
  initialParams <- c(initialCorrelationParams, initialLogTaus)

  # Run the optimization
  message("Starting likelihood optimization...")
  optimResult <- optim(
    par = initialParams,
    fn = computeNegativeLogLikelihood,
    logRrMatrix = logRrMatrix,
    seLogRrMatrix = seLogRrMatrix,
    nDatabases = nDatabases,
    # method = "BFGS",
     method = "L-BFGS-B",
    control = list(maxit = 2000)
  )

  optimResult <- nlm(
    p = initialParams,
    f = computeNegativeLogLikelihood,
    iterlim = 5000,      # Increased iteration limit
    gradtol = 1e-8,      # Stricter gradient tolerance
    steptol = 1e-8,
    logRrMatrix = logRrMatrix,
    seLogRrMatrix = seLogRrMatrix,
    nDatabases = nDatabases
  )


  if (optimResult$convergence != 0) {
    warning("Optimization may not have converged. Code: ", optimResult$convergence)
  }

  # Reconstruct the final matrices from optimized parameters
  finalParams <- optimResult$par
  reconstructedFinal <- reconstructMatrices(finalParams, nDatabases)

  estimatedCorrelation <- reconstructedFinal$correlationMatrix
  colnames(estimatedCorrelation) <- databaseIds
  rownames(estimatedCorrelation) <- databaseIds

  estimatedTaus <- reconstructedFinal$taus
  names(estimatedTaus) <- databaseIds

  message("Optimization finished.")

  return(list(
    correlationMatrix = estimatedCorrelation,
    standardDeviations = estimatedTaus
  ))
}

