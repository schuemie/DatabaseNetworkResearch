# Code for computing correlation between databases based on their estimates, using a frequentist
# approach. This doesn't really seem to work.

library(mvtnorm)
library(dplyr)
library(tidyr)

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
      dmvnorm(logRrSubset, mean = meanVectorSubset, sigma = totalCovarianceSubset, log = TRUE)
    }, error = function(e) {
      # Return a very small number if the covariance matrix is not invertible
      return(-1e10)
    })
    totalLogLikelihood <- totalLogLikelihood + logLikelihood
  }
  # print(paste(totalLogLikelihood, paste(params, collapse = ",")))
  print(totalLogLikelihood)
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

  finalCorrelationMatrix <- cov2cor(finalCovarianceMatrix)
  return(list(
    mean = finalMean,
    covarianceMatrix = finalCovarianceMatrix,
    correlationMatrix = finalCorrelationMatrix
  ))
}
