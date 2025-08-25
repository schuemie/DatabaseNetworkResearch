# Use multivariate normal distribution to generalize old systematic error model to multiple
# databases.

library(mvtnorm)
library(dplyr)

fitSystematicErrorModel <- function(data) {
  # Ensure database IDs are factors to maintain a consistent order
  databaseIds <- unique(data$databaseId)
  data$databaseId <- factor(data$databaseId, levels = databaseIds)
  nDatabases <- length(databaseIds)

  # Split data into a list, where each element contains the data for one control
  controlDataList <- split(data, data$outcomeId)

  # Filter out controls that don't have data in all databases
  controlDataList <- controlDataList[sapply(controlDataList, nrow) == nDatabases]

  if (length(controlDataList) == 0) {
    stop("No controls found with data in all specified databases.")
  }

  calculateNegativeLogLikelihood <- function(params, dataList) {
    # First nDatabases elements of params are the means
    meanVector <- params[1:nDatabases]

    # Remaining elements of params are for the Cholesky matrix
    choleskyParams <- params[(nDatabases + 1):length(params)]
    choleskyMatrix <- matrix(0, nrow = nDatabases, ncol = nDatabases)
    choleskyMatrix[lower.tri(choleskyMatrix, diag = TRUE)] <- choleskyParams

    # Reconstruct the covariance matrix: Sigma = L * L^T
    covarianceMatrix <- choleskyMatrix %*% t(choleskyMatrix)

    # Calculate the log-likelihood for each negative control and sum them up
    logLikelihoods <- sapply(dataList, function(controlData) {
      # Order data by database ID to ensure vector elements match
      controlData <- controlData[order(controlData$databaseId), ]

      # y_j is the vector of observed log RRs for control j
      y <- controlData$logRr

      # D_j is the diagonal matrix of squared standard errors (variances)
      d <- diag(controlData$seLogRr^2)

      # The covariance of the observed y_j is Sigma + D_j
      observedCovariance <- covarianceMatrix + d

      # Calculate the log-density of the multivariate normal distribution
      ll <- tryCatch({
        dmvnorm(y, mean = meanVector, sigma = observedCovariance, log = TRUE)
      }, error = function(e) {
        # Return a very small number if the covariance matrix is not invertible
        return(-1e10)
      })
      return(ll)
    })
    return(-sum(logLikelihoods))
  }
  initialMean <- rep(0, nDatabases)
  initialCholesky <- t(chol(diag(nDatabases) * 0.1))
  initialCholeskyParams <- initialCholesky[lower.tri(initialCholesky, diag = TRUE)]
  initialParams <- c(initialMean, initialCholeskyParams)
  optimizationResult <- optim(
    par = initialParams,
    fn = calculateNegativeLogLikelihood,
    dataList = controlDataList,
    method = "BFGS")
  finalParams <- optimizationResult$par
  finalMean <- finalParams[1:nDatabases]
  names(finalMean) <- databaseIds

  finalCholeskyParams <- finalParams[(nDatabases + 1):length(finalParams)]
  finalCholesky <- matrix(0, nrow = nDatabases, ncol = nDatabases)
  finalCholesky[lower.tri(finalCholesky, diag = TRUE)] <- finalCholeskyParams
  finalCovarianceMatrix <- finalCholesky %*% t(finalCholesky)

  rownames(finalCovarianceMatrix) <- databaseIds
  colnames(finalCovarianceMatrix) <- databaseIds
  return(list(
    mean = finalMean,
    covarianceMatrix = finalCovarianceMatrix,
    optimizationResult = optimizationResult
  ))
}

# Simulate data ------------------------------------------------------------------------------------
set.seed(123)
nControls <- 100
databaseIds <- c("DB1", "DB2", "DB3")
nDatabases <- length(databaseIds)
trueMean <- c(0.1, -0.05, 0.0)
trueCovariance <- matrix(c(0.1, 0.05, 0.02,
                           0.05, 0.2, 0.08,
                           0.02, 0.08, 0.15), nrow = nDatabases)
trueSystematicError <- rmvnorm(nControls, mean = rep(0, nDatabases), sigma = trueCovariance)
simulatedData <- data.frame()
for (j in 1:nControls) {
  for (i in 1:nDatabases) {
    se <- runif(1, 0.1, 0.4)
    randomError <- rnorm(1, 0, se)
    observedLogRr <- trueMean[i] + trueSystematicError[j, i] + randomError
    simulatedData <- rbind(simulatedData, data.frame(
      outcomeId = j,
      databaseId = databaseIds[i],
      logRr = observedLogRr,
      seLogRr = se
    ))
  }
}
print("True Mean Vector:")
print(trueMean)
print("True Covariance Matrix:")
print(trueCovariance)
data <- simulatedData

# Real data ----------------------------------------------------------------------------------------
ncEstimates <- readRDS("ncEstimates.rds")
targetName <- "DPP4I"
comparatorName <- "GLP1RA"
subset <- ncEstimates |>
  filter(analysisId == 7, targetName == !!targetName, comparatorName == !!comparatorName) |>
  select(databaseId, outcomeId, logRr, seLogRr)
data <- subset


# Fit model and output results ---------------------------------------------------------------------
fitResult <- fitSystematicErrorModel(data)
print("Estimated Mean Vector:")
print(fitResult$mean)
print("Estimated Covariance Matrix:")
print(fitResult$covarianceMatrix)

correlationMatrix <- cov2cor(fitResult$covarianceMatrix)
print(correlationMatrix)


# Compare database bias distributions from full model to per-database models -----------------------
groups <- data |>
  group_by(databaseId) |>
  group_split()
fun <- function(group) {
  null <- EmpiricalCalibration::fitNull(group$logRr, group$seLogRr)
  return(tibble(databaseId = group$databaseId[1],
                mean = null[1],
                variance = null[2]^2))
}
bind_rows(lapply(groups, fun))

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
      ll <- dmvnorm(yPrime, mean = rep(as.numeric(muHat), nDatabases), sigma = V, log = TRUE)
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

calibratedEstimates <- data |>
  group_by(outcomeId) |>
  group_map(~ calibrateCiRandomEffects(fitResult, .x)) |>
  bind_rows()

# Type 1 error:
calibratedEstimates |>
  mutate(significant = ciLower > 1 | ciUpper < 1) |>
  summarise(mean(significant))
# mean(significant)
# 1        0.08235294

# Geometric mean precision:
exp(mean(log(1 / calibratedEstimates$seLogRr)))
# [1] 2.742278


# Compare to current procedure ---------------------------------------------------------------------
performMetaAnalysis <- function(group) {
  meta <- meta::metagen(group$logRr, group$seLogRr, sm = "RR", iterations = 300)
  s <- summary(meta)
  rnd <- s$random
  return(
    tibble(
      estimate = exp(rnd$TE),
      ciLower = exp(rnd$lower),
      ciUpper = exp(rnd$upper),
      i2 = s$I2,
      nDatabases = nrow(group),
      logRr = rnd$TE,
      seLogRr = rnd$seTE
    )
  )
}

groups <- data |>
  group_split(outcomeId)
maEstimates <- lapply(groups, performMetaAnalysis) |>
  bind_rows()
null <- EmpiricalCalibration::fitNull(maEstimates$logRr, maEstimates$seLogRr)
calibratedEstimates2 <- EmpiricalCalibration::calibrateConfidenceInterval(
  maEstimates$logRr,
  maEstimates$seLogRr,
  EmpiricalCalibration::convertNullToErrorModel(null)
)
# Type 1 error:
calibratedEstimates2 |>
  mutate(significant = logLb95Rr > 0 | logUb95Rr < 0) |>
  summarise(mean(significant))
# mean(significant)
# 1        0.07058824

# Geometric mean precision:
exp(mean(log(1 / calibratedEstimates2$seLogRr)))
# [1] 1.858919

head(calibratedEstimates)
head(calibratedEstimates2)

ggplot(tibble(x = calibratedEstimates2$logRr, y = calibratedEstimates$logRr), aes(x = x, y = y)) +
  geom_abline(slope = 1) +
  geom_point() +
  scale_x_continuous("Old calibration") +
  scale_y_continuous("New calibration")

ggplot(tibble(x = calibratedEstimates2$seLogRr, y = calibratedEstimates$seLogRr), aes(x = x, y = y)) +
  geom_abline(slope = 1) +
  geom_point() +
  scale_x_continuous("Old calibration") +
  scale_y_continuous("New calibration")

# Compute tau distribution and 'posteriors' ----------------------------------------------------------
legendLabel <- "Htn"
minDatabases <- 6
estimates <- readRDS(sprintf("estimates_%s.rds", legendLabel))

atLeastNdbs <- estimates |>
  group_by(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl) |>
  summarise(nDatabases = n(), .groups = "drop") |>
  filter(nDatabases >= minDatabases)

groups <- estimates |>
  inner_join(
    atLeastNdbs,
    by = join_by(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl)
  ) |>
  group_by(targetId, targetName, comparatorId, comparatorName, analysisName)|>
  group_split()

# group = groups[[1]]
computeTau <- function(group) {
  ncs <- group |>
    filter(negativeControl)
  fitResult <- fitSystematicErrorModel(ncs)
  outcomeGroups <- group |>
    group_by(outcomeId) |>
    group_split()

  tauSamplesNcs <- list()
  tauSamplesHois <- list()
  rows <- list()
  # outcomeGroup = outcomeGroups[[1]]
  for (i in seq_along(outcomeGroups)) {
    outcomeGroup = outcomeGroups[[i]]
    keyRow <- outcomeGroup |>
      head(1) |>
      select(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl)
    maEstimate <- calibrateCiRandomEffects(fitResult, outcomeGroup)
    row <- keyRow |>
      mutate(tau = sqrt(maEstimate$tau2),
             tau95Lb = sqrt(maEstimate$tau2CiLower),
             tau95Ub = sqrt(maEstimate$tau2CiUpper))

    rows[[i]] <- row
    tauSample <- sqrt(exp(rnorm(10, log(maEstimate$tau2), maEstimate$seLogTau2)))
    if (keyRow$negativeControl) {
      tauSamplesNcs[[length(tauSamplesNcs) + 1]] <- tauSample
    } else {
      tauSamplesHois[[length(tauSamplesHois) + 1]] <- tauSample
    }
  }
  return(list(row = bind_rows(rows),
              tauSampleNcs = do.call(c, tauSamplesNcs),
              tauSampleHois = do.call(c, tauSamplesHois)))
}

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "mvtnorm")
ParallelLogger::clusterRequire(cluster, "dplyr")
snow::clusterExport(cluster, c("fitSystematicErrorModel", "calibrateCiRandomEffects"))
tauRowsAndSamples <- ParallelLogger::clusterApply(cluster, groups, computeTau)
ParallelLogger::stopCluster(cluster)

taus <- bind_rows(lapply(tauRowsAndSamples, function(x) x$row))
tauSamples <- list()
for (i in seq_along(tauRowsAndSamples)) {
  key <- sprintf("%s-%s", tauRowsAndSamples[[i]]$row$analysisName[1], TRUE)
  tauSamples[[key]] <- c(tauSamples[[key]], tauRowsAndSamples[[i]]$tauSampleNcs)
  key <- sprintf("%s-%s", tauRowsAndSamples[[i]]$row$analysisName[1], FALSE)
  tauSamples[[key]] <- c(tauSamples[[key]], tauRowsAndSamples[[i]]$tauSampleHois)
}





#
#
# calibratePvalueMultiDb <- function(fitResult, newData) {
#   databaseIds <- names(fitResult$mean)
#   nDatabases <- length(databaseIds)
#
#   # Order the new data to match the model's structure
#   newData <- newData[match(databaseIds, newData$databaseId), ]
#
#   if (any(is.na(newData$logRr))) {
#     stop("New data is missing for one or more databases in the model.")
#   }
#
#   # 1. Assemble vectors and matrices
#   y <- newData$logRr
#   mu <- fitResult$mean
#   sigmaSys <- fitResult$covarianceMatrix
#   d <- diag(newData$seLogRr^2)
#
#   # 2. Calculate total covariance
#   sigmaTotal <- sigmaSys + d
#
#   # 3. Calculate the squared Mahalanobis distance
#   # (y - mu) %*% solve(sigmaTotal) %*% (y - mu)
#   diff <- y - mu
#   mahalanobisSq <- t(diff) %*% solve(sigmaTotal) %*% diff
#
#   # 4. Calculate p-value from chi-squared distribution
#   # Degrees of freedom = number of databases
#   calibratedP <- pchisq(mahalanobisSq, df = nDatabases, lower.tail = FALSE)
#
#   return(as.numeric(calibratedP))
# }




