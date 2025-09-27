source("SystematicErrorModel.R")

# Simulate data ------------------------------------------------------------------------------------
set.seed(123)
nControls <- 200
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
data <- simulatedData
model <- fitSystematicErrorModel(data)

print("True Mean Vector:")
print(trueMean)
print("Estimated Mean Vector:")
print(model$mean)

print("True Covariance Matrix:")
print(trueCovariance)
print("Estimated Covariance Matrix:")
print(model$covarianceMatrix)

groups <- data |>
  group_by(outcomeId) |>
  group_split()
calibratedEstimates <- lapply(groups, function(x) calibrateCiRandomEffects(model, x))
calibratedEstimates <- bind_rows(calibratedEstimates)
calibratedEstimates |>
  summarise(mean(ciLower < 1 & ciUpper > 1))

newData <- data |>
  filter(outcomeId == 1)
# 0.96
calibrateCiRandomEffects(fitResult = model, newData = newData)


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




