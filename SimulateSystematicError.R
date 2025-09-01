# Simulate systematic error in a network of databases. We assume a generative model where each
# outcome is susceptible to a random set of biases, and in each database a random set of biases
# play out.

library(dplyr)
library(ggplot2)

# Simulation settings ------------------------------------------------------------------------------
createSimulationSettings <- function(
    nDatabases = 5,
    nNegativeControls = 50,
    trueLogRr = log(2),
    trueTau = 0.25,
    nBiasSources = 10,
    biasSourcePrevalences = runif(nBiasSources, 0, 1),
    biasSourceSd = 0.1,
    biasOutcomeSd = 0.1,
    minDatabaseSizeMultiplier = 0.5,
    maxDatabaseSizeMultiplier = 2,
    minSe = 0.05,
    maxSe = 0.5
) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(args)
}

# Various functions --------------------------------------------------------------------------------
plotSystematicErrorDistributions <- function(logRrs, seLogRrs) {
  x <- seq(log(0.1), log(10), length.out = 100)
  compute <- function(x, mcmc) {
    yMcmc <- dnorm(rep(x, nrow(mcmc$chain)), mean = mcmc$chain[, 1], sd = 1/sqrt(mcmc$chain[, 2]))
    return(quantile(yMcmc, c(0.025, 0.5, 0.975)))
  }
  plotData <- list()
  for (i in seq_len(settings$nDatabases)) {
    null <- EmpiricalCalibration::fitMcmcNull(
      logRr = logRrs[seq_len(settings$nNegativeControls), i],
      seLogRr = seLogRrs[seq_len(settings$nNegativeControls), i]
    )
    ys <- sapply(x, compute, mcmc = attr(null, "mcmc"))
    y <- ys[2, ]
    yMaxLb <- ys[1, ]
    yMaxUb <- ys[3, ]
    normFactor <- max(ys[2, ])
    y <- y / normFactor
    yMaxLb <- yMaxLb / normFactor
    yMaxUb <- yMaxUb / normFactor
    plotData[[i]] <- tibble(
      databaseId = sprintf("Database %s", i),
      x = x,
      yMax = y,
      yMaxLb = yMaxLb,
      yMaxUb =  yMaxUb,
      yMin = 0
    )
    # Are we calibrated (using leave-one-out)?
    # EmpiricalCalibration::plotCalibration(
    #   logRr = logRrs[seq_len(settings$nNegativeControls), i],
    #   seLogRr = seLogRrs[seq_len(settings$nNegativeControls), i],
    #   useMcmc = F
    # )
  }
  plotData <- bind_rows(plotData)
  breaks <- c(0.25, 1, 4, 8)
  ggplot(plotData) +
    geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.4) +
    geom_vline(xintercept = 0, size = 0.8) +
    geom_ribbon(aes(x = .data$x, ymax = .data$yMax, ymin = .data$yMin), fill = "#FF2700", alpha = 0.6, data = plotData) +
    geom_ribbon(aes(x = .data$x, ymax = .data$yMaxUb, ymin = .data$yMax), fill = "#FF2700", alpha = 0.3, data = plotData) +
    coord_cartesian(xlim = log(c(0.1, 10)), ylim = c(0, 2)) +
    scale_x_continuous("Systematic Error", breaks = log(breaks), labels = breaks) +
    facet_grid(databaseId ~ .) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 14),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          strip.text.x = element_text(size = 14),
          strip.text.y.left = element_text(size = 14, angle = 0, hjust = 0),
          strip.background = element_blank(),
          panel.spacing.y = unit(0, "lines"))
}

applyCurrentApproach <- function(logRrs, seLogRrs, settings) {
  estimates <- list()
  taus <- c()
  for (i in seq_len(nrow(logRrs))) {
    data <- tibble(logRr = logRrs[i, ], seLogRr = seLogRrs[i, ])
    estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(data, showProgressBar = FALSE)
    estimates[[i]] <- tibble(logRr = estimate$mu,
                             seLogRr = estimate$muSe,
                             logLb = estimate$mu95Lb,
                             logUb = estimate$mu95Ub)
    taus[[i]] <- estimate$tau
  }
  estimates <- bind_rows(estimates)
  null <- EmpiricalCalibration::fitMcmcNull(
    logRr = estimates$logRr[seq_len(settings$nNegativeControls)],
    seLogRr = estimates$seLogRr[seq_len(settings$nNegativeControls)]
  )
  estimateHoi <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr = tail(estimates$logRr, 1),
    seLogRr = tail(estimates$seLogRr, 1),
    model = EmpiricalCalibration::convertNullToErrorModel(null)
  )
  return(estimateHoi)
}

applyGeneralizedModel <- function(logRrs, seLogRrs, settings) {
  fitSystematicErrorModel <- function(data) {
    # Ensure database IDs are factors to maintain a consistent order
    databaseIds <- unique(data$databaseId)
    data$databaseId <- factor(data$databaseId, levels = databaseIds)
    nDatabases <- length(databaseIds)
    controlDataList <- split(data, data$outcomeId)
    controlDataList <- controlDataList[sapply(controlDataList, nrow) == nDatabases]
    if (length(controlDataList) == 0) {
      stop("No controls found with data in all specified databases.")
    }
    calculateNegativeLogLikelihood <- function(params, dataList) {
      meanVector <- params[1:nDatabases]
      choleskyParams <- params[(nDatabases + 1):length(params)]
      choleskyMatrix <- matrix(0, nrow = nDatabases, ncol = nDatabases)
      choleskyMatrix[lower.tri(choleskyMatrix, diag = TRUE)] <- choleskyParams
      covarianceMatrix <- choleskyMatrix %*% t(choleskyMatrix)
      logLikelihoods <- sapply(dataList, function(controlData) {
        controlData <- controlData[order(controlData$databaseId), ]
        y <- controlData$logRr
        d <- diag(controlData$seLogRr^2)
        observedCovariance <- covarianceMatrix + d
        ll <- tryCatch({
          mvtnorm::dmvnorm(y, mean = meanVector, sigma = observedCovariance, log = TRUE)
        }, error = function(e) {
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
      covarianceMatrix = finalCovarianceMatrix
    ))
  }
  data <- tibble(
    logRr = as.vector(logRrs[seq_len(settings$nNegativeControls), ]),
    seLogRr = as.vector(seLogRrs[seq_len(settings$nNegativeControls), ]),
    databaseId = rep(seq_len(settings$nDatabases), each = settings$nNegativeControls),
    outcomeId =rep(seq_len(settings$nNegativeControls), settings$nDatabases)
  )
  model <- fitSystematicErrorModel(data)

  calibrateCiRandomEffects <- function(fitResult, newData) {
    modelDbIds <- names(fitResult$mean)
    availableDbIds <- intersect(modelDbIds, newData$databaseId)
    if (length(availableDbIds) == 0) {
      stop("None of the databases in newData are present in the systematic error model.")
    }
    meanSys <- fitResult$mean[availableDbIds]
    covSys <- fitResult$covarianceMatrix[availableDbIds, availableDbIds, drop = FALSE]
    newData <- newData[newData$databaseId %in% availableDbIds, ]
    newData <- newData[match(availableDbIds, newData$databaseId), ] # Ensure order
    nDatabases <- length(availableDbIds)

    if (nDatabases == 1) {
      muHatFinal <- newData$logRr - meanSys
      muVarFinal <- covSys[1, 1] + newData$seLogRr^2
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
      yPrime <- newData$logRr - meanSys
      sigmaTotal <- covSys + diag(newData$seLogRr^2)
      ones <- rep(1, nDatabases)
      negLogLike <- function(logTau2) {
        tau2 <- exp(logTau2)
        if (is.infinite(tau2)) return(1e10)
        V <- sigmaTotal + diag(tau2, nDatabases)
        VInv <- try(solve(V), silent = TRUE)
        if (inherits(VInv, "try-error")) return(1e10)
        muHat <- (t(ones) %*% VInv %*% yPrime) / (t(ones) %*% VInv %*% ones)
        ll <- mvtnorm::dmvnorm(yPrime, mean = rep(as.numeric(muHat), nDatabases), sigma = V, log = TRUE)
        return(-ll)
      }
      opt <- optim(par = log(0.01), fn = negLogLike, method = "BFGS", hessian = TRUE)
      estimatedTau2 <- exp(opt$par)
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
      VFinal <- sigmaTotal + diag(estimatedTau2, nDatabases)
      VFinalInv <- solve(VFinal)
      muHatFinal <- (t(ones) %*% VFinalInv %*% yPrime) / (t(ones) %*% VFinalInv %*% ones)
      muVarFinal <- 1 / (t(ones) %*% VFinalInv %*% ones)
      ciLowerLog <- muHatFinal - 1.96 * sqrt(muVarFinal)
      ciUpperLog <- muHatFinal + 1.96 * sqrt(muVarFinal)
      return(data.frame(
        logRr = as.numeric(muHatFinal),
        logLb95Rr =  as.numeric(ciLowerLog),
        logUb95Rr =  as.numeric(ciUpperLog),
        seLogRr = sqrt(muVarFinal),
        tau2 = estimatedTau2,
        seLogTau2 = seLogTau2,
        tau2CiLower = tau2CiLower,
        tau2CiUpper = tau2CiUpper
      ))
    }
  }
  newData <- tibble(
    logRr = logRrs[settings$nNegativeControls + 1, ],
    seLogRr = seLogRrs[settings$nNegativeControls + 1, ],
    databaseId = seq_len(settings$nDatabases)
  )
  estimate <- calibrateCiRandomEffects(model, newData)
  return(estimate)
}

# Simulation function ------------------------------------------------------------------------------
# settings = createSimulationSettings()
simulateOne <- function(seed, settings) {
  set.seed(seed)

  # Compute standard error for each outcome
  # Draw the database size multiplier, reflecting some databases are bigger than others:
  dbSizeMultipliers <- runif(settings$nDatabases, settings$minDatabaseSizeMultiplier, settings$maxDatabaseSizeMultiplier)
  # Draw the base standard error for each outcome, reflecting some outcomes are more prevalent than others:
  outcomeSes <- runif(settings$nNegativeControls + 1, settings$minSe, settings$maxSe)
  # Multiply the two to get the SE per outcome in each database:
  seLogRrs <- outer(outcomeSes, dbSizeMultipliers)

  # Compute bias for each outcome
  # First, compute which database is vulnerable to which source of bias:
  biasSourcesPerDb <- matrix(rbinom(settings$nDatabases * settings$nBiasSources, 1, settings$biasSourcePrevalences),
                             nrow = settings$nBiasSources,
                             ncol = settings$nDatabases)
  # Compute the mean bias caused by each source:
  biasSourceMean <- rnorm(settings$nBiasSources, 0, settings$biasSourceSd)
  # Compute the bias caused by each source for each outcome:
  biasSourceOutcome <- matrix(rnorm(settings$nBiasSources * (settings$nNegativeControls + 1), biasSourceMean, settings$biasSourceSd),
                              nrow = settings$nBiasSources,
                              ncol = settings$nNegativeControls + 1)
  # Multiply the two to see how much bias we have for each outcome in each database:
  biasOutcomeDb <- t(biasSourceOutcome) %*% biasSourcesPerDb

  # Compute observed effect sizes
  # Draw the true effect for the outcome of interest in each database:
  trueLogRrPerDb <- rnorm(settings$nDatabases, settings$trueLogRr, settings$trueTau)
  # Draw the observed effect for each outcome. Only for the outcome of interest do we add the true effect
  # to the bias:
  logRrs <- matrix(rnorm((settings$nNegativeControls + 1) * settings$nDatabases, biasOutcomeDb + c(rep(0, settings$nNegativeControls), 1) %*% t(trueLogRrPerDb), seLogRrs),
                   nrow = settings$nNegativeControls + 1,
                   ncol = settings$nDatabases)

  # Confirmation: Fit systematic error models per database and compare to observed to see if our
  # simulation looks like the real thing:
  # plotSystematicErrorDistributions(logRrs, seLogRrs)

  # Use current approach: meta-analysis per outcome:
  # estimateCurrent <- applyCurrentApproach(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings)


  estimateGenmodel <- applyGeneralizedModel(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings)

  results <- bind_rows(
    estimateCurrent |>
      mutate(tau = NA) |>
      select(logRr, logLb95Rr, logUb95Rr, seLogRr, tau) |>
      mutate(method = "Current",
             coverage = logLb95Rr < settings$trueLogRr & logUb95Rr > settings$trueLogRr),
    estimateGenmodel |>
      mutate(tau = sqrt(tau2)) |>
      select(logRr, logLb95Rr, logUb95Rr, seLogRr, tau) |>
      mutate(method = "Generalized model",
             coverage = logLb95Rr < settings$trueLogRr & logUb95Rr > settings$trueLogRr)
  )
  return(results)
}

# Run simulation -----------------------------------------------------------------------------------
cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
snow::clusterExport(cluster, c("applyCurrentApproach", "applyGeneralizedModel"))

settings <- createSimulationSettings()
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
results |>
  group_by(method) |>
  summarise(coverage = mean(coverage),
            precision = exp(mean(log(1/seLogRr ^ 2))))
