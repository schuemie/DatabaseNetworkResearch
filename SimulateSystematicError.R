# Simulate systematic error in a network of databases. We assume a generative model where each
# outcome is susceptible to a random set of biases, and in each database a random set of biases
# play out.

source("SystematicErrorModel.R")
library(dplyr)
library(ggplot2)

# Need to have these packages installed:
# install.packages("mvtnorm")
# install.packages("EmpiricalCalibration")
# install.packages("EvidenceSynthesis")
# install.packages("ParallelLogger")

# Simulation settings ------------------------------------------------------------------------------
#' Create simulation settings
#'
#' @param nDatabases                Number of databases.
#' @param nNegativeControls         Number of negative controls.
#' @param nOutcomesOfInterest       Number of outcomes of interest.
#' @param trueLogRr                 True log effect size for the outcomes of interest.
#' @param trueTau                   True tau (heterogeneity) for the effect size of the outcomes of interest.
#' @param nBiasSources              Number of distinct bias sources.
#' @param minBiasSourcePrevalence   Minimum prevalences of sources of bias across databases.
#' @param maxBiasSourcePrevalence   Maximum prevalences of sources of bias across databases.
#' @param biasSourceSd              Standard deviation of mean bias caused by each source.
#' @param biasOutcomeSd             Standard deviation of the bias from a source around its mean (across outcomes).
#' @param minDatabaseSizeMultiplier Minimum multiplier for the standard error in a database.
#' @param maxDatabaseSizeMultiplier Maximum multiplier for the standard error in a database.
#' @param minSe                     Minimum standard error for an outcome, before applying the database multiplier.
#' @param maxSe                     Maximum standard error for an outcome, before applying the database multiplier.
#'
#' @returns
#' A settings object to be used in `simulateOne()`.
#'
createSimulationSettings <- function(
    nDatabases = 5,
    nNegativeControls = 50,
    nOutcomesOfInterest = 10,
    trueLogRr = log(2),
    trueTau = 0.20,
    nBiasSources = 10,
    minBiasSourcePrevalence = 0,
    maxBiasSourcePrevalence = 1,
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
plotSystematicErrorDistributions <- function(logRrs, seLogRrs, settings) {
  # Compute the systematic error distribution within each database and plot it. Used to check if
  # distributions look similar to what is observed in real-world

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


computeWithinDatabaseCoverage <- function(logRrs, seLogRrs, settings, trueLogRrsPerDb) {
  # Perform calibration within each database, and compute coverage of within-database estimates.
  # Used to verify calibration still works in this specific simulation scenario.

  coverage <- c()
  for (i in seq_len(settings$nDatabases)) {
    null <- EmpiricalCalibration::fitMcmcNull(
      logRr = logRrs[seq_len(settings$nNegativeControls), i],
      seLogRr = seLogRrs[seq_len(settings$nNegativeControls), i]
    )
    estimatesHois <- EmpiricalCalibration::calibrateConfidenceInterval(
      logRr = tail(logRrs[, i], settings$nOutcomesOfInterest),
      seLogRr = tail(seLogRrs[, i], settings$nOutcomesOfInterest),
      model = EmpiricalCalibration::convertNullToErrorModel(null)
    )
    trueLogRrs <- tail(trueLogRrsPerDb[, i], settings$nOutcomesOfInterest)
    coverage <- c(coverage, estimatesHois$logLb95Rr < trueLogRrs & estimatesHois$logUb95Rr > trueLogRrs)
  }
  return(mean(coverage))
}

applyCurrentApproach <- function(logRrs, seLogRrs, settings, bayesian = TRUE) {
  # Apply the current approach, where we first meta-analyse all outcomes, then use the meta-analytic
  # estimates to fit the empirical null and calibrat.

  estimates <- list()
  for (i in seq_len(nrow(logRrs))) {
    data <- tibble(logRr = logRrs[i, ], seLogRr = seLogRrs[i, ])
    if (bayesian) {
      estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(data, showProgressBar = FALSE)
      estimates[[i]] <- tibble(logRr = estimate$mu,
                               seLogRr = estimate$muSe,
                               logLb = estimate$mu95Lb,
                               logUb = estimate$mu95Ub)
    } else {
      meta <- meta::metagen(data$logRr, data$seLogRr, sm = "RR", control = list(maxiter=1000, stepadj=0.5))
      s <- summary(meta)
      rnd <- s$random
      estimates[[i]] <- tibble(logRr = rnd$TE,
                               seLogRr = rnd$seTE,
                               logLb = rnd$lower,
                               logUb = rnd$upper)
    }
  }
  estimates <- bind_rows(estimates)
  null <- EmpiricalCalibration::fitMcmcNull(
    logRr = estimates$logRr[seq_len(settings$nNegativeControls)],
    seLogRr = estimates$seLogRr[seq_len(settings$nNegativeControls)]
  )
  estimatesHois <- EmpiricalCalibration::calibrateConfidenceInterval(
    logRr = tail(estimates$logRr, settings$nOutcomesOfInterest),
    seLogRr = tail(estimates$seLogRr, settings$nOutcomesOfInterest),
    model = EmpiricalCalibration::convertNullToErrorModel(null)
  )
  return(estimatesHois)
}

applyNaiveApproach <- function(logRrs, seLogRrs, settings, bayesian = TRUE) {
  # Apply the naive approach: calibrate estimates within each database, then meta-analyse calibrated
  # estimates. This should in theory be bad because it does not take into account that systematic
  # error is calibrated between databases.

  calibratedEstimates <- list()
  for (i in seq_len(settings$nDatabases)) {
    null <- EmpiricalCalibration::fitMcmcNull(
      logRr = logRrs[seq_len(settings$nNegativeControls), i],
      seLogRr = seLogRrs[seq_len(settings$nNegativeControls), i]
    )
    calibratedEstimates[[i]] <- EmpiricalCalibration::calibrateConfidenceInterval(
      logRr = tail(logRrs[, i], settings$nOutcomesOfInterest),
      seLogRr = tail(seLogRrs[, i], settings$nOutcomesOfInterest),
      model = EmpiricalCalibration::convertNullToErrorModel(null)
    ) |>
      mutate(outcomeId = seq_len(settings$nOutcomesOfInterest))
  }
  calibratedEstimates <- bind_rows(calibratedEstimates)
  groups <- calibratedEstimates |>
    group_by(outcomeId) |>
    group_split()
  estimatesHois <- list()
  for (group in groups) {
    if (bayesian) {
      estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(group, showProgressBar = FALSE)
      estimatesHois[[group$outcomeId[1]]] <- tibble(logRr = estimate$mu,
                                                    seLogRr = estimate$muSe,
                                                    logLb95Rr = estimate$mu95Lb,
                                                    logUb95Rr = estimate$mu95Ub,
                                                    tau = estimate$tau)
    } else {
      meta <- meta::metagen(group$logRr, group$seLogRr, sm = "RR", control = list(maxiter=1000, stepadj=0.5))
      s <- summary(meta)
      rnd <- s$random
      estimatesHois[[group$outcomeId[1]]] <- tibble(logRr = rnd$TE,
                                                    seLogRr = rnd$seTE,
                                                    logLb95Rr = rnd$lower,
                                                    logUb95Rr = rnd$upper,
                                                    tau = s$tau)
    }
  }
  estimatesHois <- bind_rows(estimatesHois)
  return(estimatesHois)
}

applyGeneralizedModel <- function(logRrs, seLogRrs, settings) {
  # Apply Martijn's generalized calibration model. Currently only supports non-Bayesian approach.

  data <- tibble(
    logRr = as.vector(logRrs[seq_len(settings$nNegativeControls), ]),
    seLogRr = as.vector(seLogRrs[seq_len(settings$nNegativeControls), ]),
    databaseId = rep(seq_len(settings$nDatabases), each = settings$nNegativeControls),
    outcomeId = rep(seq_len(settings$nNegativeControls), settings$nDatabases)
  )
  model <- fitSystematicErrorModel(data)
  estimates <- list()
  for (i in seq_len(settings$nOutcomesOfInterest)) {
    newData <- tibble(
      logRr = logRrs[settings$nNegativeControls + i, ],
      seLogRr = seLogRrs[settings$nNegativeControls + i, ],
      databaseId = seq_len(settings$nDatabases)
    )
    estimates[[i]] <- calibrateCiRandomEffects(model, newData)

  }
  estimates <- bind_rows(estimates)
  return(estimates)
}

# Simulation function ------------------------------------------------------------------------------
# settings = createSimulationSettings()
simulateOne <- function(seed, settings) {
  set.seed(seed)

  # Compute standard error for each outcome
  # Draw the database size multiplier, reflecting some databases are bigger than others:
  dbSizeMultipliers <- runif(settings$nDatabases, settings$minDatabaseSizeMultiplier, settings$maxDatabaseSizeMultiplier)
  # Draw the base standard error for each outcome, reflecting some outcomes are more prevalent than others:
  outcomeSes <- runif(settings$nNegativeControls + settings$nOutcomesOfInterest, settings$minSe, settings$maxSe)
  # Multiply the two to get the SE per outcome in each database:
  seLogRrs <- outer(outcomeSes, dbSizeMultipliers)

  # Compute bias for each outcome
  # First, sample the prevalence of each bias source across databases:
  biasSourcePrevalences <- runif(settings$nBiasSources, settings$minBiasSourcePrevalence, settings$maxBiasSourcePrevalence)
  # Compute which database is vulnerable to which source of bias:
  biasSourcesPerDb <- matrix(rbinom(settings$nDatabases * settings$nBiasSources, 1, biasSourcePrevalences),
                             nrow = settings$nBiasSources,
                             ncol = settings$nDatabases)
  # Compute the mean bias caused by each source:
  biasSourceMean <- rnorm(settings$nBiasSources, 0, settings$biasSourceSd)
  # Compute the bias caused by each source for each outcome:
  biasSourceOutcome <- matrix(rnorm(settings$nBiasSources * (settings$nNegativeControls + settings$nOutcomesOfInterest), biasSourceMean, settings$biasSourceSd),
                              nrow = settings$nBiasSources,
                              ncol = settings$nNegativeControls + settings$nOutcomesOfInterest)
  # Multiply the two to see how much bias we have for each outcome in each database:
  biasOutcomeDb <- t(biasSourceOutcome) %*% biasSourcesPerDb

  # Is bias correlated?
  # cor(biasOutcomeDb)

  # Compute observed effect sizes
  # Set true effect to 0 for negative controls and draw the true effect for the outcome of interest in each database:
  trueLogRrsPerDb <- rbind(
    matrix(rep(0, settings$nNegativeControls * settings$nDatabases),
           nrow = settings$nNegativeControls,
           ncol = settings$nDatabases),
    matrix(rnorm(settings$nDatabases * settings$nOutcomesOfInterest, settings$trueLogRr, settings$trueTau),
           nrow = settings$nOutcomesOfInterest,
           ncol = settings$nDatabases)
  )
  # Draw the observed effect for each outcome. Only for the outcome of interest do we add the true effect
  # to the bias:
  logRrs <- matrix(rnorm((settings$nNegativeControls + settings$nOutcomesOfInterest) * settings$nDatabases,
                         biasOutcomeDb + trueLogRrsPerDb,
                         seLogRrs),
                   nrow = settings$nNegativeControls + settings$nOutcomesOfInterest,
                   ncol = settings$nDatabases)


  results <- tibble()
  na <- as.numeric(NA)

  # Confirmation: Fit systematic error models per database and compare to observed to see if our
  # simulation looks like the real thing. (Doesn't make sense when multi-threading)
  # plotSystematicErrorDistributions(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings)

  # Confirmation: Compute within database coverage to confirm our calibration procedure holds under these conditions:
  # coverageWithinDbs <- computeWithinDatabaseCoverage(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings, trueLogRrsPerDb = trueLogRrsPerDb)
  # results <- results |>
  #   bind_rows(
  #     tibble(logRr = na, logLb95Rr = na, logUb95Rr = na, seLogRr = na, tau = na) |>
  #       mutate(method = "Within database",
  #              seed = !!seed,
  #              outcomeId = na,
  #              coverage = coverageWithinDbs)
  #   )

  # Use current approach: Bayesian meta-analysis per outcome, then calibrate:
  # estimatesCurrent <- applyCurrentApproach(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings)
  # results <- results |>
  #   bind_rows(
  #     estimatesCurrent |>
  #       mutate(tau = na) |>
  #       select(logRr, logLb95Rr, logUb95Rr, seLogRr, tau) |>
  #       mutate(method = "Current",
  #              seed = !!seed,
  #              outcomeId = seq_len(settings$nOutcomesOfInterest),
  #              coverage = logLb95Rr < settings$trueLogRr & logUb95Rr > settings$trueLogRr)
  #   )

  # Use current approach but non-Bayesian:
  # estimatesNonBayesian <- applyCurrentApproach(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings, bayesian = FALSE)
  # results <- results |>
  #   bind_rows(
  #     estimatesNonBayesian |>
  #       mutate(tau = na) |>
  #       select(logRr, logLb95Rr, logUb95Rr, seLogRr, tau) |>
  #       mutate(method = "Current (non-Bayesian)",
  #              seed = !!seed,
  #              outcomeId = seq_len(settings$nOutcomesOfInterest),
  #              coverage = logLb95Rr < settings$trueLogRr & logUb95Rr > settings$trueLogRr)
  #   )

  # Use naive approach: first calibrate per-database estimates, then meta-analyse calibrated estimates
  # estimatesNaive <- applyNaiveApproach(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings)
  # results <- results |>
  #   bind_rows(
  #     estimatesNaive |>
  #       select(logRr, logLb95Rr, logUb95Rr, seLogRr, tau) |>
  #       mutate(method = "Naive",
  #              seed = !!seed,
  #              outcomeId = seq_len(settings$nOutcomesOfInterest),
  #              coverage = logLb95Rr < settings$trueLogRr & logUb95Rr > settings$trueLogRr)
  #   )

  # Use Martijn's frequentist generalized model:
  estimatesGenModel <- applyGeneralizedModel(logRrs = logRrs, seLogRrs = seLogRrs, settings = settings)
  results <- results |>
    bind_rows(
      estimatesGenModel |>
        mutate(tau = sqrt(tau2), logLb95Rr = log(ciLower), logUb95Rr = log(ciUpper)) |>
        select(logRr, logLb95Rr, logUb95Rr, seLogRr, tau) |>
        mutate(method = "Generalized model",
               seed = !!seed,
               outcomeId = seq_len(settings$nOutcomesOfInterest),
               coverage = logLb95Rr < settings$trueLogRr & logUb95Rr > settings$trueLogRr)
    )

  return(results)
}

# Run simulation -----------------------------------------------------------------------------------
cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "tidyr")
snow::clusterExport(cluster, c("fitSystematicErrorModel", "calibrateCiRandomEffects"))
snow::clusterExport(cluster, c("computeWithinDatabaseCoverage", "applyCurrentApproach", "applyNaiveApproach", "applyGeneralizedModel"))

settings <- createSimulationSettings(nDatabases = 5)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
results |>
  group_by(method) |>
  summarise(coverage = mean(coverage),
            precision = exp(mean(log(1/seLogRr ^ 2))))
# method                 coverage precision
# Current                   0.942      15.5
# Current (non-Bayesian)    0.939      15.7
# Generalized model         0.921      17.8
# Naive                     0.915      19.4
# Within database           0.939      NA


settings <- createSimulationSettings(nDatabases = 10)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
results |>
  group_by(method) |>
  summarise(coverage = mean(coverage),
            precision = exp(mean(log(1/seLogRr ^ 2))))
# method                 coverage precision
# Current                   0.955      21.3
# Current (non-Bayesian)    0.955      21.4
# Generalized model         0.922      25.9
# Naive                     0.82       45.8
# Within database           0.944      NA


settings <- createSimulationSettings(nDatabases = 5,
                                     minDatabaseSizeMultiplier = 1,
                                     maxDatabaseSizeMultiplier = 1,
                                     minSe = 0.2,
                                     maxSe = 0.4)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
results |>
  group_by(method) |>
  summarise(coverage = mean(coverage),
            precision = exp(mean(log(1/seLogRr ^ 2))))
# method                 coverage precision
# Current                   0.934      17.8
# Current (non-Bayesian)    0.937      17.0
# Generalized model         0.926      18.5
# Naive                     0.91       21.1
# Within database           0.938      NA

ParallelLogger::stopCluster(cluster)



# Explore single example -----------------------------------------------------------------
settings <- createSimulationSettings(nDatabases = 5,
                                     nNegativeControls = 100,
                                     nOutcomesOfInterest = 1000,
                                     minSe = 0.01,
                                     maxSe = 0.05)
set.seed(1)

dbSizeMultipliers <- runif(settings$nDatabases, settings$minDatabaseSizeMultiplier, settings$maxDatabaseSizeMultiplier)
# Draw the base standard error for each outcome, reflecting some outcomes are more prevalent than others:
outcomeSes <- runif(settings$nNegativeControls + settings$nOutcomesOfInterest, settings$minSe, settings$maxSe)
# Multiply the two to get the SE per outcome in each database:
seLogRrs <- outer(outcomeSes, dbSizeMultipliers)

# Compute bias for each outcome
# First, sample the prevalence of each bias source across databases:
biasSourcePrevalences <- runif(settings$nBiasSources, settings$minBiasSourcePrevalence, settings$maxBiasSourcePrevalence)
# Compute which database is vulnerable to which source of bias:
biasSourcesPerDb <- matrix(rbinom(settings$nDatabases * settings$nBiasSources, 1, biasSourcePrevalences),
                           nrow = settings$nBiasSources,
                           ncol = settings$nDatabases)
# Compute the mean bias caused by each source:
biasSourceMean <- rnorm(settings$nBiasSources, 0, settings$biasSourceSd)
# Compute the bias caused by each source for each outcome:
biasSourceOutcome <- matrix(rnorm(settings$nBiasSources * (settings$nNegativeControls + settings$nOutcomesOfInterest), biasSourceMean, settings$biasSourceSd),
                            nrow = settings$nBiasSources,
                            ncol = settings$nNegativeControls + settings$nOutcomesOfInterest)
# Multiply the two to see how much bias we have for each outcome in each database:
biasOutcomeDb <- t(biasSourceOutcome) %*% biasSourcesPerDb

# Compute observed effect sizes
# Set true effect to 0 for negative controls and draw the true effect for the outcome of interest in each database:
trueLogRrsPerDb <- rbind(
  matrix(rep(0, settings$nNegativeControls * settings$nDatabases),
         nrow = settings$nNegativeControls,
         ncol = settings$nDatabases),
  matrix(rnorm(settings$nDatabases * settings$nOutcomesOfInterest, settings$trueLogRr, settings$trueTau),
         nrow = settings$nOutcomesOfInterest,
         ncol = settings$nDatabases)
)
# Draw the observed effect for each outcome. Only for the outcome of interest do we add the true effect
# to the bias:
logRrs <- matrix(rnorm((settings$nNegativeControls + settings$nOutcomesOfInterest) * settings$nDatabases,
                       biasOutcomeDb + trueLogRrsPerDb,
                       seLogRrs),
                 nrow = settings$nNegativeControls + settings$nOutcomesOfInterest,
                 ncol = settings$nDatabases)

data <- tibble(
  logRr = as.vector(logRrs[seq_len(settings$nNegativeControls), ]),
  seLogRr = as.vector(seLogRrs[seq_len(settings$nNegativeControls), ]),
  databaseId = rep(seq_len(settings$nDatabases), each = settings$nNegativeControls),
  outcomeId = rep(seq_len(settings$nNegativeControls), settings$nDatabases)
)
model <- fitSystematicErrorModel(data)
model$mean
colMeans(biasOutcomeDb)

model$covarianceMatrix
cov(biasOutcomeDb)

estimates <- list()
for (i in seq_len(settings$nOutcomesOfInterest)) {
  newData <- tibble(
    logRr = logRrs[settings$nNegativeControls + i, ],
    seLogRr = seLogRrs[settings$nNegativeControls + i, ],
    databaseId = seq_len(settings$nDatabases)
  )
  estimates[[i]] <- calibrateCiRandomEffects(model, newData)

}
estimates <- bind_rows(estimates)
mean(sqrt(estimates$tau2CiLower) < settings$trueTau & sqrt(estimates$tau2CiUpper) > settings$trueTau)
mean(log(estimates$ciLower) < settings$trueLogRr & log(estimates$ciUpper) > settings$trueLogRr)
