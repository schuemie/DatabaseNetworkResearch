# Simulation studies to assess the impact of database differences (population and capture process)
# on the accuracy of database network estimates. We assume there is a target population to which we
# wish to generalize.

library(dplyr)
library(ggplot2)

# Simulation settings ------------------------------------------------------------------------------
#' Title
#'
#' @param nSubgroups              Number of population subgroups.
#' @param trueSubgroupLogRrsMean  Mean of the true log RR per subgroup.
#' @param trueSubgroupLogRrsSd    SD of the true log RR per subgroup.
#' @param nDatabases              Number of databases.
#' @param seLogRrs                Standard error per database (related to sample size).
#' @param nCaptureProcessChars    Number of data capture process characteristics.
#' @param bias0                   Baseline bias.
#' @param biasCpcSd               SD of the bias associated with each data capture characteristic.
#' @param doOvers                 How often the study is repeated until p < 0.05. (Publication bias)
#'
#' @returns
#' A setting object
#'
#' @export
createSimulationSettings <- function(
    nSubgroups = 10,
    trueSubgroupLogRrsMean = 1,
    trueSubgroupLogRrsSd = 1,
    nDatabases = 2,
    seLogRrs = runif(nDatabases, 0.05, 0.25),
    nCaptureProcessChars = 4,
    bias0 = 0.1,
    biasCpcSd = 0.5,
    doOvers = 1
) {
  args <- list()
  for (name in names(formals())) {
    args[[name]] <- get(name)
  }
  return(args)
}

# Simulation ---------------------------------------------------------------------------------------
# settings = createSimulationSettings()
simulateOne <- function(seed, settings) {
  set.seed(seed)
  # Add one database to test reproducability:
  nDatabasesPlusOne <- settings$nDatabases + 1
  seLogRrs <- c(settings$seLogRrs, tail(settings$seLogRrs, 1))

  targetMixture <- runif(settings$nSubgroups, 0, 1)
  targetMixture <- targetMixture / sum(targetMixture) # Subgroup mixture in target population to which we wish to generalize
  trueSubgroupLogRrs <- rnorm(settings$nSubgroups,
                              settings$trueSubgroupLogRrsMean,
                              settings$trueSubgroupLogRrsSd)
  trueTargetLogRr <- sum(targetMixture * trueSubgroupLogRrs)
  biasCpc <- rnorm(settings$nCaptureProcessChars, 0, settings$biasCpcSd) # Bias associated with each data capture process characteristic.

  minP <- 1
  nStudies <- 0
  for (i in seq_len(settings$doOvers)) {
    nStudies <- nStudies + 1
    databaseMixtures <- matrix(runif(nDatabasesPlusOne * settings$nSubgroups),
                               nrow = nDatabasesPlusOne,
                               ncol = settings$nSubgroups)
    databaseMixtures <- databaseMixtures / rowSums(databaseMixtures) # Subgroup mixture per database
    trueDatabaseLogRr <- databaseMixtures %*% trueSubgroupLogRrs
    databaseCpChars = matrix(rbinom(nDatabasesPlusOne * settings$nCaptureProcessChars, 1, 0.1),
                             nrow = nDatabasesPlusOne,
                             ncol = settings$nCaptureProcessChars) # Data capture process characteristics per DB (binary)
    databaseBias <- settings$bias0 + databaseCpChars %*% biasCpc

    # Observed effects:
    databaseLogRrs <- rnorm(nDatabasesPlusOne, trueDatabaseLogRr + databaseBias, seLogRrs)

    data <- data.frame(
      logRr = databaseLogRrs[seq_len(settings$nDatabases)],
      seLogRr = seLogRrs[seq_len(settings$nDatabases)]
    )
    meta <- meta::metagen(data$logRr, data$seLogRr, sm = "RR", control = list(maxiter=1000, stepadj=0.5), prediction = TRUE)
    p <- summary(meta)$random$p
    if (p < minP) {
      minP <- p
      bestStudy <- list(databaseLogRrs = databaseLogRrs, data = data, meta = meta)
      if (p < 0.05) {
        break
      }
    }
  }
  databaseLogRrs <- bestStudy$databaseLogRrs
  data <- bestStudy$data
  meta <- bestStudy$meta

  s <- summary(meta)
  reEstimate <- s$random
  feEstimate <- s$fixed
  breEstimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(data, showProgressBar = FALSE)
  traces <- attr(breEstimate, "traces")
  predictions <- rnorm(nrow(traces), traces[, 1], traces[, 2])
  predictionInterval <- HDInterval::hdi(predictions, credMass = 0.95)
  lbs <- data$logRr + qnorm(0.025) * data$seLogRr
  ubs <- data$logRr + qnorm(0.975) * data$seLogRr
  nSignificant <- sum(lbs > 0 | ubs < 0)
  rows <- tibble(
    method = c("Fixed effects", "Random effects", "Bayesian random effects"),
    seed = seed,
    trueTargetLogRr = trueTargetLogRr,
    logRr = c(feEstimate$TE, reEstimate$TE, breEstimate$mu),
    logRrLb = c(feEstimate$lower, reEstimate$lower, breEstimate$mu95Lb),
    logRrUb = c(feEstimate$upper, reEstimate$upper, breEstimate$mu95Ub),
    logPiLb = c(feEstimate$lower, meta$lower.predict, predictionInterval[1]),
    logPiUb = c(feEstimate$upper, meta$upper.predict, predictionInterval[2]),
    tau = c(0, s$tau, breEstimate$tau),
    tauLb = c(0, s$lower.tau, breEstimate$tau95Lb),
    tauUb = c(0, s$upper.tau, breEstimate$tau95Ub),
    tauSample = c(list(0), list(0), list(sample(traces[, 2], 1000))),
    replicationLogRr = tail(databaseLogRrs, 1),
    replicationLogLb = replicationLogRr + qnorm(0.025) * tail(seLogRrs, 1),
    replicationLogUb = replicationLogRr + qnorm(0.975) * tail(seLogRrs, 1),
    nSignificant = nSignificant,
    nDatabases = settings$nDatabases,
    nStudies = nStudies
  )
  return(rows)
}

computePerformance <- function(results) {
  oddsTest <- function(lb1,hr1,ub1,lb2,hr2,ub2) {
    s1 <- (ub1 - lb1)/(2*1.96)
    s2 <- (ub2 - lb2)/(2*1.96)
    se <- sqrt(s1^2 + s2^2)
    z <- (hr2 - hr1)/se
    dat <- 2*pnorm(-abs(z))
    return(dat)
  }

  metrics <- results |>
    mutate(
      precision = 1/((logRrUb - logRrLb) / qnorm(0.975) * 2)^2,
      coverage = trueTargetLogRr >= logRrLb & trueTargetLogRr <= logRrUb,
      type1 = if_else(trueTargetLogRr == 0, logRrLb > 0 | logRrUb < 0, 0),
      type2 = if_else(trueTargetLogRr != 0, logRrLb <= 0 & logRrUb >= 0, 0),
      reproEstimateInCi = replicationLogRr >= logRrLb & replicationLogRr <= logRrUb,
      reproDisagreeCi = oddsTest(logRrLb, logRr, logRrUb, replicationLogLb, replicationLogRr, replicationLogUb) < 0.05,
      sign = logRrLb > 0 | logRrUb < 0,
      signReproSign = sign & (replicationLogLb > 0 | replicationLogUb < 0),
      precisionPi = 1/((logPiUb - logPiLb) / qnorm(0.975) * 2)^2,
      coveragePi = trueTargetLogRr >= logPiLb & trueTargetLogRr <= logPiUb,
      type1Pi = if_else(trueTargetLogRr == 0, logPiLb > 0 | logPiUb < 0, 0),
      type2Pi = if_else(trueTargetLogRr != 0, logPiLb <= 0 & logPiUb >= 0, 0),
      reproEstimateInPi = replicationLogRr >= logPiLb & replicationLogRr <= logPiUb,
      reproDisagreePi = oddsTest(logPiLb, (logPiLb + logPiUb) / 2, logPiUb, replicationLogLb, replicationLogRr, replicationLogUb) < 0.05,
      signPi = logPiLb > 0 | logPiUb < 0,
      signPiReproSign = signPi & (replicationLogLb > 0 | replicationLogUb < 0),
    ) |>
    group_by(method) |>
    summarise(
      precision = exp(mean(log(precision))),
      coverage = mean(coverage),
      type1 = mean(type1),
      type2 = mean(type2),
      reproEstimateInCi = mean(reproEstimateInCi),
      reproDisagreeCi = mean(reproDisagreeCi, na.rm = TRUE),
      sign = mean(sign),
      signReproSign = mean(signReproSign),
      precisionPi = exp(mean(log(precisionPi))),
      coveragePi = mean(coveragePi),
      type1Pi = mean(type1Pi),
      type2Pi = mean(type2Pi),
      reproEstimateInPi = mean(reproEstimateInPi),
      reproDisagreePi = mean(reproDisagreePi, na.rm = TRUE),
      signPi = mean(signPi, na.rm = TRUE),
      signPiReproSign = mean(signPiReproSign, na.rm = TRUE),
      meanNstudies = mean(nStudies)
    )
  return(metrics)
}

plotTauPosterior <- function(results) {
  x <- seq(from = 0, to = 2, length.out = 100)
  priorData <- tibble(
    tau = c(x, x),
    y = c(dnorm(x, mean = 0, sd = 0.5) * 2, dnorm(x, mean = 0, sd = 0.33) * 2),
    `Half-normal` = rep(c("SD = 0.5", "SD = 0.33"), each = length(x))
  )
  vizData <- tibble(tau = results |>
                      filter(method == "Bayesian random effects") |>
                      pull(tauSample) |>
                      unlist())
  ggplot(vizData, aes(x = tau)) +
    geom_density(fill = "#3f845a", alpha = 0.75) +
    geom_line(aes(y = y, linetype = `Half-normal`), data = priorData) +
    scale_x_continuous("Tau") +
    scale_y_continuous("Density") +
    scale_linetype_manual(values = c("dashed", "dotted")) +
    coord_cartesian(xlim = c(0,2))
}

poolSes <- function(se1, se2) {
  return(sqrt(1/((1/se1^2) + (1/se2^2))))
}

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")

allRows <- list()
for (nDatabases in c(1, 2, 4, 10)) {
  for (trueEffect in c("null", "fixed", "random")) {
    for (publicationBias in c(FALSE, TRUE)) {
      message(sprintf("Simulating %d databases, true effect is %s, with%s publication bias",
                      nDatabases,
                      trueEffect,
                      if (publicationBias) "" else " no"))
      if (nDatabases == 1) {
        seLogRrs <- poolSes(0.1, 0.2)
      } else if (nDatabases == 2) {
        seLogRrs <- c(0.1, 0.2)
      } else {
        seLogRrs <- runif(nDatabases, 0.05, 0.25)
      }
      settings <- createSimulationSettings(
        nDatabases = nDatabases,
        seLogRrs = seLogRrs,
        trueSubgroupLogRrsMean = if (trueEffect == "null") 0 else log(3),
        trueSubgroupLogRrsSd = if (trueEffect == "random") 0.5 else 0,
        doOvers = if (publicationBias) 10 else 1
      )
      results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
      results <- bind_rows(results)
      metrics <- computePerformance(results)
      rows <- tibble(
        nDatabases = !!nDatabases,
        trueEffect = !!trueEffect,
        publicationBias = !!publicationBias
      ) |>
        bind_cols(metrics)
      allRows[[length(allRows) + 1]] <- rows
    }
  }
}
allRows <- bind_rows(allRows)
readr::write_csv(allRows, "SimulationResults.csv")

# Check if tau posterios resembles empirical one
settings <- createSimulationSettings(
  nDatabases = 10,
  trueSubgroupLogRrsMean = 0,
  trueSubgroupLogRrsSd = 0
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
plotTauPosterior(results)

ParallelLogger::stopCluster(cluster)


# Some specific examples ---------------------------------------------------------------------------
library(dplyr)
source("ForestPlot.R")

data <- tibble(
  logRr = c(0.9, 1.1),
  seLogRr = c(0.09, 0.11)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example1.svg", plot, width = 8, height = 3)


data <- tibble(
  logRr = c(0.9, 1.1, 1, 1.05),
  seLogRr = c(0.09, 0.11, 0.3, 0.5)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example1b.svg", plot, width = 8, height = 3)

data <- tibble(
  logRr = c(0.9),
  seLogRr = c(0.09)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example1c.svg", plot, width = 8, height = 3)

data <- tibble(
  logRr = c(.1, 1.5),
  seLogRr = c(0.75, 0.08)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example2.svg", plot, width = 8, height = 3)

data <- tibble(
  logRr = c(0.4, 2),
  seLogRr = c(0.03, 0.08)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example3.svg", plot, width = 8, height = 3)

data <- tibble(
  logRr = c(1, 1, 1, 1),
  seLogRr = c(0.08, 0.6, 0.61, 0.62)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example4.svg", plot, width = 8, height = 3)

data <- tibble(
  logRr = c(2, 1, 0.15, -0.15),
  seLogRr = c(0.04, 0.06, 0.81, 0.8)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example5.svg", plot, width = 8, height = 3)

data <- tibble(
  logRr = c(2, 2.1, 0.15, -0.15),
  seLogRr = c(0.04, 0.06, 0.08, 0.06)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example6.svg", plot, width = 8, height = 3)

