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
    biasCpcSd = 0.5
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

  databaseMixtures <- matrix(runif(nDatabasesPlusOne * settings$nSubgroups),
                             nrow = nDatabasesPlusOne,
                             ncol = settings$nSubgroups)
  databaseMixtures <- databaseMixtures / rowSums(databaseMixtures) # Subgroup mixture per database
  trueDatabaseLogRr <- databaseMixtures %*% trueSubgroupLogRrs

  databaseCpChars = matrix(rbinom(nDatabasesPlusOne * settings$nCaptureProcessChars, 1, 0.1),
                           nrow = nDatabasesPlusOne,
                           ncol = settings$nCaptureProcessChars) # Data capture process characteristics per DB (binary)
  biasCpc <- rnorm(settings$nCaptureProcessChars, 0, settings$biasCpcSd) # Bias associated with each data capture process characteristic.
  databaseBias <- settings$bias0 + databaseCpChars %*% biasCpc

  # Observed effects:
  databaseLogRrs <- rnorm(nDatabasesPlusOne, trueDatabaseLogRr + databaseBias, seLogRrs)

  data <- data.frame(
    logRr = databaseLogRrs[seq_len(settings$nDatabases)],
    seLogRr = seLogRrs[seq_len(settings$nDatabases)]
  )
  lbs <- data$logRr + qnorm(0.025) * data$seLogRr
  ubs <- data$logRr + qnorm(0.975) * data$seLogRr
  nSignificant <- sum(lbs > 0 | ubs < 0)

  meta <- meta::metagen(data$logRr, data$seLogRr, sm = "RR", control = list(maxiter=1000, stepadj=0.5), prediction = TRUE)
  s <- summary(meta)
  feEstimate <- s$fixed
  reEstimate <- s$random
  breEstimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(data, showProgressBar = FALSE)
  traces <- attr(breEstimate, "traces")
  predictions <- rnorm(nrow(traces), traces[, 1], traces[, 2])
  predictionInterval <- HDInterval::hdi(predictions, credMass = 0.95)
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
    nDatabases = settings$nDatabases
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
      signPiReproSign = signPi & (replicationLogLb > 0 | replicationLogUb < 0)
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
      signPiReproSign = mean(signPiReproSign, na.rm = TRUE)
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

# 2 DBs, non-null effect
settings <- createSimulationSettings(
  nDatabases = 2,
  seLogRrs = c(0.2, 0.4),
  trueSubgroupLogRrsMean = 1,
  trueSubgroupLogRrsSd = 0.5
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
computePerformance(results)
# method                  precision coverage type1   type2 reproEstimateInCi reproDisagreeCi  sign signReproSign precisionPi coveragePi type1Pi type2Pi reproEstimateInPi reproDisagreePi signPi signPiReproSign
# 1 Bayesian random effects     0.430     0.98     0 0.21                 0.83            0.07  0.79          0.56      0.169        1          0 0.69                 0.97            0      0.31            0.23
# 2 Fixed effects               1.95      0.85     0 0.01000              0.58            0.11  0.99          0.71      1.95         0.85       0 0.01000              0.58            0.11   0.99            0.71
# 3 Random effects              1.12      0.89     0 0.1                  0.66            0.07  0.9           0.65      0.0194       1          0 1                    1               0      0               0

# 1 DB, non-null effect
settings <- createSimulationSettings(
  nDatabases = 1,
  seLogRrs = poolSes(0.2, 0.4),
  trueSubgroupLogRrsMean = 1,
  trueSubgroupLogRrsSd = 0.5
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
computePerformance(results)
# method                  precision coverage type1 type2 reproEstimateInCi reproDisagreeCi  sign signReproSign precisionPi coveragePi type1Pi type2Pi reproEstimateInPi reproDisagreePi signPi signPiReproSign
# 1 Bayesian random effects     0.225     0.97     0  0.59              0.92            0.07  0.41          0.4        0.109       0.99       0    0.95              0.99             0     0.05            0.05
# 2 Fixed effects               1.95      0.66     0  0.08              0.42            0.3   0.92          0.89       1.95        0.66       0    0.08              0.42             0.3   0.92            0.89
# 3 Random effects              1.95      0.66     0  0.08              0.42            0.3   0.92          0.89      NA          NA          0   NA                NA              NaN   NaN               0

# 2 DBs, null effect
settings <- createSimulationSettings(
  nDatabases = 2,
  seLogRrs = c(0.2, 0.4),
  trueSubgroupLogRrsMean = 0,
  trueSubgroupLogRrsSd = 0
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
computePerformance(results)
# method                  precision coverage type1 type2 reproEstimateInCi reproDisagreeCi  sign signReproSign precisionPi coveragePi type1Pi type2Pi reproEstimateInPi reproDisagreePi signPi signPiReproSign
# 1 Bayesian random effects     0.417     0.97  0.03     0              0.79            0.06  0.03          0         0.163        1       0          0              0.96            0.02   0               0
# 2 Fixed effects               1.95      0.75  0.25     0              0.44            0.1   0.25          0.02      1.95         0.75    0.25       0              0.44            0.1    0.25            0.02
# 3 Random effects              0.950     0.89  0.11     0              0.61            0.07  0.11          0.01      0.0153       1       0          0              1               0      0               0

# 1 DB, null effect
settings <- createSimulationSettings(
  nDatabases = 1,
  seLogRrs = poolSes(0.2, 0.4),
  trueSubgroupLogRrsMean = 0,
  trueSubgroupLogRrsSd = 0
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
computePerformance(results)
# method                  precision coverage type1 type2 reproEstimateInCi reproDisagreeCi  sign signReproSign precisionPi coveragePi type1Pi type2Pi reproEstimateInPi reproDisagreePi signPi signPiReproSign
# 1 Bayesian random effects     0.229     1     0        0              0.99            0     0             0          0.111       1       0          0               1              0      0               0
# 2 Fixed effects               1.95      0.77  0.23     0              0.6             0.28  0.23          0.02       1.95        0.77    0.23       0               0.6            0.28   0.23            0.02
# 3 Random effects              1.95      0.77  0.23     0              0.6             0.28  0.23          0.02      NA          NA      NA          0              NA            NaN    NaN               0


# 10 DBs, null effect
settings <- createSimulationSettings(
  nDatabases = 10,
   trueSubgroupLogRrsMean = 0,
  trueSubgroupLogRrsSd = 0
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
plotTauPosterior(results)
computePerformance(results)
# method                  precision coverage type1 type2 reproEstimateInCi reproDisagreeCi  sign signReproSign precisionPi coveragePi type1Pi type2Pi reproEstimateInPi reproDisagreePi signPi signPiReproSign
# 1 Bayesian random effects      10.4     0.64  0.36     0              0.52            0.33  0.36          0.19        1.10       0.92    0.08       0              0.85            0.1    0.08            0.04
# 2 Fixed effects               867.      0.05  0.95     0              0.11            0.44  0.95          0.51      867.         0.05    0.95       0              0.11            0.44   0.95            0.51
# 3 Random effects               17.0     0.51  0.49     0              0.47            0.37  0.49          0.27        1.67       0.88    0.12       0              0.8             0.12   0.12            0.06



ggplot(results, aes(x = tau)) +
  geom_density()

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

