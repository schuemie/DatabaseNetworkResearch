# Simulation studies to assess the impact of database differences (population and capture process)
# on the accuracy of database network estimates. We assume there is a target population to which we
# wish to generalize.

library(dplyr)

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
    seLogRrs = runif(nDatabases, 0.01, 0.1),
    nCaptureProcessChars = 10,
    bias0 = 0.1,
    biasCpcSd = 0.1
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
  targetMixture <- runif(settings$nSubgroups, 0, 1)
  targetMixture <- targetMixture / sum(targetMixture) # Subgroup mixture in target population to which we wish to generalize
  trueSubgroupLogRrs <- rnorm(settings$nSubgroups,
                              settings$trueSubgroupLogRrsMean,
                              settings$trueSubgroupLogRrsSd)
  trueTargetLogRr <- sum(targetMixture * trueSubgroupLogRrs)

  databaseMixtures <- matrix(runif(settings$nDatabases * settings$nSubgroups),
                             nrow = settings$nDatabases,
                             ncol = settings$nSubgroups)
  databaseMixtures <- databaseMixtures / rowSums(databaseMixtures) # Subgroup mixture per database
  trueDatabaseLogRr <- databaseMixtures %*% trueSubgroupLogRrs

  databaseCpChars = matrix(rbinom(settings$nDatabases * settings$nCaptureProcessChars, 1, 0.5),
                           nrow = settings$nDatabases,
                           ncol = settings$nCaptureProcessChars) # Data capture process characteristics per DB (binary)
  biasCpc <- rnorm(settings$nCaptureProcessChars, 0, settings$biasCpcSd) # Bias associated with each data capture process characteristic.
  databaseBias <- settings$bias0 + databaseCpChars %*% biasCpc

  # Observed effects:
  databaseLogRrs <- rnorm(settings$nDatabases, trueDatabaseLogRr + databaseBias, settings$seLogRrs)

  if (settings$nDatabases == 1) {
    row <- dplyr::tibble(
      seed = seed,
      trueTargetLogRr = trueTargetLogRr,
      maLogRr = databaseLogRrs[1],
      maLogRrLb = maLogRr + qnorm(0.025) * settings$seLogRrs[1],
      maLogRrUb = maLogRr + qnorm(0.975) * settings$seLogRrs[1],
      piLogRrLb = maLogRrLb,
      piLogRrUb = maLogRrUb
    )
  } else {
    data <- data.frame(
      logRr = databaseLogRrs,
      seLogRr = settings$seLogRrs
    )
    maEstimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(data)
    traces <- attr(maEstimate, "traces")
    predictions <- rnorm(nrow(traces), traces[, 1], traces[, 2])
    predictionInterval <- HDInterval::hdi(predictions, credMass = 0.95)
    row <- data.frame(
      seed = seed,
      trueTargetLogRr = trueTargetLogRr,
      maLogRr = maEstimate$mu,
      maLogRrLb = maEstimate$mu95Lb,
      maLogRrUb = maEstimate$mu95Ub,
      piLogRrLb = predictionInterval[1],
      piLogRrUb = predictionInterval[2]
    )
  }
  return(row)
  # EvidenceSynthesis::plotMetaAnalysisForest(data, paste("Site", seq_len(settings$nDatabases)), maEstimate, showLikelihood = FALSE)
}

computePerformance <- function(results) {
  results <- bind_rows(results)
  results <- results |>
    mutate(coverageMa = trueTargetLogRr >= maLogRrLb & trueTargetLogRr <= maLogRrUb,
           coveragePi = trueTargetLogRr >= piLogRrLb & trueTargetLogRr <= piLogRrUb,
           type1Ma = if_else(trueTargetLogRr == 0, maLogRrLb > 0 | maLogRrUb < 0, 0),
           type2Ma = if_else(trueTargetLogRr != 0, maLogRrLb <= 0 & maLogRrUb >= 0, 0),
           type1Pi = if_else(trueTargetLogRr == 0, piLogRrLb > 0 | piLogRrUb < 0, 0),
           type2Pi = if_else(trueTargetLogRr != 0, piLogRrLb <= 0 & piLogRrUb >= 0, 0)) |>
    summarise(coverageMa = mean(coverageMa),
              coveragePi = mean(coveragePi),
              type1Ma = mean(type1Ma),
              type2Ma = mean(type2Ma),
              type1Pi = mean(type1Pi),
              type2Pi = mean(type2Pi)) |>
    as.data.frame()
  return(results)
}

poolSes <- function(se1, se2) {
  return(sqrt(1/((1/se1^2) + (1/se2^2))))
}

cluster <- ParallelLogger::makeCluster(10)

# 2 DBs, non-null effect
settings <- createSimulationSettings(
  nDatabases = 2,
  seLogRrs = c(0.2, 0.4),
  trueSubgroupLogRrsMean = 1,
  trueSubgroupLogRrsSd = 0.5
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
computePerformance(results)
# coverageMa coveragePi type1Ma type2Ma type1Pi type2Pi
#    0.97          1       0    0.19       0     0.6

# 1 DB, non-null effect
settings <- createSimulationSettings(
  nDatabases = 1,
  seLogRrs = poolSes(0.2, 0.4),
  trueSubgroupLogRrsMean = 1,
  trueSubgroupLogRrsSd = 0.5
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
computePerformance(results)
# coverageMa coveragePi type1Ma type2Ma type1Pi type2Pi
# 1       0.81       0.81       0       0       0       0

# 2 DBs, null effect
settings <- createSimulationSettings(
  nDatabases = 2,
  seLogRrs = c(0.2, 0.4),
  trueSubgroupLogRrsMean = 0,
  trueSubgroupLogRrsSd = 0
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
computePerformance(results)
# coverageMa coveragePi type1Ma type2Ma type1Pi type2Pi
#    0.99          1    0.01       0       0       0

# 1 DB, null effect
settings <- createSimulationSettings(
  nDatabases = 1,
  seLogRrs = poolSes(0.2, 0.4),
  trueSubgroupLogRrsMean = 0,
  trueSubgroupLogRrsSd = 0
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
computePerformance(results)
# coverageMa coveragePi type1Ma type2Ma type1Pi type2Pi
# 1       0.76       0.76    0.24       0    0.24       0


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

