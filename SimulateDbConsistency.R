# Simulation studies to assess the impact of database differences (population and capture process)
# on the accuracy of database network estimates. We assume there is a target population to which we
# wish to generalize.

source("PredictionInterval.R")
library(dplyr)
library(ggplot2)
library(ggh4x)
library(gtools)

# Simulation settings ------------------------------------------------------------------------------
#' Title
#'
#' @param nSubgroups              Number of population subgroups.
#' @param trueSubgroupLogRrsMean  Mean of the true log RR per subgroup.
#' @param trueSubgroupLogRrsSd    SD of the true log RR per subgroup.
#' @param nDatabases              Number of databases.
#' @param subgroupMixVariety      The overall mixture of subgroups for the network. Set to 0.1 for
#'                                skewed distributions, 1 for mixed, 10 for uniform across
#'                                subgroups.
#' @param subgroupConsistency     How much each database's subgroup mixture differs from the network
#'                                centroid. Set to 1 have each database differ from centroid, 500 to
#'                                have them all be alike.
#' @param targetSubgroupMixVariety The mixture of subgroups in the target population. Set to 0.1 for
#'                                 skewed distributions, 1 for mixed, 10 for uniform across
#'                                 subgroups. Note: the target mixture is also used for the held-
#'                                 out study.
#' @param cpcConsistency          How much each database's capture process characteristics will
#'                                differ. Set to 0.1 for low consistency, or 25 for high.
#' @param targetCpcMixVariety     The mixture of capture process characteristics in the replication
#'                                study Set to 0.5 for mixed characteristics, and 0.9 or 0.1
#'                                for having either all or no of the characteristics.
#' @param seLogRrs                Standard error per database (related to sample size).
#' @param nCaptureProcessChars    Number of data capture process characteristics.
#' @param bias0                   Baseline bias.
#' @param biasCpcSd               SD of the bias associated with each data capture characteristic.
#' @param doOvers                 How often the study is repeated until p < 0.05. (Publication bias)
#'
#' @details
#' Subgroup mixtures are sampled in two stages:
#'
#' 1. Sample a centroid mixture from a Dirichlet with alpha = subgroupMixVariety. A higher alpha
#' will mean the centroid will be a more equal mix of all subgroups.
#' 2. Sample mixtures for each database from a Dirichlet with alpha = centroid * subgroupConsistency
#' + epsilon. Larger values of subgroupConsistency will mean databases will be more like the
#' centroid.
#'
#' Data capture characteristics are also sampled in two stages:
#' 1. Sample the prevalence of each characteristic from a beta distribution with alpha = 1/cpcConsistency
#' and beta = 1/cpcConsistency. Higher cpcConsistency means the prevalences will be closer to 0 and 1.
#' 2. Sample the binary characteristics for each database from a binomial with 1 trial and probability
#' equal to the prevalences sampled in step 1.
#'
#' The replication study samples the subgroup mixture from a Dirichlet with alpha = targetCpcMixVariety.
#' Higher values means the replication study has a more even mixture of all subgroups.
#'
#' The replication study samples its characteristics from a binomial with 1 trial and probability
#' equal to targetCpcMixVariety.
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
    subgroupMixVariety = 1,
    subgroupConsistency = 200,
    targetSubgroupMixVariety = 10,
    seLogRrs = runif(nDatabases, 0.1, 0.3),
    nCaptureProcessChars = 20,
    cpcConsistency = 25,
    targetCpcMixVariety = 0.9,
    bias0 = 0.1,
    biasCpcSd = 0.1,
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

  trueSubgroupLogRrs <- rnorm(settings$nSubgroups,
                              settings$trueSubgroupLogRrsMean,
                              settings$trueSubgroupLogRrsSd)
  biasCpc <- rnorm(settings$nCaptureProcessChars, 0, settings$biasCpcSd) # Bias associated with each data capture process characteristic.

  # Subgroup mixture in target population to which we wish to generalize. This will also be the
  # mixture in the replication study:
  targetMixture <- rdirichlet(1, rep(settings$targetSubgroupMixVariety, settings$nSubgroups))
  trueTargetLogRr <- sum(targetMixture * trueSubgroupLogRrs)

  # Capture process characteristics in the replication study:
  targetCpChars <- rbinom(settings$nCaptureProcessChars, 1, settings$targetCpcMixVariety)
  targetBias <- settings$bias0 + targetCpChars %*% biasCpc

  # Estimate in the repication study:
  targetSeLogRr <- tail(settings$seLogRrs, 1)
  targetLogRr <- rnorm(1, trueTargetLogRr + targetBias, targetSeLogRr)

  minP <- 1
  nStudies <- 0
  for (i in seq_len(settings$doOvers)) {
    nStudies <- nStudies + 1

    # Subgroup mixture per database:
    centroid <- gtools::rdirichlet(1, rep(settings$subgroupMixVariety, settings$nSubgroups))
    databaseMixtures <- gtools::rdirichlet(settings$nDatabases, centroid * settings$subgroupConsistency + 1e-5)

    # Data capture characteristics per database (binary):
    cpcPrevalences <- rbeta(settings$nCaptureProcessChars, shape1 = 1 / settings$cpcConsistency, shape2 = 1 / settings$cpcConsistency)
    databaseCpChars = matrix(rbinom(settings$nDatabases * settings$nCaptureProcessChars, 1, cpcPrevalences),
                             byrow = TRUE,
                             nrow = settings$nDatabases,
                             ncol = settings$nCaptureProcessChars)

    # Observed effects:
    trueDatabaseLogRr <- databaseMixtures %*% trueSubgroupLogRrs
    databaseBias <- settings$bias0 + databaseCpChars %*% biasCpc
    databaseLogRrs <- rnorm(settings$nDatabases, trueDatabaseLogRr + databaseBias, settings$seLogRrs)

    data <- data.frame(
      logRr = databaseLogRrs,
      seLogRr = settings$seLogRrs
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
  predictionInterval <- computePredictionInterval(breEstimate)
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
    logPiUb = c(feEstimate$upper, meta$upper.predict, predictionInterval[3]),
    tau = c(0, s$tau, breEstimate$tau),
    tauLb = c(0, s$lower.tau, breEstimate$tau95Lb),
    tauUb = c(0, s$upper.tau, breEstimate$tau95Ub),
    tauSample = c(list(0), list(0), list(sample(attr(breEstimate, "traces")[, 2], 1000))),
    replicationLogRr = tail(databaseLogRrs, 1),
    replicationLogLb = replicationLogRr + qnorm(0.025) * targetSeLogRr,
    replicationLogUb = replicationLogRr + qnorm(0.975) * targetSeLogRr,
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
  # priorData <- tibble(
  #   tau = c(x, x),
  #   y = c(dnorm(x, mean = 0, sd = 0.5) * 2, dnorm(x, mean = 0, sd = 0.33) * 2),
  #   `Half-normal` = rep(c("SD = 0.5", "SD = 0.33"), each = length(x))
  # )
  priorData <- tibble(
    tau = x,
    y = dnorm(x, mean = 0, sd = 0.5) * 2,
    `Half-normal` = rep(c("SD = 0.5"), length(x))
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
    coord_cartesian(xlim = c(0,1))
}

poolSes <- function(se1, se2) {
  return(sqrt(1/((1/se1^2) + (1/se2^2))))
}

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
ParallelLogger::clusterRequire(cluster, "gtools")
snow::clusterExport(cluster, "computePredictionInterval")

allRows <- list()
for (nDatabases in c(1, 2, 4, 10)) {
  for (trueEffect in c("null", "fixed", "random")) {
    for (publicationBias in c(FALSE, TRUE)) {
      for (databaseHeterogeneity in c("high", "low")) {
        message(sprintf("Simulating %d databases, true effect is %s, with%s publication bias",
                        nDatabases,
                        trueEffect,
                        if (publicationBias) "" else " no"))
        if (nDatabases == 1) {
          seLogRrs <- poolSes(0.2, 0.2)
        } else if (nDatabases == 2) {
          seLogRrs <- c(0.2, 0.2)
        } else {
          seLogRrs <- runif(nDatabases, 0.05, 0.25)
        }
        settings <- createSimulationSettings(
          nDatabases = nDatabases,
          seLogRrs = seLogRrs,
          trueSubgroupLogRrsMean = if (trueEffect == "null") 0 else log(2),
          trueSubgroupLogRrsSd = if (trueEffect == "random") 1 else 0,
          subgroupMixVariety = if (databaseHeterogeneity == "high") 10 else 1,
          subgroupConsistency = if (databaseHeterogeneity == "high") 2 else 200,
          cpcConsistency = if (databaseHeterogeneity == "high") 0.1 else 25,
          doOvers = if (publicationBias) 10 else 1
        )
        results <- ParallelLogger::clusterApply(cluster, 1:1000, simulateOne, settings = settings)
        results <- bind_rows(results)
        metrics <- computePerformance(results)
        rows <- tibble(
          nDatabases = !!nDatabases,
          trueEffect = !!trueEffect,
          publicationBias = !!publicationBias,
          databaseHeterogeneity = !!databaseHeterogeneity
        ) |>
          bind_cols(metrics)
        allRows[[length(allRows) + 1]] <- rows
      }
    }
  }
}
allRows <- bind_rows(allRows)
readr::write_csv(allRows, "SimulationResults.csv")

# Check if tau posteriors resembles empirical one
settings <- createSimulationSettings(
  nDatabases = 5,
  trueSubgroupLogRrsMean = 0,
  trueSubgroupLogRrsSd = 0,
)
results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
results <- bind_rows(results)
plotTauPosterior(results)
ggsave("TauPosterior_Simulations.png", width = 5, height = 4)
ParallelLogger::stopCluster(cluster)


# settings <- createSimulationSettings(
#   nDatabases = 4,
#   # seLogRrs = c(0.2, 0.2),
#   trueSubgroupLogRrsMean = log(2),
#   trueSubgroupLogRrsSd = 0
# )
# results <- ParallelLogger::clusterApply(cluster, 1:100, simulateOne, settings = settings)
# results <- bind_rows(results)
# computePerformance(results)

# Analyse simulation results -----------------------------------------------------------------------
allRows <- readr::read_csv("SimulationResults.csv")

# Precision by nDatabases, interval
nDatabases <- allRows |>
  distinct(nDatabases) |>
  arrange(nDatabases) |>
  mutate(x = row_number())
vizData <- bind_rows(
  allRows |>
    filter(publicationBias == FALSE) |>
    select(method, precision, nDatabases) |>
    mutate(interval = "Confidence interval"),
  allRows |>
    filter(publicationBias == FALSE) |>
    select(method, precision = precisionPi, nDatabases) |>
    mutate(interval = "Prediction interval")
) |>
  inner_join(nDatabases, by = join_by(nDatabases))

ggplot(vizData, aes(x = x, y = precision, group = x)) +
  geom_violin(scale = "width", fill = "#3f845a", alpha = 0.75) +
  scale_x_continuous("Numer of databases in study", breaks = nDatabases$x, labels = nDatabases$nDatabases) +
  scale_y_log10("Precision") +
  facet_grid(interval ~ method) +
  theme(
    panel.grid.minor.x = element_blank()
  )
ggsave("SimPrecision.png", width = 5.5, height = 4)

# Coverage by nDatabases, interval
nDatabases <- allRows |>
  distinct(nDatabases) |>
  arrange(nDatabases) |>
  mutate(x = row_number())
vizData <- bind_rows(
  allRows |>
    filter(publicationBias == FALSE) |>
    select(method, coverage, nDatabases) |>
    mutate(interval = "Confidence interval"),
  allRows |>
    filter(publicationBias == FALSE) |>
    select(method, coverage = coveragePi, nDatabases) |>
    mutate(interval = "Prediction interval")
) |>
  inner_join(nDatabases, by = join_by(nDatabases))

ggplot(vizData, aes(x = x, y = coverage, group = x)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_violin(scale = "width", fill = "#3f845a", alpha = 0.75) +
  scale_x_continuous("Numer of databases in study", breaks = nDatabases$x, labels = nDatabases$nDatabases) +
  scale_y_continuous("Coverage of the 95% interval", limits = c(0, 1)) +
  facet_grid(interval ~ method) +
  theme(
    panel.grid.minor.x = element_blank()
  )
ggsave("SimCoverage.png", width = 5.5, height = 4)

# Type 1 and 2 by nDatabases, interval
nDatabases <- allRows |>
  distinct(nDatabases) |>
  arrange(nDatabases) |>
  mutate(x = row_number())
vizData <- bind_rows(
  allRows |>
    filter(publicationBias == FALSE) |>
    transmute(method, error = if_else(trueEffect == "null", type1, type2), nDatabases, type = if_else(trueEffect == "null", "Type 1", "Type 2"), trueEffect) |>
    mutate(interval = "Confidence interval"),
  allRows |>
    filter(publicationBias == FALSE) |>
    transmute(method, error = if_else(trueEffect == "null", type1Pi, type2Pi), nDatabases, type = if_else(trueEffect == "null", "Type 1", "Type 2"), trueEffect) |>
    mutate(interval = "Prediction interval")
) |>
  inner_join(nDatabases, by = join_by(nDatabases))
ref <- vizData |>
  filter(type == "Type 1") |>
  distinct(method, interval, type) |>
  mutate(reference = 0.05)

ggplot(vizData, aes(x = x, y = error, color = trueEffect)) +
  geom_hline(aes(yintercept = reference), linetype = "dashed", data = ref) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c("#d99b77", "#3f845a", "#73655d")) +
  scale_x_continuous("Numer of databases in study", breaks = nDatabases$x, labels = nDatabases$nDatabases) +
  scale_y_continuous("Error", limits = c(0, 1)) +
  guides(color=guide_legend(title="True effect")) +
  facet_nested(interval ~ method + type) +
  theme(
    panel.grid.minor.x = element_blank(),
    legend.position = "top"
  )
ggsave("SimType1And2Error.png", width = 5.5, height = 4)

# Replication per nDatabases, interval, true effect
nDatabases <- allRows |>
  distinct(nDatabases) |>
  arrange(nDatabases) |>
  mutate(x = row_number())
vizData <- bind_rows(
  allRows |>
    transmute(method, fraction = signReproSign / sign, nDatabases, trueEffect, publicationBias) |>
    mutate(interval = "Confidence interval"),
  allRows |>
    transmute(method, fraction = signPiReproSign / signPi, nDatabases, trueEffect, publicationBias) |>
    mutate(interval = "Prediction interval")
) |>
  inner_join(nDatabases, by = join_by(nDatabases)) |>
  mutate(publicationBias = if_else(publicationBias, "Yes", "No"))

ggplot(vizData, aes(x = x, y = fraction, group = x, color = publicationBias)) +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c("#3f845a", "#d99b77", "#73655d")) +
  scale_x_continuous("Numer of databases in study", breaks = nDatabases$x, labels = nDatabases$nDatabases) +
  scale_y_continuous("Fraction of significant results that replicate", limits = c(0, 1)) +
  guides(color=guide_legend(title="Publication bias")) +
  facet_nested(interval ~ method + trueEffect) +
  theme(
    panel.grid.minor.x = element_blank(),
    legend.position = "top"
  )
ggsave("SimReplication.png", width = 6, height = 4)


# Agreement per interval
nDatabases <- allRows |>
  distinct(nDatabases) |>
  arrange(nDatabases) |>
  mutate(x = row_number())
vizData <- bind_rows(
  allRows |>
    transmute(method, reproDisagreeCi, nDatabases, trueEffect, publicationBias) |>
    mutate(interval = "Confidence interval"),
  allRows |>
    transmute(method, reproDisagreeCi = reproDisagreePi, nDatabases, trueEffect, publicationBias) |>
    mutate(interval = "Prediction interval")
) |>
  inner_join(nDatabases, by = join_by(nDatabases)) |>
  mutate(publicationBias = if_else(publicationBias, "Yes", "No"))

ggplot(vizData, aes(x = x, y = reproDisagreeCi, group = x, color = publicationBias)) +
  geom_hline(yintercept = 0.05, linetype = "dashed") +
  geom_point(alpha = 0.75) +
  scale_color_manual(values = c("#3f845a", "#d99b77", "#73655d")) +
  scale_x_continuous("Numer of databases in study", breaks = nDatabases$x, labels = nDatabases$nDatabases) +
  scale_y_continuous("Fraction rejection of null", limits = c(0, 1)) +
  guides(color=guide_legend(title="Publication bias")) +
  facet_nested(interval ~ method + trueEffect) +
  theme(
    panel.grid.minor.x = element_blank(),
    legend.position = "top"
  )
ggsave("SimReplicationSignificantDifferent.png", width = 5.5, height = 4)





# Some specific examples ---------------------------------------------------------------------------
library(dplyr)
source("ForestPlot.R")

library(EvidenceSynthesis)
set.seed(123)
populations <- simulatePopulations(createSimulationSettings(n = 20000))
labels <- paste("Site", seq_along(populations))
fitModelInDatabase <- function(population, type) {
  cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId),
                                            data = population,
                                            modelType = "cox"
  )
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
  approximation <- approximateLikelihood(cyclopsFit, parameter = "x", approximation = type)
  return(approximation)
}
approximations <- lapply(populations, fitModelInDatabase, type = "grid with gradients")
# approximations <- do.call("rbind", approximations)
estimate <- computeBayesianMetaAnalysis(approximations)
plot <- plotMetaAnalysisForest(approximations, labels, estimate, xLabel = "Hazard Ratio", showLikelihood = FALSE, fileName = "example1.svg")
ggsave("example2.svg", plot = plot,  width = 7.5, height = 1 + 5 * 0.3)
plotMetaAnalysisForest(approximations, labels, estimate, xLabel = "Hazard Ratio", fileName = "example3.svg")

vizData <- tibble(
  x = 1:100,
  y = dnorm(1:100, 50, 12)
)
ggplot(vizData, aes(x = x, y = y)) +
  geom_area(fill = "#C0504D", color = "#4F1D1B", size = 3) +
  theme_void()
ggsave("exampleNormal.svg", width = 7, height = 2)


data <- lapply(populations, fitModelInDatabase, type = "normal")
plot <- plotForest(bind_rows(data))
ggsave("example4.svg", plot, width = 8, height = 4)

data <- tibble(
  logRr = c(2, 2.1, 0.15, -0.15),
  seLogRr = c(0.04, 0.06, 0.08, 0.06)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example5.svg", plot, width = 8, height = 3)

data <- tibble(
  logRr = c(0.9),
  seLogRr = c(0.09)
)
plot <- plotForest(data)
grid::grid.draw(plot)
ggsave("example6.svg", plot, width = 8, height = 2.5)




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

# Find cool example ---------------------------------------------
library(EvidenceSynthesis)
library(dplyr)
source("ForestPlot.R")

# Make up arbitrary example with 2 databases where each rejects null, meta-analysis rejects null,
# but prediction interval includes 1:
data <- tibble(
  logRr = c(0.75, 1.1),
  seLogRr = c(0.3, 0.40)
)
plotForest(data)
maTarget <- computeBayesianMetaAnalysis(data)
piTarget <- computePredictionInterval(maTarget)
maTarget
piTarget


# Constraints for corresponding example:
# specified number of DBs
# All non-significant
# Same BRFX estimate
# PI significant

cost <- function(p, nDatabases, maTarget) {
  data <- dplyr::tibble(
    logRr = p[1:nDatabases * 2 - 1],
    seLogRr = p[1:nDatabases * 2]
  )
  if (any(data$seLogRr < 0)) {
    return(Inf)
  }
  # All non-significant
  ciLb <- data$logRr + qnorm(0.025) *  data$seLogRr
  cost <- sum(pmax(0, ciLb))

  # Same BRFX estimate
  ma <- suppressMessages(EvidenceSynthesis::computeBayesianMetaAnalysis(data, showProgressBar = FALSE))
  cost <- cost + abs(ma$mu - maTarget$mu)
  cost <- cost + abs(ma$muSe - maTarget$muSe)

  # PI significant
  predictionInterval <- computePredictionInterval(ma)
  cost <- cost + pmax(0, -predictionInterval[1])

  # print(paste(cost, paste(p, collapse = ",")))
  return(cost)
}

nDatabases <- 4

# Use genetic algorithm to find a solution
cluster <- parallel::makeCluster(10)
snow::clusterExport(cluster, "computePredictionInterval")
snow::clusterExport(cluster, "cost")
snow::clusterRequire(cluster, "dplyr")
doParallel::registerDoParallel(cluster)

# cluster <- ParallelLogger::makeCluster(10)
# snow::clusterExport(cluster, "computePredictionInterval")
# ParallelLogger::clusterRequire(cluster, "dplyr")

solution <- GA::ga(
  type = "real-valued",
  fitness = function(x, nDatabases, maTarget) -cost(x, nDatabases, maTarget),
  nDatabases = nDatabases,
  maTarget = maTarget,
  maxiter = 50,
  lower = rep(0, nDatabases * 2),
  upper = rep(2, nDatabases * 2),
  parallel = cluster
)
data <- tibble(
  logRr = solution$par[1:nDatabases * 2 - 1],
  seLogRr = solution$par[1:nDatabases * 2]
)

# Simulated annealing
# solution <- optimization::optim_sa(
#   fun = function(x) cost(x, nDatabases, maTarget),
#   start = rep(1, nDatabases * 2),
#   lower = rep(0, nDatabases * 2),
#   upper = rep(2, nDatabases * 2)
# )
# data <- tibble(
#   logRr = solution$par[1:nDatabases * 2 - 1],
#   seLogRr = solution$par[1:nDatabases * 2]
# )


plotForest(data)
ma <- computeBayesianMetaAnalysis(data)
ma
maTarget

