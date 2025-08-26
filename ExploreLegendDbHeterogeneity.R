# Fetch information from LEGEND T2DM to explore how much heterogeneity in effect estimates exists
# between databases.

library(DatabaseConnector)
library(dplyr)
library(ggplot2)
library(ggh4x)

# Fetch estimates from LEGEND T2DM -----------------------------------------------------------------
connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)
schema <- "legendt2dm_class_results_v2"

# executeSql(connection, "COMMIT;")
connection <- connect(connectionDetails)
negativeControlIds <- renderTranslateQuerySql(
  connection = connection,
  sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
  schema = schema,
  snakeCaseToCamelCase = TRUE
)[, 1]
sql <- "
SELECT database_id,
    target_id,
    target.exposure_name AS target_name,
    comparator_id,
    comparator.exposure_name AS comparator_name,
    outcome_id,
    analysis_id,
    CASE WHEN analysis_id = 7 THEN 'unadjusted' ELSE 'PS matched' END AS analysis_name,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
INNER JOIN @schema.exposure_of_interest target
  ON target_id = target.exposure_id
INNER JOIN @schema.exposure_of_interest comparator
  ON comparator_id = comparator.exposure_id
WHERE se_log_rr IS NOT NULL
    AND analysis_id IN (7, 8)
    AND target.exposure_name LIKE '%main ot2'
    AND comparator.exposure_name LIKE '%main ot2'
    AND database_id NOT LIKE 'Meta-analysis%';
"
estimates <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     snakeCaseToCamelCase = TRUE)
estimates <- estimates |>
  mutate(databaseId = as.factor(databaseId),
         analysisName = as.factor(analysisName),
         targetName = as.factor(gsub(" main ot2", "", targetName)),
         comparatorName = as.factor(gsub(" main ot2", "", comparatorName)),
         negativeControl = outcomeId %in% negativeControlIds)
disconnect(connection)
legendLabel <- "T2dm"
minDatabases <- 3
saveRDS(estimates, sprintf("estimates_%s.rds", legendLabel))


# Fetch estimates from LEGEND Hypertension ---------------------------------------------------------
connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)
schema <- "legend"

connection <- connect(connectionDetails)
negativeControlIds <- renderTranslateQuerySql(
  connection = connection,
  sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
  schema = schema,
  snakeCaseToCamelCase = TRUE
)[, 1]
outcomeOfInterestIds <- renderTranslateQuerySql(
  connection = connection,
  sql = "SELECT outcome_id FROM @schema.outcome_of_interest;",
  schema = schema,
  snakeCaseToCamelCase = TRUE
)[, 1]
sql <- "
SELECT DISTINCT database_id,
    target_id,
    target.exposure_name AS target_name,
    comparator_id,
    comparator.exposure_name AS comparator_name,
    outcome_id,
    analysis_id,
    CASE WHEN analysis_id = 1 THEN 'PS stratification' ELSE 'PS matching' END AS analysis_name,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
INNER JOIN @schema.single_exposure_of_interest target
  ON target_id = target.exposure_id
INNER JOIN @schema.exposure_group target_group
  ON target_group.exposure_id = target.exposure_id
INNER JOIN @schema.single_exposure_of_interest comparator
  ON comparator_id = comparator.exposure_id
INNER JOIN @schema.exposure_group comparator_group
  ON comparator_group.exposure_id = comparator.exposure_id
WHERE se_log_rr IS NOT NULL
    AND analysis_id IN (1, 3)
    AND target_group.exposure_group = 'Drug class'
    AND comparator_group.exposure_group = 'Drug class'
    AND outcome_id IN (@outcome_ids)
    AND database_id != 'Meta-analysis';
"
estimates <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     outcome_ids = c(negativeControlIds, outcomeOfInterestIds),
                                     snakeCaseToCamelCase = TRUE)
estimates <- estimates |>
  mutate(databaseId = as.factor(databaseId),
         analysisName = as.factor(analysisName),
         targetName = as.factor(targetName),
         comparatorName = as.factor(comparatorName),
         negativeControl = outcomeId %in% negativeControlIds)
disconnect(connection)
legendLabel <- "Htn"
minDatabases <- 6
saveRDS(estimates, sprintf("estimates_%s.rds", legendLabel))


# Fetch estimates from ASD-004 SCCS --------------------------------------------
connectionDetails <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  server = keyring::key_get("ohdaResultsServer"),
  user = keyring::key_get("ohdaResultsUser"),
  password = keyring::key_get("ohdaResultsPassword")
)
schema <- "asd_004"

connection <- connect(connectionDetails)
sql <- "SELECT c.era_id,
  outcome_id,
  nesting_cohort_id,
  unblind_for_evidence_synthesis AS unblind,
  r.analysis_id,
  a.description AS analysis_name,
  d.cdm_source_abbreviation AS database_name,
  CASE WHEN true_effect_size = 1 THEN 1 ELSE 0 END AS negative_control,
  log_rr,
  se_log_rr
FROM @schema.sccs_result r
INNER JOIN @schema.sccs_diagnostics_summary ds
  ON ds.exposures_outcome_set_id = r.exposures_outcome_set_id
    AND ds.covariate_id = r.covariate_id
    AND ds.analysis_id = r.analysis_id
    AND ds.database_id = r.database_id
INNER JOIN @schema.sccs_exposures_outcome_set eos
  ON ds.exposures_outcome_set_id = eos.exposures_outcome_set_id
INNER JOIN (
  SELECT DISTINCT exposures_outcome_set_id,
    covariate_id,
    analysis_id,
    era_id
  FROM @schema.sccs_covariate
  ) c
  ON ds.exposures_outcome_set_id = c.exposures_outcome_set_id
    AND ds.covariate_id = c.covariate_id
    AND ds.analysis_id = c.analysis_id
INNER JOIN @schema.cg_cohort_definition cd
  ON eos.outcome_id = cd.cohort_definition_id
INNER JOIN @schema.sccs_analysis a
  ON r.analysis_id = a.analysis_id
INNER JOIN @schema.database_meta_data d
  ON r.database_id = d.database_id
INNER JOIN @schema.sccs_exposure e
  ON r.exposures_outcome_set_id = e.exposures_outcome_set_id
    AND c.era_id = e.era_id
WHERE se_log_rr IS NOT NULL;"
estimates <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     snakeCaseToCamelCase = TRUE)
estimates <- estimates |>
  transmute(databaseId = as.factor(databaseName),
            analysisName = as.factor(analysisName),
            targetId = eraId,
            targetName = as.factor(eraId),
            comparatorId = eraId,
            comparatorName = as.factor(eraId),
            outcomeId = outcomeId,
            negativeControl = negativeControl == 1,
            logRr,
            seLogRr)
disconnect(connection)
legendLabel <- "Sccs"
minDatabases <- 6
saveRDS(estimates, sprintf("estimates_%s.rds", legendLabel))

# Compute tau --------------------------------------------------------------------------------------
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
  group_by(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl)|>
  group_split()

# group = groups[[1]]
computeTau <- function(group) {
  keyRow <- group |>
    head(1) |>
    select(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl)
  maEstimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(group, showProgressBar = FALSE)
  traces <- attr(maEstimate, "traces")
  tauSample <- sample(traces[, 2], 100)
  row <- maEstimate |>
    select(mu, mu95Lb, mu95Ub, tau, tau95Lb, tau95Ub) |>
    bind_cols(keyRow) |>
    mutate(signficant = mu95Lb > 0 | mu95Ub < 0) |>
    mutate(type = if_else(negativeControl,
                          "Negative control",
                          if_else(signficant,
                                  "Outcome of interest and significant",
                                  "Outcome of interest")))
  return(list(row = row, tauSample = tauSample))
}

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
tauRowsAndSamples <- ParallelLogger::clusterApply(cluster, groups, computeTau)
ParallelLogger::stopCluster(cluster)

taus <- bind_rows(lapply(tauRowsAndSamples, function(x) x$row))
groups <- taus |>
  distinct(analysisName, type)
tauSamples <- list()
for (i in seq_len(nrow(groups))) {
  row <- groups[i, ]
  fun <- function(x) {
    if (x$row$analysisName == row$analysisName && x$row$type == row$type) {
      return(x$tauSample)
    } else {
      return(c())
    }
  }
  tauSample <- do.call(c, lapply(tauRowsAndSamples, fun))
  if (length(tauSample) > 10000) {
    tauSample <- sample(tauSample, 10000)
  }
  tauSamples[[sprintf("%s-%s", row$analysisName, row$type)]] <- tauSample
}

saveRDS(taus, sprintf("taus_%s.rds", legendLabel))
saveRDS(tauSamples, sprintf("tauSamples_%s.rds", legendLabel))


# Plot tau distributions ---------------------------------------------------------------------------
taus <- readRDS(sprintf("taus_%s.rds", legendLabel))

taus |>
  mutate(se = (tau95Ub - tau95Lb) / 2*qnorm(0.975)) |>
  group_by(type, analysisName) |>
  summarise(medianSe = median(se), .groups = "drop")

priorTimesX <- function(x){
  (dnorm(x, 0, 0.5) * 2) * x
}
expectedTauUnderPrior <- integrate(priorTimesX, lower = 0, upper = Inf)$value

vizData <- taus |>
  mutate(se = (tau95Ub - tau95Lb) / 2*qnorm(0.975)) |>
  filter(se < 0.4) |>
  mutate(type = gsub(" and", "\nand", type),
         analysisName = if_else(analysisName == "unadjusted", "Unadjusted", analysisName))

ggplot(vizData, aes(y = tau)) +
  geom_hline(yintercept = 0, size = 1) +
  geom_hline(yintercept = expectedTauUnderPrior, linetype = "dashed") +
  geom_boxplot(fill = "#3f845a", alpha = 0.75) +
  scale_y_continuous("Tau") +
  facet_nested(~ type + analysisName) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )
ggsave(sprintf("TauDistributions_%s.png", legendLabel), width = 6, height = 5)

# Plot tau posterior -------------------------------------------------------------------------------
tauSamples <- readRDS(sprintf("tauSamples_%s.rds", legendLabel))

x <- seq(from = 0, to = 2, length.out = 100)
priorData <- tibble(
  tau = c(x, x),
  y = c(dnorm(x, mean = 0, sd = 0.5) * 2, dnorm(x, mean = 0, sd = 0.33) * 2),
  `Half-normal` = rep(c("SD = 0.5", "SD = 0.33"), each = length(x))
)

vizData <- list()
for (i in seq_along(tauSamples)) {
  row <- strsplit(names(tauSamples[i]), "-")[[1]]
  names(row) <- c("analysisName", "type")
  row <- as_tibble(t(row)) |>
    mutate(type = gsub(" and", "\nand", type),
           analysisName = if_else(analysisName == "unadjusted", "Unadjusted", analysisName))
  vizData[[i]] <- bind_cols(row, tibble(tau = tauSamples[[i]]))
}
vizData <- bind_rows(vizData)
ggplot(vizData, aes(x = tau)) +
  geom_density(fill = "#3f845a", alpha = 0.75) +
  geom_line(aes(y = y, linetype = `Half-normal`), data = priorData) +
  scale_x_continuous("Tau") +
  scale_y_continuous("Density") +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  coord_cartesian(xlim = c(0,2)) +
  facet_grid(analysisName ~ type)
ggsave(sprintf("TauPosteriors_%s.png", legendLabel), width = 6, height = 5)
# ggsave(sprintf("TauPosteriorsCalibrated_%s.png", legendLabel), width = 6, height = 5)

# Compute correlation matrix between databases using Bayesian model --------------------------------
library(brms)
library(dplyr)
library(tidyr)


# Simulation
set.seed(42)
nOutcomes <- 50
trueMu1 <- 0.1  # True average log(RR) for Db1
trueMu2 <- 0.15 # True average log(RR) for Db2
trueSigma1 <- 0.3 # True standard deviation of log(RR) for Db1
trueSigma2 <- 0.4 # True standard deviation of log(RR) for Db2
trueRho <- 0.2   # This is the true correlation we want to recover
covMatrix <- matrix(c(trueSigma1^2, trueRho * trueSigma1 * trueSigma2,
                      trueRho * trueSigma1 * trueSigma2, trueSigma2^2),
                    nrow = 2)
trueLogRrs <- MASS::mvrnorm(nOutcomes, mu = c(trueMu1, trueMu2), Sigma = covMatrix)
simulatedDataWide <- tibble(
  outcomeId = 1:nOutcomes,
  trueLogRrDb1 = trueLogRrs[, 1],
  trueLogRrDb2 = trueLogRrs[, 2],
  seLogRrDb1 = rlnorm(nOutcomes, meanlog = -2, sdlog = 0.5),
  seLogRrDb2 = rlnorm(nOutcomes, meanlog = -1.8, sdlog = 0.6)
) %>%
  mutate(
    logRrDb1 = rnorm(nOutcomes, mean = trueLogRrDb1, sd = seLogRrDb1),
    logRrDb2 = rnorm(nOutcomes, mean = trueLogRrDb2, sd = seLogRrDb2)
  )
dataLong <- simulatedDataWide %>%
  select(outcomeId, logRrDb1, seLogRrDb1, logRrDb2, seLogRrDb2) %>%
  pivot_longer(
    cols = -outcomeId,
    names_to = c(".value", "databaseId"),
    names_pattern = "(.+)(Db[12])"
  )
print(head(dataLong))

# Real data
estimates <- readRDS(sprintf("estimates_%s.rds", legendLabel))
analysisName <- unique(estimates$analysisName)[2]
analysisName

atLeastNddbs <- estimates |>
  filter(analysisName == !!analysisName) |>
  group_by(targetId, comparatorId, outcomeId) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n >= minDatabases) |>
  arrange(targetId, comparatorId, outcomeId) |>
  mutate(dummyOutcomeId = row_number())

dataLong<- estimates |>
  filter(analysisName == !!analysisName) |>
  inner_join(atLeastNddbs, by = join_by(targetId, comparatorId, outcomeId)) |>
  select(databaseId, outcomeId = dummyOutcomeId, logRr, seLogRr)


# The formula now models logRr, accounting for measurement error via se(seLogRr).
# - `0 + databaseId`: This estimates the average logRr for each databaseId (μ1, μ2) without a global intercept.
# - `(0 + databaseId | outcomeId)`: This is the key part. It allows the effect for each
#   outcome to vary from the database average. Crucially, it estimates the
#   correlation between these outcome-specific variations across the databases.
#   This correlation is our parameter of interest, rho (ρ).
modelFormula <- bf(logRr | se(seLogRr) ~ 0 + databaseId + (0 + databaseId | outcomeId))

# Set weakly informative priors for the new model parameters.
# 'b': The fixed effects (the average logRRs for each databaseId).
# 'sd': The standard deviations of the outcome-specific effects (σ1, σ2).
# 'cor': The correlation matrix for the outcome-specific effects (contains ρ).
priors <- c(
  prior(normal(0, 2), class = "b"),
  prior(exponential(1), class = "sd"),
  prior(lkj(1), class = "cor")
)

# Fit the Bayesian model.
# This might take a few minutes to run.
correlationModelFit <- brm(
  formula = modelFormula,
  data = dataLong,
  prior = priors,
  chains = 4,
  iter = 11000,
  warmup = 1000,
  cores = parallel::detectCores(),
  control = list(adapt_delta = 0.95)
)

print(summary(correlationModelFit))
results <- summary(correlationModelFit)$random[[1]]
results$name <- rownames(results)
rownames(results) <- NULL
results <- results |>
  filter(grepl("^cor", name)) |>
  mutate(databaseId1 = stringr::str_extract(name, "cor\\(databaseId([a-zA-Z]+),databaseId([a-zA-Z]+)\\)", group = 1),
         databaseId2 = stringr::str_extract(name, "cor\\(databaseId([a-zA-Z]+),databaseId([a-zA-Z]+)\\)", group = 2)) |>
  select(databaseId1, databaseId2, rho = Estimate, lb = "l-95% CI", ub = "u-95% CI")

fullMatrix <- bind_rows(results,
                        results |>
                          mutate(temp = databaseId1,
                                 databaseId1 = databaseId2,
                                 databaseId2 = temp),
                        tibble(
                          databaseId1 = unique(c(results$databaseId1, results$databaseId2)),
                          databaseId2 = databaseId1,
                          rho = 1
                        )) |>
  select(-temp)
fullMatrix <- fullMatrix |>
  pivot_wider(id_cols = "databaseId1", names_from = "databaseId2", values_from = "rho", names_sort = TRUE)
readr::write_csv(fullMatrix, sprintf("LegendDbCorrelation_%s_%s.csv", gsub(" ", "", analysisName), legendLabel))

#
# posteriorSamples <- posterior_samples(correlationModelFit, pars = "cor_outcomeId__databaseIdDb1__databaseIdDb2")
# colnames(posteriorSamples) <- "rho"
# posteriorMedian <- median(posteriorSamples$rho)
# credibleInterval <- quantile(posteriorSamples$rho, probs = c(0.025, 0.975))
#
# cat("\n--- Final Results ---\n")
# cat("Posterior Median for Correlation (ρ):", round(posteriorMedian, 3), "\n")
# cat("95% Credible Interval for Correlation (ρ): [", round(credibleInterval[1], 3), ",", round(credibleInterval[2], 3), "]\n")
# cat("The true value used in the simulation was:", trueRho, "\n")
#
#
# # Visualize the posterior distribution of the correlation parameter.
# ggplot(posteriorSamples, aes(x = rho)) +
#   geom_density(fill = "skyblue", alpha = 0.7) +
#   geom_vline(xintercept = posteriorMedian, color = "blue", linetype = "dashed", size = 1) +
#   geom_vline(xintercept = credibleInterval, color = "red", linetype = "dotted") +
#   labs(
#     title = "Posterior Distribution of the Correlation (ρ)",
#     subtitle = paste0("Median: ", round(posteriorMedian, 2),
#                       ", 95% CrI: [", round(credibleInterval[1], 2), ", ", round(credibleInterval[2], 2), "]"),
#     x = "Correlation (ρ)",
#     y = "Density"
#   ) +
#   theme_minimal()


# Compute tau for calibrated estimates -------------------------------------------------------------
# We shouldn't really do this as systematic error will be correlated between databases, but it can
# give a rough idea of how much heterogeneity remains once we adjust for measured systematic error
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
  null <- EmpiricalCalibration::fitNull(ncs$logRr, ncs$seLogRr)
  estimates <- EmpiricalCalibration::calibrateConfidenceInterval(
    group$logRr,
    group$seLogRr,
    EmpiricalCalibration::convertNullToErrorModel(null)
  )
  estimates <- group |>
    select(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl) |>
    bind_cols(estimates)

  outcomeGroups <- estimates |>
    group_by(outcomeId) |>
    group_split()

  tauSamplesNcs <- list()
  tauSamplesHois <- list()
  rows <- list()
  # i = 1
  for (i in seq_along(outcomeGroups)) {
    outcomeGroup = outcomeGroups[[i]]
    keyRow <- outcomeGroup |>
      head(1) |>
      select(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl)
    maEstimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(outcomeGroup, showProgressBar = FALSE)
    traces <- attr(maEstimate, "traces")
    tauSample <- sample(traces[, 2], 100)
    row <- keyRow |>
      mutate(tau = maEstimate$tau,
             tau95Lb = maEstimate$tau95Lb,
             tau95Ub = maEstimate$tau95Ub)
    rows[[i]] <- row
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
ParallelLogger::clusterRequire(cluster, "dplyr")
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
for (key in names(tauSamples)) {
  tauSamples[[key]] <- sample(tauSamples[[key]], 10000)
}

saveRDS(taus, sprintf("tausCalibrated_%s.rds", legendLabel))
saveRDS(tauSamples, sprintf("tauSamplesCalibrated_%s.rds", legendLabel))


# tauSamples <- readRDS(sprintf("tauSamples_%s.rds", legendLabel))
# median(tauSamples[[1]])
# tauSamples <- readRDS(sprintf("tauSamplesCalibrated_%s.rds", legendLabel))
