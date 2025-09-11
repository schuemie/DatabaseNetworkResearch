# Fetch information from LEGEND T2DM to explore how much heterogeneity in effect estimates exists
# between databases.

library(DatabaseConnector)
library(dplyr)
library(ggplot2)
library(ggh4x)

# Fetch estimates from LEGEND T2DM -----------------------------------------------------------------
legendLabel <- "T2dm"
minDatabases <- 9

connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)
schema <- "legendt2dm_class_results"

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
    cohort_method_result.outcome_id,
    COALESCE(outcome_of_interest.outcome_name, negative_control_outcome.outcome_name) AS outcome_name,
    analysis_id,
    CASE WHEN analysis_id = 7 THEN 'unadjusted' ELSE 'PS matched' END AS analysis_name,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
INNER JOIN @schema.exposure_of_interest target
  ON target_id = target.exposure_id
INNER JOIN @schema.exposure_of_interest comparator
  ON comparator_id = comparator.exposure_id
LEFT JOIN @schema.outcome_of_interest
  ON cohort_method_result.outcome_id = outcome_of_interest.outcome_id
LEFT JOIN @schema.negative_control_outcome
  ON cohort_method_result.outcome_id = negative_control_outcome.outcome_id
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
         outcomeName = as.factor(gsub("_", " ", gsub("outcome/", "", outcomeName))),
         negativeControl = outcomeId %in% negativeControlIds)
disconnect(connection)
saveRDS(estimates, sprintf("estimates_%s.rds", legendLabel))


# Fetch estimates from LEGEND Hypertension ---------------------------------------------------------
legendLabel <- "Htn"
minDatabases <- 6
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
    cohort_method_result.outcome_id,
    COALESCE(outcome_of_interest.outcome_name, negative_control_outcome.outcome_name) AS outcome_name,
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
LEFT JOIN @schema.outcome_of_interest
  ON cohort_method_result.outcome_id = outcome_of_interest.outcome_id
LEFT JOIN @schema.negative_control_outcome
  ON cohort_method_result.outcome_id = negative_control_outcome.outcome_id
WHERE se_log_rr IS NOT NULL
    AND analysis_id IN (1, 3)
    AND target_group.exposure_group = 'Drug class'
    AND comparator_group.exposure_group = 'Drug class'
    AND cohort_method_result.outcome_id IN (@outcome_ids)
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
         outcomeName = as.factor(outcomeName),
         negativeControl = outcomeId %in% negativeControlIds)
disconnect(connection)

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


# Fetch estimates for Semanaion --------------------------------------------------------------------
connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)

schema <- "semanaion"

# executeSql(connection, "COMMIT;")
connection <- connect(connectionDetails)
negativeControlIds <- renderTranslateQuerySql(
  connection = connection,
  sql = "SELECT DISTINCT outcome_id FROM @schema.cm_target_comparator_outcome WHERE true_effect_size = 1;",
  schema = schema,
  snakeCaseToCamelCase = TRUE
)[, 1]
sql <- "
SELECT cdm_source_abbreviation AS database_id,
    target_id,
    target.cohort_name AS target_name,
    comparator_id,
    comparator.cohort_name AS comparator_name,
    outcome_id,
    log_rr,
    se_log_rr
FROM @schema.cm_result
INNER JOIN @schema.cg_cohort_definition target
  ON target_id = target.cohort_definition_id
INNER JOIN @schema.cg_cohort_definition comparator
  ON comparator_id = comparator.cohort_definition_id
INNER JOIN @schema.database_meta_data
  ON cm_result.database_id = database_meta_data.database_id
WHERE se_log_rr IS NOT NULL;
"
estimates <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     snakeCaseToCamelCase = TRUE)
estimates <- estimates |>
  mutate(databaseId = as.factor(databaseId),
         targetName = as.factor(gsub("and ", "with ", gsub("New user of |with prior T2DM |treatment ", "", targetName))),
         comparatorName = as.factor(gsub("and ", "with ", gsub("New user of |with prior T2DM |treatment ", "", comparatorName))),
         negativeControl = outcomeId %in% negativeControlIds,
         analysisId = 1,
         analysisName = "PS matching")

disconnect(connection)
legendLabel <- "Naion"
minDatabases <- 10
saveRDS(estimates, sprintf("estimates_%s.rds", legendLabel))

# Compute tau --------------------------------------------------------------------------------------
estimates <- readRDS(sprintf("estimates_%s.rds", legendLabel))
if (file.exists(sprintf("Diagnostics_%s.rds", legendLabel))) {
  analysisName <- unique(estimates$analysisName)
  analysisName <- analysisName[grepl("match", analysisName)]
  diagnostics <- readRDS(sprintf("Diagnostics_%s.rds", legendLabel))
  estimates <- estimates |>
    left_join(diagnostics |>
                mutate(analysisName = !! analysisName))
}

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

  if (all(!is.na(group$unblind))) {
    unblindedGroup <- group |>
      filter(unblind)
    if (nrow(unblindedGroup) == 0) {
      maEstimateUnblinded <- tibble(
        tauUnblinded = as.numeric(NA),
        tau95LbUnblinded = as.numeric(NA),
        tau95UbUnblinded = as.numeric(NA)
      )
      tauSampleUnblinded <- NULL
    } else {
      maEstimateUnblinded <- EvidenceSynthesis::computeBayesianMetaAnalysis(unblindedGroup, showProgressBar = FALSE) |>
        select(tauUnblinded = tau,
               tau95LbUnblinded = tau95Lb,
               tau95UbUnblinded = tau95Ub)
      traces <- attr(maEstimateUnblinded, "traces")
      tauSampleUnblinded <- sample(traces[, 2], 100)
    }
  } else {
    maEstimateUnblinded <- tibble(
      tauUnblinded = as.numeric(NA),
      tau95LbUnblinded = as.numeric(NA),
      tau95UbUnblinded = as.numeric(NA)
    )
    tauSampleUnblinded <- NULL
  }
  row <- maEstimate |>
    select(mu, mu95Lb, mu95Ub, tau, tau95Lb, tau95Ub) |>
    bind_cols(maEstimateUnblinded) |>
    bind_cols(keyRow) |>
    mutate(signficant = mu95Lb > 0 | mu95Ub < 0) |>
    mutate(type = if_else(negativeControl,
                          "Negative control",
                          if_else(signficant,
                                  "Outcome of interest and significant",
                                  "Outcome of interest")))

  return(list(row = row, tauSample = tauSample, tauSampleUnblinded = tauSampleUnblinded))
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

  fun2 <- function(x) {
    if (x$row$analysisName == row$analysisName && x$row$type == row$type) {
      return(x$tauSampleUnblinded)
    } else {
      return(c())
    }
  }
  tauSampleUnblinded <- do.call(c, lapply(tauRowsAndSamples, fun2))
  if (!is.null(tauSamplesUnblinded)) {
    if (length(tauSampleUnblinded) > 10000) {
      tauSampleUnblinded <- sample(tauSampleUnblinded, 10000)
    }
    tauSamples[[sprintf("%s-%s-Unblinded", row$analysisName, row$type)]] <- tauSampleUnblinded
  }
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
  geom_boxplot(fill = "#336B91", alpha = 0.75) +
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

vizData <- list()
for (i in seq_along(tauSamples)) {
  row <- strsplit(names(tauSamples)[i], "-")[[1]]
  if (length(row) == 2) {
    row[3] <- FALSE
  }
  names(row) <- c("analysisName", "type", "unblinded")

  row <- as_tibble(t(row)) |>
    mutate(type = gsub(" and", "\nand", type),
           analysisName = if_else(analysisName == "unadjusted", "Unadjusted", analysisName),
           unblinded = if_else(unblinded == "Unblinded", "Unblinded\nonly", "All"))
  vizData[[i]] <- bind_cols(row, tibble(tau = tauSamples[[i]]))
}
vizData <- bind_rows(vizData)
ggplot(vizData, aes(x = tau)) +
  geom_density(fill = "#336B91", alpha = 0.75) +
  geom_line(aes(y = y, linetype = `Half-normal`), data = priorData) +
  scale_x_continuous("Tau") +
  scale_y_continuous("Density") +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  coord_cartesian(xlim = c(0,1)) +
  facet_nested(analysisName ~ type + unblinded)
ggsave(sprintf("TauPosteriors_%s.png", legendLabel), width = 8, height = 5)
# ggsave(sprintf("TauPosteriorsCalibrated_%s.png", legendLabel), width = 6, height = 5)

ggplot(vizData |> filter(analysisName == "PS matched", type == "Negative control", unblinded == "All"), aes(x = tau)) +
  geom_density(fill = "#336B91", alpha = 0.75) +
  geom_line(aes(y = y, linetype = `Half-normal`), data = priorData) +
  scale_x_continuous("Tau") +
  scale_y_continuous("Density") +
  scale_linetype_manual(values = c("dashed", "dotted")) +
  coord_cartesian(xlim = c(0,1))
ggsave(sprintf("TauPosteriorsNcOnly_%s.svg", legendLabel), width = 6, height = 5)

vizData2 <- vizData |>
  filter(grepl("match", analysisName),
         type %in% c("Negative control", "Outcome of interest\nand significant"),
         unblinded != "All") |>
  mutate(type = factor(type, levels = c("Negative control",
                                        "Outcome of interest\nand significant",
                                        "Outcome of interest\nnot significant")))
ggplot(vizData2, aes(x = tau, group = type, fill = type)) +
  geom_density(alpha = 0.5) +
  scale_x_continuous("\u03C4") +
  scale_y_continuous("Density") +
  scale_fill_manual(values = c("#EB6622", "#336B91", "#11A08A")) +
  coord_cartesian(xlim = c(0, 0.5)) +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )
ggsave(sprintf("TauPosteriorsOverlay_%s.png", legendLabel), width = 4, height = 4)




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

analysisName <- unique(estimates$analysisName)[1]
analysisName

estimates |>
  group_by(databaseId, analysisName) |>
  summarise(
    p05 = quantile(seLogRr, 0.05),
    p25 = quantile(seLogRr, 0.25),
    p50 = quantile(seLogRr, 0.55),
    p75 = quantile(seLogRr, 0.75),
    p95 = quantile(seLogRr, 0.95)
  ) |>
  arrange(p50) |>
  print(n = 100)

poweredDbs <- estimates  |>
  group_by(databaseId, analysisName) |>
  summarise(medianSe = median(seLogRr), .groups = "drop") |>
  filter(medianSe < 1) |>
  group_by(databaseId) |>
  summarise(analysisCount = n()) |>
  filter(analysisCount == length(unique(estimates$analysisName)))

atLeastNddbs <- estimates |>
  filter(analysisName == !!analysisName) |>
  group_by(targetId, comparatorId, outcomeId) |>
  summarise(n = n(), .groups = "drop") |>
  # filter(n >= minDatabases) |>
  arrange(targetId, comparatorId, outcomeId) |>
  mutate(dummyOutcomeId = row_number())

dataLong<- estimates |>
  filter(analysisName == !!analysisName) |>
  inner_join(atLeastNddbs, by = join_by(targetId, comparatorId, outcomeId)) |>
  inner_join(poweredDbs, by = join_by(databaseId)) |>
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
  mutate(databaseId1 = stringr::str_extract(name, "cor\\(databaseId([a-zA-Z0-9_]+),databaseId([a-zA-Z0-9_]+)\\)", group = 1),
         databaseId2 = stringr::str_extract(name, "cor\\(databaseId([a-zA-Z0-9_]+),databaseId([a-zA-Z0-9_]+)\\)", group = 2)) |>
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


# Non-Bayesian correlation metric ------------------------------------------------------------------
# source("ComputeDatabaseCorrelation.R")
# estimates <- readRDS(sprintf("estimates_%s.rds", legendLabel))
# analysisName <- unique(estimates$analysisName)[2]
# analysisName
#
# atLeastNddbs <- estimates |>
#   filter(analysisName == !!analysisName) |>
#   group_by(targetId, comparatorId, outcomeId) |>
#   summarise(n = n(), .groups = "drop") |>
#   filter(n >= minDatabases) |>
#   arrange(targetId, comparatorId, outcomeId) |>
#   mutate(dummyOutcomeId = row_number())
#
# dataLong <- estimates |>
#   filter(analysisName == !!analysisName) |>
#   inner_join(atLeastNddbs, by = join_by(targetId, comparatorId, outcomeId)) |>
#   select(databaseId, outcomeId = dummyOutcomeId, logRr, seLogRr)
#
# model <- estimateCorrelationMatrix(dataLong)
# model

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

# Explore individual TCO examples ------------------------------------------------------------------
source("ForestPlot.R")
estimates <- readRDS(sprintf("estimates_%s.rds", legendLabel))
analysisName <- unique(estimates$analysisName)
analysisName <- analysisName[grepl("match", analysisName)]
estimates <- estimates |>
  filter(analysisName == !!analysisName)

if (file.exists(sprintf("Diagnostics_%s.rds", legendLabel)) && grepl("match", analysisName)) {
  diagnostics <- readRDS(sprintf("Diagnostics_%s.rds", legendLabel))
  estimates <- estimates |>
    left_join(diagnostics)
}

atLeastNdbs <- estimates |>
  group_by(targetId, comparatorId, outcomeId) |>
  summarise(nDatabases = n(), .groups = "drop") |>
  filter(nDatabases >= minDatabases)

highPowerTcos <- estimates |>
  inner_join(
    atLeastNdbs,
    by = join_by(targetId, comparatorId, outcomeId)
  ) |>
  group_by(targetId, targetName, comparatorId, comparatorName, outcomeId, outcomeName, negativeControl) |>
  summarise(maxSeLogRr = max(seLogRr), .groups = "drop") |>
  arrange(maxSeLogRr)
if ("JMDC" %in% estimates$databaseId) {
  hasJmdc <- estimates |>
    filter(analysisName == !!analysisName, databaseId == "JMDC") |>
    distinct(targetId, comparatorId, outcomeId)
  highPowerTcos <- highPowerTcos |>
    inner_join(hasJmdc)
}

highPowerTcos

# Individual examples
example <- highPowerTcos[35, ]
example <- highPowerTcos |>
  filter(targetId == 102100000, comparatorId == 402100000, outcomeId == 6)
example <- highPowerTcos |>
  filter(targetId == 1, comparatorId == 2, outcomeId == 2)
example
dbGroupings <- tibble(
  databaseId = c("CCAE", "CUIMC", "Germany_DA", "MDCD", "MDCR", "OptumDod", "OptumEHR", "SIDIAP", "UK_IMRD", "US_Open_Claims", "VA-OMOP", "CUMC", "IMSG", "JMDC", "NHIS_NSC", "Optum", "Panther"),
  country = c("USA", "USA", "Germany", "USA", "USA", "USA", "USA", "Spain", "UK", "USA", "USA", "USA", "Germany", "Japan", "South Korea", "USA", "USA"),
  type = c("Claims", "EHR", "EHR", "Claims", "Claims", "Claims", "EHR", "SIDIAP", "EHR", "Claims", "EHR", "EHR", "EHR", "Claims", "Claims", "Claims", "EHR"),
)
exampleEstimates <- estimates |>
  filter(analysisName == !!analysisName) |>
  inner_join(example, by = join_by(targetId, targetName, comparatorId, comparatorName, outcomeId, outcomeName, negativeControl)) |>
  inner_join(dbGroupings, by = join_by(databaseId)) |>
  arrange(databaseId)

plotForest(data = exampleEstimates,
           labels = exampleEstimates$databaseId,
           exclude = !exampleEstimates$unblind,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           fileName = "Symposium/exampleTco.png")

subset <- filter(exampleEstimates, type == "EHR")
plotForest(data = subset,
           labels = subset$databaseId,
           exclude = !exampleEstimates$unblind,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           fileName = "Symposium/exampleTco_EHRs.png")

subset <- filter(exampleEstimates, type == "Claims")
plotForest(data = subset,
           labels = subset$databaseId,
           exclude = !exampleEstimates$unblind,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           fileName = "Symposium/exampleTco_Claims.png")

subset <- filter(exampleEstimates, country == "USA")
plotForest(data = subset,
           labels = subset$databaseId,
           exclude = !exampleEstimates$unblind,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           fileName = "Symposium/exampleTco_USA.png")

# Picking example for symposium:
outcomesOfInterest <- c("Acute myocardial infarction",
                        "Heart failure",
                        "Hospitalization with heart failure",
                        "Ischemic stroke",
                        "Stroke",
                        "Venous thromboembolic events",
                        "Venous thromboembolic events ",
                        "Ingrowing nail",
                        "Impacted cerumen",
                        "Contusion of knee")
examples <- highPowerTcos |>
  filter(outcomeName %in% outcomesOfInterest) |>
  group_by(outcomeName) |>
  slice_min(order_by = maxSeLogRr, n = 5)  |>
  ungroup()

for (i in seq_len(nrow(examples))) {
  example <- examples[i, ]
  exampleEstimates <- estimates |>
    filter(analysisName == !!analysisName) |>
    inner_join(example, by = join_by(targetId, targetName, comparatorId, comparatorName, outcomeId, outcomeName)) |>
    arrange(databaseId)
  plotForest(data = exampleEstimates,
             labels = exampleEstimates$databaseId,
             exclude = !exampleEstimates$unblind,
             showFixedEffects = FALSE,
             showRandomEffects = FALSE,
             title = sprintf("%s vs %s\n%s", example$targetName, example$comparatorName, example$outcomeName),
             fileName = sprintf("Symposium/Patrick/Forest_t%d_c_%d_o%d_%s.png", example$targetId, example$comparatorId, example$outcomeId, legendLabel))

}

estimates |> group_by(databaseId) |> summarise(mean(unblind))
subset <- estimates |> filter(databaseId == "CUMC")

