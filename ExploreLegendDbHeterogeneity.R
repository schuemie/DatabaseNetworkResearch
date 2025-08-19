# Fetch information from LEGEND T2DM to explore how much heterogeneity in effect estimates exists
# between databases.

library(DatabaseConnector)
library(dplyr)
library(ggplot2)
library(ggh4x)

# Fetch negative control estimates from LEGEND T2DM ------------------------------------------------
connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)

# schema <- "legendt2dm_drug_results"
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

# Compute tau -----------------------------------------------------------------------
groups <- estimates |>
  group_by(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl)|>
  group_split()


# group = groups[[1]]
computeTau <- function(group) {
  keyRow <- group |>
    head(1) |>
    select(targetId, targetName, comparatorId, comparatorName, outcomeId, analysisName, negativeControl)
  maEstimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(group, showProgressBar = FALSE)
  row <- keyRow |>
    mutate(tau = maEstimate$tau,
           tau95Lb = maEstimate$tau95Lb,
           tau95Ub = maEstimate$tau95Ub)
  return(row)
}

cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
taus <- ParallelLogger::clusterApply(cluster, groups, computeTau)
ParallelLogger::stopCluster(cluster)

taus <- bind_rows(taus)
saveRDS(taus, "taus.rds")

# Plot tau distributions ---------------------------------------------------------------------------
taus <- readRDS("taus.rds")

taus |>
  mutate(se = (tau95Ub - tau95Lb) / 2*qnorm(0.975)) |>
  group_by(negativeControl, analysisName) |>
  summarise(medianSe = median(se), .groups = "drop")

priorTimesX <- function(x){
  (dnorm(x, 0, 0.5) * 2) * x
}
expectedTauUnderPrior <- integrate(priorTimesX, lower = 0, upper = Inf)$value

vizData <- taus |>
  mutate(se = (tau95Ub - tau95Lb) / 2*qnorm(0.975)) |>
  filter(se < 0.5) |>
  mutate(negativeControl = if_else(negativeControl, "Negative control", "Outcome of interest"),
         analysisName = if_else(analysisName == "unadjusted", "Unadjusted", analysisName))

ggplot(vizData, aes(y = tau)) +
  geom_hline(yintercept = 0, size = 1) +
  geom_hline(yintercept = expectedTauUnderPrior, linetype = "dashed") +
  geom_boxplot(fill = "#3f845a", alpha = 0.75) +
  scale_y_continuous("Tau") +
  facet_nested(~ negativeControl + analysisName) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  )
ggsave("TauDistributions.png", width = 6, height = 5)
