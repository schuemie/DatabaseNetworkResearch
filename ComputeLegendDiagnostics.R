# Compute diagnostics for all TCOs and write to table

library(dplyr)
library(DatabaseConnector)

# LEGEND T2DM --------------------------------------------------------------------------------------

connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)
schema <- "legendt2dm_class_results"

estimates <- readRDS("estimates_T2dm.rds")
exposuresOfInterest <- unique(c(estimates$targetId, estimates$comparatorId))

connection <- connect(connectionDetails)

# Balance
sql <- "
SELECT database_id,
    target_id + 1000000 AS target_id,
    comparator_id + 1000000 AS comparator_id,
    MAX(ABS(std_diff_after)) AS max_abs_std_diff_mean
FROM @schema.covariate_balance
WHERE outcome_id = 0
  AND analysis_id = 5 -- PS matching
  AND target_id IN (@exposures_of_interest)
  AND comparator_id IN (@exposures_of_interest)
GROUP BY database_id,
    target_id,
    comparator_id;
"
balance <- renderTranslateQuerySql(connection = connection,
                                   sql = sql,
                                   schema = schema,
                                   exposures_of_interest = exposuresOfInterest - 1000000,
                                   snakeCaseToCamelCase = TRUE)

# Equipoise
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    CASE WHEN
        target_equipoise < comparator_equipoise THEN target_equipoise
        ELSE comparator_equipoise
    END AS min_equipoise
FROM (
    SELECT database_id,
        target_id,
        comparator_id,
        SUM(CASE WHEN preference_score BETWEEN @min AND @max THEN target_density ELSE 0 END)/SUM(target_density) AS target_equipoise,
        SUM(CASE WHEN preference_score BETWEEN @min AND @max THEN comparator_density ELSE 0 END)/SUM(comparator_density) AS comparator_equipoise
    FROM @schema.preference_score_dist
    WHERE analysis_id = 8
      AND target_id IN (@exposures_of_interest)
      AND comparator_id IN (@exposures_of_interest)
    GROUP BY database_id,
        target_id,
        comparator_id
) tmp;
"

equipoise <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     min = 0.3,
                                     max = 0.7,
                                     exposures_of_interest = exposuresOfInterest,
                                     snakeCaseToCamelCase = TRUE)

# EASE
negativeControlIds <- renderTranslateQuerySql(connection = connection,
                                              sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
                                              schema = schema,
                                              snakeCaseToCamelCase = TRUE)[, 1]
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    outcome_id,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
WHERE outcome_id IN (@negative_control_ids)
    AND se_log_rr IS NOT NULL
    AND analysis_id = 8
    AND target_id IN (@exposures_of_interest)
    AND comparator_id IN (@exposures_of_interest)
"
ncEstimates <- renderTranslateQuerySql(connection = connection,
                                       sql = sql,
                                       schema = schema,
                                       negative_control_ids = negativeControlIds,
                                       exposures_of_interest = exposuresOfInterest,
                                       snakeCaseToCamelCase = TRUE)

groups <- ncEstimates |>
  filter(!grepl("Meta-analysis", databaseId)) |>
  group_by(databaseId, targetId, comparatorId) %>%
  group_split()

computeEase <- function(group) {
  if (nrow(group) >= 5) {
    null <- EmpiricalCalibration::fitMcmcNull(group$logRr, group$seLogRr)
    ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
    row <- group |>
      head(1) |>
      select(-outcomeId, -logRr, -seLogRr) |>
      mutate(ease = !!ease$ease)
    return(row)
  } else {
    return(NULL)
  }
}
cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
ease <- ParallelLogger::clusterApply(cluster, groups, computeEase)
ParallelLogger::stopCluster(cluster)
ease <- bind_rows(ease)

# Combine diagnostics into single table
diagnostics <- equipoise |>
  full_join(ease,
            by = join_by(databaseId, targetId, comparatorId)) |>
  full_join(balance,
            by = join_by(databaseId, targetId, comparatorId)) |>
  filter(!grepl("Meta-analysis", databaseId)) |>
  mutate(unblind = !is.na(maxAbsStdDiffMean) &
           !is.na(minEquipoise) &
           !is.na(ease) &
           maxAbsStdDiffMean < 0.15 &
           minEquipoise > 0.25 &
           ease < 0.25)
diagnostics |>
  filter(unblind) |>
  group_by(databaseId) |>
  count()

saveRDS(diagnostics, "Diagnostics_T2dm.rds")

disconnect(connection)


# LEGEND Hypertension ------------------------------------------------------------------------------

connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("ohdsiPostgresUser"),
  password = keyring::key_get("ohdsiPostgresPassword")
)
schema <- "legend"

estimates <- readRDS("estimates_Htn.rds")
exposuresOfInterest <- unique(c(estimates$targetId, estimates$comparatorId))

connection <- connect(connectionDetails)


test2 <- querySql(connection, "SELECT * FROM legend.covariate_balance WHERE analysis_id = 4 AND target_id = 2;")
unique(test2$outcome_id)
test2 |> filter(is.na(outcome_id)) |> count()
test |> filter(outcome_id == 7) |> summarise(max(abs(std_diff_after), na.rm = TRUE))

# Balance
sql <- "
SELECT database_id,
    target_id AS target_id,
    comparator_id AS comparator_id,
    MAX(ABS(std_diff_after)) AS max_abs_std_diff_mean
FROM @schema.covariate_balance
WHERE outcome_id IS NULL
  AND analysis_id = 4 -- PS matching, ITT
  AND target_id IN (@exposures_of_interest)
  AND comparator_id IN (@exposures_of_interest)
GROUP BY database_id,
    target_id,
    comparator_id;
"
balance <- renderTranslateQuerySql(connection = connection,
                                   sql = sql,
                                   schema = schema,
                                   exposures_of_interest = exposuresOfInterest,
                                   snakeCaseToCamelCase = TRUE)

# Equipoise
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    CASE WHEN
        target_equipoise < comparator_equipoise THEN target_equipoise
        ELSE comparator_equipoise
    END AS min_equipoise
FROM (
    SELECT database_id,
        target_id,
        comparator_id,
        SUM(CASE WHEN preference_score BETWEEN @min AND @max THEN target_density ELSE 0 END)/SUM(target_density) AS target_equipoise,
        SUM(CASE WHEN preference_score BETWEEN @min AND @max THEN comparator_density ELSE 0 END)/SUM(comparator_density) AS comparator_equipoise
    FROM @schema.preference_score_dist
    WHERE target_id IN (@exposures_of_interest)
      AND comparator_id IN (@exposures_of_interest)
    GROUP BY database_id,
        target_id,
        comparator_id
) tmp;
"

equipoise <- renderTranslateQuerySql(connection = connection,
                                     sql = sql,
                                     schema = schema,
                                     min = 0.3,
                                     max = 0.7,
                                     exposures_of_interest = exposuresOfInterest,
                                     snakeCaseToCamelCase = TRUE)

# EASE
negativeControlIds <- renderTranslateQuerySql(connection = connection,
                                              sql = "SELECT outcome_id FROM @schema.negative_control_outcome;",
                                              schema = schema,
                                              snakeCaseToCamelCase = TRUE)[, 1]
sql <- "
SELECT database_id,
    target_id,
    comparator_id,
    outcome_id,
    log_rr,
    se_log_rr
FROM @schema.cohort_method_result
WHERE outcome_id IN (@negative_control_ids)
    AND se_log_rr IS NOT NULL
    AND analysis_id = 3
    AND target_id IN (@exposures_of_interest)
    AND comparator_id IN (@exposures_of_interest)
"
ncEstimates <- renderTranslateQuerySql(connection = connection,
                                       sql = sql,
                                       schema = schema,
                                       negative_control_ids = negativeControlIds,
                                       exposures_of_interest = exposuresOfInterest,
                                       snakeCaseToCamelCase = TRUE)

groups <- ncEstimates |>
  filter(!grepl("Meta-analysis", databaseId)) |>
  group_by(databaseId, targetId, comparatorId) %>%
  group_split()

computeEase <- function(group) {
  if (nrow(group) >= 5) {
    null <- EmpiricalCalibration::fitMcmcNull(group$logRr, group$seLogRr)
    ease <- EmpiricalCalibration::computeExpectedAbsoluteSystematicError(null)
    row <- group |>
      head(1) |>
      select(-outcomeId, -logRr, -seLogRr) |>
      mutate(ease = !!ease$ease)
    return(row)
  } else {
    return(NULL)
  }
}
cluster <- ParallelLogger::makeCluster(10)
ParallelLogger::clusterRequire(cluster, "dplyr")
ease <- ParallelLogger::clusterApply(cluster, groups, computeEase)
ParallelLogger::stopCluster(cluster)
ease <- bind_rows(ease)

# Combine diagnostics into single table
diagnostics <- equipoise |>
  full_join(ease,
            by = join_by(databaseId, targetId, comparatorId)) |>
  full_join(balance,
            by = join_by(databaseId, targetId, comparatorId)) |>
  filter(!grepl("Meta-analysis", databaseId)) |>
  mutate(unblind = !is.na(maxAbsStdDiffMean) &
           !is.na(minEquipoise) &
           !is.na(ease) &
           maxAbsStdDiffMean < 0.15 &
           minEquipoise > 0.25 &
           ease < 0.25)
diagnostics |>
  filter(unblind) |>
  group_by(databaseId) |>
  count()

saveRDS(diagnostics, "Diagnostics_Htn.rds")

disconnect(connection)



