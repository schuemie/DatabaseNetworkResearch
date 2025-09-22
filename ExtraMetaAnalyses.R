# Rerun semaglutide naion meta-analysis to produce tau and PI
library(DatabaseConnector)
library(dplyr)
library(readxl)
source("PredictionInterval.R")
source("ForestPlot.R")

connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("ohdsiPostgresServer"),
                 keyring::key_get("ohdsiPostgresShinyDatabase"), sep = "/"),
  user = keyring::key_get("semaUser"),
  password = keyring::key_get("semaPassword")
)

schema <- "semanaion"

# executeSql(connection, "COMMIT;")
connection <- connect(connectionDetails)

# SCCS per-database estimates
sql <- "SELECT dm.cdm_source_abbreviation AS database_id,
  era_id AS exposure_id,
  e.cohort_name AS exposure_name,
  outcome_id,
  o.cohort_name AS outcome_name,
  a.analysis_id,
  a.description AS analysis_description,
  unblind_for_evidence_synthesis AS unblind,
  calibrated_rr,
  calibrated_ci_95_lb,
  calibrated_ci_95_ub
FROM @schema.sccs_result r
INNER JOIN @schema.sccs_diagnostics_summary ds
  ON ds.database_id = r.database_id
    AND ds.exposures_outcome_set_id = r.exposures_outcome_set_id
    AND ds.covariate_id = r.covariate_id
    AND ds.analysis_id = r.analysis_id
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
INNER JOIN @schema.cg_cohort_definition e
  ON era_id = e.cohort_definition_id
LEFT JOIN @schema.cg_cohort_definition o
  ON outcome_id = o.cohort_definition_id
INNER JOIN @schema.sccs_analysis a
  ON r.analysis_id = a.analysis_id
INNER JOIN @schema.database_meta_data dm
  ON r.database_id = dm.database_id;"
sccsResults <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  schema = schema,
  snakeCaseToCamelCase = TRUE
)

# SCCS likelihood profiles
sql <- "SELECT dm.cdm_source_abbreviation AS database_id,
  era_id AS exposure_id,
  e.cohort_name AS exposure_name,
  outcome_id,
  o.cohort_name AS outcome_name,
  a.analysis_id,
  a.description AS analysis_description,
  log_rr,
  log_likelihood
FROM @schema.sccs_likelihood_profile lp
INNER JOIN @schema.sccs_exposures_outcome_set eos
  ON lp.exposures_outcome_set_id = eos.exposures_outcome_set_id
INNER JOIN (
  SELECT DISTINCT exposures_outcome_set_id,
    covariate_id,
    analysis_id,
    era_id
  FROM @schema.sccs_covariate
  ) c
  ON lp.exposures_outcome_set_id = c.exposures_outcome_set_id
    AND lp.covariate_id = c.covariate_id
    AND lp.analysis_id = c.analysis_id
INNER JOIN @schema.cg_cohort_definition e
  ON era_id = e.cohort_definition_id
LEFT JOIN @schema.cg_cohort_definition o
  ON outcome_id = o.cohort_definition_id
INNER JOIN @schema.sccs_analysis a
  ON lp.analysis_id = a.analysis_id
INNER JOIN @schema.database_meta_data dm
  ON lp.database_id = dm.database_id;"
sccsProfiles <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  schema = schema,
  snakeCaseToCamelCase = TRUE
)

# SCCS meta-analysis estimates
# sql <- "SELECT r.*,
#   era_id,
#   outcome_id,
#   nesting_cohort_id,
#   unblind,
#   n_total_databases,
#   a.description AS analysis_name
# FROM @schema.es_sccs_result r
# INNER JOIN @schema.es_sccs_diagnostics_summary ds
#   ON ds.exposures_outcome_set_id = r.exposures_outcome_set_id
#     AND ds.covariate_id = r.covariate_id
#     AND ds.analysis_id = r.analysis_id
# INNER JOIN @schema.sccs_exposures_outcome_set eos
#   ON ds.exposures_outcome_set_id = eos.exposures_outcome_set_id
# INNER JOIN (
#   SELECT DISTINCT exposures_outcome_set_id,
#     covariate_id,
#     analysis_id,
#     era_id
#   FROM @schema.sccs_covariate
#   ) c
#   ON ds.exposures_outcome_set_id = c.exposures_outcome_set_id
#     AND ds.covariate_id = c.covariate_id
#     AND ds.analysis_id = c.analysis_id
# INNER JOIN (
#   SELECT exposures_outcome_set_id,
#     covariate_id,
#     analysis_id,
#     COUNT(*) AS n_total_databases
#   FROM @schema.sccs_diagnostics_summary
#   GROUP BY exposures_outcome_set_id,
#     covariate_id,
#     analysis_id
#   ) ds2
#   ON ds2.exposures_outcome_set_id = c.exposures_outcome_set_id
#     AND ds2.covariate_id = c.covariate_id
#     AND ds2.analysis_id = c.analysis_id
# INNER JOIN @schema.sccs_analysis a
#   ON r.analysis_id = a.analysis_id;"
# sccsMaResults <- renderTranslateQuerySql(
#   connection = connection,
#   sql = sql,
#   schema = schema,
#   snakeCaseToCamelCase = TRUE
# )

# Cohort method per database estimates
sql <- "
SELECT cdm_source_abbreviation AS database_id,
  cm_result.target_id,
  target.cohort_name AS target_name,
  cm_result.comparator_id,
  comparator.cohort_name AS comparator_name,
  cm_result.outcome_id,
  o.cohort_name AS outcome_name,
  a.analysis_id,
  a.description AS analysis_description,
  unblind_for_evidence_synthesis AS unblind,
  calibrated_rr,
  calibrated_ci_95_lb,
  calibrated_ci_95_ub
FROM @schema.cm_result
INNER JOIN @schema.cg_cohort_definition target
  ON cm_result.target_id = target.cohort_definition_id
INNER JOIN @schema.cg_cohort_definition comparator
  ON cm_result.comparator_id = comparator.cohort_definition_id
INNER JOIN @schema.database_meta_data
  ON cm_result.database_id = database_meta_data.database_id
INNER JOIN @schema.cm_analysis a
  ON cm_result.analysis_id = a.analysis_id
INNER JOIN @schema.cm_diagnostics_summary ds
  ON cm_result.database_id = ds.database_id
    AND cm_result.analysis_id = ds.analysis_id
    AND cm_result.target_id = ds.target_id
    AND cm_result.comparator_id = ds.comparator_id
    AND cm_result.outcome_id = ds.outcome_id
LEFT JOIN @schema.cg_cohort_definition o
  ON cm_result.outcome_id = o.cohort_definition_id;
"

cmResults <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  schema = schema,
  snakeCaseToCamelCase = TRUE
)

# Cohort method likelihood profiles
sql <- "
SELECT cdm_source_abbreviation AS database_id,
  cm_likelihood_profile.target_id,
  target.cohort_name AS target_name,
  cm_likelihood_profile.comparator_id,
  comparator.cohort_name AS comparator_name,
  cm_likelihood_profile.outcome_id,
  o.cohort_name AS outcome_name,
  a.analysis_id,
  a.description AS analysis_description,
  log_rr,
  log_likelihood
FROM @schema.cm_likelihood_profile
INNER JOIN @schema.cg_cohort_definition target
  ON cm_likelihood_profile.target_id = target.cohort_definition_id
INNER JOIN @schema.cg_cohort_definition comparator
  ON cm_likelihood_profile.comparator_id = comparator.cohort_definition_id
INNER JOIN @schema.database_meta_data
  ON cm_likelihood_profile.database_id = database_meta_data.database_id
INNER JOIN @schema.cm_analysis a
  ON cm_likelihood_profile.analysis_id = a.analysis_id
LEFT JOIN @schema.cg_cohort_definition o
  ON cm_likelihood_profile.outcome_id = o.cohort_definition_id;
"
cmProfiles <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  schema = schema,
  snakeCaseToCamelCase = TRUE
)

# Redo meta-analyses -------------------------------------------------------------------------------

# group = groups[[1]]
doSingleMetaAnalysis <- function(group) {
  outcomeId <- group$outcomeId[1]
  writeLines(sprintf("Performing meta-analysis for outcome %s", outcomeId))

  if (nrow(group) == 1) {
    stop(outcomeId)
  }
  profiles <- group |>
    rename(point = logRr, value = logLikelihood) |>
    group_by(databaseId) |>
    group_split()
  maEstimate <- suppressMessages(EvidenceSynthesis::computeBayesianMetaAnalysis(profiles, showProgressBar = FALSE))
  predictionInterval <- computePredictionInterval(maEstimate)
  maEstimate <- maEstimate |>
    mutate(
      piLb = predictionInterval[1],
      piUb = predictionInterval[2],
      outcomeId = !!outcomeId)
  return(maEstimate)
}

calibrateEstimate <- function(maEstimates) {
  ncMaEstimates <- maEstimates |>
    filter(outcomeId < 10000)

  null <- EmpiricalCalibration::fitMcmcNull(ncMaEstimates$logRr, ncMaEstimates$seLogRr)
  # EmpiricalCalibration::plotCalibrationEffect(ncMaEstimates$logRr, ncMaEstimates$seLogRr, null = null)
  hoiMaEstimate <- maEstimates |>
    filter(outcomeId > 10000)
  calibratedHoiMaEstimate <- EmpiricalCalibration::calibrateConfidenceInterval(
    hoiMaEstimate$logRr,
    hoiMaEstimate$seLogRr,
    model = EmpiricalCalibration::convertNullToErrorModel(null)
  )
  return(calibratedHoiMaEstimate)
}

dbRenames <- tibble(
  databaseId = c("Epic Legacy CUMC MERGE", "JHME", "LRx-US9-LAAD 202405", "Merative CCAE", "Merative MDCD", "Merative MDCR", "OHSU", "Optum EHR", "OPTUM Extended SES", "PharMetrics", "STARR", "USC", "VA-OMOP", "WashU"),
  databaseName = c("Epic Legacy CUMC MERGE", "JHME", "IQVIA", "CCAE", "MDCD", "MDCR", "OHSU", "Optum EHR", "Clinformatics", "PharMetrics", "STARR", "USC", "VA", "WashU")
)

outcomes <- tibble(
  outcomeName = c("Nonarteric anterior ischemic neuropathy with index date correction and 2dxGCA",
                  "Nonarteric anterior ischemic neuropathy with index date correction and 2nd dx and 2dxGCA"),
  outcomeId = c("Sensitive", "Specific")
)


# SCCS
exposureName <- "semaglutide exposures"

sccsResultsUnblinded <- sccsResults |>
  filter(exposureName == !!exposureName, unblind == 1) |>
  inner_join(dbRenames, by = join_by(databaseId))
sccsProfilesUnblinded <- sccsProfiles |>
  filter(exposureName == !!exposureName) |>
  inner_join(sccsResultsUnblinded |> select(databaseId, exposureId, outcomeId, analysisId, databaseName), by = join_by(outcomeId, databaseId, exposureId, analysisId))
sccsOutcomes <- outcomes |>
  mutate(rr = c(1.32, 1.50),
         lb = c(1.14, 1.26),
         ub = c(1.54, 1.79))


for (i in seq_len(nrow(sccsOutcomes))) {
  outcomeName <- sccsOutcomes$outcomeName[i]
  dbEstimates <- sccsResultsUnblinded |>
    filter(outcomeName == !!outcomeName, !is.na(calibratedRr)) |>
    select(label = databaseName, rr = calibratedRr, lb = calibratedCi95Lb, ub = calibratedCi95Ub) |>
    arrange(label)
  writeLines(sprintf("%s: %0.2f (%0.2f-%0.2f)", dbEstimates$label, dbEstimates$rr, dbEstimates$lb, dbEstimates$ub))
  groups <- sccsProfilesUnblinded |>
    filter(outcomeId < 10000 | outcomeName == !!outcomeName) |>
    group_by(outcomeId) |>
    group_split()
  maEstimates <- lapply(groups, doSingleMetaAnalysis)
  maEstimates <- bind_rows(maEstimates)
  calibratedHoiMaEstimate <- calibrateEstimate(maEstimates)
  writeLines(sprintf("%0.2f (%0.2f-%0.2f)", exp(calibratedHoiMaEstimate$logRr), exp(calibratedHoiMaEstimate$logLb95Rr), exp(calibratedHoiMaEstimate$logUb95Rr)))

  hoiMaEstimate <- maEstimates |>
    filter(outcomeId > 10000)
  plotForestGivenEstimates(dbEstimates,
                           labels = dbEstimates$label,
                           maEstimate = sccsOutcomes[i, ],
                           maLabel = "Bayesian random effects",
                           tau = hoiMaEstimate$tau,
                           tauLb = hoiMaEstimate$tau95Lb,
                           tauUb = hoiMaEstimate$tau95Ub,
                           predictionInterval = c(hoiMaEstimate$piLb, hoiMaEstimate$piUb),
                           xLabel = "Hazard Ratio",
                           limits = c(0.1, 10),
                           alpha = 0.05,
                           title = NULL,
                           fileName = sprintf("Symposium/SemaNaionSccs%s.png", sccsOutcomes$outcomeId[i]))

}


# Cohort method
targetName <- "New user of semaglutide as 2nd-line treatment with prior T2DM and prior metformin"

cmResultsUnblinded <- cmResults |>
  filter(targetName == !!targetName, unblind == 1) |>
  inner_join(dbRenames, by = join_by(databaseId))
cmProfilesUnblinded <- cmProfiles |>
  filter(targetName == !!targetName) |>
  inner_join(cmResultsUnblinded |> select(databaseId, targetId, comparatorId, outcomeId, analysisId, databaseName), by = join_by(outcomeId, databaseId, targetId, comparatorId, analysisId))

cmComparatorOutcomes <- outcomes |>
  cross_join(tibble(comparatorName = c("New user of dulaglutide as 2nd-line treatment with prior T2DM and prior metformin",
                                       "New user of empagliflozin as 2nd-line treatment with prior T2DM and prior metformin",
                                       "New user of sitagliptin as 2nd-line treatment with prior T2DM and prior metformin",
                                       "New user of glipizide as 2nd-line treatment with prior T2DM and prior metformin"))) |>
  mutate(rr = c(0.93, 1.44, 1.30, 1.23,     1.28, 2.27, 1.64, 1.50),
         lb = c(0.46, 0.78, 0.56, 0.66,     0.57, 1.16, 0.69, 0.70),
         ub = c(1.91, 2.68, 3.01, 2.28,     2.88, 4.46, 3.90, 3.21))

unique(cmResultsUnblinded$comparatorName)

for (i in seq_len(nrow(cmComparatorOutcomes))) {
  outcomeName <- cmComparatorOutcomes$outcomeName[i]
  comparatorName <- cmComparatorOutcomes$comparatorName[i]
  dbEstimates <- cmResultsUnblinded |>
    filter(comparatorName == !!comparatorName, outcomeName == !!outcomeName, !is.na(calibratedRr)) |>
    select(label = databaseName, rr = calibratedRr, lb = calibratedCi95Lb, ub = calibratedCi95Ub) |>
    arrange(label)
  writeLines(sprintf("%s: %0.2f (%0.2f-%0.2f)", dbEstimates$label, dbEstimates$rr, dbEstimates$lb, dbEstimates$ub))
  groups <- cmProfilesUnblinded |>
    filter(comparatorName == !!comparatorName, outcomeId < 10000 | outcomeName == !!outcomeName) |>
    group_by(outcomeId) |>
    group_split()
  maEstimates <- lapply(groups, doSingleMetaAnalysis)
  maEstimates <- bind_rows(maEstimates)
  calibratedHoiMaEstimate <- calibrateEstimate(maEstimates)
  writeLines(sprintf("%0.2f (%0.2f-%0.2f)", exp(calibratedHoiMaEstimate$logRr), exp(calibratedHoiMaEstimate$logLb95Rr), exp(calibratedHoiMaEstimate$logUb95Rr)))
  hoiMaEstimate <- maEstimates |>
    filter(outcomeId > 10000)
  comparatorName <- gsub("New user of | as 2nd-line treatment with prior T2DM and prior metformin", "", cmComparatorOutcomes$comparatorName[i])
  plotForestGivenEstimates(dbEstimates,
                           labels = dbEstimates$label,
                           maEstimate = cmComparatorOutcomes[i, ],
                           maLabel = "Bayesian random effects",
                           tau = hoiMaEstimate$tau,
                           tauLb = hoiMaEstimate$tau95Lb,
                           tauUb = hoiMaEstimate$tau95Ub,
                           predictionInterval = c(hoiMaEstimate$piLb, hoiMaEstimate$piUb),
                           xLabel = "Hazard Ratio",
                           limits = c(0.1, 10),
                           alpha = 0.05,
                           title = NULL,
                           fileName = sprintf("Symposium/SemaNaionCm%s_%s.png", cmComparatorOutcomes$outcomeId[i], comparatorName))
}


# All Sema-NAOIN studies ---------------------------------------------------------------------------
studies <- read_excel("semaglutide NAION papers.xlsx")

studies <- studies |>
  mutate(label = if_else(Which == "SCCS",
                         sprintf("%s et al.; %s", gsub(" .*", "", Authors), gsub(";.*", "", Database)),
                         sprintf("%s et al.; %s; %s", gsub(" .*", "", Authors), gsub(";.*", "", Database), gsub(" .*", "", Comparator))),
         #  label = sprintf("%s et al. %s", gsub(" .*", "", Authors), substr(Date, 1, 4)),
         logRr = log(RR),
         seLogRr = (log(UB) - log(LB)) / (2 * qnorm(0.975)))


sccs <- studies |>
  filter(Which == "SCCS")
plot <- plotForest(data = sccs,
           labels = sccs$label)
ggsave("Symposium/LiteratureSccs.png", plot, width = 12, height = 1 + 16 * 0.3, dpi = 300)
cohort <- studies |>
  filter(Which == "Cohort")
plot <- plotForest(data = cohort,
           labels = cohort$label)
ggsave("Symposium/LiteratureCohort.png", plot, width = 14, height = 1 + 23 * 0.3, dpi = 300)
