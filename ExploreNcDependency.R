# Fetch data from LEGEND T2DM negative controls (NCs) to see how consistent systematic error is
# across databases is, both at an aggregate level (across all NCs), and per NC.

library(DatabaseConnector)
library(dplyr)
library(ggplot2)

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
WHERE outcome_id IN (@negative_control_ids)
    AND se_log_rr IS NOT NULL
    AND analysis_id IN (7, 8)
    AND target.exposure_name LIKE '%main ot2'
    AND comparator.exposure_name LIKE '%main ot2';
"
ncEstimates <- renderTranslateQuerySql(connection = connection,
                                       sql = sql,
                                       schema = schema,
                                       negative_control_ids = negativeControlIds,
                                       snakeCaseToCamelCase = TRUE)
ncEstimates <- ncEstimates |>
  filter(!grepl("Meta-analysis", databaseId)) |>
  mutate(databaseId = as.factor(databaseId),
         analysisName = as.factor(analysisName),
         targetName = as.factor(gsub(" main ot2", "", targetName)),
         comparatorName = as.factor(gsub(" main ot2", "", comparatorName)))
saveRDS(ncEstimates, "ncEstimates.rds")
disconnect(connection)

# Fetch negative control estimates from Semanaion --------------------------------------------------
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
WHERE outcome_id IN (@negative_control_ids)
    AND se_log_rr IS NOT NULL;
"
ncEstimates <- renderTranslateQuerySql(connection = connection,
                                       sql = sql,
                                       schema = schema,
                                       negative_control_ids = negativeControlIds,
                                       snakeCaseToCamelCase = TRUE)
ncEstimates <- ncEstimates |>
  mutate(databaseId = as.factor(databaseId),
         targetName = as.factor(gsub("and ", "with ", gsub("New user of |with prior T2DM |treatment ", "", targetName))),
         comparatorName = as.factor(gsub("and ", "with ", gsub("New user of |with prior T2DM |treatment ", "", comparatorName))),
         analysisId = 1,
         analysisName = "PS matching")

saveRDS(ncEstimates, "ncEstimatesSemanaion.rds")
disconnect(connection)



# Compute bias distributions -----------------------------------------------------------------------
ncEstimates <- readRDS("ncEstimates.rds")

groups <- ncEstimates |>
  group_by(databaseId,
           targetId,
           targetName,
           comparatorId,
           comparatorName,
           analysisId,
           analysisName) |>
  group_split()

plotData <- list()
tableData <- list()
limits <- c(0.1, 10)

for (i in 1:length(groups)) {
  estimates <- groups[[i]]
  keyRow <- estimates |>
    head(1) |>
    mutate(label = sprintf("%s - %s", targetName, comparatorName)) |>
    select(targetId, comparatorId, analysisId, databaseId, analysisName, label)

  null <- EmpiricalCalibration::fitMcmcNull(logRr = estimates$logRr,
                                            seLogRr = estimates$seLogRr)

  mcmc <- attr(null, "mcmc")
  mean <- null[1]
  sd <- 1 / sqrt(null[2])
  lb95Mean <- quantile(mcmc$chain[, 1], 0.025)
  ub95Mean <- quantile(mcmc$chain[, 1], 0.975)
  ub95Sd <- 1 / sqrt(quantile(mcmc$chain[, 2], 0.025))
  lb95Sd <- 1 / sqrt(quantile(mcmc$chain[, 2], 0.975))
  tableData[[i]] <- keyRow |>
    mutate(mean = mean,
           sd = sd,
           xMin = mean - sd,
           xMax = mean + sd)
  x <- seq(log(limits[1]), log(limits[2]), length.out = 100)
  compute <- function(x) {
    yMcmc <- dnorm(rep(x, nrow(mcmc$chain)), mean = mcmc$chain[, 1], sd = 1/sqrt(mcmc$chain[, 2]))
    return(quantile(yMcmc, c(0.025, 0.5, 0.975)))
  }
  ys <- sapply(x, compute)
  y <- ys[2, ]
  yMaxLb <- ys[1, ]
  yMaxUb <- ys[3, ]
  normFactor <- max(ys[2, ])
  y <- y / normFactor
  yMaxLb <- yMaxLb / normFactor
  yMaxUb <- yMaxUb / normFactor
  plotData[[i]] <- keyRow |>
    bind_cols(
      tibble(x = x,
             yMax = y,
             yMaxLb = yMaxLb,
             yMaxUb =  yMaxUb,
             yMin = 0)
    )
}

tableData <- bind_rows(tableData)
plotData <- bind_rows(plotData)

# Plot systematic error distributions --------------------------------------------------------------
plotFigure <- function(tableData, plotData, title) {
  breaks <- c(0.25, 1, 4, 8)
  plot <- ggplot(tableData) +
    geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.4) +
    geom_vline(xintercept = 0, size = 0.8) +
    geom_ribbon(aes(x = .data$x, ymax = .data$yMax, ymin = .data$yMin), fill = "#FF2700", alpha = 0.6, data = plotData) +
    geom_ribbon(aes(x = .data$x, ymax = .data$yMaxUb, ymin = .data$yMax), fill = "#FF2700", alpha = 0.3, data = plotData) +
    coord_cartesian(xlim = log(limits), ylim = c(0, 2)) +
    scale_x_continuous("Systematic Error", breaks = log(breaks), labels = breaks) +
    facet_grid(databaseId ~ label, switch = "y") +
    ggtitle(title) +
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
  return(plot)
}

plotFigure(tableData |> filter(analysisId == 7), plotData |> filter(analysisId == 7), "Unadjusted")
ggsave("SysErrorDistUnadj.png", width = 14, height = 5, dpi = 300)

plotFigure(tableData |> filter(analysisId == 8), plotData |> filter(analysisId == 8), "PS matched")
ggsave("SysErrorDistMatched.png", width = 14, height = 5, dpi = 300)

# Plot individual estimates ------------------------------------------------------------------------
databases <- ncEstimates |>
  distinct(databaseId) |>
  arrange(desc(databaseId)) |>
  mutate(y = row_number())


plotFigure <- function(targetName, comparatorName) {
  # Analysis ID 7 is unadjusted (no PS adjustment)
  # Keeping only estimates with SE < 0.25 to get rid of noise due to random error
  vizData <- ncEstimates |>
    filter(analysisId == 7, targetName == !!targetName, comparatorName == !!comparatorName) |>
    filter(seLogRr < 0.25) |>
    inner_join(databases, by = join_by(databaseId)) |>
    mutate(y = y + seLogRr) |>
    arrange(outcomeId, databaseId)

  breaks <- c(0.25, 1, 4, 8)
  plot <- ggplot(vizData, aes(x = logRr, y = y, group = outcomeId, color = as.factor(outcomeId))) +
    geom_hline(yintercept = databases$y) +
    geom_vline(xintercept = log(breaks), colour = "#AAAAAA", lty = 1, size = 0.4) +
    geom_vline(xintercept = 0, size = 0.8) +
    geom_point(alpha = 0.75) +
    geom_path(, alpha = 0.35) +
    scale_x_continuous("Hazard Ratio", breaks = log(breaks), labels = breaks) +
    scale_y_continuous("Standard Error", breaks = databases$y + 0.5, labels = databases$databaseId) +
    coord_cartesian(xlim = c(log(0.25), log(4)), ylim = c(min(databases$y), max(databases$y) + 0.99)) +
    ggtitle(sprintf("%s - %s", targetName, comparatorName)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          axis.text = element_text(size = 14),
          axis.ticks.y = element_blank(),
          legend.position = "none")
  return(plot)
}
plotFigure("GLP1RA", "SU")
ggsave("NcEstimatesAcrossDbsGLP1RA_SU.png", width = 10, height = 6, dpi = 300)

plotFigure("DPP4I", "GLP1RA")
ggsave("NcEstimatesAcrossDbsDPP4I_GLP1RA.png", width = 10, height = 6, dpi = 300)
