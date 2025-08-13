library(DatabaseConnector)
library(dplyr)
library(Matrix)
library(umap)
library(ggplot2)

connectionDetails <- createConnectionDetails(
  dbms = "postgresql",
  server = paste(keyring::key_get("assureServer"), keyring::key_get("assureDatabase"), sep = "/"),
  user = keyring::key_get("assureUser"),
  password = keyring::key_get("assurePassword")
)

profilesDatabaseSchema <- "dp_temp"
profilesTable <- "db_profile_results"
vocabularyDatabaseSchema <- "vocabulary_20240229"

# Fetch data -------------------------------------------------------------------
conditionRecordCountAnalysisId <- 401
ingredientRecordCountAnalysisId <- 901
numberOfPersonsAnalysisId <- 1
lengthOfOpAnalysisId <- 108
numberOfOpAnalysisId <- 113

connection <- connect(connectionDetails)
sql <- "
SELECT cdm_source_name,
  analysis_id,
  CAST(stratum_1 AS INT) AS stratum_id,
  count_value
FROM @profiles_database_schema.@profiles_table
WHERE analysis_id IN (@analysis_ids);
"
normalizationData <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  profiles_database_schema = profilesDatabaseSchema,
  profiles_table = profilesTable,
  analysis_ids = c(numberOfPersonsAnalysisId,
                   lengthOfOpAnalysisId,
                   numberOfOpAnalysisId),
  snakeCaseToCamelCase = TRUE
)
sql <- "
SELECT cdm_source_name,
  ancestor_concept_id AS concept_id,
  concept_name,
  SUM(count_value) AS record_count
FROM @profiles_database_schema.@profiles_table
INNER JOIN @vocabulary_database_schema.concept_ancestor
  ON CAST(stratum_1 AS INT) = descendant_concept_id
INNER JOIN (
	SELECT concept_id,
	  concept_name
	FROM @vocabulary_database_schema.concept
	INNER JOIN (
	  SELECT *
	  FROM @vocabulary_database_schema.concept_ancestor
	  WHERE ancestor_concept_id = 441840 /* SNOMED clinical finding */
	  AND (min_levels_of_separation > 2
		OR descendant_concept_id IN (433736, 433595, 441408, 72404, 192671, 137977, 434621, 437312, 439847, 4171917, 438555, 4299449, 375258, 76784, 40483532, 4145627, 434157, 433778, 258449, 313878)
		)
	) temp
	  ON concept_id = descendant_concept_id
	WHERE concept_name NOT LIKE '%finding'
		AND concept_name NOT LIKE 'Disorder of%'
		AND concept_name NOT LIKE 'Finding of%'
		AND concept_name NOT LIKE 'Disease of%'
		AND concept_name NOT LIKE 'Injury of%'
		AND concept_name NOT LIKE '%by site'
		AND concept_name NOT LIKE '%by body site'
		AND concept_name NOT LIKE '%by mechanism'
		AND concept_name NOT LIKE '%of body region'
		AND concept_name NOT LIKE '%of anatomical site'
		AND concept_name NOT LIKE '%of specific body structure%'
		AND domain_id = 'Condition'
) valid_groups
	ON ancestor_concept_id = valid_groups.concept_id
WHERE analysis_id = @analysis_id
GROUP BY cdm_source_name,
  ancestor_concept_id,
  concept_name
"
conditionConceptCounts <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  profiles_database_schema = profilesDatabaseSchema,
  profiles_table = profilesTable,
  vocabulary_database_schema = vocabularyDatabaseSchema,
  analysis_id = conditionRecordCountAnalysisId,
  snakeCaseToCamelCase = TRUE
)

sql <- "
SELECT cdm_source_name,
  concept_id,
  concept_name,
  SUM(count_value) AS record_count
FROM @profiles_database_schema.@profiles_table
INNER JOIN @vocabulary_database_schema.concept
  ON CAST(stratum_1 AS INT) = concept_id
WHERE analysis_id = @analysis_id
GROUP BY cdm_source_name,
  concept_id,
  concept_name
"
drugConceptCounts <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  profiles_database_schema = profilesDatabaseSchema,
  profiles_table = profilesTable,
  vocabulary_database_schema = vocabularyDatabaseSchema,
  analysis_id = ingredientRecordCountAnalysisId,
  snakeCaseToCamelCase = TRUE
)

disconnect(connection)

profileData <- list(
  normalizationData = normalizationData,
  conditionConceptCounts = conditionConceptCounts,
  drugConceptCounts = drugConceptCounts
)
saveRDS(profileData, "e:/temp/profileData.rds")

# Normalize data ---------------------------------------------------------------
profileData <- readRDS("e:/temp/profileData.rds")
nDatabases <- profileData$normalizationData |>
  summarise(n_distinct(cdmSourceName)) |>
  pull()

# Restrict to concept IDs found in all databases:
commonConditionConceptIds <- profileData$conditionConceptCounts |>
  group_by(conceptId) |>
  summarise(n = n_distinct(cdmSourceName)) |>
  filter(n == nDatabases) |>
  pull(conceptId)
conditionConceptCounts <- profileData$conditionConceptCounts |>
  filter(conceptId %in% commonConditionConceptIds)

# Compute total days observed (approximate)
meanDaysPerObservationPeriod <- profileData$normalizationData |>
  filter(analysisId == lengthOfOpAnalysisId) |>
  mutate(daysObserved = (15 + stratumId * 30) * countValue) |>
  group_by(cdmSourceName) |>
  summarise(daysObserved = sum(daysObserved) / sum(countValue))
nObservationPeriod <- profileData$normalizationData |>
  filter(analysisId == numberOfOpAnalysisId)  |>
  mutate(nObservationPeriod = stratumId * countValue) |>
  group_by(cdmSourceName) |>
  summarise(nObservationPeriod = sum(nObservationPeriod))
daysObserved <- inner_join(
  meanDaysPerObservationPeriod,
  nObservationPeriod,
  by = join_by(cdmSourceName)
) |>
  mutate(daysObserved = daysObserved * nObservationPeriod) |>
  select("cdmSourceName", "daysObserved")

# Combine and normalize
vectors <- bind_rows(
  profileData$drugConceptCounts,
  conditionConceptCounts
) |>
  inner_join(daysObserved, by = join_by(cdmSourceName)) |>
  mutate(value = recordCount / daysObserved)

# Reindex
covariateIdToConceptId <- vectors |>
  distinct(conceptId) |>
  arrange(conceptId) |>
  mutate(covariateId = row_number())
databaseIdToCdmSourceName <- vectors |>
  distinct(cdmSourceName) |>
  arrange(cdmSourceName) |>
  mutate(databaseId = row_number())

vectors <- vectors |>
  inner_join(covariateIdToConceptId, by = join_by(conceptId)) |>
  inner_join(databaseIdToCdmSourceName, by = join_by(cdmSourceName)) |>
  select("databaseId", databaseName = "cdmSourceName", "covariateId", covariateName = "conceptName", covariateValue = "value")

# Compute distance -------------------------------------------------------------
sMatrix <- sparseMatrix(
  i = vectors$databaseId,
  j = vectors$covariateId,
  x = vectors$covariateValue,
  dims = c(max(vectors$databaseId), max(vectors$covariateId)),
  dimnames = list(
    databaseIdToCdmSourceName$cdmSourceName,#paste0("db", 1:max(vectors$databaseId)), # Row names
    covariateIdToConceptId$covariateId#paste0("cov", 1:max(vectors$covariateId)) # Column names
  )
)
rowNorms <- sqrt(rowSums(sMatrix^2))
sMatrix <- sMatrix / rowNorms
cosineSimilarity <- crossprod(t(sMatrix))

readr::write_csv(as.data.frame(as.matrix(cosineSimilarity)), "e:/temp/cs.csv")

# Plot -------------------------------------------------------------------------
umapSettings <- umap::umap.defaults
umapSettings$metric <- "cosine"
umapSettings$n_neighbors <- 2
umapSettings$min_dist <- 0.001
map <- umap::umap(cosineSimilarity, input = "dist", config = umapSettings)

vizData <- tibble(
  databaseName = rownames(map$layout),
  x = map$layout[, 1],
  y = map$layout[, 2]
)

ggplot(vizData, aes(x = x, y = y)) +
  # geom_point(alpha = 0.6, size = 2, color = "#4285F4") +
  geom_label(aes(label = databaseName), hjust = 0.5, vjust = 0.5) +
  theme_minimal()



