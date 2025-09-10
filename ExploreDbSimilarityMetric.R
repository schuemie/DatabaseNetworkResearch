# Compute similarity between databases based on aggregate statistcs (already collected for database
# diagnostics).
# Currently restricting to drugs (ingredients) and conditions, because
# 1. Other characteristics like demographics will likely get swamped by these anyway.
# 2. Drugs and ingredients will be correlated with everything else, and are therefore good proxies.
# 3. We can make drugs and conditions very similar in terms of standard concepts, so focusing on the
#    content differences, not the coding differences.

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
# Needed for normalization of counts
numberOfPersonsAnalysisId <- 1
numberOfOpAnalysisId <- 113

# Demographics
genderAnalysisId <- 2
ageAnalysisId  <- 101

# Observation periods
lengthOfOpAnalysisId <- 108

# Visits
visitAnalysisId <- 200

# Conditions
conditionRecordCountAnalysisId <- 401

# Drugs
ingredientRecordCountAnalysisId <- 901


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
  analysis_id,
  CAST(stratum_1 AS INT) AS stratum_id,
  count_value
FROM @profiles_database_schema.@profiles_table
WHERE analysis_id IN (@analysis_ids);
"
otherFeatures <- renderTranslateQuerySql(
  connection = connection,
  sql = sql,
  profiles_database_schema = profilesDatabaseSchema,
  profiles_table = profilesTable,
  analysis_ids = c(genderAnalysisId,
                   ageAnalysisId,
                   visitAnalysisId),
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
  otherFeatures = otherFeatures,
  conditionConceptCounts = conditionConceptCounts,
  drugConceptCounts = drugConceptCounts
)
saveRDS(profileData, "e:/temp/profileData.rds")

# Normalize data ---------------------------------------------------------------
profileData <- readRDS("e:/temp/profileData.rds")
nDatabases <- profileData$normalizationData |>
  summarise(n_distinct(cdmSourceName)) |>
  pull()
nPersons <- profileData$normalizationData |>
  filter(analysisId == numberOfPersonsAnalysisId) |>
  select(cdmSourceName, nPersons = countValue)

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
drugAndConditionvectors <- bind_rows(
  profileData$drugConceptCounts |>
    mutate(type = "Drug"),
  conditionConceptCounts |>
    mutate(type = "Condition")
) |>
  inner_join(daysObserved, by = join_by(cdmSourceName)) |>
  mutate(covariateValue = recordCount / daysObserved) |>
  select(databaseId = cdmSourceName, covariateId = conceptId, covariateValue, type) |>
  as_tibble()

otherVectors <- bind_rows(
  otherFeatures |>
    mutate(type = if_else(analysisId == visitAnalysisId, "Visit", "Demographic")),
  normalizationData |>
    filter(analysisId == lengthOfOpAnalysisId) |>
    mutate(type = "Observation period")
) |>
  inner_join(nPersons, by = join_by(cdmSourceName  )) |>
  mutate(covariateId = stratumId * 1000 + analysisId,
         covariateValue = countValue / nPersons) |>
  select(databaseId = cdmSourceName, covariateId, covariateValue, type) |>
  as_tibble()

vectors <- bind_rows(
  drugAndConditionvectors,
  otherVectors
)


# Compute distance -------------------------------------------------------------
library(dplyr)
library(tidyr)
library(text2vec)

wideDf <- vectors |>
  pivot_wider(
    id_cols = c(databaseId, type),
    names_from = covariateId,
    values_from = covariateValue,
    values_fill = list(covariateValue = 0)
  )
computeCosineSimilarityByType <- function(data) {
  # Remove databaseId and type columns for similarity computation
  mat <- as.matrix(data |> select(-databaseId, -type))
  rownames(mat) <- data$databaseId
  simMatrix <- sim2(mat, method = "cosine", norm = "l2")

  # Convert similarity matrix to tidy format
  simDf <- as.data.frame(as.table(simMatrix))
  colnames(simDf) <- c("databaseId1", "databaseId2", "cosineSimilarity")
  simDf$type <- data$type[1]
  return(simDf)
}
results <- wideDf |>
  group_split(type) |>
  lapply(computeCosineSimilarityByType) |>
  bind_rows()

results <- results |>
  group_by(databaseId1, databaseId2) |>
  summarise(similarity = mean(cosineSimilarity), .groups = "drop")
fullMatrix <- results |>
  pivot_wider(id_cols = "databaseId1", names_from = "databaseId2", values_from = "similarity", names_sort = TRUE)
readr::write_csv(fullMatrix, "e:/temp/DatabaseCharacteristicsSimilarity.csv")

# Heat map with hierarchical clustering ----------------------------------------
matrix <- as.matrix(fullMatrix |> select(-databaseId1))
rownames(matrix) <- fullMatrix$databaseId1
matrix = 1 - matrix # Turn similarity into distance
hc <- hclust(as.dist(matrix), method = "average")

library(pheatmap)
pheatmap(
  matrix,                 # similarity matrix
  clustering_method = "average", # match hclust method
  clustering_distance_rows = as.dist(matrix),
  clustering_distance_cols = as.dist(matrix),
  display_numbers = FALSE,    # show numbers if you want
  treeheight_row = 20,        # height of row dendrogram
  treeheight_col = 20,         # height of col dendrogram
  filename = "e:/temp/DatabaseCharacteristicsSimilarity.png"
)
