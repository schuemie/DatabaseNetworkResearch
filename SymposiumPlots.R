# Plots for the global symposium plenary
source("ForestPlot.R")

# Forest plots showing examples of 2 and 4 databases, both with (almost) same meta-analytic estimate,
# but 4 DBs has narrower prediction interval
data <- tibble(
  logRr =   c(1.011, 1.075),
  seLogRr = c(0.175, 0.402)
)
plotForest(data = data,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           showPredictionInterval = FALSE,
           fileName = "Symposium/ExampleSameCi2Dbs.png")
plotForest(data = data,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           showPredictionInterval = TRUE,
           fileName = "Symposium/ExampleSameCi2Dbs_Pi.png")

data <- tibble(
  logRr =   c(0.877, 1.07, 1.05, 1.24),
  seLogRr = c(0.526, 0.6, 0.608, 0.8)
)
plotForest(data = data,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           showPredictionInterval = FALSE,
           fileName = "Symposium/ExampleSameCi4Dbs.png")
plotForest(data = data,
           showFixedEffects = FALSE,
           showRandomEffects = FALSE,
           showPredictionInterval = TRUE,
           fileName = "Symposium/ExampleSameCi4Dbs_Pi.png")
