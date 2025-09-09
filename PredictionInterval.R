library(dplyr)

computePredictionInterval <- function(estimate) {
  traces <- attr(estimate, "traces")
  # Truncate extremely small tau estimates to improve stability:
  traces[traces[, 2] < 0.01, 2] <- 0.01
  gridMin <- min(qnorm(0.025, traces[, 1], traces[, 2]))
  gridMax <- max(qnorm(0.975, traces[, 1], traces[, 2]))
  grid <- seq(gridMin, gridMax, length.out = 4000)
  predictiveDensity <- sapply(grid, function(x) mean(dnorm(x, mean = traces[, 1], sd = traces[, 2])))
  # plot(grid, predictiveDensity)
  findHdiFromGrid <- function(grid, density, credMass = 0.95) {
    probMass <- density / sum(density)
    sortedIndices <- order(probMass, decreasing = TRUE)
    sortedProbMass <- probMass[sortedIndices]
    cumulativeProb <- cumsum(sortedProbMass)
    hdiPointCount <- which(cumulativeProb >= credMass)[1]
    hdiIndices <- sortedIndices[1:hdiPointCount]
    hdiInterval <- range(grid[hdiIndices])
    return(hdiInterval)
  }
  predictionInterval <- findHdiFromGrid(grid, predictiveDensity, credMass = 0.95)
  predictionEstimate <- weighted.mean(grid, predictiveDensity)
  return(c(predictionInterval[1], predictionEstimate, predictionInterval[2]))

  # To verify: use very large sample:
  # predictionInterval
  # predictions <- do.call(c, lapply(seq_len(nrow(traces)), function(i) rnorm(100000, traces[i, 1], traces[i, 2])))
  # predictionInterval <- HDInterval::hdi(predictions, credMass = 0.95)
  # predictionInterval
}
