library(dplyr)

.qmixnorm <- function(p, means, sds) {
  pmix <- function(x) {
    mean(pnorm(x, mean = means, sd = sds))
  }
  sapply(p, function(pVal) {
    if (pVal <= 0) return(-Inf)
    if (pVal >= 1) return(Inf)
    objective_function <- function(x) {
      pmix(x) - pVal
    }
    search_interval <- c(min(means) - 10 * max(sds), max(means) + 10 * max(sds))
    uniroot(objective_function, interval = search_interval)$root
  })
}

computePredictionInterval <- function(estimate) {
  traces <- attr(estimate, "traces")
  predictionInterval <- HDInterval::hdi(.qmixnorm, credMass = 0.95, means = traces[, 1], sds = traces[, 2])
  return(predictionInterval)

  # To verify: use very large sample:
  # predictionInterval
  # predictions <- do.call(c, lapply(seq_len(nrow(traces)), function(i) rnorm(100000, traces[i, 1], traces[i, 2])))
  # predictionInterval <- HDInterval::hdi(predictions, credMass = 0.95)
  # predictionInterval
  # -0.06760276  2.06989313
}
