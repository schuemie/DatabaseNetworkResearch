library(dplyr)
library(tidyr)
library(brms)

fitBayesianSystematicErrorModel <- function(data) {
  data <- data |>
    mutate(databaseId = as.factor(databaseId),
           outcomeId = as.factor(outcomeId))

  # The formula now models logRr, accounting for measurement error via se(seLogRr).
  # - `0 + databaseId`: This estimates the average logRr for each databaseId (μ1, μ2) without a global intercept.
  # - `(0 + databaseId | outcomeId)`: This is the key part. It allows the effect for each
  #   outcome to vary from the database average. Crucially, it estimates the
  #   correlation between these outcome-specific variations across the databases.
  #   This correlation is our parameter of interest, rho (ρ).
  modelFormula <- brms::bf(logRr | se(seLogRr) ~ 0 + databaseId + (0 + databaseId | outcomeId))

  # Set weakly informative priors for the new model parameters.
  # 'b': The fixed effects (the average logRRs for each databaseId).
  # 'sd': The standard deviations of the outcome-specific effects (σ1, σ2).
  # 'cor': The correlation matrix for the outcome-specific effects (contains ρ).
  priors <- c(
    brms::prior(normal(0, 2), class = "b"),
    brms::prior(exponential(1), class = "sd"),
    brms::prior(lkj(1), class = "cor")
  )

  # Fit the Bayesian model.
  # This might take a few minutes to run.
  correlationModelFit <- brms::brm(
    formula = modelFormula,
    data = data,
    prior = priors,
    chains = 4,
    iter = 11000,
    warmup = 1000,
    cores = parallel::detectCores(),
    control = list(adapt_delta = 0.95)
  )

  summary(correlationModelFit) #(look for Rhat values equal to 1.00)
  # plot(correlationModelFit)
  bayesianFit <- correlationModelFit
  return (correlationModelFit)
}


calibrateCiBayesianRandomEffects <- function(bayesianFit,
                                             newData,
                                             nSamples = 4000,
                                             tauPriorScale = 1) {
  # Extract posterior samples from the brms fit
  posteriorSamples <- as.data.frame(bayesianFit)
  x <- tibble(name = colnames(posteriorSamples))

  # Get the database IDs from the model's fixed effect names
  # (e.g., from "b_databaseIdCCAE", extract "CCAE")
  modelDbIds <- unique(gsub("b_databaseId", "", names(posteriorSamples)[startsWith(names(posteriorSamples), "b_")]))

  # Find common databases and subset/reorder the new data
  availableDbIds <- intersect(modelDbIds, newData$databaseId)
  if (length(availableDbIds) == 0) {
    stop("None of the databases in newData are present in the systematic error model.")
  }

  newData <- newData[newData$databaseId %in% availableDbIds, ]
  newData <- newData[match(availableDbIds, newData$databaseId), ] # Ensure order matches model

  nDatabases <- length(availableDbIds)

  # Select a random subset of posterior samples to work with for efficiency
  if (nrow(posteriorSamples) > nSamples) {
    sampleIndices <- sample(1:nrow(posteriorSamples), nSamples)
    posteriorSamples <- posteriorSamples[sampleIndices, ]
  } else {
    nSamples <- nrow(posteriorSamples)
  }

  # Vector to store the final calibrated log-RR samples
  calibratedMuSamples <- numeric(nSamples)
  calibratedTauSamples <- numeric(nSamples)
  ones <- rep(1, nDatabases)

  # Define a grid for tau values to compute its posterior
  tauGrid <- seq(0, 2, length.out = 1000)

  # --- Main Loop: Iterate through each posterior sample ---
  for (i in 1:nSamples) {
    postSample <- posteriorSamples[i, ]

    # --- 1. Reconstruct Systematic Error Mean & Covariance ---
    # (This section remains the same)
    meanSysNames <- paste0("b_databaseId", availableDbIds)
    meanSys <- as.numeric(postSample[meanSysNames])
    sdNames <- paste0("sd_outcomeId__databaseId", availableDbIds)
    sds <- as.numeric(postSample[sdNames])
    corSys <- matrix(0, nDatabases, nDatabases)
    diag(corSys) <- 1
    if (nDatabases > 1) {
      for (row_idx in 2:nDatabases) {
        for (col_idx in 1:(row_idx - 1)) {
          db1 <- availableDbIds[col_idx]
          db2 <- availableDbIds[row_idx]
          paramName <- paste0("cor_outcomeId__databaseId", db1, "__databaseId", db2)
          if (!paramName %in% names(postSample)) {
            paramName <- paste0("cor_outcomeId__databaseId", db2, "__databaseId", db1)
          }
          correlation <- as.numeric(postSample[paramName])
          corSys[row_idx, col_idx] <- correlation
          corSys[col_idx, row_idx] <- correlation
        }
      }
    }
    covSys <- diag(sds) %*% corSys %*% diag(sds)

    # --- 2. De-bias the new data ---
    yPrime <- newData$logRr - meanSys

    if (nDatabases == 1) {
      totalVar <- covSys[1, 1] + newData$seLogRr^2
      calibratedMuSamples[i] <- rnorm(1, mean = yPrime, sd = sqrt(totalVar))
      sampledTau2[i] <- NA
    } else {
      # --- 3. Compute the posterior for tau and sample from it ---

      sigmaTotal <- covSys + diag(newData$seLogRr^2)
      logTauPosterior <- sapply(tauGrid, function(tau) {
        # For each tau value on the grid, calculate the marginal log-likelihood of the data
        # after integrating out the pooled effect 'mu'.
        V <- sigmaTotal + diag(tau^2, nDatabases)
        VInv <- try(solve(V), silent = TRUE)
        if (inherits(VInv, "try-error")) return(-Inf)

        # This is the log of the marginal likelihood p(y' | tau)
        logLikelihood <- mvtnorm::dmvnorm(yPrime, mean = rep(0, nDatabases), sigma = V, log = TRUE)

        # Get the log of the prior p(tau)
        # Using half-Normal prior
        logPrior <- dnorm(tau, mean = 0, sd = tauPriorScale, log = TRUE) + log(2)

        # Log posterior is proportional to log-likelihood + log-prior
        return(logLikelihood + logPrior)
      })

      # Normalize and sample one value for tau from its posterior
      logTauPosterior <- logTauPosterior - max(logTauPosterior) # For numerical stability
      tauProbs <- exp(logTauPosterior)
      tauProbs <- tauProbs / sum(tauProbs)
      sampledTau <- sample(tauGrid, size = 1, prob = tauProbs)
      sampledTau2 <- sampledTau^2
      calibratedTauSamples[i] <- sampledTau2

      # --- 4. Use the sampled tau to find the posterior for mu ---
      VFinal <- sigmaTotal + diag(sampledTau2, nDatabases)
      VFinalInv <- solve(VFinal)

      muVarPost <- 1 / (t(ones) %*% VFinalInv %*% ones)
      muMeanPost <- muVarPost * (t(ones) %*% VFinalInv %*% yPrime)

      # --- 5. Draw one sample from mu's posterior ---
      calibratedMuSamples[i] <- rnorm(1, mean = muMeanPost, sd = sqrt(muVarPost))
    }
  }

  # --- Summarize the final posterior distribution ---
  # (This section remains the same)
  posteriorMedian <- median(calibratedMuSamples)
  posteriorCi <- quantile(calibratedMuSamples, probs = c(0.025, 0.975))
  posteriorTau2Median <- median(calibratedTauSamples)
  posteriorTau2Ci <- quantile(calibratedTauSamples, probs = c(0.025, 0.975))

  return(data.frame(
    estimate = exp(posteriorMedian),
    ciLower = exp(posteriorCi[1]),
    ciUpper = exp(posteriorCi[2]),
    tau2 = posteriorTau2Median,
    tau2CiLower = posteriorTau2Ci[1],
    tau2CiUpper = posteriorTau2Ci[2],
    nDatabases = nDatabases,
    logRr = posteriorMedian,
    seLogRr = sd(calibratedMuSamples)
  ))
}



fitUnifiedModel <- function(data) {
  data <- data |>
    mutate(databaseId = as.factor(databaseId),
           outcomeId = as.factor(outcomeId),
           outcomeGroup = as.factor(if_else(negativeControl == 1, "Control", outcomeId)))
  modelFormula <- brms::bf(
    logRr | se(seLogRr) ~ 0 + outcomeGroup +
      (0 + databaseId | outcomeId) +     # SHARED systematic error (as before)
      (0 + outcomeGroup | databaseId)    # UNIQUE mu and tau per group
  )
  # x <- brms::get_prior(
  #   formula = modelFormula,
  #   data = data
  # )

  priors <- c(
    brms::prior(normal(0, 0.001), class = "b", coef = "outcomeGroupControl"),
    brms::prior(normal(0, 2), class = "b"),
    brms::prior(student_t(3, 0, 0.01), class = "sd", group = "databaseId", coef = "outcomeGroupControl"),
    brms::prior(student_t(3, 0, 2.5), class = "sd"),
    brms::prior(lkj(1), class = "cor")
  )



  # Fit the single, powerful model
  fit <- brms::brm(
    formula = modelFormula,
    data = data,
    prior = priors,
    chains = 4,
    iter = 4000,
    warmup = 1000,
    cores = parallel::detectCores(),
    control = list(adapt_delta = 0.95)
  )
  summary(correlationModelFit)

}

