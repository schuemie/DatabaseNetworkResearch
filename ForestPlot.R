# Code to create for plots including fixed-effects, random-effects (DerSimonian & Laird), Bayesian
# random effects, and Bayesian random-effects prediction interval.

library(EvidenceSynthesis)
library(ggplot2)

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


plotForest <- function(data,
                       labels = paste("Site", seq_len(nrow(data))),
                       xLabel = "Hazard Ratio",
                       limits = c(0.1, 10),
                       alpha = 0.05) {
  d1 <- data.frame(
    logRr = -100,
    logLb95Ci = -100,
    logUb95Ci = -100,
    type = "header",
    label = "Source"
  )
  getEstimate <- function(approximation) {
    ci <- suppressMessages(computeConfidenceInterval(
      approximation = approximation,
      alpha = alpha
    ))
    return(data.frame(
      logRr = ci$logRr,
      logLb95Ci = log(ci$lb),
      logUb95Ci = log(ci$ub),
      type = "db"
    ))
  }
  d2 <- lapply(split(data, seq_len(nrow(data))), getEstimate)
  d2 <- do.call(rbind, d2)
  d2$label <- labels

  ma <- summary(meta::metagen(data$logRr, data$seLogRr))

  estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(data)
  predictionInterval <- computePredictionInterval(estimate)
  d3 <- data.frame(
    logRr = c(ma$fixed$TE, ma$random$TE,  estimate$logRr, predictionInterval[2]),
    logLb95Ci = c(ma$fixed$lower, ma$random$lower,  estimate$mu95Lb, predictionInterval[1]),
    logUb95Ci = c(ma$fixed$upper, ma$random$upper,  estimate$mu95Ub, predictionInterval[3]),
    type = c("ma1", "ma2", "ma3", "ma4"),
    label = c("Fixed FX", sprintf("Random FX (tau = %0.2f)", ma$tau), sprintf("Bayesian RFX (tau = %.2f)", estimate$tau), "Prediction interval")
  )

  d <- rbind(d1, d2, d3)
  d$y <- seq(nrow(d), 1)
  # d$name <- factor(d$name, levels = c(d3$name, rev(as.character(labels)), "Source"))

  # ggplot puts whisker for infinite values, but not large values:
  plotD <- d
  plotD$logLb95Ci[is.infinite(plotD$logLb95Ci)] <- -10
  plotD$logUb95Ci[is.infinite(plotD$logUb95Ci)] <- 10

  rowHeight <- 0.8
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  yLimits <- c(min(d$y) - rowHeight / 2, max(d$y) + rowHeight / 2)
  p <- ggplot2::ggplot(plotD, ggplot2::aes(x = exp(.data$logRr), y = .data$y)) +
    ggplot2::geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
    ggplot2::geom_vline(xintercept = 1, size = 0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(
      xmin = exp(.data$logLb95Ci),
      xmax = exp(.data$logUb95Ci)
    ), height = 0.15) +
    ggplot2::geom_point(size = 3, shape = 23, ggplot2::aes(fill = .data$type)) +
    ggplot2::scale_fill_manual(values = c("#000000", "#000000", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")) +
    ggplot2::scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) +
    ggplot2::coord_cartesian(xlim = limits, ylim = yLimits) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines")
    )

  # p
  d$logLb95Ci[is.infinite(d$logLb95Ci)] <- NA
  d$logUb95Ci[is.infinite(d$logUb95Ci)] <- NA
  d$logRr[exp(d$logRr) < limits[1] | exp(d$logRr) > limits[2]] <- NA
  labels <- sprintf("%0.2f (%0.2f - %0.2f)", exp(d$logRr), exp(d$logLb95Ci), exp(d$logUb95Ci))
  labels <- gsub("NA", "", labels)
  labels <- gsub(" \\( - \\)", "-", labels)
  labels <- data.frame(
    y = rep(d$y, 2),
    x = rep(1:2, each = nrow(d)),
    label = c(as.character(d$label), labels),
    stringsAsFactors = FALSE
  )
  labels$label[nrow(d) + 1] <- paste(xLabel, "(95% CI)")
  data_table <- ggplot2::ggplot(labels, ggplot2::aes(
    x = .data$x,
    y = .data$y,
    label = .data$label
  )) +
    ggplot2::geom_text(size = 4, hjust = 0, vjust = 0.5) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = nrow(d) - 0.5)) +
    ggplot2::scale_y_continuous(limits = yLimits) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none",
      panel.border = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = "white"),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = "white"),
      plot.margin = grid::unit(c(0, 0, 0.1, 0), "lines")
    ) +
    ggplot2::labs(x = "", y = "") +
    ggplot2::coord_cartesian(xlim = c(1, 3))

  plot <- gridExtra::grid.arrange(data_table, p, ncol = 2, widths = c(1.5,1))
  invisible(plot)
}
