# Code to create for plots including fixed-effects, random-effects (DerSimonian & Laird), Bayesian
# random effects, and Bayesian random-effects prediction interval.

source("PredictionInterval.R")
library(EvidenceSynthesis)
library(ggplot2)
library(dplyr)

plotForest <- function(data,
                       labels = paste("Site", seq_len(nrow(data))),
                       exclude = NULL,
                       xLabel = "Hazard Ratio",
                       limits = c(0.1, 10),
                       alpha = 0.05,
                       showFixedEffects = TRUE,
                       showRandomEffects = TRUE,
                       showBayesianRandomEffects = TRUE,
                       showPredictionInterval = TRUE,
                       title = NULL,
                       fileName = NULL) {
  d1 <- data.frame(
    logRr = -100,
    logLb95Ci = -100,
    logUb95Ci = -100,
    type = "header",
    label = "Source",
    exclude = FALSE
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
  d2 <- bind_rows(d2) |>
    mutate(label = labels)

  if (is.null(exclude)) {
    d2 <- d2 |>
      mutate(exclude = FALSE)
    filteredData <- data
  } else {
    d2 <- d2 |>
      mutate(exclude = !!exclude)
    filteredData <- data[!exclude, ]
  }

  d3 <- tibble()
  if (showFixedEffects || showRandomEffects) {
    ma <- summary(meta::metagen(filteredData$logRr, filteredData$seLogRr, prediction = TRUE, level.ci = 1 - alpha))
    if (showFixedEffects) {
      d3 <- bind_rows(
        d3,
        tibble(
          logRr = ma$TE.fixed,
          logLb95Ci = ma$lower.fixed,
          logUb95Ci = ma$upper.fixed,
          type = "ma",
          label = "Fixed effects",
          exclude = FALSE
        )
      )
    }
    if (showRandomEffects) {
      d3 <- bind_rows(
        d3,
        tibble(
          logRr = c(ma$TE.random, NA),
          logLb95Ci = c(ma$lower.random, NA),
          logUb95Ci = c(ma$upper.random, NA),
          type = c("ma", "maSub"),
          label = c("Random effects", sprintf("\u03C4 = %.2f (%.2f - %.2f)", ma$tau, ma$lower.tau, ma$upper.tau)),
          exclude = FALSE
        )
      )
    }
  }
  if (showBayesianRandomEffects) {
    estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(filteredData, alpha = alpha)
    d3 <- bind_rows(
      d3,
      tibble(
        logRr = c(estimate$logRr, NA),
        logLb95Ci = c(estimate$mu95Lb, NA),
        logUb95Ci = c(estimate$mu95Ub, NA),
        type = c("ma", "maSub"),
        label = c("Bayesian random effects", sprintf("\u03C4 = %.2f (%.2f - %.2f)", estimate$tau, estimate$tau95Lb, estimate$tau95Ub)),
        exclude = FALSE
      )
    )
  }
  if (showPredictionInterval) {
    if (showBayesianRandomEffects) {
      predictionInterval <- computePredictionInterval(estimate, alpha)
    } else if (showRandomEffects) {
      predictionInterval <- c(ma$predict$lower, ma$predict$upper)
    } else {
      predictionInterval <- c(ma$lower.fixed, ma$upper.fixed)
    }
    d3 <- bind_rows(
      d3,
      tibble(
        logRr = NA,
        logLb95Ci = predictionInterval[1],
        logUb95Ci = predictionInterval[2],
        type = "pi",
        label = "Prediction interval",
        exclude = FALSE
      )
    )
  }
  d <- rbind(d1, d2, d3)
  d <- d |>
    mutate(y = if_else(type == "maSub", 0.5, 1)) |>
    mutate(y = sum(y) - cumsum(y) + 1)

  maEstimates <- d |>
    filter(type == "ma", !is.na(logRr))
  diamondData <- tibble(
    x = exp(c(maEstimates$logLb95Ci, maEstimates$logRr, maEstimates$logUb95Ci, maEstimates$logRr)),
    y = c(maEstimates$y, maEstimates$y + 0.2, maEstimates$y, maEstimates$y - 0.2),
    group = rep(maEstimates$y, 4)
  )

  maBoundaryY <- d |>
    filter(type == "ma") |>
    summarise(max(y)) |>
    pull() + 0.5

  # ggplot puts whisker for infinite values, but not large values:
  plotD <- d
  plotD$logLb95Ci[is.infinite(plotD$logLb95Ci)] <- -10
  plotD$logUb95Ci[is.infinite(plotD$logUb95Ci)] <- 10
  plotD <- plotD |>
    filter(type != "ma")

  rowHeight <- 0.5
  breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
  yLimits <- c(min(d$y) - rowHeight / 2, max(d$y) + rowHeight / 2)
  rightPlot <- ggplot(plotD, aes(x = exp(logRr), y = y)) +
    geom_rect(xmin = -10, xmax = 10, ymin = 0, ymax = maBoundaryY, size = 0, fill = "#69AED5", alpha = 0.25, data = tibble(logRr = 1, y = 1)) +
    geom_segment(aes(x = x, y = y, xend = x, yend = yend), color = "#AAAAAA",  size = 0.2, data = data.frame(x = breaks, y = 0, yend = max(d$y) - 0.5)) +
    geom_segment(aes(x = x, y = y, xend = x, yend = yend), size = 0.5, data = data.frame(x = 1, y = 0, yend = max(d$y) - 0.5)) +
    geom_hline(aes(yintercept = max(d$y) - 0.5)) +
    geom_errorbarh(aes(xmin = exp(logLb95Ci), xmax = exp(logUb95Ci), color = exclude), height = 0.15) +
    geom_point(size = 3, shape = 16, aes(color = exclude)) +
    geom_polygon(aes(x = x, y = y, group = group), data = diamondData) +
    scale_color_manual(values = c("black", "#BBBBBB")) +
    scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) +
    coord_cartesian(xlim = limits, ylim = yLimits) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      axis.line.x.bottom = element_line(),
      plot.margin = grid::unit(c(0, 0, 0, -0.19), "lines")
    )
   rightPlot
  d$logLb95Ci[is.infinite(d$logLb95Ci)] <- NA
  d$logUb95Ci[is.infinite(d$logUb95Ci)] <- NA
  d$logRr[exp(d$logRr) < limits[1] | exp(d$logRr) > limits[2]] <- NA
  estimateLabels <- sprintf("%0.2f", exp(d$logRr))
  estimateLabels <- gsub("NA", "-", estimateLabels)
  estimateLabels[d$type %in% c("maSub", "pi", "header")] <- ""
  intervalLabels <- sprintf("(%0.2f - %0.2f)", exp(d$logLb95Ci), exp(d$logUb95Ci))
  intervalLabels <- gsub("NA", "", intervalLabels)
  intervalLabels <- gsub("\\( - \\)", "", intervalLabels)
  intervalLabels[d$type %in% c("maSub", "header")] <- ""

  textTable <- data.frame(
    y = rep(d$y, 3),
    x = rep(c(1, 2, 2.2), each = nrow(d)),
    label = c(as.character(d$label), estimateLabels, intervalLabels),
    fontface = rep(if_else(d$type %in% c("ma", "header", "pi"), "bold", "plain"), 3),
    color = rep(if_else(d$exclude, "#BBBBBB", "black"), 3)
  )
  textTable$label[nrow(d) + 1] <- paste(xLabel, "(95% CI)")
  leftPlot <- ggplot(textTable, aes(x = x, y = y, label = label)) +
    geom_rect(xmin = -10, xmax = 10, ymin = 0, ymax = maBoundaryY, size = 0, fill = "#69AED5", alpha = 0.25, data = tibble(x = 1, y = 1, label = "NA")) +
    geom_text(aes(fontface = fontface, color = color), size = 4, hjust = 0, vjust = 0.5) +
    geom_hline(aes(yintercept = max(d$y) - 0.5)) +
    labs(x = "", y = "") +
    scale_color_identity() +
    coord_cartesian(xlim = c(1, 2.75), ylim = yLimits) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      axis.text.x = element_text(color = "white"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      axis.line.x.bottom = element_line(),
      plot.margin = grid::unit(c(0, 0, 0, 0), "lines")
    )
  plot <- gridExtra::grid.arrange(leftPlot, rightPlot, ncol = 2, widths = c(1.5, 1), padding = unit(0, "line"), top = title)
  if (!is.null(fileName)) {
    ggsave(fileName, plot, width = 7, height = 1 + max(d$y) * 0.3, dpi = 300)
  }
  invisible(plot)
}
