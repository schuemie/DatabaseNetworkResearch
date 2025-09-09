# Code to create for plots including fixed-effects, random-effects (DerSimonian & Laird), Bayesian
# random effects, and Bayesian random-effects prediction interval.

source("PredictionInterval.R")
library(EvidenceSynthesis)
library(ggplot2)
library(dplyr)

plotForest <- function(data,
                       labels = paste("Site", seq_len(nrow(data))),
                       xLabel = "Hazard Ratio",
                       limits = c(0.1, 10),
                       alpha = 0.05,
                       showFixedEffects = TRUE,
                       showRandomEffects = TRUE,
                       showBayesianRandomEffects = TRUE,
                       showPredictionInterval = TRUE,
                       fileName = NULL) {
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
  d2 <- bind_rows(d2) |>
    mutate(label = labels)

  d3 <- tibble()
  diamondData <- tibble()
  if (showFixedEffects || showRandomEffects) {
    ma <- summary(meta::metagen(data$logRr, data$seLogRr, prediction = TRUE))
    if (showFixedEffects) {
      d3 <- bind_rows(
        d3,
        tibble(
          logRr = ma$TE.fixed,
          logLb95Ci = ma$lower.fixed,
          logUb95Ci = ma$upper.fixed,
          type = "ma",
          label = "Fixed FX"
        )
      )
    }
    if (showRandomEffects) {
      d3 <- bind_rows(
        d3,
        tibble(
          logRr = ma$TE.random,
          logLb95Ci = ma$lower.random,
          logUb95Ci = ma$upper.random,
          type = "ma",
          label = sprintf("Random FX (tau = %0.2f)", ma$tau)
        )
      )
    }
  }
  if (showBayesianRandomEffects) {
    estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(data)
    d3 <- bind_rows(
      d3,
      tibble(
        logRr = estimate$logRr,
        logLb95Ci = estimate$mu95Lb,
        logUb95Ci = estimate$mu95Ub,
        type = "ma",
        label = sprintf("Bayesian RFX (tau = %.2f)", estimate$tau)
      )
    )
  }
  if (showPredictionInterval) {
    if (showBayesianRandomEffects) {
      predictionInterval <- computePredictionInterval(estimate)
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
        label = "Prediction interval"
      )
    )
  }
  d <- rbind(d1, d2, d3)
  d <- d |>
    mutate(y = seq(from = nrow(d), to = 1, by = -1))

  maEstimates <- d |>
    filter(type == "ma")
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
  p <- ggplot(plotD, aes(x = exp(logRr), y = y)) +
    geom_rect(xmin = -10, xmax = 10, ymin = 0, ymax = maBoundaryY, size = 0, fill = "#69AED5", alpha = 0.25, data = tibble(logRr = 1, y = 1)) +
    geom_segment(aes(x = x, y = y, xend = x, yend = yend), color = "#AAAAAA",  size = 0.2, data = data.frame(x = breaks, y = 0, yend = max(d$y) - 0.5)) +
    geom_segment(aes(x = x, y = y, xend = x, yend = yend), size = 0.5, data = data.frame(x = 1, y = 0, yend = max(d$y) - 0.5)) +
    geom_errorbarh(aes(
      xmin = exp(logLb95Ci),
      xmax = exp(logUb95Ci)
    ), height = 0.15) +
    geom_point(size = 3, shape = 23, aes(fill = type)) +
    geom_polygon(aes(x = x, y = y, group = group), data = diamondData) +
    scale_fill_manual(values = c("#000000", "#000000", "#FFFFFF", "#FFFFFF", "#FFFFFF", "#FFFFFF")) +
    scale_x_continuous(xLabel, trans = "log10", breaks = breaks, labels = breaks) +
    # scale_y_continuous(limits = yLimits) +
    coord_cartesian(xlim = limits, ylim = yLimits) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      axis.text.x = element_text(colour = "black"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = grid::unit(c(0, 0, 0, -0.19), "lines")
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
    x = rep(c(1, 2), each = nrow(d)),
    label = c(as.character(d$label), labels),
    stringsAsFactors = FALSE
  )
  labels$label[nrow(d) + 1] <- paste(xLabel, "(95% CI)")
  data_table <- ggplot(labels, aes(x = x, y = y, label = label)) +
    geom_rect(xmin = -10, xmax = 10, ymin = 0, ymax = maBoundaryY, size = 0, fill = "#69AED5", alpha = 0.25, data = tibble(x = 1, y = 1, label = "NA")) +
    geom_text(size = 4, hjust = 0, vjust = 0.5) +
    geom_hline(aes(yintercept = nrow(d) - 0.5)) +
    labs(x = "", y = "") +
    coord_cartesian(xlim = c(1, 2.75), ylim = yLimits) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      legend.position = "none",
      panel.border = element_blank(),
      axis.text.x = element_text(colour = "white"),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.y = element_blank(),
      plot.margin = grid::unit(c(0, 0, 0, 0), "lines")
    )


  # p <- p +
  #   geom_text(aes(x = x, y = y, label = label), size = 4, hjust = 0, vjust = 0.5, data = labels |> mutate(logRr = log(x)))
  plot <- gridExtra::grid.arrange(data_table, p, ncol = 2, widths = c(1.5, 1), padding = unit(0, "line"))
  if (!is.null(fileName)) {
    ggsave(fileName, plot, width = 7, height = 1 + nrow(d) * 0.35, dpi = 300)
  }
  invisible(plot)
}
