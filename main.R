library(tidyverse)
library(lubridate)
library(gridExtra)    # to stack different plots

source("forecast_plot.R")
source("coverage_plot.R")
source("reliability_plot.R")
source("murphy_plot.R")

RELIABILITY_QUANTIL <- 0.75
MURPHY_QUANTIL <- 0.75

B_COVERAGE <- 99
N_RES_RELIABILITY <- 99

SOURCES <- c(
    "IDR cond" = "predictions/wind_IDR_S100_condWindAngle/",
    "NNQF + MLP" = "predictions/wind_NNQF_UVS_h100_50/",
    "QRF" = "predictions/wind_QRF_lag12_noAngle/"
)


get_forecasts <- function(tasks=1:12, zones=NA) {
  if (is.na(zones)) {
    zones <- 1:10
  }
  collect_data <- data.frame()

  files <- character()    # store all files that should be read in
  for (z in zones) {      # product from tasks and zones
    files <- c(files, paste0("Task", tasks, "_", z, ".csv"))
  }

  for (name in names(SOURCES)) {
    for (f in files) {
      csv_path <- paste0(SOURCES[[name]], f)
      df <- read.csv(csv_path, row.names=1) %>%
        rename(target_end_date = time, truth=y, location = zoneid) %>%
        pivot_longer(cols = paste0("X", 1:99 * 0.01), names_to = "quantile") %>%
        mutate(quantile = as.numeric(substring(quantile, 2)),
              target_end_date = ymd_hms(target_end_date),
              model = name) %>%
        select(-score)

      # wind track had NA values, filter them
      df <- filter(df, !is.na(truth))

      collect_data <- rbind(collect_data, df)
    }
  }
  return(collect_data)
}

calc_scores <- function(df) {
  scores <- df %>% group_by(model) %>%
    summarise(mean_pinball_score = mean((1 * (value > truth) - quantile) * (value - truth)),
              .groups = "drop")
  print(scores)
}

assemble_plot <- function(name, task=1, zone=NA, rel_add="hist1") {
  print("Starting...")
  df <- get_forecasts(tasks = task, zone = zone)
  print("get_forecast done")

  calc_scores(df)

  forecast_plots <- get_forecast_plots(df)
  print("forecast_plot done")

  coverage_plots <- get_coverage_plots(df, B_COVERAGE)
  print("coverage_plot done")

  reliability_plots <- get_reliability_plots(df, RELIABILITY_QUANTIL, N_RES_RELIABILITY,
                                             add_layer=rel_add, load_interrim=TRUE)
  print("reliability_plot done")

  murphy_plots <- get_murphy_plots(df, MURPHY_QUANTIL)
  # ggtitle(expr(paste("Murphy Diagram (", alpha, " = ", !!MURPHY_QUANTIL, ")")))
  print("murphy_plot done")

  assembled_plot <- grid.arrange(forecast_plots, coverage_plots, reliability_plots,
                                 murphy_plots, ncol=1,
                                 heights=c(0.258, 0.25, 0.25, 0.242))
  # for fix aspect ratios use 135 as width
  ggsave(paste0("figures/Wind/", name, ".pdf"), plot=assembled_plot,
         width=160, height=200, unit="mm", device = "pdf", dpi=300)
  print("Finished.")
  return(assembled_plot)
}

pl <- assemble_plot(name="Figure10_33", task=1:12, zone=1, rel_add="points")
