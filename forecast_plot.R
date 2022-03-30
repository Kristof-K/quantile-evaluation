library(tidyverse)
library(lubridate)
Sys.setlocale("LC_ALL", "C")

get_forecast_plots <- function(df, time_interval=NA) {
  if (is.na(time_interval)) {
    # filter for each model and quantile first 24 * 10 entries
    df <- group_by(df, model, quantile) %>%
      group_modify(function(df, model, quantile) df[17:(24 * 7 + 7),]) %>%
      ungroup()
  } else {
    df <- filter(df, target_end_date >= time_interval[1],
                 target_end_date <= time_interval[2],
                 location == 1)
  }

  median_and_truth <- filter(df, quantile == 0.5)

  levels <- c(90, 80, 50, 20)
  abs_alpha <- setNames(1 - levels / 100, paste0(levels, "%"))
  # since we draw bands over each other we need to know what are the alpha differences
  rel_alpha <- abs_alpha - lag(abs_alpha, default=0)
  lower_quantiles <-(100 - levels) / 2
  upper_quantiles <- 100 - lower_quantiles

  df_lower <- filter(df, (100 * quantile) %in% lower_quantiles) %>%
    mutate(level = paste0((1 - quantile * 2) * 100, "%")) %>%
    rename(lower = value) %>%
    select(model, target_end_date, level, lower)

  df_upper <- filter(df, (100 * quantile) %in% upper_quantiles) %>%
    mutate(level = paste0((1 - (1 - quantile) * 2) * 100, "%")) %>%
    rename(upper = value) %>%
    select(model, target_end_date, level, upper)

  bands <- left_join(df_lower, df_upper, by=c("model", "target_end_date", "level"))

  g <- ggplot(mapping=aes(x=target_end_date)) +
    facet_wrap(~model, nrow=1) +
    geom_ribbon(data=bands,
                aes(ymin=lower, ymax=upper, group=level, alpha=level), fill='deepskyblue4') +
    geom_line(data=median_and_truth, aes(y=value, group=1, color="Median")) +
    geom_line(data=median_and_truth, aes(y=truth, group=1, color="Truth")) +
    scale_x_datetime(date_breaks = "3 days") +
    scale_color_manual(name=NULL, values=c("Truth"="darkred", "Median"="deepskyblue3"),
                       guide=guide_legend(order=1, direction="horizontal")) +
    scale_alpha_manual(name="Prediction interval", values=rel_alpha,
                       guide=guide_legend(order=2, override.aes=list(alpha=abs_alpha),
                                          direction="horizontal", title.position="top", title.hjust=1)) +
    xlab(NULL) +
    ylab('Wind power output') +
    expand_limits(y = 1.15) +
    scale_y_continuous(breaks=0:4 / 4, labels=function(x) ifelse(x == 0, "0", x)) +
    theme_bw(base_size = 11)  +
    theme(legend.justification=c(1,1), legend.position=c(0.99,0.99), legend.box="horizontal",
          legend.background=element_blank(), legend.box.background=element_blank(),
          legend.box.just="right", legend.title=element_text(size=6, face = "bold"),
          legend.text=element_text(size=6), legend.margin=margin(1,2,1,70),
          legend.key.size = unit(0.4, "lines"))
  return(g)
}
