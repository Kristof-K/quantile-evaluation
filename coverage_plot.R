library(tidyverse)
Sys.setlocale("LC_ALL", "C")

coverage <- function(df) {
  exceedance <- group_by(df, model, quantile) %>%
    summarize(l = mean(truth < value), u=mean(truth <= value), .groups = "drop")
  return(exceedance)
}

# get critical value of binomial test (used for consistency bands)
get_consistency_values <- function(df) {
  # for probability p, number of trials n and level alpha, determine [k_low, k_up],
  # so that a Bin(n, p) random variable is with probability nom_level in this interval
  get_interval <- function(p, n, nom_level) {
    k_low <- qbinom((1 - nom_level) / 2, n, p) - 1
    k_low <- pmax(0, k_low)
    k_up <- qbinom(1 - (1 - nom_level) / 2, n, p)
    k_interval <- data.frame(k_low=k_low, k_up=k_up)
    colnames(k_interval) <- paste0(c("lower", "upper"), nom_level * 100)
    return(k_interval)
  }

  consistency_bands <- df %>%
    group_by(model, quantile) %>%
    summarize(count = n(), .groups = "drop") %>%
    mutate(get_interval(quantile, count, 0.9) / count,
           get_interval(quantile, count, 0.5) / count) %>%
    select(-count)
  return(consistency_bands)
}

# use resampling to get confidence values
get_confidence_values <- function(df, B) {
  coverage_df <- data.frame()
  sep_colnames <- "__"  # cannot appear in normal column names

  # get data in wide format to resample more easily
  df_wide <- pivot_wider(df, id_cols=c(-model, -quantile), names_from=c(model, quantile),
                         values_from=value, names_sep=sep_colnames)
  # store column names that were created in pivot_wider
  quantile_fcst_cols <- apply(crossing(unique(df$model), as.character(unique(df$quantile))),
                              1, function(row) paste(row, collapse=sep_colnames))

  for(i in 1:B){
    df_resampled <- df_wide %>%
      slice_sample(n = nrow(df_wide), replace = TRUE) %>%
      pivot_longer(cols=all_of(quantile_fcst_cols), names_to=c("model", "quantile"),
                   names_sep=sep_colnames) %>%
      mutate(quantile = as.numeric(quantile))
    coverage_df <- bind_rows(coverage_df, coverage(df_resampled))
  }

  # compute CIs from bootstrapped coverage
  confidence_bands <- coverage_df %>%
    group_by(model, quantile) %>%
    summarize(l_5 = quantile(l, 0.05),
              l_95 = quantile(l, 0.95),
              l_25 = quantile(l, 0.25),
              l_75 = quantile(l, 0.75),
              u_5 = quantile(u, 0.05),
              u_95 = quantile(u, 0.95),
              u_25 = quantile(u, 0.25),
              u_75 = quantile(u, 0.75),
              .groups = "drop")
  return(confidence_bands)
}

plot_coverage <- function(df, B = 1000, type = "confidence", difference = FALSE) {
  margins <- as.numeric(theme_bw()$plot.margin)
  margins[2] <- margins[2] + 20    # manipulate right margin

  # some customizations used in all plots
  my_theme <- list(
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = function(x) ifelse(x == 0, "0", x)),
    scale_y_continuous(labels = function(y) ifelse(y == 0, "0", y)),
    xlab("Quantile level"),
    ylab("Coverage"),
    theme_bw(base_size = 11),
    theme(panel.grid.major = element_line(size = 0.05),
          panel.grid.minor = element_line(size = 0.05),
          plot.margin = unit(margins, "points"),
          strip.background.x = element_blank(),  # no facet boxes in x direction
          strip.text.x = element_blank())        # no facet texts in x direction

  )

  # compute coverage on full sample
  coverage_full <- coverage(df)

  if (type == "confidence") {
    bands <- get_confidence_values(df, B)

    results <- coverage_full %>% left_join(bands, by = c("model", "quantile"))

    if (!difference) {
      g <- ggplot(results) +
        facet_wrap("model", ncol = 3) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70")+
        geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkblue", alpha = 0.2) +
        geom_ribbon(aes(x = quantile, ymin = u_5, ymax = u_95), fill = "darkred", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme
    } else {
      results <- mutate_at(results, vars("l", "u", "l_5", "l_25", "l_75", "l_95", "u_5", "u_25","u_75", "u_95"),
                           list(~ . - quantile))
      g <- ggplot(results) +
        facet_wrap("model") +
        geom_hline(yintercept = 0, size = 0.3, linetype = "solid", color = "darkgray") +
        geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkblue", alpha = 0.2) +
        geom_ribbon(aes(x = quantile, ymin = u_5, ymax = u_95), fill = "darkred", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme
    }

  } else if (type == "consistency") {
    bands <- get_consistency_values(df)

    results <- coverage_full %>% left_join(bands, by = c("model", "quantile"))

    if (!difference) {
      g <- ggplot(results) +
        facet_wrap("model", ncol = 3) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70")+
        geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
        geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme
    } else {
      results <- mutate_at(results, vars("l", "u", "lower50", "upper50", "lower90", "upper90"),
                           list(~ . - quantile))
      g <- ggplot(results) +
        facet_wrap("model") +
        geom_hline(yintercept = 0, size = 0.3, linetype = "solid", color = "darkgray") +
        geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
        geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme
    }
  }

  return(g)
}

get_coverage_plots <- function(df, B) {
  g <- plot_coverage(df, B=B, type="confidence", difference=FALSE)
  # g <- plot_coverage(df, B=B, type="consistency", difference=TRUE)
  return(g)
}
