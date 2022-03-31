library(tidyverse)
library(gridExtra)
Sys.setlocale("LC_ALL", "C")


elementary_quantile_score <- function(y_true, y_pred, theta, alpha){
  ((y_true < y_pred) - alpha) * ((theta < y_pred) - (theta < y_true))
}

get_thetas <- function(value, truth, n=1001){
  tmp <- c(value, truth)
  thetas <- seq(from = min(tmp) - 0.1, to = max(tmp) + 0.1, length.out = n)
  return(thetas)
}

get_elementary_scores <- function(truth, value, quantile, n) {
  thetas <- get_thetas(value, truth, n)

  scores <- sapply(thetas,
                   function(t) mean(elementary_quantile_score(truth, value, t, quantile)))

  return(data.frame(theta=thetas, mean_score=scores))
}

quantile_score <- function(y_true, y_pred, alpha){
  return((1 * (y_true < y_pred) - alpha) * (y_pred - y_true))
}

get_murphy_plots <- function(data, alpha) {
  df <- filter(data, quantile == alpha)

  score_label <- df %>%
    mutate(qs = quantile_score(truth, value, quantile)) %>% 
    group_by(model, quantile) %>%
    summarize(mean_qs = mean(qs),
              label = paste0(unique(model), " (", format(round(mean_qs, digits = 4), scientific=FALSE, nsmall=4),
                             ")", collapse = "\n"),
              .groups = "drop")
  
  df <- df %>%
    group_by(model, quantile) %>%
    summarize(get_elementary_scores(truth, value, quantile, n=500), .groups = "drop")
  
  df <- df %>% 
    left_join(score_label, by = c("model", "quantile"))
  
  
  ymax = max(df$mean_score)
  xmax = max(df$theta)

  g <- ggplot(df) +
    facet_wrap(~quantile, strip.position="right") +
    geom_line(aes(x=theta, y=mean_score, color=label), size=0.5) +
    # facet_wrap("quantile", scales = "free", nrow=1) +
    xlab(expression(paste("Threshold ", theta))) +
    ylab("Elementary score") +
    theme_bw(base_size = 11) +
    theme(legend.justification=c(1,1), legend.position=c(0.99,0.99),
          legend.title=element_text(size=6, face = "bold"),
          legend.text=element_text(size=6),
          legend.title.align = 0, 
          legend.text.align = 0,
          # aspect.ratio = 1,
          legend.key.size = unit(0.4, "lines"),
          legend.background = element_blank(),
          panel.grid.major = element_line(size = 0.05), 
          panel.grid.minor = element_line(size = 0.05)) + 
    scale_color_brewer(palette="Set1") +
    scale_x_continuous(breaks = 0:4 / 4, labels=function(x) ifelse(x == 0, "0", x)) +
    labs(color = "Model (quantile score)") +
    expand_limits(x = xmax, y = ymax)

  return(g)
}
