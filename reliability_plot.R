library(tidyverse)
library(latex2exp)
library(isotone)

DIGITS <- 3

# Reliability Diagrams
# Code by Daniel with adaptions
reldiag = function(x, y, alpha = 0.5, n_resamples = 999, digits = 3, region_level = 0.9) {
  if (length(alpha) > 1) {
    alpha = alpha[1]
  }

  require(isotone)
  pava = function(x,y){
    # In case of ties, isotone::gpava uses the conditional mean instead of quantile, try e.g.,
    # gpava(c(-1,-1,-1),c(-1,0,0),solver = weighted.median,ties = "secondary")

    # Wrong fix: The following step replaces y values with the respective quantile in case of ties
    # y = unlist(lapply(split(y,x),function(y) rep(quantile(y,alpha,type = 1),length(y))),use.names = FALSE)

    # New fix: Use ranking of predictor values and break ties by ordering the corresponding instances in order of decreasing observations
    ranking = match(1:length(x),order(x,y,decreasing = c(FALSE,TRUE)))

    return(gpava(ranking,y,solver = weighted.fractile,p = alpha,ties = "secondary")$x)
  }
  score = function(x, y) mean((as.numeric(x >= y) - alpha)*(x-y))
  marg = function(x) quantile(x, alpha, type = 1)
  identif = function(x, y) as.numeric(x > y) - alpha
  score_label = "QS "

  ord_x = order(x)
  x = x[ord_x]
  y = y[ord_x]

  x_rc = pava(x,y)

  res = y - x

  s = score(x,y)
  c_rc_ucond = optim(par = 0,fn = function(c) score(x+c,y),method = "Brent",lower = min(res),upper = max(res))$par
  s_rc_ucond = score(x + c_rc_ucond,y)
  s_rc = score(x_rc,y)
  s_mg = score(marg(y),y)

  mcb = s - s_rc
  umcb = s - s_rc_ucond
  cmcb = s_rc_ucond - s_rc
  dsc = s_mg - s_rc
  unc = s_mg

  # The Score is exactly equal to uMCB + cMCB - DSC + UNC.
  # However, when rounding the values there may be slight discrepancies between the rounded values.
  # We avoid this for didactic reasons by computing the score from the rounded values.
  s = sum(round(c(umcb,cmcb,-dsc,unc), digits))


  # test: mean identification zero? (t-test)
  # v = identif(x,y)
  # t = sqrt(length(v)) * mean(v)/sd(v)
  # pval_ucond = 1 - abs(pt(t,length(v)-1) - 0.5)*2

  # Unconditional calibration test
  # Coverage test: One-sided Binomial tests with Bonferroni correction
  eps = 10^-10 # avoid numerical artifacts by assuming that values with an absolute difference of less than eps are identical
  hard_cov <- sum(y < x - eps)
  soft_cov <- sum(y < x + eps)

  pval_hard = dbinom(hard_cov,length(y),alpha) + pbinom(hard_cov,length(y),alpha,FALSE)
  pval_soft = pbinom(soft_cov,size = length(y),prob = alpha)
  pval_ucond = min(pval_hard,pval_soft,0.5)*2
  # print(paste0("p-Values: hard ",pval_hard,", soft ",pval_soft))

  n_samples = n_resamples + 1 # total number of samples including observed sample
  low = floor(n_samples * (1-region_level)/2)
  up = n_samples - low
  pval_digits = ceiling(log(n_samples,10))

  resamples = sapply(1:n_resamples,function(i) x + sample(res,length(y),replace = TRUE))

  x_rc_resamples = apply(resamples, 2, function(y) pava(x,y))
  x_rc_resamples_sorted = apply(cbind(x_rc,x_rc_resamples),1,sort) - marg(res) # includes observed values + bias corrected (shifted by mean residual)

  ran_x = range(x)

  mcb_resamples = sapply(1:n_resamples,function(i) score(x,resamples[,i]) - score(x_rc_resamples[,i],resamples[,i]))
  mcb_bounds = sort(c(mcb,mcb_resamples))[c(low,up)]

  rank_obs = tail(rank(c(mcb_resamples,mcb)),1)
  pval = 1 - (rank_obs - 1)/(n_resamples + 1)

  results <- data.frame(quantile = alpha, x = x, y = y, x_rc = x_rc,
                        lower = x_rc_resamples_sorted[low,],
                        upper = x_rc_resamples_sorted[up,],
                        score = s,
                        umcb = umcb, cmcb = cmcb, mcb = mcb, dsc = dsc, unc = unc,
                        pval_cond = pval, pval_ucond = pval_ucond)
}

get_inset <- function(df, xmin=0, xmax=0, ...){
  ggplot(df, aes(x)) +
    geom_histogram(fill="gray", col="black", size=0.2, bins = 8) +
    theme_classic( base_size=5.5) +
    theme(axis.line.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())+
    expand_limits(x = c(xmin, xmax)) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE))
}

annotation_custom2 <- function (grob, xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, data) {
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = GeomCustomAnn,
        inherit.aes = TRUE, params = list(grob = grob,
                                          xmin = xmin, xmax = xmax,
                                          ymin = ymin, ymax = ymax))
}

get_annotation <- function(df, model, xmax){
  inset_plot <- get_inset(df, xmax=xmax)
  annotation_custom2(grob=ggplotGrob(inset_plot),
                     data = subset(df, model == unique(df$model)),
                     ymin = min(facet_lims$mn), ymax=max(facet_lims$mx)/4, xmin=max(facet_lims$mx)/1.5, xmax=0.975*max(facet_lims$mx))
  #ymin = 0, ymax=750, xmin=1500, xmax=2750)
  #ymin = -2000, ymax=4000, xmin=10000, xmax=15000)
}


get_reliability_plots <- function(df, quantile_level, n_resamples, add_layer="hist1",
                                  load_interrim=FALSE) {
  df <- filter(df, quantile == quantile_level)

  if (!load_interrim) {
    recal_and_bands <- df %>% group_by(model, quantile) %>%
      summarise(reldiag(value, truth, alpha=quantile, n_resamples=n_resamples), .groups="drop")

    # target variable is normalized
    recal_and_bands$lower <- pmin(1, pmax(0, recal_and_bands$lower))
    recal_and_bands$upper <- pmin(1, pmax(0, recal_and_bands$upper))

    write.csv(recal_and_bands, "interrim/wind_pava_6.csv")
  } else {
    recal_and_bands <- read.csv("interrim/wind_pava_6.csv", row.names=1)
  }

  recal_and_bands <- recal_and_bands %>%
    mutate_at(c("x_rc", "lower", "upper"), ~ replace(., .<0, 0))

  scores <- recal_and_bands %>%
    group_by(model, quantile) %>%
    distinct(across(score:pval_ucond)) %>%
    mutate(label = paste0(c("S ", "uMCB ","cMCB ","DSC ","UNC "),
                          format(round(c(score, umcb, cmcb, dsc, unc), digits=DIGITS), scientific=FALSE, nsmall=DIGITS),
                          c("", paste0(" [p = ", format(round(pval_ucond, digits = 2), nsmall=2),"]"), "", "", ""),
                          c("", "", paste0(" [p = ", format(round(pval_cond, digits = 2), nsmall=2),"]"), "", ""),
                          collapse = "\n"))

  facet_lims <- recal_and_bands %>%
    group_by(model, quantile) %>%
    summarize(mn = min(c(x, x_rc, lower, upper)), mx = max(c(x, x_rc, lower, upper)),
              .groups = "drop")

  main_plot <- ggplot(recal_and_bands, aes(x, x_rc, group=model)) +
    facet_grid(quantile ~ model) +
  {if (add_layer == "points") geom_point(aes(x, y), alpha=0.05, size=0.05)} +
  {if (add_layer == "bin2d") geom_bin2d(aes(x=x, y=y, alpha=log(..count..)), fill="#746b6b", bins=30)} +
  {if (add_layer == "bin2d") scale_alpha_continuous("Count", labels=function(x) floor(exp(x)))} +
  {if (add_layer == "hist2") geom_histogram(mapping = aes(x = x,y = 0.2*max(facet_lims$mx)*after_stat(count/max(count))),
                                            bins = 60, colour = "grey", fill = NA, size=0.3)} +
    geom_abline(intercept = 0 , slope = 1, colour="grey70") +
    #geom_point(color = "red", size=0.5) +
    # geom_step(color = "red", direction = "vh") +
    geom_smooth(aes(ymin = lower, ymax = upper), linetype = 0, stat = "identity", fill = "skyblue3") +
  {if (add_layer == "points") geom_line(aes(x,lower), color = "deepskyblue2")} +
  {if (add_layer == "points") geom_line(aes(x,upper), color = "deepskyblue2")} +
    geom_line(color = "firebrick3") +
  {if (add_layer == "hist1") geom_blank(data = facet_lims, aes(x = mx, y = mx))} +
  {if (add_layer == "hist1") geom_blank(data = facet_lims, aes(x = mn, y = mn))} +
    # geom_rug(sides = "b", alpha = 0.2, size = 0.25) +
    xlab("Forecast value") +
    ylab("Conditional quantile") +
    #labs(title = paste0(model, ":\n", target))  +
    geom_label(data = scores, mapping = aes(x = -Inf, y = Inf, label = label),
               size = 6*0.36, hjust = 0, vjust = 1, label.size = NA, alpha=0, label.padding = unit(1, "lines")) +
    scale_x_continuous(guide = guide_axis(check.overlap = TRUE), breaks=0:4 / 4,
                       labels=function(x) ifelse(x == 0, "0", x)) +
    scale_y_continuous(breaks=0:4 / 4, labels=function(x) ifelse(x == 0, "0", x)) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major = element_line(size = 0.05),
          panel.grid.minor = element_line(size = 0.05),
          legend.justification=c(1,0), legend.position=c(0.99,0.01),
          legend.background=element_blank(), legend.box.background=element_blank(),
          legend.title=element_text(size=6, face = "bold"),
          legend.text=element_text(size=6), legend.key.size = unit(0.4, "lines"),
          strip.background.x = element_blank(),  # no facet boxes in x direction
          strip.text.x = element_blank())        # no facet texts in x direction)

  if (add_layer == "hist1") {
    insets <- recal_and_bands %>%
      group_by(model, quantile) %>%
      group_map(get_annotation, facet_lims=facet_lims, .keep=TRUE)
    main_plot <- main_plot + insets
  }
  return(main_plot)
}
