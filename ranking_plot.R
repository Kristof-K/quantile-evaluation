library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)

MODELS <- c("IDR cond"="predictions/wind_IDR_S100_condWindAngle/",
            "NNQF + MLP"="predictions/wind_NNQF_UVS_50/",
            "QRF"="predictions/wind_QRF_lag12_noAngle/")
COL_NAMES <- c("IDR.cond"="IDR cond", "NNQF...MLP"="NNQF + MLP", "QRF"="QRF")

get_scores_by_task <- function(csv_folder) {
  collect_dfs <- data.frame()

  for (task in 1:12) {
    for (zone in 1:10) {
      csv_file <- paste0(csv_folder, "Task", task, "_", zone, ".csv")
      collect_dfs <- rbind(collect_dfs, cbind(Task=task, read.csv(csv_file, row.names=1)))
    }
  }

  mean_score_by_task <- collect_dfs %>%
    pivot_longer(cols = all_of(paste0("X", 1:99 * 0.01)), names_to = "Quantile") %>%
    mutate(Quantile = as.numeric(substring(Quantile, 2)),
           score = (1 * (value >= y) - Quantile) * (value - y)) %>%
    group_by(Task) %>%
    summarise(mean_score = mean(score, na.rm=TRUE), .groups = "drop")
  return(mean_score_by_task)
}


update_task_scores <- function() {
  file_name <- "predictions/BestWind.csv"
  df <- read.csv2(file_name)

  for (name in names(MODELS)) {
    new_scores <- get_scores_by_task(MODELS[[name]]) %>%
      rename(!!name := mean_score)
    df <- left_join(df, new_scores, by="Task")
  }
  write.csv2(df, file_name)
}


plot_results <- function() {
  label_pos <- 0.13
  benchmark_pos <- 0.38
  font_size <- 5

  raw_data <- read.csv2("predictions/BestWind.csv", row.names=1) %>%
    select(-X)
  # minus task, minus benchmark, minus idr nodels = nr of participants
  p <- ncol(raw_data) - 2 - 3
  participants <- raw_data %>% select(Task, starts_with("X")) %>%
    pivot_longer(cols = -Task, names_to = "Rank") %>%
    mutate(Rank = factor(substring(Rank, 2), ordered=TRUE,
                           levels=paste(1:p)))
  benchmark <- select(raw_data, Task, Benchmark) %>%
    rename(value = Benchmark)
  models <- select(raw_data, Task, all_of(names(COL_NAMES))) %>%
    pivot_longer(cols = -Task, names_to = "Model") %>%
    mutate(Model = factor(COL_NAMES[Model], ordered=TRUE, levels=unname(COL_NAMES)))

  part_colors <- scales::seq_gradient_pal("#132B43", "#56B1F7",
                                          "Lab")(seq(0,1,length.out=p))
  model_colors <- c("#0e9e62", "#0fe428", "#83d119")
  benchmark_c <- "#e6382c"

  curves <- ggplot(mapping = aes(x=Task, y=value)) +
    geom_line(data=participants, aes(color = Rank)) +
    scale_color_manual(values=part_colors) +
    ggnewscale::new_scale_color() +
    geom_line(data=benchmark, color = benchmark_c) +
    geom_point(data=benchmark, color = benchmark_c, shape = 2) +
    geom_line(data=models, aes(color = Model), size=1.1) +
    geom_point(data=models, aes(color = Model), size=1.4, shape=1) +
    scale_x_continuous(breaks = 1:12) +
    scale_color_manual(values=model_colors) +
    xlab("Task") +
    ylab("Mean pinball score") +
    ggtitle(paste("Final Scores of Wind Track")) +
    theme_bw() +
    theme(text = element_text(size = 16), axis.text = element_text(size = 13),
          legend.position = "none")

  raw_data <- read.csv2("predictions/BestWindFinal.csv") %>%
    mutate(Text = paste0(format(round(Rating * 100, 1), nsmall=1), "%")) %>%
    arrange(desc(Rank)) %>%
    mutate(Rank = factor(paste(Rank), ordered=TRUE, levels=paste(Rank)))

  bar_colors <- c(setNames(part_colors, paste(1:p)), "Benchmark"="red",
                  setNames(model_colors, names(MODELS)))
  draw_c <- unname(bar_colors[raw_data$Label])

  bars <- ggplot(raw_data, aes(x=Rank, y=Rating, fill=Rank)) +
    geom_col(show.legend=FALSE) +
    geom_text(aes(y=label_pos, label=Text),
              color=c(benchmark_c, rep("black", p - 5 + 3), rep("white", 5)),
              angle=90, size=font_size) +
    geom_text(data=filter(raw_data, Label %in% names(MODELS)), aes(y=Rating * 0.75,
                                                              label=Label),
              angle=90,size=font_size) +
    annotate("text", x=paste(p+1), y=benchmark_pos, label="Benchmark",
             color=benchmark_c, angle=90, size=font_size) +
    scale_fill_manual(values = draw_c) +
    scale_color_manual(values = draw_c) +
    scale_x_discrete(breaks = paste(1:p)) +
    ggtitle("Linear Weighted Skill Score") +
    theme_bw() +
    theme(text = element_text(size = 16), axis.text = element_text(size = 13),
          legend.position = "none")
  # for load: ifelse(Rating <=0.3, Rating*1.33, Rating * 0.75)

  g <- grid.arrange(curves, bars, nrow=2)
  # save with height 7.5
  return(g)
}

g <- plot_results()