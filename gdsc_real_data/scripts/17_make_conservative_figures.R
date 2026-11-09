options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

results_root <- file.path(root, "results", "conservative_alpha_experiment")
figure_dir <- file.path(results_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

metrics <- read.csv(file.path(results_root, "repeat_metrics.csv"))
method_order <- c(
  "Pretrained target-min", "Pretrained target-1se", "Overall", "Individual"
)
plot_data <- metrics[metrics$method %in% method_order, ]
plot_data$method <- factor(plot_data$method, levels = method_order)
method_labels <- c(
  `Pretrained target-min` = "Pretrained\ntarget minimum",
  `Pretrained target-1se` = "Pretrained\ntarget one-SE",
  Overall = "Overall",
  Individual = "Individual"
)

p <- ggplot(plot_data, aes(method, mse, group = repeat_id)) +
  geom_line(color = "grey78", linewidth = 0.35, alpha = 0.65) +
  geom_boxplot(aes(group = method), width = 0.55, outlier.shape = NA, fill = "white", color = "grey35") +
  geom_point(aes(color = method), size = 1.7, alpha = 0.82, position = position_jitter(width = 0.06, height = 0)) +
  scale_x_discrete(labels = method_labels) +
  scale_color_manual(
    values = c(
      `Pretrained target-min` = "#3F7BD9",
      `Pretrained target-1se` = "#D55E00",
      Overall = "#009E73",
      Individual = "#999999"
    ),
    guide = "none"
  ) +
  labs(
    title = "Target-domain alpha selection",
    subtitle = "Two thousand external-selected genes; both lasso stages use lambda.1se",
    x = NULL,
    y = "Pancreatic target-domain MSE"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 9)
  )

ggsave(
  file.path(figure_dir, "conservative_alpha_mse.png"),
  p,
  width = 8.2,
  height = 5.2,
  dpi = 300
)
capture.output(sessionInfo(), file = file.path(figure_dir, "figure_sessionInfo.txt"))

