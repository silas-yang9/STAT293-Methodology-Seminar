options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

results_root <- file.path(root, "results", "full_experiment")
figure_dir <- file.path(results_root, "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
primary_dir <- file.path(results_root, "primary_1000_ln_ic50_1se")
method_order <- c("Target mean", "Individual", "Overall", "Pretrained")
method_colors <- c(
  "Target mean" = "#9E9E9E",
  "Individual" = "#D95F02",
  "Overall" = "#1B9E77",
  "Pretrained" = "#3B6FB6"
)
base_theme <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    legend.position = "bottom"
  )

primary_metrics <- read.csv(file.path(primary_dir, "repeat_metrics.csv"))
primary_metrics$method <- factor(primary_metrics$method, levels = method_order)
p_mse <- ggplot(primary_metrics, aes(x = method, y = mse, color = method)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, linewidth = 0.6) +
  geom_jitter(width = 0.10, height = 0, alpha = 0.7, size = 1.7) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  labs(
    title = "Pancreatic-cancer prediction error across repeated splits",
    x = NULL,
    y = "Target-domain MSE"
  ) +
  base_theme + theme(legend.position = "none")
ggsave(file.path(figure_dir, "primary_mse_repeated_cv.png"), p_mse, width = 7.0, height = 4.5, dpi = 300)

alpha <- read.csv(file.path(primary_dir, "alpha_distribution.csv"))
alpha$alpha <- as.numeric(as.character(alpha$alpha))
p_alpha <- ggplot(alpha, aes(x = factor(alpha), y = proportion)) +
  geom_col(fill = "#3B6FB6", width = 0.72) +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, max(alpha$proportion) * 1.12)) +
  labs(
    title = "Cross-validated alpha choices in the primary analysis",
    x = expression(alpha),
    y = "Share of outer fits"
  ) +
  base_theme
ggsave(file.path(figure_dir, "primary_alpha_distribution.png"), p_alpha, width = 6.4, height = 4.2, dpi = 300)

predictions <- read.csv(file.path(primary_dir, "predictions.csv"))
average_predictions <- aggregate(
  predicted ~ model_id + method,
  data = predictions,
  FUN = mean
)
observed <- aggregate(observed ~ model_id, data = predictions, FUN = function(x) unique(x)[1])
average_predictions <- merge(average_predictions, observed, by = "model_id")
average_predictions$method <- factor(average_predictions$method, levels = method_order)
p_prediction <- ggplot(average_predictions, aes(x = observed, y = predicted, color = method)) +
  geom_abline(slope = 1, intercept = 0, color = "grey65", linetype = 2) +
  geom_point(size = 1.8, alpha = 0.85) +
  facet_wrap(~method, ncol = 2) +
  scale_color_manual(values = method_colors, drop = FALSE) +
  labs(
    title = "Observed and cross-validated gemcitabine response",
    subtitle = "Predictions are averaged across 20 repeated five-fold splits",
    x = "Observed GDSC2 LN_IC50",
    y = "Mean out-of-fold prediction"
  ) +
  base_theme + theme(legend.position = "none")
ggsave(file.path(figure_dir, "primary_observed_vs_predicted.png"), p_prediction, width = 7.2, height = 6.1, dpi = 300)

coefficient_stability <- read.csv(file.path(primary_dir, "coefficient_stability.csv"))
feature_map <- read.csv(file.path(root, "processed", "selected_features_top2000_annotated.csv"))
coefficient_stability <- merge(coefficient_stability, feature_map[, c("ensembl_gene_id", "gene_symbol")], by = "ensembl_gene_id", all.x = TRUE)
top_genes <- coefficient_stability[coefficient_stability$method == "Pretrained", ]
top_genes <- top_genes[order(-top_genes$selection_frequency, -top_genes$mean_abs_coefficient), ]
top_genes <- head(top_genes, 15)
top_genes$label <- ifelse(is.na(top_genes$gene_symbol) | top_genes$gene_symbol == "", top_genes$ensembl_gene_id, top_genes$gene_symbol)
top_genes$direction <- ifelse(top_genes$positive_fraction_when_selected >= 0.8, "Positive", ifelse(top_genes$positive_fraction_when_selected <= 0.2, "Negative", "Mixed"))
top_genes$label <- factor(top_genes$label, levels = rev(top_genes$label))
p_genes <- ggplot(top_genes, aes(x = label, y = selection_frequency, fill = direction)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  scale_fill_manual(values = c("Negative" = "#C44E52", "Mixed" = "#9E9E9E", "Positive" = "#3B6FB6")) +
  labs(
    title = "Most stable Pretrained Lasso predictors",
    subtitle = "Direction is the sign consistency when the gene was selected",
    x = NULL,
    y = "Selection frequency",
    fill = "Coefficient direction"
  ) +
  base_theme
ggsave(file.path(figure_dir, "primary_pretrained_gene_stability.png"), p_genes, width = 7.0, height = 5.3, dpi = 300)

settings <- read.csv(file.path(results_root, "scenario_settings.csv"))
scenario_labels <- c(
  "primary_1000_ln_ic50_1se" = "Primary: 1,000 genes",
  "features_500_ln_ic50_1se" = "500 genes",
  "features_2000_ln_ic50_1se" = "2,000 genes",
  "auc_1000_1se" = "AUC endpoint",
  "lambda_min_1000_ln_ic50" = "lambda.min"
)
relative_rows <- list()
for (scenario_name in settings$scenario) {
  metrics <- read.csv(file.path(results_root, scenario_name, "repeat_metrics.csv"))
  wide <- reshape(metrics[, c("repeat_id", "method", "mse")], idvar = "repeat_id", timevar = "method", direction = "wide")
  relative_rows[[scenario_name]] <- rbind(
    data.frame(scenario = scenario_name, repeat_id = wide$repeat_id, method = "Overall", percent_mse_reduction = 100 * (wide$mse.Individual - wide$mse.Overall) / wide$mse.Individual),
    data.frame(scenario = scenario_name, repeat_id = wide$repeat_id, method = "Pretrained", percent_mse_reduction = 100 * (wide$mse.Individual - wide$mse.Pretrained) / wide$mse.Individual)
  )
}
relative <- do.call(rbind, relative_rows)
relative$scenario_label <- factor(scenario_labels[relative$scenario], levels = unname(scenario_labels[settings$scenario]))
relative$method <- factor(relative$method, levels = c("Overall", "Pretrained"))
p_sensitivity <- ggplot(relative, aes(x = scenario_label, y = percent_mse_reduction, color = method)) +
  geom_hline(yintercept = 0, color = "grey65", linetype = 2) +
  geom_boxplot(aes(group = interaction(scenario_label, method)), position = position_dodge(width = 0.65), width = 0.54, outlier.shape = NA) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.65), alpha = 0.65, size = 1.3) +
  scale_color_manual(values = method_colors[c("Overall", "Pretrained")]) +
  labs(
    title = "Sensitivity of external-data gains",
    subtitle = "Positive values indicate lower MSE than Individual Lasso",
    x = NULL,
    y = "MSE reduction relative to Individual Lasso",
    color = NULL
  ) +
  base_theme + theme(axis.text.x = element_text(angle = 22, hjust = 1))
ggsave(file.path(figure_dir, "sensitivity_mse_reduction.png"), p_sensitivity, width = 8.0, height = 4.8, dpi = 300)

writeLines(capture.output(sessionInfo()), file.path(figure_dir, "figure_sessionInfo.txt"))
