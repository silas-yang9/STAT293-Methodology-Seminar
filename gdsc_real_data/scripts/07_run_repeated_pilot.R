options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ptLasso))

processed_dir <- file.path(root, "processed")
results_dir <- file.path(root, "results", "repeated_cv")
checkpoint_dir <- file.path(results_dir, "checkpoints")
assignment_dir <- file.path(results_dir, "assignments")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(assignment_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(file.path(processed_dir, "gdsc_gemcitabine_top1000.csv.gz"), check.names = FALSE)
metadata_names <- c(
  "model_id", "model_name", "gdsc_cell_line_name", "cancer_group",
  "cancer_type_detail", "role", "DRUG_ID", "DRUG_NAME", "LN_IC50",
  "AUC", "Z_SCORE"
)
gene_names <- setdiff(names(dat), metadata_names)
x <- as.matrix(dat[, gene_names, drop = FALSE])
storage.mode(x) <- "double"
y <- dat$LN_IC50
groups <- dat$cancer_group
target_rows <- which(dat$role == "target")
source_rows <- which(dat$role == "external")

base_seed <- 2932026L
n_repeats <- 10L
n_outer <- 5L
n_inner <- 3L
alpha_grid <- seq(0, 1, by = 0.25)

make_inner_folds <- function(group_labels, fold_seed) {
  set.seed(fold_seed)
  fold_id <- integer(length(group_labels))
  for (group_name in unique(group_labels)) {
    idx <- which(group_labels == group_name)
    fold_id[idx] <- sample(rep(seq_len(n_inner), length.out = length(idx)))
  }
  fold_id
}

extract_support <- function(fit) {
  coefficient_list <- coef(fit, alpha = fit$alphahat, s = "lambda.1se")
  group_legend <- fit$fit[[1]]$group.legend
  target_index <- match("Pancreatic Carcinoma", group_legend)
  nonzero_genes <- function(coefficient_matrix) {
    coefficient_names <- rownames(coefficient_matrix)[as.vector(coefficient_matrix) != 0]
    intersect(coefficient_names, gene_names)
  }
  overall_genes <- nonzero_genes(coefficient_list$overall)
  individual_genes <- nonzero_genes(coefficient_list$individual[[target_index]])
  target_adjustment_genes <- nonzero_genes(coefficient_list$pretrain[[target_index]])
  list(
    Individual = individual_genes,
    Overall = overall_genes,
    Pretrained = union(overall_genes, target_adjustment_genes)
  )
}

for (repeat_id in seq_len(n_repeats)) {
  repeat_seed <- base_seed + (repeat_id - 1L) * 10000L
  set.seed(repeat_seed)
  target_order <- sample(target_rows)
  outer_fold <- integer(nrow(dat))
  outer_fold[target_order] <- rep(seq_len(n_outer), length.out = length(target_order))
  write.csv(
    data.frame(
      repeat_id = repeat_id,
      model_id = dat$model_id[target_rows],
      outer_fold = outer_fold[target_rows]
    ),
    file.path(assignment_dir, sprintf("repeat_%02d_assignment.csv", repeat_id)),
    row.names = FALSE
  )

  for (fold in seq_len(n_outer)) {
    prefix <- sprintf("repeat_%02d_fold_%02d", repeat_id, fold)
    prediction_path <- file.path(checkpoint_dir, paste0(prefix, "_predictions.csv"))
    support_path <- file.path(checkpoint_dir, paste0(prefix, "_support.csv"))
    diagnostic_path <- file.path(checkpoint_dir, paste0(prefix, "_diagnostics.csv"))
    if (all(file.exists(c(prediction_path, support_path, diagnostic_path)))) {
      message(sprintf("Repeat %d fold %d already complete; skipping", repeat_id, fold))
      next
    }
    test_rows <- target_rows[outer_fold[target_rows] == fold]
    target_train_rows <- setdiff(target_rows, test_rows)
    train_rows <- c(source_rows, target_train_rows)
    inner_seed <- repeat_seed + 100L * fold
    inner_fold <- make_inner_folds(groups[train_rows], inner_seed)
    message(sprintf("Fitting repeat %d/%d, fold %d/%d", repeat_id, n_repeats, fold, n_outer))

    fit <- cv.ptLasso(
      x = x[train_rows, , drop = FALSE],
      y = y[train_rows],
      groups = groups[train_rows],
      alphalist = alpha_grid,
      family = "gaussian",
      use.case = "inputGroups",
      type.measure = "mse",
      nfolds = n_inner,
      foldid = inner_fold,
      s = "lambda.1se",
      overall.lambda = "lambda.1se",
      alphahat.choice = "overall",
      group.intercepts = TRUE
    )
    prediction <- predict(
      fit,
      xtest = x[test_rows, , drop = FALSE],
      groupstest = groups[test_rows],
      ytest = y[test_rows],
      alphatype = "fixed",
      s = "lambda.1se",
      return.link = TRUE
    )
    target_mean <- mean(y[target_train_rows])
    prediction_frame <- rbind(
      data.frame(repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Individual", predicted = as.numeric(prediction$yhatind)),
      data.frame(repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Overall", predicted = as.numeric(prediction$yhatoverall)),
      data.frame(repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Pretrained", predicted = as.numeric(prediction$yhatpre)),
      data.frame(repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Target mean", predicted = target_mean)
    )
    write.csv(prediction_frame, prediction_path, row.names = FALSE)

    support_genes <- extract_support(fit)
    support_frame <- do.call(rbind, lapply(names(support_genes), function(method) {
      genes <- support_genes[[method]]
      if (length(genes) == 0L) {
        return(data.frame(repeat_id = repeat_id, fold = fold, method = method, ensembl_gene_id = NA_character_))
      }
      data.frame(repeat_id = repeat_id, fold = fold, method = method, ensembl_gene_id = genes)
    }))
    write.csv(support_frame, support_path, row.names = FALSE)
    write.csv(
      data.frame(
        repeat_id = repeat_id,
        fold = fold,
        repeat_seed = repeat_seed,
        inner_seed = inner_seed,
        selected_alpha = fit$alphahat,
        n_train_total = length(train_rows),
        n_train_target = length(target_train_rows),
        n_train_external = length(source_rows),
        n_test_target = length(test_rows),
        features = ncol(x),
        inner_folds = n_inner,
        lasso_choice = "lambda.1se"
      ),
      diagnostic_path,
      row.names = FALSE
    )
  }
}

prediction_files <- list.files(checkpoint_dir, pattern = "_predictions\\.csv$", full.names = TRUE)
support_files <- list.files(checkpoint_dir, pattern = "_support\\.csv$", full.names = TRUE)
diagnostic_files <- list.files(checkpoint_dir, pattern = "_diagnostics\\.csv$", full.names = TRUE)
stopifnot(length(prediction_files) == n_repeats * n_outer)
predictions <- do.call(rbind, lapply(prediction_files, read.csv))
supports <- do.call(rbind, lapply(support_files, read.csv))
diagnostics <- do.call(rbind, lapply(diagnostic_files, read.csv))
write.csv(predictions, file.path(results_dir, "predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(results_dir, "diagnostics.csv"), row.names = FALSE)

metric_function <- function(frame) {
  residual <- frame$observed - frame$predicted
  data.frame(
    mse = mean(residual^2),
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    correlation = if (sd(frame$predicted) == 0) NA_real_ else cor(frame$observed, frame$predicted),
    r_squared = 1 - sum(residual^2) / sum((frame$observed - mean(frame$observed))^2)
  )
}
repeat_metrics <- do.call(rbind, lapply(split(predictions, list(predictions$repeat_id, predictions$method), drop = TRUE), function(frame) {
  data.frame(repeat_id = unique(frame$repeat_id), method = unique(frame$method), metric_function(frame))
}))
repeat_metrics <- repeat_metrics[order(repeat_metrics$repeat_id, repeat_metrics$method), ]
write.csv(repeat_metrics, file.path(results_dir, "repeat_metrics.csv"), row.names = FALSE)

metric_summary <- do.call(rbind, lapply(split(repeat_metrics, repeat_metrics$method), function(frame) {
  data.frame(
    method = unique(frame$method),
    repeats = nrow(frame),
    mean_mse = mean(frame$mse),
    sd_mse = sd(frame$mse),
    min_mse = min(frame$mse),
    max_mse = max(frame$mse),
    mean_mae = mean(frame$mae),
    mean_correlation = mean(frame$correlation, na.rm = TRUE),
    mean_r_squared = mean(frame$r_squared)
  )
}))
metric_summary <- metric_summary[order(metric_summary$mean_mse), ]
write.csv(metric_summary, file.path(results_dir, "metric_summary.csv"), row.names = FALSE)

metric_wide <- reshape(repeat_metrics[, c("repeat_id", "method", "mse")], idvar = "repeat_id", timevar = "method", direction = "wide")
paired_differences <- data.frame(
  repeat_id = metric_wide$repeat_id,
  individual_minus_pretrained = metric_wide$mse.Individual - metric_wide$mse.Pretrained,
  overall_minus_pretrained = metric_wide$mse.Overall - metric_wide$mse.Pretrained,
  target_mean_minus_pretrained = metric_wide$`mse.Target mean` - metric_wide$mse.Pretrained,
  individual_minus_overall = metric_wide$mse.Individual - metric_wide$mse.Overall
)
write.csv(paired_differences, file.path(results_dir, "paired_repeat_differences.csv"), row.names = FALSE)

alpha_distribution <- as.data.frame(table(diagnostics$selected_alpha))
names(alpha_distribution) <- c("alpha", "fold_count")
alpha_distribution$proportion <- alpha_distribution$fold_count / sum(alpha_distribution$fold_count)
write.csv(alpha_distribution, file.path(results_dir, "alpha_distribution.csv"), row.names = FALSE)

supports <- supports[!is.na(supports$ensembl_gene_id), ]
if (nrow(supports) > 0L) {
  support_stability <- aggregate(
    fold ~ method + ensembl_gene_id,
    data = supports,
    FUN = length
  )
  names(support_stability)[names(support_stability) == "fold"] <- "selected_fits"
  support_stability$selection_frequency <- support_stability$selected_fits / (n_repeats * n_outer)
  support_stability <- support_stability[order(support_stability$method, -support_stability$selection_frequency), ]
} else {
  support_stability <- data.frame(method = character(), ensembl_gene_id = character(), selected_fits = integer(), selection_frequency = numeric())
}
write.csv(support_stability, file.path(results_dir, "selection_stability.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(results_dir, "sessionInfo.txt"))
writeLines(
  c(
    "Repeated GDSC gemcitabine real-data pilot",
    sprintf("Completed: %s", format(Sys.time(), tz = "Asia/Shanghai")),
    sprintf("Base seed: %d", base_seed),
    sprintf("Repeats x outer folds: %d x %d", n_repeats, n_outer),
    sprintf("Inner folds: %d", n_inner),
    sprintf("Features: %d", ncol(x)),
    sprintf("Alpha grid: %s", paste(alpha_grid, collapse = ", ")),
    "Lasso choice: lambda.1se",
    "",
    paste(capture.output(print(metric_summary, row.names = FALSE)), collapse = "\n"),
    "",
    "Mean paired MSE differences:",
    paste(capture.output(print(colMeans(paired_differences[, -1]), row.names = FALSE)), collapse = "\n")
  ),
  file.path(results_dir, "run_summary.txt")
)
print(metric_summary, row.names = FALSE)
print(colMeans(paired_differences[, -1]))
