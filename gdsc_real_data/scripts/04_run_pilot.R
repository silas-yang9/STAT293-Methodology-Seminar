options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))

suppressPackageStartupMessages({
  library(ptLasso)
})

processed_dir <- file.path(root, "processed")
results_dir <- file.path(root, "results")
fold_dir <- file.path(results_dir, "folds")
fit_dir <- file.path(results_dir, "fits")
dir.create(fold_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fit_dir, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(processed_dir, "gdsc_gemcitabine_top1000.csv.gz")
dat <- read.csv(data_path, check.names = FALSE)
metadata_names <- c(
  "model_id", "model_name", "gdsc_cell_line_name", "cancer_group",
  "cancer_type_detail", "role", "DRUG_ID", "DRUG_NAME", "LN_IC50",
  "AUC", "Z_SCORE"
)
gene_names <- setdiff(names(dat), metadata_names)
stopifnot(length(gene_names) == 1000L)
stopifnot(sum(dat$role == "target") == 30L)

x <- as.matrix(dat[, gene_names, drop = FALSE])
storage.mode(x) <- "double"
y <- dat$LN_IC50
groups <- dat$cancer_group
target_rows <- which(dat$role == "target")
source_rows <- which(dat$role == "external")

seed <- 2932026L
n_outer <- 5L
n_inner <- 3L
alpha_grid <- seq(0, 1, by = 0.25)

# Balanced random target folds, fixed before fitting.
set.seed(seed)
target_order <- sample(target_rows)
outer_fold <- integer(nrow(dat))
outer_fold[target_order] <- rep(seq_len(n_outer), length.out = length(target_order))
fold_assignment <- data.frame(
  model_id = dat$model_id[target_rows],
  outer_fold = outer_fold[target_rows]
)
write.csv(fold_assignment, file.path(processed_dir, "outer_fold_assignment.csv"), row.names = FALSE)

make_inner_folds <- function(group_labels, fold_seed) {
  set.seed(fold_seed)
  fold_id <- integer(length(group_labels))
  for (group_name in unique(group_labels)) {
    idx <- which(group_labels == group_name)
    fold_id[idx] <- sample(rep(seq_len(n_inner), length.out = length(idx)))
  }
  fold_id
}

for (fold in seq_len(n_outer)) {
  result_path <- file.path(fold_dir, sprintf("fold_%02d_predictions.csv", fold))
  support_path <- file.path(fold_dir, sprintf("fold_%02d_support.csv", fold))
  diagnostic_path <- file.path(fold_dir, sprintf("fold_%02d_diagnostics.csv", fold))
  fit_path <- file.path(fit_dir, sprintf("fold_%02d_fit.rds", fold))
  if (all(file.exists(c(result_path, support_path, diagnostic_path, fit_path)))) {
    message(sprintf("Fold %d already complete; skipping", fold))
    next
  }

  test_rows <- target_rows[outer_fold[target_rows] == fold]
  train_rows <- c(source_rows, setdiff(target_rows, test_rows))
  inner_fold <- make_inner_folds(groups[train_rows], seed + 100L * fold)
  message(sprintf("Fitting outer fold %d: train=%d, target test=%d", fold, length(train_rows), length(test_rows)))

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

  prediction_frame <- rbind(
    data.frame(model_id = dat$model_id[test_rows], fold = fold, observed = y[test_rows], method = "Individual", predicted = as.numeric(prediction$yhatind)),
    data.frame(model_id = dat$model_id[test_rows], fold = fold, observed = y[test_rows], method = "Overall", predicted = as.numeric(prediction$yhatoverall)),
    data.frame(model_id = dat$model_id[test_rows], fold = fold, observed = y[test_rows], method = "Pretrained", predicted = as.numeric(prediction$yhatpre))
  )
  write.csv(prediction_frame, result_path, row.names = FALSE)

  # Extract support from the fitted coefficients rather than predict()'s
  # support fields. With a single test group, ptLasso 1.0.2 can return support
  # fields associated with the first group in its internal list.
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
  supports <- list(
    Individual = match(individual_genes, gene_names),
    Overall = match(overall_genes, gene_names),
    Pretrained = match(union(overall_genes, target_adjustment_genes), gene_names)
  )
  support_frame <- do.call(rbind, lapply(names(supports), function(method) {
    indices <- supports[[method]]
    if (length(indices) == 0L) {
      return(data.frame(fold = fold, method = method, feature_index = NA_integer_, ensembl_gene_id = NA_character_))
    }
    data.frame(fold = fold, method = method, feature_index = indices, ensembl_gene_id = gene_names[indices])
  }))
  write.csv(support_frame, support_path, row.names = FALSE)

  diagnostics <- data.frame(
    fold = fold,
    selected_alpha = fit$alphahat,
    n_train_total = length(train_rows),
    n_train_target = sum(dat$role[train_rows] == "target"),
    n_train_external = sum(dat$role[train_rows] == "external"),
    n_test_target = length(test_rows),
    features = ncol(x),
    inner_folds = n_inner,
    lasso_choice = "lambda.1se"
  )
  write.csv(diagnostics, diagnostic_path, row.names = FALSE)
  saveRDS(fit, fit_path)
  message(sprintf("Completed outer fold %d; alpha=%.2f", fold, fit$alphahat))
}

prediction_files <- file.path(fold_dir, sprintf("fold_%02d_predictions.csv", seq_len(n_outer)))
if (!all(file.exists(prediction_files))) {
  stop("Not all folds are complete")
}
predictions <- do.call(rbind, lapply(prediction_files, read.csv))
write.csv(predictions, file.path(results_dir, "predictions.csv"), row.names = FALSE)

metric_rows <- do.call(rbind, lapply(split(predictions, predictions$method), function(frame) {
  residual <- frame$observed - frame$predicted
  data.frame(
    method = unique(frame$method),
    n = nrow(frame),
    mse = mean(residual^2),
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    correlation = cor(frame$observed, frame$predicted),
    r_squared = 1 - sum(residual^2) / sum((frame$observed - mean(frame$observed))^2)
  )
}))
metric_rows <- metric_rows[order(metric_rows$mse), ]
write.csv(metric_rows, file.path(results_dir, "metrics.csv"), row.names = FALSE)

support_files <- file.path(fold_dir, sprintf("fold_%02d_support.csv", seq_len(n_outer)))
supports <- do.call(rbind, lapply(support_files, read.csv))
supports <- supports[!is.na(supports$ensembl_gene_id), ]
stability <- aggregate(
  fold ~ method + ensembl_gene_id,
  data = supports,
  FUN = function(values) length(unique(values)) / n_outer
)
names(stability)[names(stability) == "fold"] <- "selection_frequency"
stability <- stability[order(stability$method, -stability$selection_frequency, stability$ensembl_gene_id), ]
write.csv(stability, file.path(results_dir, "selection_stability.csv"), row.names = FALSE)

diagnostic_files <- file.path(fold_dir, sprintf("fold_%02d_diagnostics.csv", seq_len(n_outer)))
diagnostics <- do.call(rbind, lapply(diagnostic_files, read.csv))
write.csv(diagnostics, file.path(results_dir, "fold_diagnostics.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(results_dir, "sessionInfo.txt"))
writeLines(
  c(
    "GDSC gemcitabine real-data pilot",
    sprintf("Completed: %s", format(Sys.time(), tz = "Asia/Shanghai")),
    sprintf("Seed: %d", seed),
    sprintf("Outer folds: %d", n_outer),
    sprintf("Inner folds: %d", n_inner),
    sprintf("Features: %d (selected using external expression only)", ncol(x)),
    sprintf("Alpha grid: %s", paste(alpha_grid, collapse = ", ")),
    "Lasso choice: lambda.1se",
    "",
    paste(capture.output(print(metric_rows, row.names = FALSE)), collapse = "\n")
  ),
  file.path(results_dir, "run_summary.txt")
)
print(metric_rows, row.names = FALSE)
