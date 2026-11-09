options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ptLasso))

processed_dir <- file.path(root, "processed")
results_root <- file.path(root, "results", "conservative_alpha_experiment")
checkpoint_dir <- file.path(results_root, "checkpoints")
assignment_dir <- file.path(results_root, "assignments")
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(assignment_dir, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(file.path(processed_dir, "gdsc_gemcitabine_top2000.csv.gz"), check.names = FALSE)
metadata_names <- c(
  "model_id", "model_name", "gdsc_cell_line_name", "cancer_group",
  "cancer_type_detail", "role", "DRUG_ID", "DRUG_NAME", "LN_IC50",
  "AUC", "Z_SCORE"
)
gene_names <- setdiff(names(dat), metadata_names)
stopifnot(length(gene_names) == 2000L)
stopifnot(sum(dat$role == "target") == 30L)
x <- as.matrix(dat[, gene_names, drop = FALSE])
storage.mode(x) <- "double"
y <- dat$LN_IC50
groups <- dat$cancer_group
target_label <- "Pancreatic Carcinoma"
target_rows <- which(dat$role == "target")
source_rows <- which(dat$role == "external")

n_repeats <- as.integer(Sys.getenv("CONSERVATIVE_REPEATS", "20"))
if (!is.finite(n_repeats) || n_repeats < 1L || n_repeats > 20L) {
  stop("CONSERVATIVE_REPEATS must be an integer from 1 to 20.")
}
base_seed <- 902293L
n_outer <- 5L
n_inner <- 3L
alpha_grid <- seq(0, 1, by = 0.1)
lambda_choice <- "lambda.1se"

get_lamhat_internal <- get("get.lamhat", envir = asNamespace("ptLasso"))

make_inner_folds <- function(group_labels, fold_seed) {
  set.seed(fold_seed)
  fold_id <- integer(length(group_labels))
  for (group_name in unique(group_labels)) {
    idx <- which(group_labels == group_name)
    fold_id[idx] <- sample(rep(seq_len(n_inner), length.out = length(idx)))
  }
  fold_id
}

gene_vector <- function(coefficient_matrix) {
  output <- rep(0, length(gene_names))
  names(output) <- gene_names
  rows <- intersect(rownames(coefficient_matrix), gene_names)
  output[rows] <- as.numeric(coefficient_matrix[rows, 1])
  output
}

pretrained_coefficients <- function(fit, target_index, alpha_value) {
  overall <- gene_vector(coef(fit, model = "overall", s = lambda_choice))
  adjustment <- gene_vector(
    coef(fit, model = "pretrain", alpha = alpha_value, s = lambda_choice)[[target_index]]
  )
  (1 - alpha_value) * overall + adjustment
}

metric_function <- function(frame) {
  residual <- frame$observed - frame$predicted
  prediction_sd <- sd(frame$predicted)
  data.frame(
    mse = mean(residual^2),
    rmse = sqrt(mean(residual^2)),
    mae = mean(abs(residual)),
    correlation = if (is.na(prediction_sd) || prediction_sd == 0) NA_real_ else cor(frame$observed, frame$predicted),
    r_squared = 1 - sum(residual^2) / sum((frame$observed - mean(frame$observed))^2)
  )
}

for (repeat_id in seq_len(n_repeats)) {
  repeat_seed <- base_seed + (repeat_id - 1L) * 10000L
  set.seed(repeat_seed)
  target_order <- sample(target_rows)
  outer_fold <- integer(nrow(dat))
  outer_fold[target_order] <- rep(seq_len(n_outer), length.out = length(target_order))
  assignment_path <- file.path(assignment_dir, sprintf("repeat_%02d_assignment.csv", repeat_id))
  if (!file.exists(assignment_path)) {
    write.csv(
      data.frame(repeat_id = repeat_id, model_id = dat$model_id[target_rows], outer_fold = outer_fold[target_rows]),
      assignment_path,
      row.names = FALSE
    )
  }

  for (fold in seq_len(n_outer)) {
    prefix <- sprintf("repeat_%02d_fold_%02d", repeat_id, fold)
    output_paths <- c(
      predictions = file.path(checkpoint_dir, paste0(prefix, "_predictions.csv")),
      diagnostics = file.path(checkpoint_dir, paste0(prefix, "_diagnostics.csv")),
      alpha_errors = file.path(checkpoint_dir, paste0(prefix, "_alpha_errors.csv")),
      coefficients = file.path(checkpoint_dir, paste0(prefix, "_coefficients.csv"))
    )
    if (all(file.exists(output_paths))) {
      message(sprintf("repeat %d fold %d complete; skipping", repeat_id, fold))
      next
    }

    test_rows <- target_rows[outer_fold[target_rows] == fold]
    target_train_rows <- setdiff(target_rows, test_rows)
    train_rows <- c(source_rows, target_train_rows)
    inner_seed <- repeat_seed + 100L * fold
    inner_fold <- make_inner_folds(groups[train_rows], inner_seed)
    message(sprintf("Conservative alpha: repeat %d/%d, fold %d/%d", repeat_id, n_repeats, fold, n_outer))

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
      s = lambda_choice,
      overall.lambda = lambda_choice,
      alphahat.choice = "overall",
      group.intercepts = TRUE
    )

    legend <- fit$fit[[1]]$group.legend
    target_index <- match(target_label, legend)
    stopifnot(!is.na(target_index))
    target_column <- paste0("group_", target_label)
    target_errors <- as.numeric(fit$errpre[, target_column])
    target_se <- vapply(seq_along(alpha_grid), function(alpha_index) {
      target_cv <- fit$fit[[alpha_index]]$fitpre[[target_index]]
      selected_lambda <- get_lamhat_internal(target_cv, lambda_choice)
      lambda_index <- which.min(abs(target_cv$lambda - selected_lambda))
      if (is.null(target_cv$cvsd)) NA_real_ else as.numeric(target_cv$cvsd[lambda_index])
    }, numeric(1))
    stopifnot(all(is.finite(target_se)))

    alpha_global <- fit$alphahat
    alpha_target_min <- unname(fit$varying.alphahat[target_index])
    best_target_index <- which.min(target_errors)
    one_se_threshold <- target_errors[best_target_index] + target_se[best_target_index]
    eligible_alpha <- which(target_errors <= one_se_threshold + 1e-12)
    alpha_target_1se <- min(alpha_grid[eligible_alpha])

    alpha_error_table <- data.frame(
      repeat_id = repeat_id,
      fold = fold,
      alpha = alpha_grid,
      target_inner_mse = target_errors,
      target_inner_se = target_se,
      eligible_under_1se = alpha_grid %in% alpha_grid[eligible_alpha]
    )
    write.csv(alpha_error_table, output_paths[["alpha_errors"]], row.names = FALSE)

    alpha_rules <- c(
      `Pretrained global-alpha` = alpha_global,
      `Pretrained target-min` = alpha_target_min,
      `Pretrained target-1se` = alpha_target_1se
    )
    prediction_list <- lapply(alpha_rules, function(alpha_value) {
      predict(
        fit,
        xtest = x[test_rows, , drop = FALSE],
        groupstest = groups[test_rows],
        ytest = y[test_rows],
        alpha = alpha_value,
        s = lambda_choice,
        return.link = TRUE
      )
    })
    reference_prediction <- prediction_list[[1]]
    predictions <- do.call(rbind, c(
      lapply(names(alpha_rules), function(method) {
        data.frame(
          repeat_id = repeat_id,
          fold = fold,
          model_id = dat$model_id[test_rows],
          observed = y[test_rows],
          method = method,
          predicted = as.numeric(prediction_list[[method]]$yhatpre),
          selected_alpha = alpha_rules[[method]]
        )
      }),
      list(
        data.frame(
          repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows],
          observed = y[test_rows], method = "Overall",
          predicted = as.numeric(reference_prediction$yhatoverall), selected_alpha = NA_real_
        ),
        data.frame(
          repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows],
          observed = y[test_rows], method = "Individual",
          predicted = as.numeric(reference_prediction$yhatind), selected_alpha = NA_real_
        ),
        data.frame(
          repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows],
          observed = y[test_rows], method = "Target mean",
          predicted = mean(y[target_train_rows]), selected_alpha = NA_real_
        )
      )
    ))
    write.csv(predictions, output_paths[["predictions"]], row.names = FALSE)

    coefficient_vectors <- list(
      `Pretrained global-alpha` = pretrained_coefficients(fit, target_index, alpha_global),
      `Pretrained target-min` = pretrained_coefficients(fit, target_index, alpha_target_min),
      `Pretrained target-1se` = pretrained_coefficients(fit, target_index, alpha_target_1se),
      Overall = gene_vector(coef(fit, model = "overall", s = lambda_choice)),
      Individual = gene_vector(coef(fit, model = "individual", s = lambda_choice)[[target_index]])
    )
    coefficients <- do.call(rbind, lapply(names(coefficient_vectors), function(method) {
      values <- coefficient_vectors[[method]]
      selected <- which(abs(values) > 1e-12)
      if (length(selected) == 0L) {
        return(data.frame(
          repeat_id = repeat_id, fold = fold, method = method,
          ensembl_gene_id = NA_character_, coefficient = NA_real_
        ))
      }
      data.frame(
        repeat_id = repeat_id,
        fold = fold,
        method = method,
        ensembl_gene_id = names(values)[selected],
        coefficient = values[selected]
      )
    }))
    write.csv(coefficients, output_paths[["coefficients"]], row.names = FALSE)

    diagnostics <- data.frame(
      repeat_id = repeat_id,
      fold = fold,
      repeat_seed = repeat_seed,
      inner_seed = inner_seed,
      alpha_global = alpha_global,
      alpha_target_min = alpha_target_min,
      alpha_target_1se = alpha_target_1se,
      target_min_inner_mse = target_errors[best_target_index],
      target_min_inner_se = target_se[best_target_index],
      one_se_threshold = one_se_threshold,
      n_train_total = length(train_rows),
      n_train_target = length(target_train_rows),
      n_train_external = length(source_rows),
      n_test_target = length(test_rows),
      features = length(gene_names),
      lambda_choice = lambda_choice
    )
    write.csv(diagnostics, output_paths[["diagnostics"]], row.names = FALSE)
    rm(fit, prediction_list, coefficient_vectors)
    invisible(gc())
  }
}

expected_fits <- n_repeats * n_outer
file_sets <- list(
  predictions = list.files(checkpoint_dir, pattern = "_predictions\\.csv$", full.names = TRUE),
  diagnostics = list.files(checkpoint_dir, pattern = "_diagnostics\\.csv$", full.names = TRUE),
  alpha_errors = list.files(checkpoint_dir, pattern = "_alpha_errors\\.csv$", full.names = TRUE),
  coefficients = list.files(checkpoint_dir, pattern = "_coefficients\\.csv$", full.names = TRUE)
)
stopifnot(all(vapply(file_sets, length, integer(1)) >= expected_fits))
repeat_ids_to_use <- seq_len(n_repeats)
read_selected_files <- function(paths) {
  frames <- lapply(paths, read.csv)
  frame <- do.call(rbind, frames)
  frame[frame$repeat_id %in% repeat_ids_to_use, , drop = FALSE]
}
predictions <- read_selected_files(file_sets$predictions)
diagnostics <- read_selected_files(file_sets$diagnostics)
alpha_errors <- read_selected_files(file_sets$alpha_errors)
coefficients <- read_selected_files(file_sets$coefficients)
write.csv(predictions, file.path(results_root, "predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(results_root, "diagnostics.csv"), row.names = FALSE)
write.csv(alpha_errors, file.path(results_root, "alpha_errors.csv"), row.names = FALSE)

repeat_metrics <- do.call(rbind, lapply(split(predictions, list(predictions$repeat_id, predictions$method), drop = TRUE), function(frame) {
  data.frame(repeat_id = unique(frame$repeat_id), method = unique(frame$method), metric_function(frame))
}))
repeat_metrics <- repeat_metrics[order(repeat_metrics$repeat_id, repeat_metrics$method), ]
write.csv(repeat_metrics, file.path(results_root, "repeat_metrics.csv"), row.names = FALSE)

method_summary <- do.call(rbind, lapply(split(repeat_metrics, repeat_metrics$method), function(frame) {
  data.frame(
    method = unique(frame$method),
    repeats = nrow(frame),
    mean_mse = mean(frame$mse),
    sd_mse = sd(frame$mse),
    median_mse = median(frame$mse),
    q025_mse = as.numeric(quantile(frame$mse, 0.025)),
    q975_mse = as.numeric(quantile(frame$mse, 0.975)),
    mean_mae = mean(frame$mae),
    mean_correlation = mean(frame$correlation, na.rm = TRUE),
    mean_r_squared = mean(frame$r_squared)
  )
}))
method_summary <- method_summary[order(method_summary$mean_mse), ]
write.csv(method_summary, file.path(results_root, "method_summary.csv"), row.names = FALSE)

metric_wide <- reshape(repeat_metrics[, c("repeat_id", "method", "mse")], idvar = "repeat_id", timevar = "method", direction = "wide")
paired <- data.frame(
  repeat_id = metric_wide$repeat_id,
  individual_minus_target_1se = metric_wide$mse.Individual - metric_wide$`mse.Pretrained target-1se`,
  overall_minus_target_1se = metric_wide$mse.Overall - metric_wide$`mse.Pretrained target-1se`,
  global_minus_target_1se = metric_wide$`mse.Pretrained global-alpha` - metric_wide$`mse.Pretrained target-1se`,
  target_min_minus_target_1se = metric_wide$`mse.Pretrained target-min` - metric_wide$`mse.Pretrained target-1se`
)
write.csv(paired, file.path(results_root, "paired_differences.csv"), row.names = FALSE)

selection_frequency <- do.call(rbind, lapply(c("alpha_global", "alpha_target_min", "alpha_target_1se"), function(column) {
  frame <- as.data.frame(table(diagnostics[[column]]))
  names(frame) <- c("alpha", "selected_fits")
  frame$rule <- column
  frame$selection_frequency <- frame$selected_fits / expected_fits
  frame
}))
write.csv(selection_frequency, file.path(results_root, "alpha_selection_frequency.csv"), row.names = FALSE)

selected_coefficients <- coefficients[!is.na(coefficients$ensembl_gene_id), ]
coefficient_summary <- do.call(rbind, lapply(split(selected_coefficients, list(selected_coefficients$method, selected_coefficients$ensembl_gene_id), drop = TRUE), function(frame) {
  data.frame(
    method = unique(frame$method),
    ensembl_gene_id = unique(frame$ensembl_gene_id),
    selected_fits = nrow(frame),
    selection_frequency = nrow(frame) / expected_fits,
    median_coefficient = median(frame$coefficient),
    mean_abs_coefficient = mean(abs(frame$coefficient)),
    positive_fraction_when_selected = mean(frame$coefficient > 0)
  )
}))
coefficient_summary <- coefficient_summary[order(coefficient_summary$method, -coefficient_summary$selection_frequency), ]
write.csv(coefficient_summary, file.path(results_root, "coefficient_stability.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(results_root, "sessionInfo.txt"))
writeLines(
  c(
    "Conservative alpha-selection experiment",
    sprintf("Completed: %s", format(Sys.time(), tz = "Asia/Shanghai")),
    sprintf("Repeats: %d", n_repeats),
    "Features: 2000",
    "First-stage and second-stage lambda: lambda.1se",
    "Target-1se rule: smallest alpha whose pancreatic inner-CV error is within one SE of the target minimum.",
    "",
    paste(capture.output(print(method_summary, row.names = FALSE)), collapse = "\n")
  ),
  file.path(results_root, "run_summary.txt")
)
print(method_summary, row.names = FALSE)

