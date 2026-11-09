options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ptLasso))

processed_dir <- file.path(root, "processed")
results_root <- file.path(root, "results", "nested_tuning_experiment")
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
all_gene_names <- setdiff(names(dat), metadata_names)
stopifnot(length(all_gene_names) == 2000L)
stopifnot(sum(dat$role == "target") == 30L)

candidate_grid <- data.frame(
  candidate_id = c("f1000_1se", "f2000_1se", "f1000_min", "f2000_min"),
  features = c(1000L, 2000L, 1000L, 2000L),
  second_stage_lambda = c("lambda.1se", "lambda.1se", "lambda.min", "lambda.min")
)
candidate_grid$first_stage_lambda <- "lambda.1se"
write.csv(candidate_grid, file.path(results_root, "candidate_grid.csv"), row.names = FALSE)

n_repeats <- as.integer(Sys.getenv("NESTED_REPEATS", "20"))
if (!is.finite(n_repeats) || n_repeats < 1L || n_repeats > 20L) {
  stop("NESTED_REPEATS must be an integer from 1 to 20.")
}
base_seed <- 731293L
n_outer <- 5L
n_inner <- 3L
alpha_grid <- seq(0, 1, by = 0.1)
groups <- dat$cancer_group
target_label <- "Pancreatic Carcinoma"
target_rows <- which(dat$role == "target")
source_rows <- which(dat$role == "external")
y <- dat$LN_IC50
x_all <- as.matrix(dat[, all_gene_names, drop = FALSE])
storage.mode(x_all) <- "double"

make_inner_folds <- function(group_labels, fold_seed) {
  set.seed(fold_seed)
  fold_id <- integer(length(group_labels))
  for (group_name in unique(group_labels)) {
    idx <- which(group_labels == group_name)
    fold_id[idx] <- sample(rep(seq_len(n_inner), length.out = length(idx)))
  }
  fold_id
}

gene_vector <- function(coefficient_matrix, gene_names) {
  output <- rep(0, length(gene_names))
  names(output) <- gene_names
  rows <- intersect(rownames(coefficient_matrix), gene_names)
  output[rows] <- as.numeric(coefficient_matrix[rows, 1])
  output
}

extract_candidate_coefficients <- function(fit, gene_names, target_index, target_alpha, second_stage_lambda) {
  overall_prediction <- gene_vector(
    coef(fit, model = "overall", s = second_stage_lambda),
    gene_names
  )
  overall_offset <- gene_vector(
    coef(fit, model = "overall", s = "lambda.1se"),
    gene_names
  )
  individual <- gene_vector(
    coef(fit, model = "individual", s = second_stage_lambda)[[target_index]],
    gene_names
  )
  target_adjustment <- gene_vector(
    coef(fit, model = "pretrain", alpha = target_alpha, s = second_stage_lambda)[[target_index]],
    gene_names
  )
  pretrained <- (1 - target_alpha) * overall_offset + target_adjustment
  list(Overall = overall_prediction, Individual = individual, Pretrained = pretrained)
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
      coefficients = file.path(checkpoint_dir, paste0(prefix, "_coefficients.csv")),
      diagnostics = file.path(checkpoint_dir, paste0(prefix, "_diagnostics.csv")),
      candidates = file.path(checkpoint_dir, paste0(prefix, "_candidates.csv"))
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
    message(sprintf("Nested tuning: repeat %d/%d, fold %d/%d", repeat_id, n_repeats, fold, n_outer))

    candidate_records <- vector("list", nrow(candidate_grid))
    candidate_predictions <- vector("list", nrow(candidate_grid))
    candidate_coefficients <- vector("list", nrow(candidate_grid))

    for (candidate_index in seq_len(nrow(candidate_grid))) {
      candidate <- candidate_grid[candidate_index, ]
      gene_names <- all_gene_names[seq_len(candidate$features)]
      x <- x_all[, gene_names, drop = FALSE]
      message(sprintf("  fitting %s", candidate$candidate_id))

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
        s = candidate$second_stage_lambda,
        overall.lambda = candidate$first_stage_lambda,
        alphahat.choice = "overall",
        group.intercepts = TRUE
      )

      legend <- fit$fit[[1]]$group.legend
      target_index <- match(target_label, legend)
      stopifnot(!is.na(target_index))
      target_column <- paste0("group_", target_label)
      stopifnot(target_column %in% colnames(fit$errpre))
      target_alpha <- unname(fit$varying.alphahat[target_index])
      alpha_row <- which(abs(fit$errpre[, "alpha"] - target_alpha) < 1e-8)[1]

      prediction <- predict(
        fit,
        xtest = x[test_rows, , drop = FALSE],
        groupstest = groups[test_rows],
        ytest = y[test_rows],
        alpha = target_alpha,
        s = candidate$second_stage_lambda,
        return.link = TRUE
      )

      candidate_records[[candidate_index]] <- data.frame(
        repeat_id = repeat_id,
        fold = fold,
        candidate_id = candidate$candidate_id,
        features = candidate$features,
        first_stage_lambda = candidate$first_stage_lambda,
        second_stage_lambda = candidate$second_stage_lambda,
        target_alpha = target_alpha,
        pretrained_inner_mse = fit$errpre[alpha_row, target_column],
        individual_inner_mse = unname(fit$errind[target_column]),
        overall_inner_mse = unname(fit$erroverall[target_column])
      )
      candidate_predictions[[candidate_index]] <- list(
        Pretrained = as.numeric(prediction$yhatpre),
        Individual = as.numeric(prediction$yhatind),
        Overall = as.numeric(prediction$yhatoverall)
      )
      candidate_coefficients[[candidate_index]] <- extract_candidate_coefficients(
        fit,
        gene_names,
        target_index,
        target_alpha,
        candidate$second_stage_lambda
      )
      rm(fit, x)
      invisible(gc())
    }

    candidate_table <- do.call(rbind, candidate_records)
    write.csv(candidate_table, output_paths[["candidates"]], row.names = FALSE)

    selection <- c(
      Pretrained = which.min(candidate_table$pretrained_inner_mse),
      Individual = which.min(candidate_table$individual_inner_mse),
      Overall = which.min(candidate_table$overall_inner_mse)
    )
    predictions <- do.call(rbind, lapply(names(selection), function(method) {
      selected_index <- selection[[method]]
      selected <- candidate_table[selected_index, ]
      data.frame(
        repeat_id = repeat_id,
        fold = fold,
        model_id = dat$model_id[test_rows],
        observed = y[test_rows],
        method = method,
        predicted = candidate_predictions[[selected_index]][[method]],
        selected_candidate = selected$candidate_id,
        selected_features = selected$features,
        selected_second_stage_lambda = selected$second_stage_lambda,
        selected_target_alpha = if (method == "Pretrained") selected$target_alpha else NA_real_
      )
    }))
    predictions <- rbind(
      predictions,
      data.frame(
        repeat_id = repeat_id,
        fold = fold,
        model_id = dat$model_id[test_rows],
        observed = y[test_rows],
        method = "Target mean",
        predicted = mean(y[target_train_rows]),
        selected_candidate = "target_mean",
        selected_features = NA_integer_,
        selected_second_stage_lambda = NA_character_,
        selected_target_alpha = NA_real_
      )
    )
    write.csv(predictions, output_paths[["predictions"]], row.names = FALSE)

    coefficients <- do.call(rbind, lapply(names(selection), function(method) {
      selected_index <- selection[[method]]
      values <- candidate_coefficients[[selected_index]][[method]]
      selected_genes <- which(abs(values) > 1e-12)
      if (length(selected_genes) == 0L) {
        return(data.frame(
          repeat_id = repeat_id, fold = fold, method = method,
          selected_candidate = candidate_table$candidate_id[selected_index],
          ensembl_gene_id = NA_character_, coefficient = NA_real_
        ))
      }
      data.frame(
        repeat_id = repeat_id,
        fold = fold,
        method = method,
        selected_candidate = candidate_table$candidate_id[selected_index],
        ensembl_gene_id = names(values)[selected_genes],
        coefficient = values[selected_genes]
      )
    }))
    write.csv(coefficients, output_paths[["coefficients"]], row.names = FALSE)

    diagnostics <- do.call(rbind, lapply(names(selection), function(method) {
      selected_index <- selection[[method]]
      selected <- candidate_table[selected_index, ]
      metric_name <- paste0(tolower(method), "_inner_mse")
      data.frame(
        repeat_id = repeat_id,
        fold = fold,
        method = method,
        selected_candidate = selected$candidate_id,
        selected_features = selected$features,
        selected_second_stage_lambda = selected$second_stage_lambda,
        selected_target_alpha = if (method == "Pretrained") selected$target_alpha else NA_real_,
        selected_inner_mse = selected[[metric_name]],
        repeat_seed = repeat_seed,
        inner_seed = inner_seed,
        n_train_total = length(train_rows),
        n_train_target = length(target_train_rows),
        n_train_external = length(source_rows),
        n_test_target = length(test_rows)
      )
    }))
    write.csv(diagnostics, output_paths[["diagnostics"]], row.names = FALSE)
  }
}

expected_fits <- n_repeats * n_outer
file_sets <- list(
  predictions = list.files(checkpoint_dir, pattern = "_predictions\\.csv$", full.names = TRUE),
  coefficients = list.files(checkpoint_dir, pattern = "_coefficients\\.csv$", full.names = TRUE),
  diagnostics = list.files(checkpoint_dir, pattern = "_diagnostics\\.csv$", full.names = TRUE),
  candidates = list.files(checkpoint_dir, pattern = "_candidates\\.csv$", full.names = TRUE)
)
stopifnot(all(vapply(file_sets, length, integer(1)) >= expected_fits))
repeat_ids_to_use <- seq_len(n_repeats)
read_selected_files <- function(paths) {
  frames <- lapply(paths, read.csv)
  frame <- do.call(rbind, frames)
  frame[frame$repeat_id %in% repeat_ids_to_use, , drop = FALSE]
}
predictions <- read_selected_files(file_sets$predictions)
coefficients <- read_selected_files(file_sets$coefficients)
diagnostics <- read_selected_files(file_sets$diagnostics)
candidates <- read_selected_files(file_sets$candidates)
write.csv(predictions, file.path(results_root, "predictions.csv"), row.names = FALSE)
write.csv(diagnostics, file.path(results_root, "diagnostics.csv"), row.names = FALSE)
write.csv(candidates, file.path(results_root, "candidate_inner_metrics.csv"), row.names = FALSE)

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
  individual_minus_pretrained = metric_wide$mse.Individual - metric_wide$mse.Pretrained,
  overall_minus_pretrained = metric_wide$mse.Overall - metric_wide$mse.Pretrained,
  target_mean_minus_pretrained = metric_wide$`mse.Target mean` - metric_wide$mse.Pretrained
)
write.csv(paired, file.path(results_root, "paired_differences.csv"), row.names = FALSE)

selection_frequency <- as.data.frame(table(diagnostics$method, diagnostics$selected_candidate))
names(selection_frequency) <- c("method", "candidate_id", "selected_fits")
selection_frequency <- selection_frequency[selection_frequency$selected_fits > 0, ]
selection_frequency$selection_frequency <- selection_frequency$selected_fits / expected_fits
write.csv(selection_frequency, file.path(results_root, "candidate_selection_frequency.csv"), row.names = FALSE)

selected_coefficients <- coefficients[!is.na(coefficients$ensembl_gene_id), ]
if (nrow(selected_coefficients) > 0L) {
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
} else {
  coefficient_summary <- data.frame()
}
write.csv(coefficient_summary, file.path(results_root, "coefficient_stability.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(results_root, "sessionInfo.txt"))
writeLines(
  c(
    "Nested target-domain tuning experiment",
    sprintf("Completed: %s", format(Sys.time(), tz = "Asia/Shanghai")),
    sprintf("Repeats: %d", n_repeats),
    sprintf("Base seed: %d", base_seed),
    sprintf("Outer folds: %d", n_outer),
    sprintf("Inner folds: %d", n_inner),
    sprintf("Alpha grid: %s", paste(alpha_grid, collapse = ", ")),
    "First-stage lambda fixed at lambda.1se; feature count and second-stage lambda selected within each outer training set.",
    "Pretrained alpha selected from pancreatic-group inner-CV error.",
    "",
    paste(capture.output(print(method_summary, row.names = FALSE)), collapse = "\n")
  ),
  file.path(results_root, "run_summary.txt")
)
print(method_summary, row.names = FALSE)

