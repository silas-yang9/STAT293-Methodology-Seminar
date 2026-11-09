options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ptLasso))

processed_dir <- file.path(root, "processed")
results_root <- file.path(root, "results", "full_experiment")
dir.create(results_root, recursive = TRUE, showWarnings = FALSE)

dat <- read.csv(file.path(processed_dir, "gdsc_gemcitabine_top2000.csv.gz"), check.names = FALSE)
metadata_names <- c(
  "model_id", "model_name", "gdsc_cell_line_name", "cancer_group",
  "cancer_type_detail", "role", "DRUG_ID", "DRUG_NAME", "LN_IC50",
  "AUC", "Z_SCORE"
)
all_gene_names <- setdiff(names(dat), metadata_names)
stopifnot(length(all_gene_names) == 2000L)
stopifnot(sum(dat$role == "target") == 30L)

scenarios <- data.frame(
  scenario = c(
    "primary_1000_ln_ic50_1se",
    "features_500_ln_ic50_1se",
    "features_2000_ln_ic50_1se",
    "auc_1000_1se",
    "lambda_min_1000_ln_ic50"
  ),
  outcome = c("LN_IC50", "LN_IC50", "LN_IC50", "AUC", "LN_IC50"),
  features = c(1000L, 500L, 2000L, 1000L, 1000L),
  lambda_choice = c("lambda.1se", "lambda.1se", "lambda.1se", "lambda.1se", "lambda.min"),
  repeats = c(20L, 10L, 10L, 10L, 10L)
)
write.csv(scenarios, file.path(results_root, "scenario_settings.csv"), row.names = FALSE)

base_seed <- 2932026L
n_outer <- 5L
n_inner <- 3L
alpha_grid <- seq(0, 1, by = 0.1)
groups <- dat$cancer_group
target_rows <- which(dat$role == "target")
source_rows <- which(dat$role == "external")

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

extract_coefficients <- function(fit, gene_names, lambda_choice) {
  coefficient_list <- coef(fit, alpha = fit$alphahat, s = lambda_choice)
  group_legend <- fit$fit[[1]]$group.legend
  target_index <- match("Pancreatic Carcinoma", group_legend)
  stopifnot(!is.na(target_index))
  overall <- gene_vector(coefficient_list$overall, gene_names)
  individual <- gene_vector(coefficient_list$individual[[target_index]], gene_names)
  target_adjustment <- gene_vector(coefficient_list$pretrain[[target_index]], gene_names)
  pretrained <- (1 - fit$alphahat) * overall + target_adjustment
  list(Individual = individual, Overall = overall, Pretrained = pretrained)
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

scenario_summaries <- list()
paired_summaries <- list()

for (scenario_index in seq_len(nrow(scenarios))) {
  setting <- scenarios[scenario_index, ]
  scenario_name <- setting$scenario
  scenario_dir <- file.path(results_root, scenario_name)
  checkpoint_dir <- file.path(scenario_dir, "checkpoints")
  assignment_dir <- file.path(scenario_dir, "assignments")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(assignment_dir, recursive = TRUE, showWarnings = FALSE)

  gene_names <- all_gene_names[seq_len(setting$features)]
  x <- as.matrix(dat[, gene_names, drop = FALSE])
  storage.mode(x) <- "double"
  y <- dat[[setting$outcome]]
  if (any(!is.finite(y))) stop(sprintf("Non-finite %s values", setting$outcome))

  for (repeat_id in seq_len(setting$repeats)) {
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
      prediction_path <- file.path(checkpoint_dir, paste0(prefix, "_predictions.csv"))
      coefficient_path <- file.path(checkpoint_dir, paste0(prefix, "_coefficients.csv"))
      diagnostic_path <- file.path(checkpoint_dir, paste0(prefix, "_diagnostics.csv"))
      if (all(file.exists(c(prediction_path, coefficient_path, diagnostic_path)))) {
        message(sprintf("%s repeat %d fold %d complete; skipping", scenario_name, repeat_id, fold))
        next
      }

      test_rows <- target_rows[outer_fold[target_rows] == fold]
      target_train_rows <- setdiff(target_rows, test_rows)
      train_rows <- c(source_rows, target_train_rows)
      inner_seed <- repeat_seed + 100L * fold
      inner_fold <- make_inner_folds(groups[train_rows], inner_seed)
      message(sprintf("Fitting %s: repeat %d/%d, fold %d/%d", scenario_name, repeat_id, setting$repeats, fold, n_outer))

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
        s = setting$lambda_choice,
        overall.lambda = setting$lambda_choice,
        alphahat.choice = "overall",
        group.intercepts = TRUE
      )
      prediction <- predict(
        fit,
        xtest = x[test_rows, , drop = FALSE],
        groupstest = groups[test_rows],
        ytest = y[test_rows],
        alphatype = "fixed",
        s = setting$lambda_choice,
        return.link = TRUE
      )
      target_mean <- mean(y[target_train_rows])
      predictions <- rbind(
        data.frame(scenario = scenario_name, repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Individual", predicted = as.numeric(prediction$yhatind)),
        data.frame(scenario = scenario_name, repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Overall", predicted = as.numeric(prediction$yhatoverall)),
        data.frame(scenario = scenario_name, repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Pretrained", predicted = as.numeric(prediction$yhatpre)),
        data.frame(scenario = scenario_name, repeat_id = repeat_id, fold = fold, model_id = dat$model_id[test_rows], observed = y[test_rows], method = "Target mean", predicted = target_mean)
      )
      write.csv(predictions, prediction_path, row.names = FALSE)

      coefficient_vectors <- extract_coefficients(fit, gene_names, setting$lambda_choice)
      coefficients <- do.call(rbind, lapply(names(coefficient_vectors), function(method) {
        values <- coefficient_vectors[[method]]
        selected <- which(abs(values) > 1e-12)
        if (length(selected) == 0L) {
          return(data.frame(scenario = scenario_name, repeat_id = repeat_id, fold = fold, method = method, ensembl_gene_id = NA_character_, coefficient = NA_real_))
        }
        data.frame(scenario = scenario_name, repeat_id = repeat_id, fold = fold, method = method, ensembl_gene_id = names(values)[selected], coefficient = values[selected])
      }))
      write.csv(coefficients, coefficient_path, row.names = FALSE)
      write.csv(
        data.frame(
          scenario = scenario_name,
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
          outcome = setting$outcome,
          inner_folds = n_inner,
          lambda_choice = setting$lambda_choice
        ),
        diagnostic_path,
        row.names = FALSE
      )
    }
  }

  expected_fits <- setting$repeats * n_outer
  prediction_files <- list.files(checkpoint_dir, pattern = "_predictions\\.csv$", full.names = TRUE)
  coefficient_files <- list.files(checkpoint_dir, pattern = "_coefficients\\.csv$", full.names = TRUE)
  diagnostic_files <- list.files(checkpoint_dir, pattern = "_diagnostics\\.csv$", full.names = TRUE)
  stopifnot(length(prediction_files) == expected_fits)
  stopifnot(length(coefficient_files) == expected_fits)
  stopifnot(length(diagnostic_files) == expected_fits)
  predictions <- do.call(rbind, lapply(prediction_files, read.csv))
  coefficients <- do.call(rbind, lapply(coefficient_files, read.csv))
  diagnostics <- do.call(rbind, lapply(diagnostic_files, read.csv))
  write.csv(predictions, file.path(scenario_dir, "predictions.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(scenario_dir, "diagnostics.csv"), row.names = FALSE)

  repeat_metrics <- do.call(rbind, lapply(split(predictions, list(predictions$repeat_id, predictions$method), drop = TRUE), function(frame) {
    data.frame(scenario = scenario_name, repeat_id = unique(frame$repeat_id), method = unique(frame$method), metric_function(frame))
  }))
  repeat_metrics <- repeat_metrics[order(repeat_metrics$repeat_id, repeat_metrics$method), ]
  write.csv(repeat_metrics, file.path(scenario_dir, "repeat_metrics.csv"), row.names = FALSE)

  method_summary <- do.call(rbind, lapply(split(repeat_metrics, repeat_metrics$method), function(frame) {
    data.frame(
      scenario = scenario_name,
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
  write.csv(method_summary, file.path(scenario_dir, "method_summary.csv"), row.names = FALSE)
  scenario_summaries[[scenario_name]] <- method_summary

  metric_wide <- reshape(repeat_metrics[, c("repeat_id", "method", "mse")], idvar = "repeat_id", timevar = "method", direction = "wide")
  paired <- data.frame(
    scenario = scenario_name,
    repeat_id = metric_wide$repeat_id,
    individual_minus_pretrained = metric_wide$mse.Individual - metric_wide$mse.Pretrained,
    overall_minus_pretrained = metric_wide$mse.Overall - metric_wide$mse.Pretrained,
    target_mean_minus_pretrained = metric_wide$`mse.Target mean` - metric_wide$mse.Pretrained,
    individual_minus_overall = metric_wide$mse.Individual - metric_wide$mse.Overall
  )
  write.csv(paired, file.path(scenario_dir, "paired_differences.csv"), row.names = FALSE)
  paired_summaries[[scenario_name]] <- paired

  alpha_distribution <- as.data.frame(table(diagnostics$selected_alpha))
  names(alpha_distribution) <- c("alpha", "fold_count")
  alpha_distribution$scenario <- scenario_name
  alpha_distribution$proportion <- alpha_distribution$fold_count / sum(alpha_distribution$fold_count)
  write.csv(alpha_distribution, file.path(scenario_dir, "alpha_distribution.csv"), row.names = FALSE)

  selected_coefficients <- coefficients[!is.na(coefficients$ensembl_gene_id), ]
  if (nrow(selected_coefficients) > 0L) {
    coefficient_summary <- do.call(rbind, lapply(split(selected_coefficients, list(selected_coefficients$method, selected_coefficients$ensembl_gene_id), drop = TRUE), function(frame) {
      data.frame(
        scenario = scenario_name,
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
  write.csv(coefficient_summary, file.path(scenario_dir, "coefficient_stability.csv"), row.names = FALSE)
}

all_method_summaries <- do.call(rbind, scenario_summaries)
all_paired <- do.call(rbind, paired_summaries)
write.csv(all_method_summaries, file.path(results_root, "method_summary_all.csv"), row.names = FALSE)
write.csv(all_paired, file.path(results_root, "paired_differences_all.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(results_root, "sessionInfo.txt"))
writeLines(
  c(
    "Full GDSC gemcitabine real-data experiment",
    sprintf("Completed: %s", format(Sys.time(), tz = "Asia/Shanghai")),
    sprintf("Base seed: %d", base_seed),
    sprintf("Outer folds: %d", n_outer),
    sprintf("Inner folds: %d", n_inner),
    sprintf("Alpha grid: %s", paste(alpha_grid, collapse = ", ")),
    "",
    paste(capture.output(print(all_method_summaries, row.names = FALSE)), collapse = "\n")
  ),
  file.path(results_root, "run_summary.txt")
)
print(all_method_summaries, row.names = FALSE)
