# Reusable functions for simulation experiments based on Figures 4 and 6 in Craig et al. (2025),
# "Pretraining and the lasso".

balanced_folds <- function(groups, nfolds, seed) {
  set.seed(seed)
  fold_id <- integer(length(groups))
  for (g in sort(unique(groups))) {
    index <- which(groups == g)
    fold_id[index] <- sample(rep(seq_len(nfolds), length.out = length(index)))
  }
  fold_id
}

safe_seed <- function(base_seed, replication, component = 0L) {
  as.integer(base_seed + 10000L * as.integer(replication) + as.integer(component))
}

simulate_figure4_data <- function(seed,
                                  k = 9L,
                                  n_per_group = 55L,
                                  p = 1000L,
                                  common_size = 10L,
                                  individual_size = 10L,
                                  coefficient = 2.5,
                                  noise_sd = 10,
                                  test_multiplier = 10L) {
  stopifnot(
    p >= common_size + k * individual_size,
    k >= 2L,
    n_per_group >= 3L,
    test_multiplier >= 1L
  )

  set.seed(seed)
  groups <- rep(seq_len(k), each = n_per_group)
  groups_test <- rep(seq_len(k), each = n_per_group * test_multiplier)

  beta_common <- numeric(p)
  beta_common[seq_len(common_size)] <- coefficient

  beta_by_group <- matrix(0, nrow = p, ncol = k)
  individual_support <- vector("list", k)
  for (g in seq_len(k)) {
    first <- common_size + (g - 1L) * individual_size + 1L
    last <- first + individual_size - 1L
    individual_support[[g]] <- first:last
    beta_by_group[, g] <- beta_common
    beta_by_group[first:last, g] <- coefficient
  }

  make_sample <- function(sample_groups) {
    x <- matrix(stats::rnorm(length(sample_groups) * p),
                nrow = length(sample_groups), ncol = p)
    mu <- numeric(length(sample_groups))
    for (g in seq_len(k)) {
      index <- which(sample_groups == g)
      mu[index] <- as.numeric(x[index, , drop = FALSE] %*% beta_by_group[, g])
    }
    error <- stats::rnorm(length(sample_groups), sd = noise_sd)
    list(x = x, y = mu + error, mu = mu, error = error)
  }

  train <- make_sample(groups)
  test <- make_sample(groups_test)

  list(
    x = train$x,
    y = train$y,
    mu = train$mu,
    error = train$error,
    groups = groups,
    xtest = test$x,
    ytest = test$y,
    mutest = test$mu,
    errortest = test$error,
    groupstest = groups_test,
    beta_common = beta_common,
    beta_by_group = beta_by_group,
    common_support = seq_len(common_size),
    individual_support = individual_support,
    empirical_snr = stats::var(train$mu) / stats::var(train$error),
    theoretical_snr = sum(beta_by_group[, 1L]^2) / noise_sd^2
  )
}

fit_figure4_setting <- function(data,
                                replication,
                                nfolds,
                                lambda_rule,
                                base_seed,
                                alpha_values = c(0, 0.5, 0.95)) {
  stopifnot(lambda_rule %in% c("lambda.1se", "lambda.min"))

  fold_seed <- safe_seed(
    base_seed,
    replication,
    component = 100L + 10L * nfolds
  )
  fold_id <- balanced_folds(data$groups, nfolds, fold_seed)

  set.seed(safe_seed(
    base_seed,
    replication,
    component = if (lambda_rule == "lambda.1se") 301L else 302L
  ))

  elapsed <- system.time({
    fit <- ptLasso::cv.ptLasso(
      x = data$x,
      y = data$y,
      groups = data$groups,
      alphalist = alpha_values,
      family = "gaussian",
      type.measure = "mse",
      nfolds = nfolds,
      foldid = fold_id,
      overall.lambda = lambda_rule,
      s = lambda_rule,
      alphahat.choice = "overall",
      group.intercepts = TRUE
    )
  })[["elapsed"]]

  model_labels <- c(
    "Overall",
    sprintf("Pretrained (alpha = %s)", as.character(alpha_values)),
    "Individual"
  )
  cv_error <- c(
    unname(fit$erroverall["overall"]),
    as.numeric(fit$errpre[, "overall"]),
    unname(fit$errind["overall"])
  )

  test_error <- numeric(length(model_labels))
  prediction_first <- predict(
    fit,
    xtest = data$xtest,
    groupstest = data$groupstest,
    ytest = data$ytest,
    alpha = alpha_values[1L],
    alphatype = "fixed",
    s = lambda_rule
  )
  test_error[1L] <- unname(prediction_first$erroverall["allGroups"])
  test_error[length(test_error)] <- unname(prediction_first$errind["allGroups"])
  test_error[2L] <- unname(prediction_first$errpre["allGroups"])

  if (length(alpha_values) > 1L) {
    for (j in 2:length(alpha_values)) {
      prediction <- predict(
        fit,
        xtest = data$xtest,
        groupstest = data$groupstest,
        ytest = data$ytest,
        alpha = alpha_values[j],
        alphatype = "fixed",
        s = lambda_rule
      )
      test_error[j + 1L] <- unname(prediction$errpre["allGroups"])
    }
  }

  data.frame(
    replication = replication,
    simulation_seed = safe_seed(base_seed, replication, 1L),
    nfolds = nfolds,
    lambda_rule = lambda_rule,
    model = model_labels,
    alpha = c(NA_real_, alpha_values, NA_real_),
    cv_error = cv_error,
    test_error = test_error,
    relative_cv_bias = 100 * (test_error - cv_error) / test_error,
    empirical_snr = data$empirical_snr,
    theoretical_snr = data$theoretical_snr,
    elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )
}

summarize_figure4 <- function(results, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  split_results <- split(
    results,
    interaction(results$nfolds, results$lambda_rule, results$model,
                drop = TRUE, lex.order = TRUE)
  )
  summary_rows <- lapply(split_results, function(x) {
    data.frame(
      nfolds = x$nfolds[1L],
      lambda_rule = x$lambda_rule[1L],
      model = x$model[1L],
      replications = nrow(x),
      mean_cv_error = mean(x$cv_error),
      mean_test_error = mean(x$test_error),
      mean_relative_bias = mean(x$relative_cv_bias),
      median_relative_bias = stats::median(x$relative_cv_bias),
      q025_relative_bias = unname(stats::quantile(x$relative_cv_bias, 0.025)),
      q975_relative_bias = unname(stats::quantile(x$relative_cv_bias, 0.975)),
      stringsAsFactors = FALSE
    )
  })
  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL
  summary_table <- summary_table[
    order(summary_table$nfolds, summary_table$model,
          summary_table$lambda_rule),
  ]
  utils::write.csv(
    summary_table,
    file.path(output_dir, "figure4_summary.csv"),
    row.names = FALSE
  )

  models <- c(
    "Overall",
    "Pretrained (alpha = 0)",
    "Pretrained (alpha = 0.5)",
    "Pretrained (alpha = 0.95)",
    "Individual"
  )
  plot_labels <- c("Overall", "alpha = 0", "alpha = 0.5",
                   "alpha = 0.95", "Individual")
  lambda_rules <- c("lambda.1se", "lambda.min")
  colors <- c("#65C57A", "#3C8DCA")

  grDevices::png(
    file.path(output_dir, "figure4_cv_bias.png"),
    width = 1800,
    height = 1100,
    res = 170
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(2, 1), mar = c(2, 4, 3, 1),
                oma = c(5.5, 5, 3, 0))

  y_limit <- range(results$relative_cv_bias, finite = TRUE)
  y_padding <- diff(y_limit) * 0.06
  y_limit <- y_limit + c(-y_padding, y_padding)

  for (fold_count in c(3L, 10L)) {
    panel <- results[results$nfolds == fold_count, , drop = FALSE]
    box_data <- list()
    positions <- numeric()
    box_colors <- character()
    counter <- 0L
    for (i in seq_along(models)) {
      for (j in seq_along(lambda_rules)) {
        counter <- counter + 1L
        keep <- panel$model == models[i] &
          panel$lambda_rule == lambda_rules[j]
        box_data[[counter]] <- panel$relative_cv_bias[keep]
        positions[counter] <- i + c(-0.18, 0.18)[j]
        box_colors[counter] <- colors[j]
      }
    }
    graphics::boxplot(
      box_data,
      at = positions,
      col = box_colors,
      border = "#333333",
      outline = TRUE,
      xaxt = "n",
      ylim = y_limit,
      ylab = "",
      main = sprintf("Number of CV folds: %d", fold_count)
    )
    graphics::abline(h = 0, lty = 2, col = "#555555")
    graphics::axis(1, at = seq_along(models), labels = FALSE)
    if (fold_count == 10L) {
      graphics::mtext(plot_labels, side = 1, at = seq_along(models), line = 1,
                      las = 2, cex = 0.8)
    }
    if (fold_count == 3L) {
      graphics::legend(
        "topright",
        legend = lambda_rules,
        fill = colors,
        border = "#333333",
        bty = "n",
        cex = 0.85
      )
    }
  }
  graphics::mtext(
    "Experiment 1: cross-validation bias",
    side = 3,
    outer = TRUE,
    line = 1,
    font = 2,
    cex = 1.1
  )
  graphics::mtext(
    "Relative CV bias: 100 x (test error - CV error) / test error",
    side = 2,
    outer = TRUE,
    line = 2.5,
    cex = 0.9
  )

  invisible(summary_table)
}

generate_figure6_components <- function(seed,
                                        k = 10L,
                                        p = 500L,
                                        active_size = 10L) {
  stopifnot(p >= active_size * (k + 1L))
  set.seed(seed)

  beta_shared <- numeric(p)
  beta_shared[seq_len(active_size)] <- abs(stats::rnorm(active_size))

  beta_individual <- matrix(0, nrow = p, ncol = k)
  for (g in seq_len(k)) {
    first <- active_size + (g - 1L) * active_size + 1L
    last <- first + active_size - 1L
    beta_individual[first:last, g] <- stats::rnorm(active_size)
  }

  list(beta_shared = beta_shared, beta_individual = beta_individual)
}

simulate_figure6_data <- function(seed,
                                  rho,
                                  components,
                                  k = 10L,
                                  n_d1_per_group = 100L,
                                  n_d2 = 100L,
                                  n_test = 100L,
                                  p = 500L,
                                  target_snr = 1.8) {
  stopifnot(rho >= 0, rho <= 1)
  set.seed(seed)

  beta_by_group <- rho * components$beta_shared +
    (1 - rho) * components$beta_individual

  average_signal_variance <- mean(colSums(beta_by_group^2))
  noise_sd <- sqrt(average_signal_variance / target_snr)

  make_group_sample <- function(group, n) {
    x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
    mu <- as.numeric(x %*% beta_by_group[, group])
    error <- stats::rnorm(n, sd = noise_sd)
    list(x = x, y = mu + error, mu = mu, error = error)
  }

  d1_parts <- lapply(seq_len(k), make_group_sample, n = n_d1_per_group)
  d1 <- list(
    x = do.call(rbind, lapply(d1_parts, `[[`, "x")),
    y = unlist(lapply(d1_parts, `[[`, "y"), use.names = FALSE),
    mu = unlist(lapply(d1_parts, `[[`, "mu"), use.names = FALSE),
    error = unlist(lapply(d1_parts, `[[`, "error"), use.names = FALSE),
    groups = rep(seq_len(k), each = n_d1_per_group)
  )
  d2 <- make_group_sample(1L, n_d2)
  test <- make_group_sample(1L, n_test)

  list(
    d1 = d1,
    d2 = d2,
    test = test,
    beta_by_group = beta_by_group,
    noise_sd = noise_sd,
    target_snr = target_snr,
    empirical_snr = stats::var(d1$mu) / stats::var(d1$error)
  )
}

fit_external_pretraining <- function(data,
                                     fold_seed,
                                     alpha_grid = seq(0, 1, by = 0.1),
                                     nfolds = 3L) {
  k <- length(unique(data$d1$groups))
  p <- ncol(data$d1$x)

  group_factor <- factor(data$d1$groups, levels = seq_len(k))
  dummies <- stats::model.matrix(~ group_factor - 1L)[, -1L, drop = FALSE]
  x_overall <- cbind(dummies, data$d1$x)
  penalty_overall <- c(rep(0, ncol(dummies)), rep(1, p))
  fold_d1 <- balanced_folds(data$d1$groups, nfolds, fold_seed)

  set.seed(fold_seed + 1L)
  fit_overall <- glmnet::cv.glmnet(
    x_overall,
    data$d1$y,
    family = "gaussian",
    type.measure = "mse",
    foldid = fold_d1,
    penalty.factor = penalty_overall,
    keep = TRUE
  )

  first_stage_lambda <- fit_overall$lambda.1se
  coefficient_vector <- as.numeric(
    stats::coef(fit_overall, s = first_stage_lambda)
  )
  feature_coefficients <- coefficient_vector[-seq_len(k)]
  support <- which(feature_coefficients != 0)

  zero_dummies_d2 <- matrix(0, nrow(data$d2$x), k - 1L)
  zero_dummies_test <- matrix(0, nrow(data$test$x), k - 1L)
  offset_d2 <- as.numeric(stats::predict(
    fit_overall,
    newx = cbind(zero_dummies_d2, data$d2$x),
    s = first_stage_lambda
  ))
  offset_test <- as.numeric(stats::predict(
    fit_overall,
    newx = cbind(zero_dummies_test, data$test$x),
    s = first_stage_lambda
  ))

  fold_d2 <- balanced_folds(rep(1L, length(data$d2$y)), nfolds,
                            fold_seed + 2L)
  alpha_fits <- vector("list", length(alpha_grid))
  cv_errors <- numeric(length(alpha_grid))

  for (i in seq_along(alpha_grid)) {
    alpha_value <- alpha_grid[i]
    if (alpha_value == 1) {
      penalty_factor <- rep(1, p)
    } else if (alpha_value == 0) {
      penalty_factor <- rep(1e9, p)
      penalty_factor[support] <- 1
    } else {
      penalty_factor <- rep(1 / alpha_value, p)
      penalty_factor[support] <- 1
    }

    set.seed(fold_seed + 100L + i)
    alpha_fits[[i]] <- glmnet::cv.glmnet(
      data$d2$x,
      data$d2$y,
      family = "gaussian",
      type.measure = "mse",
      foldid = fold_d2,
      offset = (1 - alpha_value) * offset_d2,
      penalty.factor = penalty_factor,
      keep = TRUE
    )
    cv_errors[i] <- min(alpha_fits[[i]]$cvm)
  }

  selected_index <- which.min(cv_errors)
  selected_alpha <- alpha_grid[selected_index]
  selected_fit <- alpha_fits[[selected_index]]
  prediction <- as.numeric(stats::predict(
    selected_fit,
    newx = data$test$x,
    s = "lambda.min",
    newoffset = (1 - selected_alpha) * offset_test
  ))
  second_stage_coefficients <- as.numeric(
    stats::coef(selected_fit, s = "lambda.min")
  )[-1L]
  effective_coefficients <-
    (1 - selected_alpha) * feature_coefficients +
    second_stage_coefficients

  list(
    prediction = prediction,
    selected_alpha = selected_alpha,
    cv_errors = data.frame(alpha = alpha_grid, cv_mse = cv_errors),
    overall_support_size = length(support),
    first_stage_coefficients = feature_coefficients,
    second_stage_coefficients = second_stage_coefficients,
    effective_coefficients = effective_coefficients
  )
}

fit_figure6_scenario <- function(data,
                                 replication,
                                 rho,
                                 base_seed,
                                 nfolds = 3L,
                                 alpha_grid = seq(0, 1, by = 0.1)) {
  fold_seed <- safe_seed(
    base_seed,
    replication,
    component = 500L + as.integer(100 * rho)
  )

  elapsed <- system.time({
    fold_d2 <- balanced_folds(
      rep(1L, length(data$d2$y)), nfolds, fold_seed + 1L
    )
    set.seed(fold_seed + 2L)
    fit_d2 <- glmnet::cv.glmnet(
      data$d2$x,
      data$d2$y,
      family = "gaussian",
      type.measure = "mse",
      foldid = fold_d2
    )
    prediction_d2 <- as.numeric(stats::predict(
      fit_d2, newx = data$test$x, s = "lambda.min"
    ))

    pretrained <- fit_external_pretraining(
      data,
      fold_seed = fold_seed + 20L,
      alpha_grid = alpha_grid,
      nfolds = nfolds
    )

    group1_index <- data$d1$groups == 1L
    x_group1_combined <- rbind(data$d1$x[group1_index, , drop = FALSE],
                               data$d2$x)
    y_group1_combined <- c(data$d1$y[group1_index], data$d2$y)
    fold_group1 <- balanced_folds(
      rep(1L, length(y_group1_combined)), nfolds, fold_seed + 3L
    )
    set.seed(fold_seed + 4L)
    fit_group1 <- glmnet::cv.glmnet(
      x_group1_combined,
      y_group1_combined,
      family = "gaussian",
      type.measure = "mse",
      foldid = fold_group1
    )
    prediction_group1 <- as.numeric(stats::predict(
      fit_group1, newx = data$test$x, s = "lambda.min"
    ))

    groups_all <- c(data$d1$groups, rep(1L, length(data$d2$y)))
    group_factor <- factor(groups_all, levels = seq_len(max(groups_all)))
    dummies_all <- stats::model.matrix(~ group_factor - 1L)[, -1L, drop = FALSE]
    x_all <- cbind(dummies_all, rbind(data$d1$x, data$d2$x))
    y_all <- c(data$d1$y, data$d2$y)
    penalty_all <- c(rep(0, ncol(dummies_all)), rep(1, ncol(data$d1$x)))
    fold_all <- balanced_folds(groups_all, nfolds, fold_seed + 5L)
    set.seed(fold_seed + 6L)
    fit_all <- glmnet::cv.glmnet(
      x_all,
      y_all,
      family = "gaussian",
      type.measure = "mse",
      foldid = fold_all,
      penalty.factor = penalty_all
    )
    test_dummies <- matrix(0, nrow(data$test$x), ncol(dummies_all))
    prediction_all <- as.numeric(stats::predict(
      fit_all,
      newx = cbind(test_dummies, data$test$x),
      s = "lambda.min"
    ))
  })[["elapsed"]]

  predictions <- list(
    "D2 only" = prediction_d2,
    "Pretraining" = pretrained$prediction,
    "D1 group 1 + D2" = prediction_group1,
    "D1 all + D2" = prediction_all
  )
  pse <- vapply(
    predictions,
    function(prediction) mean((data$test$y - prediction)^2),
    numeric(1)
  )

  data.frame(
    replication = replication,
    simulation_seed = safe_seed(base_seed, replication,
                                400L + as.integer(100 * rho)),
    rho = rho,
    method = names(pse),
    pse = as.numeric(pse),
    log_pse = log(as.numeric(pse)),
    selected_alpha = c(NA_real_, pretrained$selected_alpha, NA_real_, NA_real_),
    first_stage_support_size = c(
      NA_integer_, pretrained$overall_support_size, NA_integer_, NA_integer_
    ),
    target_snr = data$target_snr,
    empirical_snr = data$empirical_snr,
    elapsed_seconds = elapsed,
    stringsAsFactors = FALSE
  )
}

summarize_figure6 <- function(results, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  split_results <- split(
    results,
    interaction(results$rho, results$method, drop = TRUE, lex.order = TRUE)
  )
  summary_rows <- lapply(split_results, function(x) {
    data.frame(
      rho = x$rho[1L],
      method = x$method[1L],
      replications = nrow(x),
      mean_pse = mean(x$pse),
      median_pse = stats::median(x$pse),
      q025_pse = unname(stats::quantile(x$pse, 0.025)),
      q975_pse = unname(stats::quantile(x$pse, 0.975)),
      stringsAsFactors = FALSE
    )
  })
  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL
  summary_table <- summary_table[order(summary_table$rho, summary_table$method), ]
  utils::write.csv(
    summary_table,
    file.path(output_dir, "figure6_summary.csv"),
    row.names = FALSE
  )

  alpha_rows <- results[
    results$method == "Pretraining",
    c("replication", "rho", "selected_alpha", "first_stage_support_size",
      "empirical_snr")
  ]
  utils::write.csv(
    alpha_rows,
    file.path(output_dir, "figure6_alpha_selection.csv"),
    row.names = FALSE
  )

  methods <- c("D2 only", "Pretraining", "D1 group 1 + D2", "D1 all + D2")
  colors <- c("#E7C64A", "#4FC27A", "#F05A59", "#4B8FE2")
  panel_titles <- c(
    `0` = "No shared support",
    `0.5` = "Blended support",
    `1` = "No individual support"
  )

  grDevices::png(
    file.path(output_dir, "figure6_pse_comparison.png"),
    width = 1900,
    height = 850,
    res = 170
  )
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(7, 4.5, 4, 1), oma = c(0, 0, 3, 0))

  y_limit <- range(results$log_pse, finite = TRUE)
  y_padding <- diff(y_limit) * 0.08
  y_limit <- y_limit + c(-y_padding, y_padding)

  for (rho_value in c(0, 0.5, 1)) {
    panel <- results[results$rho == rho_value, , drop = FALSE]
    box_data <- lapply(methods, function(method) {
      panel$log_pse[panel$method == method]
    })
    pretraining_median <- stats::median(
      panel$log_pse[panel$method == "Pretraining"]
    )
    graphics::boxplot(
      box_data,
      col = colors,
      border = "#333333",
      names = FALSE,
      ylim = y_limit,
      ylab = if (rho_value == 0) "log(PSE)" else "",
      main = sprintf("%s\nrho = %s", panel_titles[as.character(rho_value)],
                     format(rho_value, trim = TRUE))
    )
    graphics::abline(h = pretraining_median, lwd = 1.5, col = "#222222")
    graphics::axis(1, at = seq_along(methods), labels = methods,
                   las = 2, cex.axis = 0.75)
  }
  graphics::mtext(
    "Experiment 2: prediction squared error",
    side = 3,
    outer = TRUE,
    line = 1,
    font = 2,
    cex = 1.1
  )

  invisible(summary_table)
}

write_session_info <- function(output_dir) {
  writeLines(
    capture.output(utils::sessionInfo()),
    file.path(output_dir, "sessionInfo.txt")
  )
}
