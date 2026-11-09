# Functions for the one-factor-at-a-time parameter sensitivity experiment.

support_metrics <- function(estimated, truth, tolerance = 1e-8) {
  stopifnot(length(estimated) == length(truth))
  selected <- abs(estimated) > tolerance
  active <- abs(truth) > tolerance

  true_positive <- sum(selected & active)
  false_positive <- sum(selected & !active)
  selected_count <- sum(selected)
  active_count <- sum(active)
  tpr <- if (active_count > 0L) true_positive / active_count else NA_real_
  fdp <- if (selected_count > 0L) false_positive / selected_count else 0
  precision <- if (selected_count > 0L) {
    true_positive / selected_count
  } else {
    0
  }
  f1 <- if (is.finite(tpr) && precision + tpr > 0) {
    2 * precision * tpr / (precision + tpr)
  } else {
    0
  }

  c(
    selected_size = selected_count,
    true_support_size = active_count,
    tpr = tpr,
    fdp = fdp,
    f1 = f1
  )
}

fit_parameter_sensitivity_scenario <- function(
    data,
    replication,
    scenario_id,
    scenario_index,
    base_seed,
    nfolds = 3L,
    alpha_grid = seq(0, 1, by = 0.1)) {
  k <- length(unique(data$d1$groups))
  p <- ncol(data$d1$x)
  fold_seed <- safe_seed(
    base_seed,
    replication,
    component = 700L + 20L * scenario_index
  )

  elapsed <- system.time({
    fold_d2 <- balanced_folds(
      rep(1L, length(data$d2$y)), nfolds, fold_seed + 1L
    )
    set.seed(fold_seed + 2L)
    fit_individual <- glmnet::cv.glmnet(
      data$d2$x,
      data$d2$y,
      family = "gaussian",
      type.measure = "mse",
      foldid = fold_d2
    )
    prediction_individual <- as.numeric(stats::predict(
      fit_individual, newx = data$test$x, s = "lambda.min"
    ))
    coefficients_individual <- as.numeric(
      stats::coef(fit_individual, s = "lambda.min")
    )[-1L]

    pretrained <- fit_external_pretraining(
      data,
      fold_seed = fold_seed + 20L,
      alpha_grid = alpha_grid,
      nfolds = nfolds
    )

    groups_all <- c(data$d1$groups, rep(1L, length(data$d2$y)))
    group_factor <- factor(groups_all, levels = seq_len(k))
    dummies_all <- stats::model.matrix(
      ~ group_factor - 1L
    )[, -1L, drop = FALSE]
    x_all <- cbind(dummies_all, rbind(data$d1$x, data$d2$x))
    y_all <- c(data$d1$y, data$d2$y)
    penalty_all <- c(rep(0, ncol(dummies_all)), rep(1, p))
    fold_all <- balanced_folds(groups_all, nfolds, fold_seed + 3L)
    set.seed(fold_seed + 4L)
    fit_overall <- glmnet::cv.glmnet(
      x_all,
      y_all,
      family = "gaussian",
      type.measure = "mse",
      foldid = fold_all,
      penalty.factor = penalty_all
    )
    test_dummies <- matrix(0, nrow(data$test$x), ncol(dummies_all))
    prediction_overall <- as.numeric(stats::predict(
      fit_overall,
      newx = cbind(test_dummies, data$test$x),
      s = "lambda.min"
    ))
    coefficients_overall_full <- as.numeric(
      stats::coef(fit_overall, s = "lambda.min")
    )
    coefficients_overall <- coefficients_overall_full[-seq_len(k)]
  })[["elapsed"]]

  predictions <- list(
    Individual = prediction_individual,
    Pretrained = pretrained$prediction,
    Overall = prediction_overall
  )
  coefficients <- list(
    Individual = coefficients_individual,
    Pretrained = pretrained$effective_coefficients,
    Overall = coefficients_overall
  )
  truth <- data$beta_by_group[, 1L]

  rows <- lapply(names(predictions), function(method_name) {
    errors <- data$test$y - predictions[[method_name]]
    support <- support_metrics(coefficients[[method_name]], truth)
    data.frame(
      replication = replication,
      scenario_id = scenario_id,
      method = method_name,
      mse = mean(errors^2),
      mae = mean(abs(errors)),
      selected_size = unname(support["selected_size"]),
      true_support_size = unname(support["true_support_size"]),
      tpr = unname(support["tpr"]),
      fdp = unname(support["fdp"]),
      f1 = unname(support["f1"]),
      selected_alpha = if (method_name == "Pretrained") {
        pretrained$selected_alpha
      } else {
        NA_real_
      },
      first_stage_support_size = if (method_name == "Pretrained") {
        pretrained$overall_support_size
      } else {
        NA_real_
      },
      target_snr = data$target_snr,
      empirical_snr = data$empirical_snr,
      elapsed_seconds = elapsed,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

parameter_sensitivity_summary <- function(results, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  split_results <- split(
    results,
    interaction(
      results$scenario_id,
      results$method,
      drop = TRUE,
      lex.order = TRUE
    )
  )
  summary_rows <- lapply(split_results, function(x) {
    data.frame(
      scenario_id = x$scenario_id[1L],
      factor = x$factor[1L],
      level = x$level[1L],
      p = x$p[1L],
      active_size = x$active_size[1L],
      target_snr = x$target_snr[1L],
      method = x$method[1L],
      replications = nrow(x),
      mean_mse = mean(x$mse),
      se_mse = stats::sd(x$mse) / sqrt(nrow(x)),
      mean_mae = mean(x$mae),
      mean_tpr = mean(x$tpr),
      mean_fdp = mean(x$fdp),
      mean_f1 = mean(x$f1),
      se_f1 = stats::sd(x$f1) / sqrt(nrow(x)),
      mean_selected_size = mean(x$selected_size),
      mean_selected_alpha = if (all(is.na(x$selected_alpha))) {
        NA_real_
      } else {
        mean(x$selected_alpha, na.rm = TRUE)
      },
      stringsAsFactors = FALSE
    )
  })
  summary_table <- do.call(rbind, summary_rows)
  rownames(summary_table) <- NULL
  summary_table <- summary_table[
    order(summary_table$scenario_id, summary_table$method),
  ]
  utils::write.csv(
    summary_table,
    file.path(output_dir, "parameter_sensitivity_summary.csv"),
    row.names = FALSE
  )

  wide_mse <- reshape(
    results[, c("replication", "scenario_id", "method", "mse")],
    idvar = c("replication", "scenario_id"),
    timevar = "method",
    direction = "wide"
  )
  scenario_metadata <- unique(results[, c(
    "scenario_id", "factor", "level", "p", "active_size", "target_snr"
  )])
  wide_mse <- merge(wide_mse, scenario_metadata, by = "scenario_id")
  pairwise_rows <- lapply(split(wide_mse, wide_mse$scenario_id), function(x) {
    data.frame(
      scenario_id = x$scenario_id[1L],
      factor = x$factor[1L],
      level = x$level[1L],
      mean_pretrained_minus_individual = mean(
        x$mse.Pretrained - x$mse.Individual
      ),
      pretrained_win_rate_vs_individual = mean(
        x$mse.Pretrained < x$mse.Individual
      ),
      mean_pretrained_minus_overall = mean(
        x$mse.Pretrained - x$mse.Overall
      ),
      pretrained_win_rate_vs_overall = mean(
        x$mse.Pretrained < x$mse.Overall
      ),
      stringsAsFactors = FALSE
    )
  })
  pairwise_table <- do.call(rbind, pairwise_rows)
  rownames(pairwise_table) <- NULL
  utils::write.csv(
    pairwise_table,
    file.path(output_dir, "parameter_sensitivity_pairwise.csv"),
    row.names = FALSE
  )

  required_scenarios <- c(
    "baseline", "p_250", "p_1000", "snr_0.8", "snr_3.0", "s_5", "s_20"
  )
  if (all(required_scenarios %in% summary_table$scenario_id)) {
    plot_parameter_sensitivity(
      summary_table,
      metric = "mse",
      output_file = file.path(output_dir, "parameter_sensitivity_mse.png")
    )
    plot_parameter_sensitivity(
      summary_table,
      metric = "f1",
      output_file = file.path(output_dir, "parameter_sensitivity_f1.png")
    )
  }

  invisible(summary_table)
}

plot_parameter_sensitivity <- function(summary_table, metric, output_file) {
  stopifnot(metric %in% c("mse", "f1"))
  methods <- c("Individual", "Overall", "Pretrained")
  colors <- c(Individual = "#E67E22", Overall = "#42B883", Pretrained = "#3F7BD9")
  point_characters <- c(Individual = 16, Overall = 17, Pretrained = 15)
  panels <- list(
    Dimension = list(ids = c("p_250", "baseline", "p_1000"), x = c(250, 500, 1000)),
    SNR = list(ids = c("snr_0.8", "baseline", "snr_3.0"), x = c(0.8, 1.8, 3.0)),
    Sparsity = list(ids = c("s_5", "baseline", "s_20"), x = c(5, 10, 20))
  )
  mean_column <- if (metric == "mse") "mean_mse" else "mean_f1"
  se_column <- if (metric == "mse") "se_mse" else "se_f1"
  y_label <- if (metric == "mse") "Mean test MSE" else "Mean support F1"

  grDevices::png(output_file, width = 1900, height = 720, res = 170)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mfrow = c(1, 3), mar = c(4.4, 4.8, 3.2, 1.2), oma = c(0, 0, 2.2, 0))

  for (panel_name in names(panels)) {
    panel_spec <- panels[[panel_name]]
    panel <- summary_table[
      summary_table$scenario_id %in% panel_spec$ids &
        summary_table$method %in% methods,
    ]
    panel$panel_order <- match(panel$scenario_id, panel_spec$ids)
    panel$x_value <- panel_spec$x[panel$panel_order]
    panel_se <- panel[[se_column]]
    panel_se[!is.finite(panel_se)] <- 0
    y_lower <- panel[[mean_column]] - panel_se
    y_upper <- panel[[mean_column]] + panel_se
    y_range <- range(c(y_lower, y_upper), finite = TRUE)
    padding <- max(diff(y_range) * 0.12, 0.02)

    graphics::plot(
      range(panel_spec$x),
      y_range + c(-padding, padding),
      type = "n",
      xlab = panel_name,
      ylab = y_label,
      main = panel_name,
      xaxt = "n"
    )
    graphics::axis(1, at = panel_spec$x, labels = panel_spec$x)
    graphics::grid(col = "#E6E6E6")

    for (method_name in methods) {
      rows <- panel[panel$method == method_name, ]
      rows <- rows[order(rows$x_value), ]
      graphics::lines(
        rows$x_value,
        rows[[mean_column]],
        col = colors[method_name],
        lwd = 2,
        type = "b",
        pch = point_characters[method_name]
      )
      graphics::arrows(
        rows$x_value,
        rows[[mean_column]] - ifelse(is.finite(rows[[se_column]]),
                                     rows[[se_column]], 0),
        rows$x_value,
        rows[[mean_column]] + ifelse(is.finite(rows[[se_column]]),
                                     rows[[se_column]], 0),
        angle = 90,
        code = 3,
        length = 0.04,
        col = colors[method_name]
      )
    }
    if (panel_name == "Dimension") {
      graphics::legend(
        "topright",
        legend = methods,
        col = colors[methods],
        pch = point_characters[methods],
        lty = 1,
        lwd = 2,
        bty = "n",
        cex = 0.82
      )
    }
  }
  graphics::mtext(
    if (metric == "mse") {
      "Parameter sensitivity: prediction error"
    } else {
      "Parameter sensitivity: support recovery"
    },
    side = 3,
    outer = TRUE,
    line = 0.5,
    font = 2,
    cex = 1.08
  )
}
