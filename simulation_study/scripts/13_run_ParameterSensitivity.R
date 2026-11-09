# Systematic one-factor-at-a-time sensitivity study for p, SNR, and sparsity.
#
# Usage:
#   Rscript simulation_study/scripts/13_run_ParameterSensitivity.R --reps=100
#
# Results are checkpointed after every replication-scenario pair and can be
# resumed safely with the same output directory.

full_args <- commandArgs(trailingOnly = FALSE)
file_args <- sub("^--file=", "", full_args[startsWith(full_args, "--file=")])
script_file <- if (length(file_args) > 0L) {
  normalizePath(file_args[1L], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path("simulation_study", "scripts", "13_run_ParameterSensitivity.R"),
    winslash = "/", mustWork = TRUE
  )
}
script_dir <- dirname(script_file)
study_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
project_root <- normalizePath(file.path(study_root, ".."), winslash = "/")

source(file.path(script_dir, "00_requirements.R"))
source(file.path(script_dir, "01_funct_PretrainedLassoSimulations.R"))
source(file.path(script_dir, "02_funct_ParameterSensitivity.R"))

parse_argument <- function(name, default) {
  prefix <- paste0("--", name, "=")
  matches <- commandArgs(trailingOnly = TRUE)
  matches <- matches[startsWith(matches, prefix)]
  if (length(matches) == 0L) return(default)
  sub(prefix, "", matches[1L], fixed = TRUE)
}

repetitions <- as.integer(parse_argument("reps", "100"))
base_seed <- as.integer(parse_argument("seed", "293700"))
output_dir <- parse_argument(
  "out", file.path(study_root, "results", "parameter_sensitivity")
)
stopifnot(repetitions >= 1L)

scenarios <- data.frame(
  scenario_id = c(
    "baseline", "p_250", "p_1000", "snr_0.8", "snr_3.0", "s_5", "s_20"
  ),
  factor = c(
    "baseline", "dimension", "dimension", "noise", "noise",
    "sparsity", "sparsity"
  ),
  level = c("baseline", "250", "1000", "0.8", "3.0", "5", "20"),
  p = c(500L, 250L, 1000L, 500L, 500L, 500L, 500L),
  active_size = c(10L, 10L, 10L, 10L, 10L, 5L, 20L),
  target_snr = c(1.8, 1.8, 1.8, 0.8, 3.0, 1.8, 1.8),
  stringsAsFactors = FALSE
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_file <- file.path(output_dir, "parameter_sensitivity_raw.csv")
failure_file <- file.path(output_dir, "parameter_sensitivity_failures.csv")

results <- if (file.exists(raw_file)) {
  utils::read.csv(raw_file, stringsAsFactors = FALSE)
} else {
  NULL
}
failures <- if (file.exists(failure_file)) {
  utils::read.csv(failure_file, stringsAsFactors = FALSE)
} else {
  NULL
}

completed_keys <- character()
if (!is.null(results) && nrow(results) > 0L) {
  counts <- table(paste(results$replication, results$scenario_id, sep = "::"))
  completed_keys <- names(counts)[counts == 3L]
}

settings <- c(
  "Parameter sensitivity study settings",
  "====================================",
  sprintf("Requested replications: %d", repetitions),
  sprintf("Base seed: %d", base_seed),
  "One-factor-at-a-time design around p=500, SNR=1.8, s=10",
  "Dimension p: 250, 500, 1000",
  "Target SNR: 0.8, 1.8, 3.0",
  "Active coefficients per shared/individual block s: 5, 10, 20",
  "D1: K=10 groups; n_k=100 per group",
  "D2: n=100 from target group 1",
  "Independent target-group test set: n=100",
  "rho=0.5; Gaussian predictors with within-block correlation",
  "Alpha grid=0,0.1,...,1; three-fold cross-validation",
  "Individual and Overall use lambda.min; Pretrained follows the paper's",
  "lambda.1se first-stage support/offset and lambda.min second stage.",
  "Metrics: test MSE, test MAE, TPR, FDP, support F1, selected size."
)
writeLines(settings, file.path(output_dir, "settings.txt"))
utils::write.csv(scenarios, file.path(output_dir, "scenarios.csv"), row.names = FALSE)

for (replication in seq_len(repetitions)) {
  for (scenario_index in seq_len(nrow(scenarios))) {
    scenario <- scenarios[scenario_index, ]
    task_key <- paste(replication, scenario$scenario_id, sep = "::")
    if (task_key %in% completed_keys) {
      message(sprintf(
        "Sensitivity replication %d/%d, %s already complete",
        replication, repetitions, scenario$scenario_id
      ))
      next
    }

    message(sprintf(
      "Sensitivity replication %d/%d, scenario %s",
      replication, repetitions, scenario$scenario_id
    ))

    scenario_rows <- tryCatch({
      component_seed <- safe_seed(
        base_seed, replication, component = 800L + scenario_index
      )
      simulation_seed <- safe_seed(
        base_seed, replication, component = 900L + scenario_index
      )
      components <- generate_figure6_components(
        seed = component_seed,
        k = 10L,
        p = scenario$p,
        active_size = scenario$active_size
      )
      data <- simulate_figure6_data(
        seed = simulation_seed,
        rho = 0.5,
        components = components,
        k = 10L,
        n_d1_per_group = 100L,
        n_d2 = 100L,
        n_test = 100L,
        p = scenario$p,
        target_snr = scenario$target_snr
      )
      rows <- fit_parameter_sensitivity_scenario(
        data = data,
        replication = replication,
        scenario_id = scenario$scenario_id,
        scenario_index = scenario_index,
        base_seed = base_seed
      )
      rows$factor <- scenario$factor
      rows$level <- scenario$level
      rows$p <- scenario$p
      rows$active_size <- scenario$active_size
      rows$rho <- 0.5
      rows$simulation_seed <- simulation_seed
      rows
    }, error = function(error) {
      failure <- data.frame(
        replication = replication,
        scenario_id = scenario$scenario_id,
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        error = conditionMessage(error),
        stringsAsFactors = FALSE
      )
      failures <<- rbind(failures, failure)
      utils::write.csv(failures, failure_file, row.names = FALSE)
      message("  FAILED: ", conditionMessage(error))
      NULL
    })

    if (!is.null(scenario_rows)) {
      if (!is.null(results)) {
        keep <- !(
          results$replication == replication &
            results$scenario_id == scenario$scenario_id
        )
        results <- results[keep, , drop = FALSE]
      }
      results <- rbind(results, scenario_rows)
      results <- results[
        order(results$replication, results$scenario_id, results$method),
      ]
      utils::write.csv(results, raw_file, row.names = FALSE)
      parameter_sensitivity_summary(results, output_dir)
      completed_keys <- union(completed_keys, task_key)
      message(sprintf("  completed in %.1f seconds", scenario_rows$elapsed_seconds[1L]))
    }
    flush.console()
  }
}

if (is.null(results) || nrow(results) == 0L) {
  stop("No parameter-sensitivity scenarios completed successfully.")
}

parameter_sensitivity_summary(results, output_dir)
write_session_info(output_dir)
message(
  "Parameter sensitivity results written to ",
  normalizePath(output_dir, winslash = "/")
)
