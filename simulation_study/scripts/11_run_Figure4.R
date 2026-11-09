# Reproduce the cross-validation bias experiment in Figure 4.
#
# Usage:
#   Rscript simulation_study/scripts/11_run_Figure4.R --reps=100
#
# Results are checkpointed after each replication and can be resumed safely.

full_args <- commandArgs(trailingOnly = FALSE)
file_args <- sub("^--file=", "", full_args[startsWith(full_args, "--file=")])
script_file <- if (length(file_args) > 0L) {
  normalizePath(file_args[1L], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path("simulation_study", "scripts", "11_run_Figure4.R"),
    winslash = "/", mustWork = TRUE
  )
}
script_dir <- dirname(script_file)
study_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
project_root <- normalizePath(file.path(study_root, ".."), winslash = "/")

source(file.path(script_dir, "00_requirements.R"))
source(file.path(script_dir, "01_funct_PretrainedLassoSimulations.R"))

parse_argument <- function(name, default) {
  prefix <- paste0("--", name, "=")
  match <- commandArgs(trailingOnly = TRUE)
  match <- match[startsWith(match, prefix)]
  if (length(match) == 0L) return(default)
  sub(prefix, "", match[1L], fixed = TRUE)
}

repetitions <- as.integer(parse_argument("reps", "100"))
base_seed <- as.integer(parse_argument("seed", "293400"))
output_dir <- parse_argument("out", file.path(study_root, "results", "figure4"))
stopifnot(repetitions >= 1L)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_file <- file.path(output_dir, "figure4_raw.csv")
failure_file <- file.path(output_dir, "figure4_failures.csv")

if (file.exists(raw_file)) {
  results <- utils::read.csv(raw_file, stringsAsFactors = FALSE)
} else {
  results <- NULL
}
if (file.exists(failure_file)) {
  failures <- utils::read.csv(failure_file, stringsAsFactors = FALSE)
} else {
  failures <- NULL
}

expected_rows_per_replication <- 2L * 2L * 5L
completed <- integer()
if (!is.null(results) && nrow(results) > 0L) {
  counts <- table(results$replication)
  completed <- as.integer(names(counts)[counts == expected_rows_per_replication])
}

settings <- c(
  "Figure 4 reproduction settings",
  "================================",
  sprintf("Requested replications: %d", repetitions),
  sprintf("Base seed: %d", base_seed),
  "K = 9; n_k = 55; n = 495; p = 1000",
  "10 common coefficients = 2.5",
  "10 disjoint individual coefficients per group = 2.5",
  "Gaussian noise standard deviation = 10",
  "Independent test size = 4950",
  "Alpha values = 0, 0.5, 0.95",
  "CV folds = 3, 10",
  "Lambda rules = lambda.1se, lambda.min",
  "The selected lambda rule is used for the first-stage offset/support and",
  "for evaluating the overall, pretrained, and individual models."
)
writeLines(settings, file.path(output_dir, "settings.txt"))

for (replication in seq_len(repetitions)) {
  if (replication %in% completed) {
    message(sprintf("Figure 4 replication %d/%d already complete", replication,
                    repetitions))
    next
  }

  message(sprintf("Figure 4 replication %d/%d", replication, repetitions))
  simulation_seed <- safe_seed(base_seed, replication, 1L)
  data <- simulate_figure4_data(simulation_seed)

  replication_rows <- tryCatch({
    rows <- list()
    index <- 0L
    for (fold_count in c(3L, 10L)) {
      for (lambda_rule in c("lambda.1se", "lambda.min")) {
        index <- index + 1L
        message(sprintf("  folds=%d; rule=%s", fold_count, lambda_rule))
        rows[[index]] <- fit_figure4_setting(
          data = data,
          replication = replication,
          nfolds = fold_count,
          lambda_rule = lambda_rule,
          base_seed = base_seed
        )
      }
    }
    do.call(rbind, rows)
  }, error = function(error) {
    failure <- data.frame(
      replication = replication,
      simulation_seed = simulation_seed,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      error = conditionMessage(error),
      stringsAsFactors = FALSE
    )
    failures <<- rbind(failures, failure)
    utils::write.csv(failures, failure_file, row.names = FALSE)
    message("  FAILED: ", conditionMessage(error))
    NULL
  })

  if (!is.null(replication_rows)) {
    if (!is.null(results)) {
      results <- results[results$replication != replication, , drop = FALSE]
    }
    results <- rbind(results, replication_rows)
    results <- results[order(results$replication, results$nfolds,
                             results$lambda_rule, results$model), ]
    utils::write.csv(results, raw_file, row.names = FALSE)
    summarize_figure4(results, output_dir)
    message(sprintf("  completed in %.1f seconds",
                    sum(unique(replication_rows[, c("nfolds", "lambda_rule",
                                                    "elapsed_seconds")])$elapsed_seconds)))
  }
  flush.console()
}

if (is.null(results) || nrow(results) == 0L) {
  stop("No Figure 4 replications completed successfully.")
}

summarize_figure4(results, output_dir)
write_session_info(output_dir)
message("Figure 4 results written to ", normalizePath(output_dir, winslash = "/"))
