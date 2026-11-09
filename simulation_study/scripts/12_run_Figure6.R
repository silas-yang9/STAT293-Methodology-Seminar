# Reproduce the external-dataset experiment in Figure 6.
#
# Usage:
#   Rscript simulation_study/scripts/12_run_Figure6.R --reps=100
#
# Results are checkpointed after each replication and can be resumed safely.

full_args <- commandArgs(trailingOnly = FALSE)
file_args <- sub("^--file=", "", full_args[startsWith(full_args, "--file=")])
script_file <- if (length(file_args) > 0L) {
  normalizePath(file_args[1L], winslash = "/", mustWork = TRUE)
} else {
  normalizePath(
    file.path("simulation_study", "scripts", "12_run_Figure6.R"),
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
base_seed <- as.integer(parse_argument("seed", "293600"))
output_dir <- parse_argument("out", file.path(study_root, "results", "figure6"))
stopifnot(repetitions >= 1L)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
raw_file <- file.path(output_dir, "figure6_raw.csv")
failure_file <- file.path(output_dir, "figure6_failures.csv")

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

expected_rows_per_replication <- 3L * 4L
completed <- integer()
if (!is.null(results) && nrow(results) > 0L) {
  counts <- table(results$replication)
  completed <- as.integer(names(counts)[counts == expected_rows_per_replication])
}

settings <- c(
  "Figure 6 reproduction settings",
  "================================",
  sprintf("Requested replications: %d", repetitions),
  sprintf("Base seed: %d", base_seed),
  "D1: K = 10; n_k = 100; p = 500",
  "D2: n = 100 from group 1",
  "Independent group-1 test set: n = 100",
  "10 shared and 10 disjoint individual coefficients per group",
  "rho = 0, 0.5, 1",
  "Gaussian noise scaled to target average SNR = 1.8",
  "Alpha grid = 0, 0.1, ..., 1",
  "Three-fold CV; first-stage support/offset uses lambda.1se;",
  "second-stage and comparison predictions use lambda.min.",
  "All-data fits include unpenalized group intercepts."
)
writeLines(settings, file.path(output_dir, "settings.txt"))

for (replication in seq_len(repetitions)) {
  if (replication %in% completed) {
    message(sprintf("Figure 6 replication %d/%d already complete", replication,
                    repetitions))
    next
  }

  message(sprintf("Figure 6 replication %d/%d", replication, repetitions))
  components <- generate_figure6_components(
    safe_seed(base_seed, replication, 300L)
  )

  replication_rows <- tryCatch({
    rows <- list()
    for (i in seq_along(c(0, 0.5, 1))) {
      rho_value <- c(0, 0.5, 1)[i]
      message(sprintf("  rho=%s", format(rho_value, trim = TRUE)))
      simulation_seed <- safe_seed(
        base_seed, replication, 400L + as.integer(100 * rho_value)
      )
      data <- simulate_figure6_data(
        seed = simulation_seed,
        rho = rho_value,
        components = components
      )
      rows[[i]] <- fit_figure6_scenario(
        data = data,
        replication = replication,
        rho = rho_value,
        base_seed = base_seed
      )
    }
    do.call(rbind, rows)
  }, error = function(error) {
    failure <- data.frame(
      replication = replication,
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
    results <- results[order(results$replication, results$rho,
                             results$method), ]
    utils::write.csv(results, raw_file, row.names = FALSE)
    summarize_figure6(results, output_dir)
    message(sprintf("  completed in %.1f seconds",
                    sum(unique(replication_rows[, c("rho", "elapsed_seconds")])$elapsed_seconds)))
  }
  flush.console()
}

if (is.null(results) || nrow(results) == 0L) {
  stop("No Figure 6 replications completed successfully.")
}

summarize_figure6(results, output_dir)
write_session_info(output_dir)
message("Figure 6 results written to ", normalizePath(output_dir, winslash = "/"))
