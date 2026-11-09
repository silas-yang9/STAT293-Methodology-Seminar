# Install and load the R packages used by the STAT 293 project.
# Run from the repository root with: Rscript 00_requirements.R

script_file <- tryCatch(
  normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
  error = function(e) NA_character_
)
project_root <- if (is.na(script_file)) {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
} else {
  dirname(script_file)
}

project_library <- file.path(project_root, ".Rlib")
dir.create(project_library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(normalizePath(project_library, winslash = "/"), .libPaths()))

cran_packages <- c(
  "glmnet", "ggplot2", "gridExtra", "knitr", "rmarkdown"
)
missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_cran) > 0L) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("ptLasso", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  remotes::install_github("erincr/ptLasso", upgrade = "never")
}

required_packages <- c("glmnet", "ptLasso", "ggplot2", "knitr", "rmarkdown")
still_missing <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(still_missing) > 0L) {
  stop("Packages could not be installed: ", paste(still_missing, collapse = ", "))
}

suppressPackageStartupMessages({
  library(glmnet)
  library(ptLasso)
  library(ggplot2)
  library(knitr)
  library(rmarkdown)
})

versions <- vapply(
  required_packages,
  function(package) as.character(utils::packageVersion(package)),
  character(1)
)
message(
  "R package requirements are available: ",
  paste(sprintf("%s %s", names(versions), versions), collapse = "; ")
)
