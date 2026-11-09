# Package requirements for the STAT 293 Pretrained Lasso project.
#
# This file is normally sourced by a runner in simulation_study/scripts.
# Packages are loaded from the project-local .Rlib first so that the simulation
# scripts use the same versions across sessions.

if (!exists("project_root", inherits = FALSE)) {
  requirements_file <- normalizePath(
    sys.frame(1)$ofile, winslash = "/", mustWork = TRUE
  )
  project_root <- normalizePath(
    file.path(dirname(requirements_file), "..", ".."),
    winslash = "/", mustWork = TRUE
  )
}
project_lib <- normalizePath(
  file.path(project_root, ".Rlib"), winslash = "/", mustWork = TRUE
)
.libPaths(c(project_lib, .libPaths()))

required_packages <- c("glmnet", "ptLasso")
report_packages <- c("knitr", "rmarkdown")
all_packages <- c(required_packages, report_packages)

missing_packages <- all_packages[
  !vapply(all_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ",
    paste(missing_packages, collapse = ", "),
    ". Install them into .Rlib before continuing."
  )
}

suppressPackageStartupMessages({
  library(glmnet)
  library(ptLasso)
})

package_versions <- vapply(
  all_packages,
  function(pkg) as.character(utils::packageVersion(pkg)),
  character(1)
)

message(
  "Loaded project packages: ",
  paste(sprintf("%s %s", names(package_versions), package_versions),
        collapse = "; ")
)
