options(stringsAsFactors = FALSE)

script_path <- normalizePath(sub("^--file=", "", commandArgs()[grep("^--file=", commandArgs())]), winslash = "/")
root <- dirname(dirname(script_path))
course_root <- dirname(root)
.libPaths(c(file.path(course_root, ".Rlib"), .libPaths()))
suppressPackageStartupMessages(library(ptLasso))

processed_dir <- file.path(root, "processed")
results_dir <- file.path(root, "results")
fold_dir <- file.path(results_dir, "folds")
fit_dir <- file.path(results_dir, "fits")
dat <- read.csv(file.path(processed_dir, "gdsc_gemcitabine_top1000.csv.gz"), check.names = FALSE)
metadata_names <- c(
  "model_id", "model_name", "gdsc_cell_line_name", "cancer_group",
  "cancer_type_detail", "role", "DRUG_ID", "DRUG_NAME", "LN_IC50",
  "AUC", "Z_SCORE"
)
gene_names <- setdiff(names(dat), metadata_names)

for (fold in 1:5) {
  fit <- readRDS(file.path(fit_dir, sprintf("fold_%02d_fit.rds", fold)))
  coefficient_list <- coef(fit, alpha = fit$alphahat, s = "lambda.1se")
  group_legend <- fit$fit[[1]]$group.legend
  target_index <- match("Pancreatic Carcinoma", group_legend)
  stopifnot(!is.na(target_index))
  nonzero_genes <- function(coefficient_matrix) {
    coefficient_names <- rownames(coefficient_matrix)[as.vector(coefficient_matrix) != 0]
    intersect(coefficient_names, gene_names)
  }
  overall_genes <- nonzero_genes(coefficient_list$overall)
  individual_genes <- nonzero_genes(coefficient_list$individual[[target_index]])
  target_adjustment_genes <- nonzero_genes(coefficient_list$pretrain[[target_index]])
  support_genes <- list(
    Individual = individual_genes,
    Overall = overall_genes,
    Pretrained = union(overall_genes, target_adjustment_genes)
  )
  support_frame <- do.call(rbind, lapply(names(support_genes), function(method) {
    genes <- support_genes[[method]]
    if (length(genes) == 0L) {
      return(data.frame(fold = fold, method = method, feature_index = NA_integer_, ensembl_gene_id = NA_character_))
    }
    data.frame(fold = fold, method = method, feature_index = match(genes, gene_names), ensembl_gene_id = genes)
  }))
  write.csv(
    support_frame,
    file.path(fold_dir, sprintf("fold_%02d_support.csv", fold)),
    row.names = FALSE
  )
  message(sprintf("Rebuilt support for fold %d", fold))
}

support_files <- file.path(fold_dir, sprintf("fold_%02d_support.csv", 1:5))
supports <- do.call(rbind, lapply(support_files, read.csv))
supports <- supports[!is.na(supports$ensembl_gene_id), ]
if (nrow(supports) > 0L) {
  stability <- aggregate(
    fold ~ method + ensembl_gene_id,
    data = supports,
    FUN = function(values) length(unique(values)) / 5
  )
  names(stability)[names(stability) == "fold"] <- "selection_frequency"
  stability <- stability[order(stability$method, -stability$selection_frequency, stability$ensembl_gene_id), ]
} else {
  stability <- data.frame(method = character(), ensembl_gene_id = character(), selection_frequency = numeric())
}
write.csv(stability, file.path(results_dir, "selection_stability.csv"), row.names = FALSE)
