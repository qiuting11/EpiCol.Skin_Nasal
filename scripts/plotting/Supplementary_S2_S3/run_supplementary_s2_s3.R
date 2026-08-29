#!/usr/bin/env Rscript

repo_root <- Sys.getenv("PAPER_SUBMISSION_CODE", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
script_dir <- file.path(repo_root, "scripts", "plotting", "Supplementary_S2_S3")

run_script <- function(file) {
  source(file.path(script_dir, file), local = new.env(parent = globalenv()))
}

run_script("Fig_S2B_genus_eCF_eCM_deployment.R")
run_script("Fig_S2C_species_eCF_eCM_architecture.R")
run_script("Fig_S2D_to_S2G_conservative_filtering_delta_refined.R")
run_script("Fig_S3A_controlled_parameter_eCF_bubble.R")
run_script("Fig_S3B_controlled_parameter_mcnemar_delta.R")
run_script("Fig_S3C_to_S3F_strict_sensitive_presence_taxonomy.R")
