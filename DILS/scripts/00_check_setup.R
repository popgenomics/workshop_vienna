# Script 00 - Check that the workshop can run on this computer
#
# This script installs nothing. It only reports whether the project,
# packages and pedagogical tables are ready.
# It does not train a forest and it does not reveal scientific results.
#
# Run it from the DILS workshop project root.

required_packages = c("tidyverse", "abcrf")
missing_packages = required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0) {
  stop(
    "Please install these R packages before the workshop: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(tidyverse)
library(abcrf)

source("R/helpers.R")
source("R/teaching_helpers.R")

stop_if_not_project_root()

message("Working directory: ", getwd())
message("R version: ", paste(R.version$major, R.version$minor, sep = "."))
message("tidyverse version: ", as.character(utils::packageVersion("tidyverse")))
message("ggplot2 version: ", as.character(utils::packageVersion("ggplot2")))
message("abcrf version: ", as.character(utils::packageVersion("abcrf")))
message("ranger version: ", as.character(utils::packageVersion("ranger")))

pedagogical_tables = c(
  "data/reference_summary_statistics.tsv",
  "data/reference_parameters.tsv",
  "data/observed_summary_statistics.tsv",
  "data/goodness_of_fit_data.tsv",
  "data/locus_model_choice_data.tsv"
)
missing_tables = pedagogical_tables[!file.exists(pedagogical_tables)]
if (length(missing_tables) > 0) {
  stop(
    "Missing pedagogical tables: ",
    paste(missing_tables, collapse = ", "),
    call. = FALSE
  )
}
message("All five pedagogical tables are present.")

reference_ss = load_reference_summary_statistics()
reference_params = load_reference_parameters()
observed = load_observed_summary_statistics()
gof_data = load_goodness_of_fit_data()
locus_data = load_locus_model_choice_data()

message(
  "reference_summary_statistics: ",
  nrow(reference_ss), " rows x ", ncol(reference_ss), " columns"
)
message(
  "reference_parameters: ",
  nrow(reference_params), " rows x ", ncol(reference_params), " columns"
)
message(
  "observed_summary_statistics: ",
  nrow(observed), " row x ", ncol(observed), " columns"
)
message(
  "goodness_of_fit_data: ",
  nrow(gof_data), " rows x ", ncol(gof_data), " columns"
)
message(
  "locus_model_choice_data: ",
  nrow(locus_data), " rows x ", ncol(locus_data), " columns"
)

required_reference_meta = c(
  "migration_status", "scenario", "model", "simulation_id", "genomic_Ne"
)
missing_meta = setdiff(required_reference_meta, names(reference_ss))
if (length(missing_meta) > 0) {
  stop("Reference statistics table is missing: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}
if (!"simulation_id" %in% names(reference_params)) {
  stop("The parameter table has no simulation_id column.", call. = FALSE)
}
if (!"FST_avg" %in% names(observed)) {
  stop("The observed table has no FST_avg column.", call. = FALSE)
}
if (!"origin" %in% names(gof_data) || !"origin" %in% names(locus_data)) {
  stop("A goodness-of-fit or locus table is missing the origin column.", call. = FALSE)
}
message("Essential columns are present.")

if (any(c("model", "scenario", "migration_status") %in% names(observed))) {
  stop("The observed table should not contain model-class labels.", call. = FALSE)
}

ensure_results_dir()
probe_file = file.path("results", "00_setup_write_test.txt")
writeLines("ok", probe_file)
if (!file.exists(probe_file)) {
  stop("Cannot write into results/.", call. = FALSE)
}
unlink(probe_file)
message("The results/ directory is writable.")

message(
  "This workshop uses only the prepared tables in data/. ",
  "You do not need the original DILS run directory."
)
message("Setup check finished. Continue with teaching/index.qmd and script 01.")
