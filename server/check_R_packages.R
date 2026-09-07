# Check DILS R packages and teaching tables (base R only).
# Intended to be called from server/check_prerequisites.sh or:
#   Rscript server/check_R_packages.R
#
# Exit status 1 if a blocking dependency is missing.

tested_r = "4.3.3"
tested = list(
  tidyverse = "2.0.0",
  ggplot2 = "3.5.2",
  dplyr = "1.1.4",
  readr = "2.1.5",
  abcrf = "1.9",
  ranger = "0.17.0"
)
required_packages = c("tidyverse", "abcrf", "ranger")
strict_versions = list(abcrf = "1.9")

dils_tables = c(
  "data/reference_summary_statistics.tsv",
  "data/reference_parameters.tsv",
  "data/observed_summary_statistics.tsv",
  "data/goodness_of_fit_data.tsv",
  "data/locus_model_choice_data.tsv"
)

args = commandArgs(trailingOnly = FALSE)
file_arg = sub("^--file=", "", grep("^--file=", args, value = TRUE))
if (length(file_arg) == 1 && nzchar(file_arg)) {
  script_path = normalizePath(file_arg, winslash = "/", mustWork = TRUE)
  repo_root = dirname(dirname(script_path))
} else {
  repo_root = getwd()
}
dils_root = file.path(repo_root, "DILS")

blocking = 0L
warn = 0L

note = function(...) {
  cat(paste0(...), "\n", sep = "")
}

fail = function(...) {
  blocking <<- blocking + 1L
  cat("ERROR: ", paste0(...), "\n", sep = "")
}

warn_msg = function(...) {
  warn <<- warn + 1L
  cat("WARNING: ", paste0(...), "\n", sep = "")
}

r_version = paste(R.version$major, R.version$minor, sep = ".")
note("R version: ", r_version, " (tested: ", tested_r, ")")
if (!identical(r_version, tested_r)) {
  warn_msg("R ", r_version, " differs from tested ", tested_r, ".")
}

if (!dir.exists(dils_root)) {
  fail("DILS directory not found at ", dils_root)
} else {
  note("DILS directory: ", dils_root)
}

installed_versions = installed.packages()[, "Version"]
pkg_ok = function(pkg) {
  pkg %in% names(installed_versions)
}

for (pkg in required_packages) {
  if (!pkg_ok(pkg)) {
    fail("R package missing: ", pkg)
  }
}

for (pkg in names(tested)) {
  if (!pkg_ok(pkg)) {
    if (!pkg %in% required_packages) {
      warn_msg("Optional/indirect package not installed: ", pkg)
    }
    next
  }
  got = unname(installed_versions[[pkg]])
  expected = tested[[pkg]]
  note(pkg, " version: ", got, " (tested: ", expected, ")")
  if (pkg %in% names(strict_versions) && !identical(got, strict_versions[[pkg]])) {
    fail(
      pkg, " ", got, " is installed; workshop requires ",
      strict_versions[[pkg]],
      ". Do not silently replace this version; random-forest results can change."
    )
  } else if (!identical(got, expected)) {
    warn_msg(pkg, " ", got, " differs from tested ", expected, ".")
  }
}

if (dir.exists(dils_root)) {
  for (rel in dils_tables) {
    path = file.path(dils_root, rel)
    if (!file.exists(path)) {
      fail("Missing DILS table: ", path)
    } else {
      info = file.info(path)
      note("Found ", rel, " (", info$size, " bytes)")
    }
  }
  dictionary = file.path(dils_root, "data", "data_dictionary.tsv")
  if (!file.exists(dictionary)) {
    warn_msg("data_dictionary.tsv is missing (not required to run the forests).")
  }
}

note("Blocking errors: ", blocking, "; version/presence warnings: ", warn)
if (blocking > 0L) {
  quit(save = "no", status = 1L)
}
quit(save = "no", status = 0L)
