# Script 05 - Genomic heterogeneity in migration
#
# Biological question
# Given the demographic scenario from script 03, is migration
# homogeneous (1M) or heterogeneous (2M) across the genome?
#
# Before running the code
# Distinguish 1M, 2M, modeBarrier: bimodal, ancient migration under AM,
# ongoing migration under IM/SC, and possibly different barrier classes
# for M12 and M21.
#
# Data used here
# Simulations from selected_scenario only. SI has no 1M/2M comparison.
#
# Run from the project root, section by section.

library(tidyverse)
library(abcrf)

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()
ensure_results_dir()

set.seed(20260915)

# ---------------------------------------------------------------------------
# Task 1 - TODO-05-1
# Copy the scenario allocation from script 03.
# ---------------------------------------------------------------------------

selected_scenario = student_todo(
  "TODO-05-1",
  "Set this to IM, SC, AM or SI from script 03."
)
check_student_choice(
  "TODO-05-1",
  selected_scenario,
  allowed = c("IM", "SC", "AM", "SI"),
  label = "scenario"
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-05-2
# Define 1M and 2M in this workshop.
# ---------------------------------------------------------------------------

meaning_1M = student_todo(
  "TODO-05-2",
  "One sentence: what 1M means in a given direction."
)
meaning_2M = student_todo(
  "TODO-05-2",
  "One sentence: what 2M means under the bimodal barrier model."
)
if (!is.character(meaning_1M) || nchar(trimws(meaning_1M)) < 20) {
  stop("TODO-05-2: write a sentence for 1M.", call. = FALSE)
}
if (!is.character(meaning_2M) || nchar(trimws(meaning_2M)) < 20) {
  stop("TODO-05-2: write a sentence for 2M.", call. = FALSE)
}

if (selected_scenario == "SI") {
  message(
    "SI has no migration parameter. There is no 1M versus 2M ",
    "comparison under this scenario. Skip the forest and continue ",
    "with script 06."
  )
}

if (selected_scenario != "SI") {
  reference_ss = load_reference_summary_statistics()
  observed = load_observed_summary_statistics()
  result_prefix = paste0("05_migration_", selected_scenario)
  expected_models = switch(
    selected_scenario,
    IM = c("IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N"),
    SC = c("SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N"),
    AM = c("AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N")
  )

  training_ss = reference_ss %>%
    filter(.data$scenario == selected_scenario)
  if (!setequal(unique(training_ss$model), expected_models)) {
    stop("The selected scenario does not contain exactly the four expected models.", call. = FALSE)
  }

  print(count(training_ss, genomic_migration, model))
  count_classes(training_ss$genomic_migration, label = "genomic_migration")

  student_checkpoint(
    "CP-05-1",
    c(
      "Under AM, does 2M describe current or ancient migration?",
      "Can barrier loci differ between M12 and M21?"
    )
  )

  pred_cols = predictor_names(training_ss)
  cleaned = drop_uninformative_predictors(
    training_predictors = training_ss[, pred_cols],
    observed_predictors = observed[, pred_cols]
  )
  training_data = data.frame(
    genomic_migration = factor(
      training_ss$genomic_migration,
      levels = c("homogeneous", "heterogeneous")
    ),
    cleaned$training,
    check.names = FALSE
  )
  observed_data = cleaned$observed

  # ---------------------------------------------------------------------------
  # Task 3 - TODO-05-3
  # Train the 1M versus 2M forest.
  # ---------------------------------------------------------------------------

  rf_model = student_todo(
    "TODO-05-3",
    "Call abcrf(genomic_migration ~ ., data = training_data, ntree = n_trees, paral = FALSE)."
  )
  if (!inherits(rf_model, "abcrf")) {
    stop("TODO-05-3: rf_model must be an abcrf object.", call. = FALSE)
  }

  print(rf_model)
  save_abcrf_object(rf_model, paste0(result_prefix, "_abcrf.rds"))
  message("Out-of-bag prior error rate (prior.err): ", rf_model$prior.err)
  print(rf_model$model.rf$confusion.matrix)

  prediction = predict(
    rf_model,
    obs = observed_data,
    training = training_data,
    ntree = n_trees,
    paral = FALSE
  )
  print_abcrf_prediction(prediction)
  save_abcrf_object(prediction, paste0(result_prefix, "_prediction.rds"))

  plot_abcrf_diagnostics(
    rf_model = rf_model,
    training_data = training_data,
    observed_predictors = observed_data,
    response_name = "genomic_migration",
    file_prefix = result_prefix
  )

  message("Use this migration allocation to fill selected_genomic_migration in script 06.")
}

# ---------------------------------------------------------------------------
# Task 4 - TODO-05-4
# Reject a genome-scan overclaim.
# ---------------------------------------------------------------------------

why_2M_is_not_a_locus_scan = student_todo(
  "TODO-05-4",
  "Explain why selecting 2M does not already identify barrier loci."
)
if (!is.character(why_2M_is_not_a_locus_scan) || nchar(trimws(why_2M_is_not_a_locus_scan)) < 40) {
  stop("TODO-05-4: write a short paragraph.", call. = FALSE)
}

# Interpret the output
# heterogeneous means 2M under the bimodal barrier model used here.

# What this result does not show
# Selecting 2M already identifies the barrier loci. That sentence is false.
# Script 10 asks which observed loci resemble isolation simulations.

# Take-home message
# 2M is a genome-wide statement: loci do not all share the same migration
# rate. It is not a list of barrier genes.

# Optional extension
# If the scenario is AM, rewrite your conclusion so that it refers to
# ancient migration only.
