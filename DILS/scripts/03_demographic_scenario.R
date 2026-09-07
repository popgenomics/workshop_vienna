# Script 03 - Demographic scenario
#
# Biological question
# Within the migration-status group chosen in script 02, which
# demographic history is better supported?
#   If ongoing migration: IM versus SC
#   If isolation:         AM versus SI
#
# Before running the code
# Why would comparing IM and SC while retaining AM and SI simulations
# answer a different and poorly defined question?
#
# Data used here
# Only the simulations that belong to the selected_migration_status
# group, plus the observed summary-statistic row.
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
# Task 1 - TODO-03-1
# Copy the allocation from script 02.
# ---------------------------------------------------------------------------

selected_migration_status = student_todo(
  "TODO-03-1",
  "Set this to \"migration\" or \"isolation\" from script 02."
)
check_student_choice(
  "TODO-03-1",
  selected_migration_status,
  allowed = c("migration", "isolation"),
  label = "migration-status allocation"
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-03-2
# Answer the poorly defined comparison in your own words.
# ---------------------------------------------------------------------------

why_not_all_four_scenarios = student_todo(
  "TODO-03-2",
  "Write why IM versus SC must not keep AM and SI simulations in the training table."
)
if (!is.character(why_not_all_four_scenarios) || nchar(trimws(why_not_all_four_scenarios)) < 40) {
  stop("TODO-03-2: write at least two sentences.", call. = FALSE)
}

reference_ss = load_reference_summary_statistics()
observed = load_observed_summary_statistics()
result_prefix = paste0("03_demographic_scenario_", selected_migration_status)

if (selected_migration_status == "migration") {
  message("Branch: ongoing migration. Compare all IM models versus all SC models.")
  expected_scenarios = c("IM", "SC")
  expected_models = c(
    "IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N",
    "SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N"
  )
} else {
  message("Branch: isolation. Compare all AM models versus all SI models.")
  expected_scenarios = c("AM", "SI")
  expected_models = c(
    "AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N",
    "SI_1N", "SI_2N"
  )
}

# ---------------------------------------------------------------------------
# Task 3 - TODO-03-3
# Build the conditional training subset.
# ---------------------------------------------------------------------------

training_ss = student_todo(
  "TODO-03-3",
  "Filter reference_ss so that scenario is in expected_scenarios."
)
if (!is.data.frame(training_ss) || !"scenario" %in% names(training_ss)) {
  stop("TODO-03-3: training_ss must be a filtered data frame.", call. = FALSE)
}
if (!setequal(unique(training_ss$scenario), expected_scenarios)) {
  stop("TODO-03-3: the subset must contain exactly the expected scenarios.", call. = FALSE)
}
if (!setequal(unique(training_ss$model), expected_models)) {
  stop("TODO-03-3: the subset does not contain exactly the expected models.", call. = FALSE)
}

message("Number of simulations retained: ", nrow(training_ss))
print(count(training_ss, scenario, model))
count_classes(training_ss$scenario, label = "scenario")

student_checkpoint(
  "CP-03-1",
  c(
    "What is the response column now?",
    "Why are class sizes allowed to differ for AM versus SI?"
  )
)

pred_cols = predictor_names(training_ss)
cleaned = drop_uninformative_predictors(
  training_predictors = training_ss[, pred_cols],
  observed_predictors = observed[, pred_cols]
)
training_data = data.frame(
  scenario = factor(training_ss$scenario, levels = expected_scenarios),
  cleaned$training,
  check.names = FALSE
)
observed_data = cleaned$observed
if (any(names(training_data)[-1] %in% metadata_columns)) {
  stop("Metadata leaked into the predictors.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Task 4 - TODO-03-4
# Train the forest for this conditional comparison.
# ---------------------------------------------------------------------------

rf_model = student_todo(
  "TODO-03-4",
  "Call abcrf(scenario ~ ., data = training_data, ntree = n_trees, paral = FALSE)."
)
if (!inherits(rf_model, "abcrf")) {
  stop("TODO-03-4: rf_model must be an abcrf object.", call. = FALSE)
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
  response_name = "scenario",
  file_prefix = result_prefix
)

# ---------------------------------------------------------------------------
# Task 5 - TODO-03-5
# Write a conditional conclusion, not a general one.
# ---------------------------------------------------------------------------

conditional_conclusion = student_todo(
  "TODO-03-5",
  "Write one sentence that starts from the group selected in script 02."
)
if (!is.character(conditional_conclusion) || nchar(trimws(conditional_conclusion)) < 30) {
  stop("TODO-03-5: write a conditional sentence.", call. = FALSE)
}

# Interpret the output
# The allocation is IM versus SC, or AM versus SI, inside one group only.

# What this result does not show
# It does not rank all four scenarios at once. It does not prove the
# selected scenario is true.

# Take-home message
# Each later comparison is conditioned on this choice.

# Optional extension
# Inspect the confusion matrix and say which two scenarios are mixed.

message("Use this scenario allocation to fill selected_scenario in scripts 04, 05 and 06.")
