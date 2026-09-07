# Script 04 - Genomic heterogeneity in Ne
#
# Biological question
# Given the demographic scenario selected in script 03, is Ne
# homogeneous (1N) or heterogeneous (2N) across the genome?
#
# Before running the code
# Link 1N, 2N, linked selection and variation of Ne among regions.
# Predict whether 2N identifies selected loci directly.
#
# Data used here
# Only simulations from selected_scenario, plus the observed row.
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
# Task 1 - TODO-04-1
# Copy the scenario allocation from script 03.
# ---------------------------------------------------------------------------

selected_scenario = student_todo(
  "TODO-04-1",
  "Set this to IM, SC, AM or SI from script 03."
)
check_student_choice(
  "TODO-04-1",
  selected_scenario,
  allowed = c("IM", "SC", "AM", "SI"),
  label = "scenario"
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-04-2
# Match each Ne class to a biological meaning.
# ---------------------------------------------------------------------------

meaning_1N = student_todo(
  "TODO-04-2",
  "One sentence: what homogeneous Ne (1N) assumes."
)
meaning_2N = student_todo(
  "TODO-04-2",
  "One sentence: what heterogeneous Ne (2N) represents, including linked selection."
)
if (!is.character(meaning_1N) || nchar(trimws(meaning_1N)) < 20) {
  stop("TODO-04-2: write a sentence for 1N.", call. = FALSE)
}
if (!is.character(meaning_2N) || nchar(trimws(meaning_2N)) < 20) {
  stop("TODO-04-2: write a sentence for 2N.", call. = FALSE)
}

reference_ss = load_reference_summary_statistics()
observed = load_observed_summary_statistics()
result_prefix = paste0("04_Ne_", selected_scenario)

expected_models = switch(
  selected_scenario,
  IM = c("IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N"),
  SC = c("SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N"),
  AM = c("AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N"),
  SI = c("SI_1N", "SI_2N")
)

training_ss = reference_ss %>%
  filter(.data$scenario == selected_scenario)
if (!setequal(unique(training_ss$model), expected_models)) {
  stop("The selected scenario does not contain exactly the expected models.", call. = FALSE)
}

print(count(training_ss, genomic_Ne, model))
count_classes(training_ss$genomic_Ne, label = "genomic_Ne")

student_checkpoint(
  "CP-04-1",
  c(
    "What is the response column?",
    "Why must this forest stay inside one demographic scenario?"
  )
)

pred_cols = predictor_names(training_ss)
cleaned = drop_uninformative_predictors(
  training_predictors = training_ss[, pred_cols],
  observed_predictors = observed[, pred_cols]
)
training_data = data.frame(
  genomic_Ne = factor(
    training_ss$genomic_Ne,
    levels = c("homogeneous", "heterogeneous")
  ),
  cleaned$training,
  check.names = FALSE
)
observed_data = cleaned$observed

# ---------------------------------------------------------------------------
# Task 3 - TODO-04-3
# Train the 1N versus 2N forest.
# ---------------------------------------------------------------------------

rf_model = student_todo(
  "TODO-04-3",
  "Call abcrf(genomic_Ne ~ ., data = training_data, ntree = n_trees, paral = FALSE)."
)
if (!inherits(rf_model, "abcrf")) {
  stop("TODO-04-3: rf_model must be an abcrf object.", call. = FALSE)
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
  response_name = "genomic_Ne",
  file_prefix = result_prefix
)

# ---------------------------------------------------------------------------
# Task 4 - TODO-04-4
# Reject a causal overclaim.
# ---------------------------------------------------------------------------

why_2N_is_not_selected_loci = student_todo(
  "TODO-04-4",
  "Explain why choosing 2N is not the same as identifying selected loci."
)
if (!is.character(why_2N_is_not_selected_loci) || nchar(trimws(why_2N_is_not_selected_loci)) < 40) {
  stop("TODO-04-4: write a short paragraph.", call. = FALSE)
}

# Interpret the output
# homogeneous means 1N. heterogeneous means 2N.

# What this result does not show
# 2N is a statistical representation of genomic variation in Ne. It does
# not name the loci under selection.

# Take-home message
# Linked selection can be modelled as heterogeneous Ne. That is not a
# genome scan for selected genes.

# Optional extension
# Compare piA_avg in 1N versus 2N simulations inside this scenario.

message("Use this Ne allocation to fill selected_genomic_Ne in script 06.")
