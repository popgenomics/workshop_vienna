# Script 02 - Ongoing migration versus current isolation
#
# Biological question
# After the split, do the two populations still exchange genes today
# (IM or SC), or is there no ongoing migration (AM or SI)?
#
# Before running the code
# Predict which scenarios belong to migration and which belong to
# isolation. Then predict why AM is grouped with SI.
#
# Data used here
# All rows of the reference summary-statistic table, plus the observed
# summary-statistic row.
#
# Run from the project root, section by section.

library(tidyverse)
library(abcrf)

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()
ensure_results_dir()

set.seed(20260915)

reference_ss = load_reference_summary_statistics()
observed = load_observed_summary_statistics()

print(count(reference_ss, migration_status, scenario))
print(unique(reference_ss$model))

# ---------------------------------------------------------------------------
# Task 1 - TODO-02-1
# Name the two classes of this comparison.
# ---------------------------------------------------------------------------

migration_scenarios = student_todo(
  "TODO-02-1",
  "Set migration_scenarios to the two scenario names with ongoing gene flow."
)
isolation_scenarios = student_todo(
  "TODO-02-1",
  "Set isolation_scenarios to the two scenario names with no ongoing gene flow."
)

if (!setequal(migration_scenarios, c("IM", "SC"))) {
  stop("TODO-02-1: ongoing migration is IM and SC.", call. = FALSE)
}
if (!setequal(isolation_scenarios, c("AM", "SI"))) {
  stop("TODO-02-1: current isolation is AM and SI.", call. = FALSE)
}

student_checkpoint(
  "CP-02-1",
  c(
    "Why does ancestral migration (AM) belong to current isolation?",
    "Which metadata column stores this grouping?"
  )
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-02-2
# Identify the response and convert it to a factor.
# ---------------------------------------------------------------------------

response_name = student_todo(
  "TODO-02-2",
  "Set response_name to the metadata column that codes migration versus isolation."
)
check_student_choice(
  "TODO-02-2",
  response_name,
  allowed = "migration_status",
  label = "response column"
)

training_ss = reference_ss
count_classes(training_ss[[response_name]], label = response_name)

pred_cols = predictor_names(training_ss)
if (!identical(pred_cols, names(observed))) {
  stop("Observed columns do not match the reference predictors.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Task 3 - TODO-02-3
# Name columns that would leak the answer if used as predictors.
# ---------------------------------------------------------------------------

leaky_columns = student_todo(
  "TODO-02-3",
  "Give a character vector of metadata names that must not enter the forest."
)
if (!is.character(leaky_columns) || length(leaky_columns) < 3) {
  stop("TODO-02-3: list at least three metadata columns.", call. = FALSE)
}
if (!all(leaky_columns %in% metadata_columns)) {
  stop("TODO-02-3: every leaky column must be one of metadata_columns.", call. = FALSE)
}
if (!"model" %in% leaky_columns) {
  stop("TODO-02-3: model would leak the answer. Include it.", call. = FALSE)
}

cleaned = drop_uninformative_predictors(
  training_predictors = training_ss[, pred_cols],
  observed_predictors = observed[, pred_cols]
)

training_data = data.frame(
  migration_status = factor(training_ss$migration_status),
  cleaned$training,
  check.names = FALSE
)
observed_data = cleaned$observed

if (!is.factor(training_data$migration_status)) {
  stop("The response must be a factor before abcrf().", call. = FALSE)
}
if (any(names(training_data)[-1] %in% metadata_columns)) {
  stop("Metadata leaked into the predictors.", call. = FALSE)
}

student_checkpoint(
  "CP-02-2",
  c(
    "Why does abcrf() refuse a character response?",
    "What would happen if model were included among predictors?"
  )
)

# ---------------------------------------------------------------------------
# Task 4 - TODO-02-4
# Complete the classification forest.
# ---------------------------------------------------------------------------

rf_model = student_todo(
  "TODO-02-4",
  "Call abcrf(migration_status ~ ., data = training_data, ntree = n_trees, paral = FALSE)."
)
if (!inherits(rf_model, "abcrf")) {
  stop("TODO-02-4: rf_model must be an abcrf object.", call. = FALSE)
}

print(rf_model)
save_abcrf_object(rf_model, "02_migration_vs_isolation_abcrf.rds")
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
save_abcrf_object(prediction, "02_migration_vs_isolation_prediction.rds")

# ---------------------------------------------------------------------------
# Task 5 - TODO-02-5
# Name the four quantities printed above.
# ---------------------------------------------------------------------------

meaning_of_allocation = student_todo(
  "TODO-02-5",
  "One sentence: what allocation is."
)
meaning_of_post_prob = student_todo(
  "TODO-02-5",
  "One sentence: what post.prob is, and why it is not the vote share."
)
if (!is.character(meaning_of_allocation) || nchar(trimws(meaning_of_allocation)) < 20) {
  stop("TODO-02-5: write a full sentence for allocation.", call. = FALSE)
}
if (!is.character(meaning_of_post_prob) || nchar(trimws(meaning_of_post_prob)) < 20) {
  stop("TODO-02-5: write a full sentence for post.prob.", call. = FALSE)
}

plot_abcrf_diagnostics(
  rf_model = rf_model,
  training_data = training_data,
  observed_predictors = observed_data,
  response_name = "migration_status",
  file_prefix = "02_migration_vs_isolation"
)

# Interpret the output
# Copy allocation, the votes and post.prob into the answer sheet.
# Also copy the OOB error. These numbers describe this comparison only.

# What this result does not show
# Winning IM+SC versus AM+SI does not identify IM versus SC, and it does
# not prove that the winning group is true.

# Take-home message
# The first forest answers one biological question: ongoing gene flow or
# not. The next script conditions on your allocation.

# Optional extension
# Lower n_trees to 50 only for a timing test. Do not trust that forest.

message("Use this allocation to fill selected_migration_status in script 03.")
