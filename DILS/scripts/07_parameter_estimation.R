# Script 07 - Parameter estimation
#
# Biological question
# Under the elementary model assembled in script 06, what values of the
# demographic parameters are compatible with the observed summaries?
#
# Before running the code
# Distinguish classification ABC-RF, regression ABC-RF, a factor
# response, a numeric response, prior, posterior, a central credibility
# interval, and an HPD interval (not computed here).
#
# Data used here
# Simulations of selected_complete_model only, paired by simulation_id
# with their parameters, plus the observed row.
#
# Run from the project root, section by section.

library(tidyverse)
library(abcrf)

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()
ensure_results_dir("parameter_models")

set.seed(20260915)

# ---------------------------------------------------------------------------
# Task 1 - TODO-07-1
# Copy the complete model name from script 06.
# ---------------------------------------------------------------------------

selected_complete_model = student_todo(
  "TODO-07-1",
  "Set this to the name assembled in script 06."
)
check_student_choice(
  "TODO-07-1",
  selected_complete_model,
  allowed = allowed_complete_models,
  label = "complete model"
)

reference_ss = load_reference_summary_statistics()
reference_params = load_reference_parameters()
observed = load_observed_summary_statistics()

stats_model = reference_ss %>%
  filter(.data$model == selected_complete_model)
params_model = reference_params %>%
  filter(.data$model == selected_complete_model)

if (!identical(stats_model$simulation_id, params_model$simulation_id)) {
  stop("simulation_id does not match between statistics and parameters.", call. = FALSE)
}

parameter_columns = setdiff(names(params_model), metadata_columns)
defined_parameters = parameter_columns[
  vapply(parameter_columns, function(p) !all(is.na(params_model[[p]])), logical(1))
]
message("Parameters defined under this model:")
print(defined_parameters)

student_checkpoint(
  "CP-07-1",
  c(
    "Which parameters are structurally absent (all NA) under this model?",
    "Is the response of regAbcrf() a factor or a number?"
  )
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-07-2
# Choose one defined parameter to inspect first.
# ---------------------------------------------------------------------------

chosen_parameter = student_todo(
  "TODO-07-2",
  "Set this to one name from defined_parameters, for example N1 or Tsplit."
)
check_student_choice(
  "TODO-07-2",
  chosen_parameter,
  allowed = defined_parameters,
  label = "parameter"
)

print(summary(params_model[[chosen_parameter]]))
message("This summary is the prior distribution used in the reference table.")

# ---------------------------------------------------------------------------
# Task 3 - TODO-07-3
# State the difference between the two ABC-RF tasks.
# ---------------------------------------------------------------------------

classification_versus_regression = student_todo(
  "TODO-07-3",
  "Write how abcrf() and regAbcrf() differ in response type and in what predict() returns."
)
if (!is.character(classification_versus_regression) || nchar(trimws(classification_versus_regression)) < 40) {
  stop("TODO-07-3: write a short paragraph.", call. = FALSE)
}

pred_cols = predictor_names(stats_model)
cleaned = drop_uninformative_predictors(
  training_predictors = stats_model[, pred_cols],
  observed_predictors = observed[, pred_cols]
)
observed_data = cleaned$observed

training_data = data.frame(
  parameter = params_model[[chosen_parameter]],
  cleaned$training,
  check.names = FALSE
)
names(training_data)[1] = chosen_parameter

# ---------------------------------------------------------------------------
# Task 4 - TODO-07-4
# Fit the regression forest for the chosen parameter.
# ---------------------------------------------------------------------------

model_rf = student_todo(
  "TODO-07-4",
  paste0(
    "Call regAbcrf(", chosen_parameter, " ~ ., data = training_data, ",
    "ntree = n_trees, paral = FALSE)."
  )
)
if (!inherits(model_rf, "regAbcrf")) {
  stop("TODO-07-4: model_rf must be a regAbcrf object.", call. = FALSE)
}

prediction = predict(
  model_rf,
  obs = observed_data,
  training = training_data,
  quantiles = c(0.025, 0.975)
)

quantile_matrix = as.matrix(prediction$quantiles)
posterior_one = tibble(
  parameter = chosen_parameter,
  posterior_expectation = as.numeric(prediction$expectation),
  posterior_median = as.numeric(prediction$med),
  posterior_variance = as.numeric(prediction$variance),
  quantile_0.025 = as.numeric(quantile_matrix[, 1]),
  quantile_0.975 = as.numeric(quantile_matrix[, 2])
)
print(posterior_one)
saveRDS(
  model_rf,
  file.path("results", "parameter_models", paste0(chosen_parameter, ".rds"))
)

# ---------------------------------------------------------------------------
# Task 5 - TODO-07-5
# Interpret the posterior without calling it the true value.
# ---------------------------------------------------------------------------

posterior_interpretation = student_todo(
  "TODO-07-5",
  "Interpret the median and the 2.5%-97.5% interval without saying true value or HPD."
)
if (!is.character(posterior_interpretation) || nchar(trimws(posterior_interpretation)) < 40) {
  stop("TODO-07-5: write a short paragraph.", call. = FALSE)
}

# Interpret the output
# quantile_0.025 and quantile_0.975 form a 95% central credibility interval.
# They are posterior quantiles, not an HPD interval.

# What this result does not show
# The posterior median is not the true parameter of the real populations.
# Estimating one model does not prove that the model is correct.

# Take-home message
# Classification chooses a class. Regression estimates a number under
# that class, using the same observed summaries.

# Optional extension
# Loop over every name in defined_parameters, as in solutions/07.

message("The remaining parameters can be estimated later from solutions/07.")
