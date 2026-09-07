# Script 10 - Locus-specific classification
#
# Biological question
# Under the selected IM/SC 2M model, which observed loci look more like
# complete two-way isolation than like bidirectional migration?
#
# Before running the code
# Reconstruct the four-step argument: genome-wide 2M, monolocus training
# classes, classification forest, then classification of each locus.
#
# Data used here
# data/locus_model_choice_data.tsv
#
# Run from the project root, section by section.

library(tidyverse)
library(abcrf)

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()
ensure_results_dir("locus_model_choice")

workshop_seed = 20260915
set.seed(workshop_seed)

# ---------------------------------------------------------------------------
# Task 1 - TODO-10-1
# Copy the complete model from script 06. The analysis is relevant only
# for IM or SC with heterogeneous migration (2M).
# ---------------------------------------------------------------------------

selected_complete_model = student_todo(
  "TODO-10-1",
  "Set this to the IM_2M_* or SC_2M_* name assembled in script 06."
)
check_student_choice(
  "TODO-10-1",
  selected_complete_model,
  allowed = allowed_complete_models,
  label = "complete model"
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-10-2
# Assign the two forest sizes. Both default to n_trees, but they have
# different statistical roles.
# ---------------------------------------------------------------------------

n_classification_trees = student_todo(
  "TODO-10-2",
  "Set n_classification_trees to n_trees. This forest chooses allocation."
)
n_posterior_probability_trees = student_todo(
  "TODO-10-2",
  "Set n_posterior_probability_trees to n_trees. This forest estimates post.prob."
)
if (!identical(as.integer(n_classification_trees), n_trees)) {
  stop("TODO-10-2: use n_trees for the classification forest.", call. = FALSE)
}
if (!identical(as.integer(n_posterior_probability_trees), n_trees)) {
  stop("TODO-10-2: use n_trees for the posterior-probability forest.", call. = FALSE)
}

ongoing_two_rate = grepl("^(IM|SC)_2M_", selected_complete_model)
if (!ongoing_two_rate) {
  stop(
    "A locus-specific migration-versus-isolation analysis is relevant here\n",
    "only after selecting ongoing migration (IM or SC) with heterogeneous\n",
    "migration among loci (2M).",
    call. = FALSE
  )
}

locus_data = load_locus_model_choice_data()
table_model = unique(locus_data$fitted_model)
if (!identical(table_model, selected_complete_model)) {
  stop(
    "The locus table was built under ", paste(table_model, collapse = ", "),
    " but selected_complete_model is ", selected_complete_model, ".",
    call. = FALSE
  )
}

message("Selected complete model: ", selected_complete_model)
message("Locus table: ", nrow(locus_data), " rows x ", ncol(locus_data), " columns")
message("Classification trees: ", n_classification_trees)
message("Posterior-probability trees: ", n_posterior_probability_trees)

# ---------------------------------------------------------------------------
# Split training simulations and observed loci
# ---------------------------------------------------------------------------

training_rows = locus_data %>%
  filter(.data$origin %in% c("training_migration", "training_isolation"))
observed_loci = locus_data %>%
  filter(.data$origin == "observed_locus")

count_classes(training_rows$true_class, label = "true_class")
message("Number of observed loci: ", nrow(observed_loci))

# ---------------------------------------------------------------------------
# Select predictors, following model_comp_2pop_locus.R by column name
# ---------------------------------------------------------------------------
# Excluded by name:
#   dataset, all metadata, any std column, any pearson column,
#   minDivAB_avg, minDivAB_std, maxDivAB_avg, maxDivAB_std,
#   Gmin_avg, Gmin_std, Gmax_avg, Gmax_std.
# Also drop a statistic if it is not numeric, contains a non-finite
# value, is entirely missing, or has sd < 1e-4 in either training class
# or in all training simulations. Variance is measured on simulations
# only.

candidate_statistics = locus_predictor_names(locus_data)
excluded_by_name = c(
  "dataset",
  locus_metadata_columns,
  "minDivAB_avg",
  "minDivAB_std",
  "maxDivAB_avg",
  "maxDivAB_std",
  "Gmin_avg",
  "Gmin_std",
  "Gmax_avg",
  "Gmax_std"
)

predictor_selection = lapply(candidate_statistics, function(nm) {
  x_all = training_rows[[nm]]
  x_mig = training_rows[[nm]][training_rows$true_class == "migration"]
  x_iso = training_rows[[nm]][training_rows$true_class == "isolation"]
  sd_migration = if (is.numeric(x_mig)) stats::sd(x_mig) else NA_real_
  sd_isolation = if (is.numeric(x_iso)) stats::sd(x_iso) else NA_real_
  sd_all_training = if (is.numeric(x_all)) stats::sd(x_all) else NA_real_

  reason = "kept"
  kept = TRUE
  if (nm %in% excluded_by_name) {
    reason = "excluded_by_name"
    kept = FALSE
  } else if (grepl("std", nm, ignore.case = TRUE)) {
    reason = "std_column"
    kept = FALSE
  } else if (grepl("pearson", nm, ignore.case = TRUE)) {
    reason = "pearson_column"
    kept = FALSE
  } else if (!is.numeric(x_all)) {
    reason = "non_numeric"
    kept = FALSE
  } else if (all(is.na(x_all))) {
    reason = "all_missing"
    kept = FALSE
  } else if (any(!is.finite(x_all))) {
    reason = "non_finite"
    kept = FALSE
  } else if (is.na(sd_migration) || sd_migration < 1e-4) {
    reason = "low_sd_migration"
    kept = FALSE
  } else if (is.na(sd_isolation) || sd_isolation < 1e-4) {
    reason = "low_sd_isolation"
    kept = FALSE
  } else if (is.na(sd_all_training) || sd_all_training < 1e-4) {
    reason = "low_sd_all_training"
    kept = FALSE
  }

  tibble(
    variable = nm,
    kept = kept,
    reason = reason,
    sd_migration = sd_migration,
    sd_isolation = sd_isolation,
    sd_all_training = sd_all_training
  )
}) %>%
  bind_rows()

print(predictor_selection)
readr::write_tsv(
  predictor_selection,
  "results/locus_model_choice/10_predictor_selection.tsv"
)

kept_predictors = predictor_selection$variable[predictor_selection$kept]
if (length(kept_predictors) == 0) {
  stop("No predictor remained after filtering.", call. = FALSE)
}
message("Predictors before filtering: ", length(candidate_statistics))
message("Predictors after filtering: ", length(kept_predictors))
message("Kept: ", paste(kept_predictors, collapse = ", "))

training_predictors = training_rows[, kept_predictors, drop = FALSE]
observed_predictors = observed_loci[, kept_predictors, drop = FALSE]

observed_not_numeric = names(observed_predictors)[
  !vapply(observed_predictors, is.numeric, logical(1))
]
if (length(observed_not_numeric) > 0) {
  stop(
    "Observed loci have non-numeric values for: ",
    paste(observed_not_numeric, collapse = ", "),
    call. = FALSE
  )
}
observed_not_finite = names(observed_predictors)[
  vapply(observed_predictors, function(x) any(!is.finite(x)), logical(1))
]
if (length(observed_not_finite) > 0) {
  stop(
    "Observed loci have non-finite values for: ",
    paste(observed_not_finite, collapse = ", "),
    call. = FALSE
  )
}

if (!identical(names(training_predictors), names(observed_predictors))) {
  stop("Training and observed predictors differ.", call. = FALSE)
}

leaked_metadata = intersect(names(training_predictors), locus_metadata_columns)
if (length(leaked_metadata) > 0) {
  stop(
    "Metadata leaked into predictors: ",
    paste(leaked_metadata, collapse = ", "),
    call. = FALSE
  )
}

# ---------------------------------------------------------------------------
# Response factor and ABC-RF forest
# ---------------------------------------------------------------------------
# Predictors are not standardised. The forest is trained on raw values.

locus_class = factor(
  training_rows$true_class,
  levels = c("migration", "isolation")
)

if (!is.factor(locus_class)) {
  stop("locus_class must be a factor.", call. = FALSE)
}
if (!setequal(as.character(unique(locus_class)), c("migration", "isolation"))) {
  stop("Both migration and isolation must be present in the training classes.", call. = FALSE)
}
message("Class counts in the response factor:")
print(table(locus_class, useNA = "ifany"))

training_data = data.frame(
  locus_class,
  training_predictors,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

student_checkpoint(
  "CP-10-1",
  c(
    "What are the two training classes?",
    "Does this forest distinguish M12 = 0 from M21 = 0?"
  )
)

# ---------------------------------------------------------------------------
# Task 3 - TODO-10-3
# Train the classification forest.
# ---------------------------------------------------------------------------

rf_model = student_todo(
  "TODO-10-3",
  "Call abcrf(locus_class ~ ., data = training_data, ntree = n_classification_trees, paral = FALSE)."
)
if (!inherits(rf_model, "abcrf")) {
  stop("TODO-10-3: rf_model must be an abcrf object.", call. = FALSE)
}

oob_error = rf_model$model.rf$prediction.error
message("Out-of-bag error: ", oob_error)
message("Confusion matrix (ranger OOB):")
confusion_object = rf_model$model.rf$confusion.matrix
if (is.null(confusion_object) || length(confusion_object) == 0) {
  stop("The ranger forest did not return a confusion matrix.", call. = FALSE)
}
print(confusion_object)

confusion_table = tibble::rownames_to_column(
  as.data.frame(confusion_object),
  var = "true_class"
)
readr::write_tsv(
  confusion_table,
  "results/locus_model_choice/10_confusion_matrix.tsv"
)

class_counts = table(locus_class, useNA = "ifany")
oob_diagnostics = tibble(
  oob_error_global = as.numeric(oob_error),
  n_true_migration = as.integer(class_counts[["migration"]]),
  n_true_isolation = as.integer(class_counts[["isolation"]]),
  oob_error_migration = as.numeric(confusion_object["migration", "class.error"]),
  oob_error_isolation = as.numeric(confusion_object["isolation", "class.error"]),
  n_predictors_before_filtering = length(candidate_statistics),
  n_predictors_after_filtering = length(kept_predictors)
)
print(oob_diagnostics)
readr::write_tsv(
  oob_diagnostics,
  "results/locus_model_choice/10_oob_diagnostics.tsv"
)

importance_values = rf_model$model.rf$variable.importance
if (is.null(importance_values) || length(importance_values) == 0) {
  stop("The ranger forest did not return variable importance.", call. = FALSE)
}
importance_table = tibble(
  variable = names(importance_values),
  importance = as.numeric(importance_values)
) %>%
  arrange(dplyr::desc(.data$importance)) %>%
  mutate(
    rank = dplyr::row_number(),
    variable_type = dplyr::if_else(
      grepl("^LD[0-9]+$", .data$variable),
      "derived_LDA_axis",
      "summary_statistic"
    )
  )
print(importance_table)
readr::write_tsv(
  importance_table,
  "results/locus_model_choice/10_variable_importance.tsv"
)

importance_file = "results/locus_model_choice/10_variable_importance.png"
grDevices::png(importance_file, width = 1100, height = 800, res = 120)
abcrf::variableImpPlot(rf_model)
grDevices::dev.off()
message("Saved ", importance_file)

save_abcrf_object(rf_model, "locus_model_choice/10_locus_abcrf.rds")

# ---------------------------------------------------------------------------
# Predict each observed locus
# ---------------------------------------------------------------------------

prediction = predict(
  rf_model,
  observed_predictors,
  training = training_data,
  ntree = n_posterior_probability_trees,
  paral = FALSE
)
print_abcrf_prediction(prediction)
save_abcrf_object(prediction, "locus_model_choice/10_locus_prediction.rds")

n_observed = nrow(observed_loci)
if (length(prediction$allocation) != n_observed) {
  stop(
    "predict() returned ", length(prediction$allocation),
    " allocations for ", n_observed, " observed loci.",
    call. = FALSE
  )
}
if (length(prediction$post.prob) != n_observed) {
  stop(
    "predict() returned ", length(prediction$post.prob),
    " posterior probabilities for ", n_observed, " observed loci.",
    call. = FALSE
  )
}

allocation = as.character(prediction$allocation)
if (anyNA(allocation) || any(!nzchar(allocation))) {
  stop("Some allocations are missing.", call. = FALSE)
}
unexpected_allocation = setdiff(unique(allocation), c("migration", "isolation"))
if (length(unexpected_allocation) > 0) {
  stop(
    "Allocations other than migration and isolation were produced: ",
    paste(unexpected_allocation, collapse = ", "),
    call. = FALSE
  )
}

posterior_probability_of_allocation = as.numeric(prediction$post.prob)
if (any(!is.finite(posterior_probability_of_allocation))) {
  stop("Some post.prob values are missing or non-finite.", call. = FALSE)
}
if (any(posterior_probability_of_allocation < 0 | posterior_probability_of_allocation > 1)) {
  stop("Some post.prob values fall outside [0, 1].", call. = FALSE)
}

# post.prob is the estimated posterior probability of the allocated class.
# Taking its complement is meaningful here only because migration and
# isolation are the two mutually exclusive and exhaustive training classes.
# This conversion must not be generalized to a comparison with three or
# more classes.
estimated_p_isolation = if_else(
  allocation == "isolation",
  posterior_probability_of_allocation,
  1 - posterior_probability_of_allocation
)
estimated_p_migration = 1 - estimated_p_isolation

# allocation comes from the classification-forest votes.
# post.prob comes from a second, regression forest.
# They are therefore not required to lie on the same side of 0.5.
# ---------------------------------------------------------------------------
# Task 4 - TODO-10-4
# Define candidate_barrier from allocation, not from the 0.5 line.
# ---------------------------------------------------------------------------

candidate_barrier = student_todo(
  "TODO-10-4",
  "Set candidate_barrier to allocation == \"isolation\". Do not use estimated_p_isolation > 0.5."
)
if (!identical(as.vector(candidate_barrier), allocation == "isolation")) {
  stop(
    "TODO-10-4: candidate_barrier must be allocation == \"isolation\". ",
    "The 0.5 line is not the allocation rule.",
    call. = FALSE
  )
}

if (any(estimated_p_isolation < 0 | estimated_p_isolation > 1) ||
    any(estimated_p_migration < 0 | estimated_p_migration > 1)) {
  stop("Some estimated probabilities fall outside [0, 1].", call. = FALSE)
}
if (any(abs(estimated_p_isolation + estimated_p_migration - 1) > 1e-10)) {
  stop("estimated_p_isolation + estimated_p_migration is not 1.", call. = FALSE)
}

locus_probabilities = tibble(
  locus_id = observed_loci$locus_id,
  allocation = allocation,
  posterior_probability_of_allocation = posterior_probability_of_allocation,
  estimated_p_isolation = estimated_p_isolation,
  estimated_p_migration = estimated_p_migration,
  candidate_barrier = candidate_barrier
)

# vote is the number of trees supporting each class. The vote fraction is
# not the posterior probability estimated by the second ABC-RF forest.
if (!is.null(prediction$vote)) {
  vote = as.data.frame(prediction$vote, stringsAsFactors = FALSE)
  if ("migration" %in% names(vote) && "isolation" %in% names(vote)) {
    locus_probabilities = locus_probabilities %>%
      mutate(
        vote_migration = vote$migration,
        vote_isolation = vote$isolation,
        vote_fraction_isolation = vote_isolation / (vote_migration + vote_isolation)
      )
  }
}

interpretation_stats = intersect(
  c(
    "FST_avg",
    "netdivAB_avg",
    "divAB_avg",
    "sf_avg",
    "ss_avg",
    "piA_avg",
    "piB_avg"
  ),
  names(observed_loci)
)
locus_probabilities = bind_cols(
  locus_probabilities,
  observed_loci[, interpretation_stats, drop = FALSE]
)

print(head(locus_probabilities, 8))
readr::write_tsv(
  locus_probabilities,
  "results/locus_model_choice/10_locus_probabilities.tsv"
)

# ---------------------------------------------------------------------------
# Figures and class summary
# ---------------------------------------------------------------------------
# Loci are ordered by estimated_p_isolation, not by a genomic coordinate.
# No physical map is available in this table.
# The line at 0.5 is a posterior-equiprobability reference. It is not the
# rule used to produce allocation or candidate_barrier.

plot_data = locus_probabilities %>%
  arrange(.data$estimated_p_isolation) %>%
  mutate(
    locus_rank = dplyr::row_number(),
    allocation = factor(.data$allocation, levels = c("migration", "isolation"))
  )

equiprobability_caption = paste(
  "The line at 0.5 marks posterior equiprobability of isolation and migration.",
  "Allocation is decided by the classification forest, not by this 0.5 line."
)

p_loci = ggplot(
  plot_data,
  aes(x = .data$locus_rank, y = .data$estimated_p_isolation, colour = .data$allocation)
) +
  geom_hline(yintercept = 0.5, colour = "grey40", linewidth = 0.4) +
  geom_point(size = 1.1, alpha = 0.85) +
  scale_colour_manual(
    values = c(migration = "#2A9D8F", isolation = "#B85C38")
  ) +
  labs(
    title = "Locus-specific classification under the selected 2M model",
    subtitle = "Loci ordered by estimated P(isolation); colours are classification-forest allocations",
    x = "Loci ordered by estimated P(isolation)",
    y = "Estimated P(isolation)",
    colour = "Allocation",
    caption = equiprobability_caption
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    plot.caption = element_text(hjust = 0, size = 8)
  )

print(p_loci)
ggsave(
  "results/locus_model_choice/10_locus_probabilities.png",
  p_loci,
  width = 9,
  height = 5.5,
  dpi = 120
)

p_dist = ggplot(locus_probabilities, aes(x = .data$estimated_p_isolation)) +
  geom_histogram(binwidth = 0.05, fill = "grey55", colour = "white", boundary = 0) +
  geom_vline(xintercept = 0.5, colour = "black", linewidth = 0.4) +
  labs(
    title = "Distribution of estimated P(isolation)",
    x = "Estimated P(isolation)",
    y = "Number of loci",
    caption = equiprobability_caption
  ) +
  theme_bw() +
  theme(plot.caption = element_text(hjust = 0, size = 8))

print(p_dist)
ggsave(
  "results/locus_model_choice/10_probability_distribution.png",
  p_dist,
  width = 7,
  height = 4.5,
  dpi = 120
)

n_allocation_probability_disagreements = sum(
  (locus_probabilities$allocation == "isolation") !=
    (locus_probabilities$estimated_p_isolation > 0.5)
)

class_summary = tibble(
  n_loci = nrow(locus_probabilities),
  n_allocated_migration = sum(locus_probabilities$allocation == "migration"),
  n_allocated_isolation = sum(locus_probabilities$allocation == "isolation"),
  proportion_allocated_isolation = mean(locus_probabilities$allocation == "isolation"),
  median_p_isolation = median(locus_probabilities$estimated_p_isolation),
  min_p_isolation = min(locus_probabilities$estimated_p_isolation),
  max_p_isolation = max(locus_probabilities$estimated_p_isolation),
  n_allocation_probability_disagreements = n_allocation_probability_disagreements
)
print(class_summary)
readr::write_tsv(
  class_summary,
  "results/locus_model_choice/10_class_summary.tsv"
)

# These quantities are classification probabilities, not p-values.
# Loci are not declared "significant". No FDR is applied.
# No extra probability threshold is applied.
# candidate_barrier is exactly allocation == "isolation".

run_metadata = tibble(
  seed = workshop_seed,
  n_classification_trees = n_classification_trees,
  n_posterior_probability_trees = n_posterior_probability_trees,
  R_version = paste(R.version$major, R.version$minor, sep = "."),
  abcrf_version = as.character(utils::packageVersion("abcrf")),
  ranger_version = as.character(utils::packageVersion("ranger"))
)
print(run_metadata)
readr::write_tsv(
  run_metadata,
  "results/locus_model_choice/10_run_metadata.tsv"
)

message(
  "How to read an isolation allocation:\n",
  "- the locus summaries look more like simulations with M12 = M21 = 0;\n",
  "- the locus is a candidate associated with a barrier to gene flow;\n",
  "- this is not proof that the locus itself is causal;\n",
  "- linked selection, low recombination, data error or model misspecification can look similar;\n",
  "- the analysis is conditional on the multilocus model ", selected_complete_model,
  " and its posterior;\n",
  "- the probability expresses classification uncertainty, not a p-value;\n",
  "- a high FST is neither necessary nor sufficient to conclude that a barrier is present;\n",
  "- this classification cannot tell whether the barrier acts on M12, M21 or both."
)

# ---------------------------------------------------------------------------
# Task 5 - TODO-10-5
# List conclusions you must not draw.
# ---------------------------------------------------------------------------

forbidden_conclusions = student_todo(
  "TODO-10-5",
  "Give a character vector of overclaims you reject (causal locus, significance, FDR, both directions, FST alone)."
)
if (!is.character(forbidden_conclusions) || length(forbidden_conclusions) < 4) {
  stop("TODO-10-5: list at least four forbidden conclusions.", call. = FALSE)
}

# Interpret the output
# allocation comes from the classification forest.
# post.prob comes from the second forest.
# The 0.5 line is an equiprobability reference only.

# What this result does not show
# The locus is not shown to be causal, significant, or directional.

# Take-home message
# An isolation allocation is a candidate, conditional on the selected
# 2M model. It is not a proof of reproductive isolation.

# Optional extension
# Compare FST_avg between allocated classes without treating FST as the
# classifier.

message("Locus-specific classification finished.")
