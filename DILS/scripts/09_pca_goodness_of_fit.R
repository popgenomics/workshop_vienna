# Script 09 - PCA goodness of fit
#
# Biological question
# Can the selected and fitted model reproduce the observed summaries
# jointly, not one statistic at a time?
#
# Before running the code
# Should the observed dataset be used to define the PCA axes?
#
# Data used here
# data/goodness_of_fit_data.tsv. Axes are fitted on simulations only.
#
# Run from the project root, section by section.

library(tidyverse)

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()
ensure_results_dir("goodness_of_fit")

gof_data = load_goodness_of_fit_data()
predictor_cols = gof_predictor_names(gof_data)
fitted_label = paste(unique(gof_data$fitted_model), collapse = ", ")

# ---------------------------------------------------------------------------
# Task 1 - TODO-09-1
# Decide whether the observation may define the axes.
# ---------------------------------------------------------------------------

observation_defines_axes = student_todo(
  "TODO-09-1",
  "Set this to FALSE or TRUE, then be ready to justify the choice."
)
if (!identical(observation_defines_axes, FALSE)) {
  stop(
    "TODO-09-1: the observation must not define the PCA axes. Set this to FALSE.",
    call. = FALSE
  )
}

simulation_data = gof_data %>%
  filter(.data$origin != "observed")
observed_data = gof_data %>%
  filter(.data$origin == "observed")
if (nrow(observed_data) != 1) {
  stop("Expected exactly one observed row.", call. = FALSE)
}

simulation_predictors = simulation_data %>%
  select(all_of(predictor_cols))
observed_predictors = observed_data %>%
  select(all_of(predictor_cols))

# The next block is provided. It drops non-informative predictors using
# simulations only. You do not rewrite it.

diagnose_predictor = function(name) {
  sim_x = simulation_predictors[[name]]
  obs_x = observed_predictors[[name]][1]
  simulation_mean = NA_real_
  simulation_median = NA_real_
  simulation_sd = NA_real_
  if (is.numeric(sim_x)) {
    simulation_mean = mean(sim_x, na.rm = TRUE)
    simulation_median = median(sim_x, na.rm = TRUE)
    simulation_sd = stats::sd(sim_x)
  }
  observed_value = if (is.numeric(obs_x)) as.numeric(obs_x) else NA_real_
  absolute_difference = if (is.finite(simulation_mean) && is.finite(observed_value)) {
    abs(observed_value - simulation_mean)
  } else {
    NA_real_
  }
  reason = NA_character_
  if (!is.numeric(sim_x)) {
    reason = "non_numeric"
  } else if (all(is.na(sim_x))) {
    reason = "all_missing"
  } else if (any(!is.finite(sim_x))) {
    reason = "non_finite"
  } else if (is.na(simulation_sd) || simulation_sd == 0) {
    reason = "zero_variance"
  } else if (simulation_sd < 1e-5) {
    reason = "near_zero_variance"
  }
  list(
    keep = is.na(reason),
    row = tibble(
      variable = name,
      reason = reason,
      simulation_mean = simulation_mean,
      simulation_median = simulation_median,
      simulation_sd = simulation_sd,
      observed_value = observed_value,
      absolute_difference = absolute_difference
    )
  )
}

diagnostics = lapply(predictor_cols, diagnose_predictor)
keep = vapply(diagnostics, function(x) x$keep, logical(1))
dropped_predictors = bind_rows(lapply(diagnostics[!keep], function(x) x$row))
if (nrow(dropped_predictors) > 0) {
  machine_tol = .Machine$double.eps * 100
  constant_mismatch = dropped_predictors %>%
    filter(.data$reason == "zero_variance") %>%
    filter(is.finite(.data$observed_value), is.finite(.data$simulation_mean)) %>%
    filter(
      abs(.data$observed_value - .data$simulation_mean) >
        pmax(machine_tol, abs(.data$simulation_mean) * machine_tol)
    )
  if (nrow(constant_mismatch) > 0) {
    stop("A constant simulated predictor disagrees with the observation.", call. = FALSE)
  }
  readr::write_tsv(dropped_predictors, "results/goodness_of_fit/09_pca_dropped_predictors.tsv")
  message("Removed ", nrow(dropped_predictors), " non-informative predictor(s).")
}
simulation_predictors = simulation_predictors[, keep, drop = FALSE]
observed_predictors = observed_predictors[, keep, drop = FALSE]
message("Predictors after filtering: ", ncol(simulation_predictors))

student_checkpoint(
  "CP-09-1",
  c(
    "Were dropped variables chosen using the observation?",
    "Where does standardisation happen: in the table, or in prcomp()?"
  )
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-09-2
# Fit PCA on simulations only.
# ---------------------------------------------------------------------------

pca_model = student_todo(
  "TODO-09-2",
  "Call prcomp(simulation_predictors, center = TRUE, scale. = TRUE)."
)
if (!inherits(pca_model, "prcomp")) {
  stop("TODO-09-2: pca_model must be a prcomp object.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Task 3 - TODO-09-3
# Project the observation with predict(), not by refitting.
# ---------------------------------------------------------------------------

observed_coordinates = student_todo(
  "TODO-09-3",
  "Project the observation with predict(pca_model, newdata = observed_predictors)."
)
observed_coordinates = as.data.frame(observed_coordinates)
simulation_coordinates = as.data.frame(pca_model$x)

pca_coordinates = bind_rows(
  bind_cols(simulation_data %>% select(origin, dataset_id), simulation_coordinates),
  bind_cols(observed_data %>% select(origin, dataset_id), observed_coordinates)
)
eigenvalues = pca_model$sdev^2
pca_variance = tibble(
  component = paste0("PC", seq_along(eigenvalues)),
  standard_deviation = pca_model$sdev,
  proportion_of_variance = eigenvalues / sum(eigenvalues),
  cumulative_proportion = cumsum(eigenvalues / sum(eigenvalues))
)
print(head(pca_variance, 8))
readr::write_tsv(pca_coordinates, "results/goodness_of_fit/09_pca_coordinates.tsv")
readr::write_tsv(pca_variance, "results/goodness_of_fit/09_pca_variance.tsv")

pc1_percent = 100 * pca_variance$proportion_of_variance[[1]]
pc2_percent = 100 * pca_variance$proportion_of_variance[[2]]
plot_data = pca_coordinates %>%
  select(origin, dataset_id, PC1, PC2)
prior_points = plot_data %>%
  filter(.data$origin == "prior_predictive")
posterior_points = plot_data %>%
  filter(.data$origin == "posterior_predictive")
observed_point = plot_data %>%
  filter(.data$origin == "observed")

p_pca = ggplot() +
  geom_point(
    data = prior_points,
    aes(x = PC1, y = PC2),
    colour = "grey70",
    alpha = 0.15,
    size = 0.5,
    stroke = 0
  ) +
  geom_point(
    data = posterior_points,
    aes(x = PC1, y = PC2),
    colour = "#2A9D8F",
    alpha = 0.20,
    size = 0.5,
    stroke = 0
  ) +
  geom_point(
    data = observed_point,
    aes(x = PC1, y = PC2),
    colour = "black",
    shape = 8,
    size = 1.5
  ) +
  labs(
    title = paste("PCA goodness of fit under", fitted_label),
    subtitle = "Grey = prior predictive; turquoise = posterior predictive; black star = observed",
    x = sprintf("PC1 (%.1f%% of variance)", pc1_percent),
    y = sprintf("PC2 (%.1f%% of variance)", pc2_percent)
  ) +
  theme_bw()
print(p_pca)
ggsave(
  "results/goodness_of_fit/09_pca_goodness_of_fit.png",
  p_pca,
  width = 8,
  height = 6.5,
  dpi = 120
)

# ---------------------------------------------------------------------------
# Task 4 - TODO-09-4
# Reject a global-validation overclaim.
# ---------------------------------------------------------------------------

why_pc1_pc2_is_not_proof = student_todo(
  "TODO-09-4",
  "Explain why an observation near the centre of PC1-PC2 does not prove that the model is true."
)
if (!is.character(why_pc1_pc2_is_not_proof) || nchar(trimws(why_pc1_pc2_is_not_proof)) < 40) {
  stop("TODO-09-4: write a short paragraph.", call. = FALSE)
}

# Interpret the output
# Near the centre: no obvious mismatch on these two axes.
# On the periphery: imperfect fit.
# Outside the cloud: the model cannot jointly reproduce the summaries.

# What this result does not show
# Proximity on PC1 and PC2 is not a proof that the model is true.

# Take-home message
# PCA is a multivariate posterior predictive check. It complements, and
# does not replace, the univariate plots.

# Optional extension
# Read the dropped-predictor table and say which jSFS cells were empty.
