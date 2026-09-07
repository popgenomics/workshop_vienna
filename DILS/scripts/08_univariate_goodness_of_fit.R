# Script 08 - Univariate goodness of fit
#
# Biological question
# Can the selected and fitted model reproduce each observed summary
# statistic, one at a time?
#
# Before running the code
# Model choice is relative. Goodness of fit asks whether the fitted
# model can reproduce the data.
#
# Data used here
# data/goodness_of_fit_data.tsv
# Origins: prior_predictive, posterior_predictive, observed.
#
# Run from the project root, section by section.

library(tidyverse)

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()
ensure_results_dir("goodness_of_fit")

gof_data = load_goodness_of_fit_data()
print(table(gof_data$origin, useNA = "ifany"))
print(classic_statistics)

# ---------------------------------------------------------------------------
# Task 1 - TODO-08-1
# State the two questions in your own words.
# ---------------------------------------------------------------------------

relative_versus_absolute = student_todo(
  "TODO-08-1",
  "Write two sentences: model choice is relative; goodness of fit is absolute."
)
if (!is.character(relative_versus_absolute) || nchar(trimws(relative_versus_absolute)) < 40) {
  stop("TODO-08-1: write two sentences.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Task 2 - TODO-08-2
# Choose which classic statistics to plot.
# Use names from classic_statistics. You may choose a short list.
# ---------------------------------------------------------------------------

statistics_to_plot = student_todo(
  "TODO-08-2",
  "Give a character vector of names from classic_statistics."
)
if (!is.character(statistics_to_plot) || length(statistics_to_plot) < 2) {
  stop("TODO-08-2: choose at least two classic statistics.", call. = FALSE)
}
missing_chosen = setdiff(statistics_to_plot, classic_statistics)
if (length(missing_chosen) > 0) {
  stop("TODO-08-2: unknown names: ", paste(missing_chosen, collapse = ", "), call. = FALSE)
}

student_checkpoint(
  "CP-08-1",
  c(
    "Which origin is prior predictive?",
    "How is the observation drawn: density or vertical line?"
  )
)

gof_colours = c(
  prior_predictive = "grey55",
  posterior_predictive = "#2A9D8F",
  observed = "black"
)
fitted_label = paste(unique(gof_data$fitted_model), collapse = ", ")

simulated_long = gof_data %>%
  filter(.data$origin != "observed") %>%
  select(origin, all_of(statistics_to_plot)) %>%
  pivot_longer(
    cols = all_of(statistics_to_plot),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(statistic = factor(.data$statistic, levels = statistics_to_plot))

observed_long = gof_data %>%
  filter(.data$origin == "observed") %>%
  select(all_of(statistics_to_plot)) %>%
  pivot_longer(
    cols = all_of(statistics_to_plot),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  mutate(statistic = factor(.data$statistic, levels = statistics_to_plot))

p_univariate = ggplot(
  simulated_long,
  aes(x = value, colour = origin, fill = origin)
) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  geom_vline(
    data = observed_long,
    aes(xintercept = value),
    colour = "black",
    linewidth = 0.7
  ) +
  facet_wrap(~ statistic, scales = "free") +
  scale_colour_manual(values = gof_colours) +
  scale_fill_manual(values = gof_colours) +
  labs(
    title = paste("Univariate goodness of fit under", fitted_label),
    subtitle = "Grey = prior predictive; turquoise = posterior predictive; black line = observed",
    x = "Summary statistic (raw value)",
    y = "Density",
    colour = "Origin",
    fill = "Origin"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")
print(p_univariate)

n_plot = length(statistics_to_plot)
n_rows = ceiling(n_plot / 2)
ggsave(
  "results/goodness_of_fit/08_univariate_distributions.png",
  p_univariate,
  width = 10,
  height = max(4, 3 * n_rows),
  dpi = 120,
  limitsize = FALSE
)

# ---------------------------------------------------------------------------
# Task 3 - TODO-08-3
# Identify the three visual elements.
# ---------------------------------------------------------------------------

what_grey_density_is = student_todo(
  "TODO-08-3",
  "One sentence: what the grey density represents."
)
what_black_line_is = student_todo(
  "TODO-08-3",
  "One sentence: why the observation is a line, never a density."
)
if (!is.character(what_grey_density_is) || nchar(trimws(what_grey_density_is)) < 15) {
  stop("TODO-08-3: describe the grey density.", call. = FALSE)
}
if (!is.character(what_black_line_is) || nchar(trimws(what_black_line_is)) < 15) {
  stop("TODO-08-3: describe the black line.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# Task 4 - TODO-08-4
# Rebuild dils_tail_probability.
# If observation equals the simulated median, return 0.5.
# If observation is above the median, return the proportion of simulations
# strictly above the observation.
# Otherwise return the proportion strictly below it.
# ---------------------------------------------------------------------------

dils_tail_probability = student_todo(
  "TODO-08-4",
  "Replace this with a function(distribution, observation) that follows the DILS rule above."
)
if (!is.function(dils_tail_probability)) {
  stop("TODO-08-4: provide a function.", call. = FALSE)
}

posterior_data = gof_data %>%
  filter(.data$origin == "posterior_predictive")
observed_row = gof_data %>%
  filter(.data$origin == "observed")

univariate_gof = lapply(classic_statistics, function(stat) {
  distribution = posterior_data[[stat]]
  observation = observed_row[[stat]]
  tibble(
    statistic = stat,
    observed = observation,
    posterior_mean = mean(distribution),
    posterior_median = median(distribution),
    quantile_2.5 = quantile(distribution, 0.025, names = FALSE),
    quantile_97.5 = quantile(distribution, 0.975, names = FALSE),
    p_lower_tail = mean(distribution < observation),
    p_upper_tail = mean(distribution > observation),
    p_tail = dils_tail_probability(distribution, observation)
  )
}) %>%
  bind_rows()

univariate_gof = univariate_gof %>%
  mutate(p_fdr = p.adjust(p_tail, method = "fdr")) %>%
  mutate(
    across(
      c(p_lower_tail, p_upper_tail, p_tail, p_fdr),
      function(x) round(x, 5)
    )
  )
print(univariate_gof)
readr::write_tsv(univariate_gof, "results/goodness_of_fit/08_univariate_gof.tsv")

# ---------------------------------------------------------------------------
# Task 5 - TODO-08-5
# Interpret p_tail without treating it as a classical p-value.
# ---------------------------------------------------------------------------

p_tail_limits = student_todo(
  "TODO-08-5",
  "State that p_tail follows DILS, is not a two-sided p-value, never validates the model when large, and is FDR-corrected on the 30 classic statistics."
)
if (!is.character(p_tail_limits) || nchar(trimws(p_tail_limits)) < 40) {
  stop("TODO-08-5: write a short paragraph.", call. = FALSE)
}

# Interpret the output
# Look for posterior clouds that move toward the black line, and for
# statistics that remain in a tail.

# What this result does not show
# A large p_tail does not prove that the model is true.
# FDR here is not a locus-significance test.

# Take-home message
# Model choice picks a winner among candidates. Goodness of fit asks
# whether that winner can reproduce the observed summaries.

# Optional extension
# Set statistics_to_plot = classic_statistics and inspect every panel.
