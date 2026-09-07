# Script 01 - Explore the reference tables
#
# Biological question
# What is stored in a DILS reference table, and how does it differ from
# the observed dataset?
#
# Before running the code
# Predict: why the observed table has no true Ne, split time or migration
# rate. Write your prediction in teaching/student_answer_sheet.md.
#
# Data used here
# data/reference_summary_statistics.tsv
# data/reference_parameters.tsv
# data/observed_summary_statistics.tsv
#
# Run from the project root, section by section.

library(tidyverse)

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()
ensure_results_dir()

# ---------------------------------------------------------------------------
# Data used here
# ---------------------------------------------------------------------------

reference_ss = load_reference_summary_statistics()
reference_params = load_reference_parameters()
observed = load_observed_summary_statistics()

message("Use names() and glimpse() if you need to rediscover a column.")
message("Metadata columns provided by the workshop:")
print(metadata_columns)

glimpse(reference_ss[, seq_len(12)])
glimpse(reference_params[, seq_len(12)])
glimpse(observed[, seq_len(8)])

# ---------------------------------------------------------------------------
# Task 1 - TODO-01-1
# Classify these columns by role: metadata, parameter or summary_statistic.
# ---------------------------------------------------------------------------

column_roles = student_todo(
  "TODO-01-1",
  "Create a named character vector for simulation_id, model, FST_avg, N1 and fA2_fB1."
)

check_named_roles(
  "TODO-01-1",
  column_roles,
  expected = c(
    simulation_id = "metadata",
    model = "metadata",
    FST_avg = "summary_statistic",
    N1 = "parameter",
    fA2_fB1 = "summary_statistic"
  )
)

student_checkpoint(
  "CP-01-1",
  c(
    "Which of these five names would be allowed as an ABC-RF predictor?",
    "Why is N1 absent from the observed table?"
  )
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-01-2
# Count simulations and complete models.
# Use nrow() and n_distinct(). Do not type a remembered number.
# ---------------------------------------------------------------------------

n_simulations = student_todo(
  "TODO-01-2",
  "Set n_simulations to nrow(reference_ss)."
)
n_models = student_todo(
  "TODO-01-2",
  "Set n_models to n_distinct(reference_ss$model)."
)

if (!identical(as.integer(n_simulations), nrow(reference_ss))) {
  stop("TODO-01-2: n_simulations must equal nrow(reference_ss).", call. = FALSE)
}
if (!identical(as.integer(n_models), dplyr::n_distinct(reference_ss$model))) {
  stop("TODO-01-2: n_models must equal n_distinct(reference_ss$model).", call. = FALSE)
}

message("Simulations per complete model:")
print(count(reference_ss, model) %>% arrange(model))
print(table(reference_ss$migration_status, reference_ss$scenario, useNA = "ifany"))

student_checkpoint(
  "CP-01-2",
  c(
    "Is one row one locus or one multilocus simulation?",
    "Why does SI appear with NA in genomic_migration?"
  )
)

# ---------------------------------------------------------------------------
# Task 3 - TODO-01-3
# Compare simulated predictors with the observed table.
# ---------------------------------------------------------------------------

simulated_predictors = student_todo(
  "TODO-01-3",
  "Set simulated_predictors to predictor_names(reference_ss)."
)
observed_predictors = student_todo(
  "TODO-01-3",
  "Set observed_predictors to names(observed)."
)

if (!identical(simulated_predictors, predictor_names(reference_ss))) {
  stop("TODO-01-3: use predictor_names(reference_ss) for the simulated predictors.", call. = FALSE)
}
if (!identical(observed_predictors, names(observed))) {
  stop("TODO-01-3: use names(observed) for the observed predictors.", call. = FALSE)
}
if (!identical(simulated_predictors, observed_predictors)) {
  stop("The simulated predictors and the observed columns should be identical.", call. = FALSE)
}
message("Simulated and observed predictor names match, in the same order.")
message("Classic statistics: ", length(classic_statistics))
message("jSFS cells: ", length(jsfs_names(reference_ss)))

if (!identical(reference_ss$simulation_id, reference_params$simulation_id)) {
  stop("simulation_id does not match between the two reference tables.", call. = FALSE)
}
message("simulation_id pairs each statistic row with its parameter row.")

student_checkpoint(
  "CP-01-3",
  c(
    "What does fA2_fB1 mean in a folded jSFS with nMin = 6?",
    "Why were fA0_fB0 and fA6_fB6 removed?"
  )
)

# ---------------------------------------------------------------------------
# Task 4 - TODO-01-4
# Write your own definition of a reference table.
# ---------------------------------------------------------------------------

reference_table_definition = student_todo(
  "TODO-01-4",
  "Replace this call with one character string of your own definition."
)
if (!is.character(reference_table_definition) || length(reference_table_definition) != 1) {
  stop("TODO-01-4: provide a single character string.", call. = FALSE)
}
if (nchar(trimws(reference_table_definition)) < 40) {
  stop("TODO-01-4: write at least one full sentence.", call. = FALSE)
}
message("Recorded definition:")
message(reference_table_definition)

# ---------------------------------------------------------------------------
# Interpret the output
# ---------------------------------------------------------------------------
# Each row is one multilocus simulation, not one locus.
# Metadata label the simulation. Parameters were used to simulate it.
# Summary statistics are computed from the simulated sequences.
# The observed table has statistics only, because the true parameters of
# the real populations are unknown.

# What this result does not show
# Exploring the tables does not identify the best demographic model.

# Take-home message
# ABC compares observed summaries to simulated summaries. Parameters stay
# on the simulation side until a later estimation step.

# Optional extension
# Draw netdivAB_avg against FST_avg and colour points by scenario.
# Keep axes labelled in words, not only by colour.

p_div = ggplot(reference_ss, aes(x = netdivAB_avg, y = FST_avg, colour = scenario)) +
  geom_point(alpha = 0.08, size = 0.4) +
  labs(
    title = "Net divergence and FST in the reference table",
    subtitle = "Each point is one multilocus simulation",
    x = "Mean net divergence (Da)",
    y = "Mean FST"
  ) +
  theme_bw() +
  guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2)))
print(p_div)
ggsave("results/01_netdiv_vs_fst.png", p_div, width = 7.5, height = 5.5, dpi = 120)

message("Exploration finished. No ABC-RF model was trained.")
