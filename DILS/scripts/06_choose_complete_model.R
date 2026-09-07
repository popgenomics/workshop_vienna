# Script 06 - Assemble the complete elementary model
#
# Biological question
# Which of the 14 DILS elementary models is implied by the previous
# binary decisions?
#
# Before running the code
# This script trains no forest. Predict why assembling four binary
# choices is not the same as training one 14-class forest.
#
# Data used here
# No summary-statistic table. Only the allocations from scripts 02-05.
#
# Run from the project root, section by section.

source("R/helpers.R")
source("R/teaching_helpers.R")
stop_if_not_project_root()

# ---------------------------------------------------------------------------
# Task 1 - TODO-06-1
# Fill the four previous allocations.
# ---------------------------------------------------------------------------

selected_migration_status = student_todo(
  "TODO-06-1",
  "Copy \"migration\" or \"isolation\" from script 02."
)
selected_scenario = student_todo(
  "TODO-06-1",
  "Copy IM, SC, AM or SI from script 03."
)
selected_genomic_Ne = student_todo(
  "TODO-06-1",
  "Copy \"homogeneous\" or \"heterogeneous\" from script 04."
)
selected_genomic_migration = student_todo(
  "TODO-06-1",
  "Copy \"homogeneous\" or \"heterogeneous\" from script 05. Use NA_character_ if the scenario is SI."
)

check_student_choice(
  "TODO-06-1",
  selected_migration_status,
  allowed = c("migration", "isolation"),
  label = "migration status"
)
check_student_choice(
  "TODO-06-1",
  selected_scenario,
  allowed = c("IM", "SC", "AM", "SI"),
  label = "scenario"
)
check_student_choice(
  "TODO-06-1",
  selected_genomic_Ne,
  allowed = c("homogeneous", "heterogeneous"),
  label = "Ne class"
)

if (selected_migration_status == "migration" && !selected_scenario %in% c("IM", "SC")) {
  stop("Under ongoing migration, selected_scenario must be IM or SC.", call. = FALSE)
}
if (selected_migration_status == "isolation" && !selected_scenario %in% c("AM", "SI")) {
  stop("Under isolation, selected_scenario must be AM or SI.", call. = FALSE)
}

if (selected_scenario == "SI") {
  selected_genomic_migration = NA_character_
} else {
  check_student_choice(
    "TODO-06-1",
    selected_genomic_migration,
    allowed = c("homogeneous", "heterogeneous"),
    label = "migration class"
  )
}

student_checkpoint(
  "CP-06-1",
  c(
    "Are these four values mutually compatible?",
    "Which label is unused if the scenario is SI?"
  )
)

# ---------------------------------------------------------------------------
# Task 2 - TODO-06-2
# Build the elementary model name.
# ---------------------------------------------------------------------------
# Use genomic_label_to_code() if you want a helper:
#   homogeneous -> 1N or 1M
#   heterogeneous -> 2N or 2M
# SI has the form SI_1N or SI_2N.
# Other scenarios have the form SC_2M_2N, IM_1M_1N, ...

selected_complete_model = student_todo(
  "TODO-06-2",
  "Paste scenario, migration code and Ne code into one of the 14 allowed names."
)
check_student_choice(
  "TODO-06-2",
  selected_complete_model,
  allowed = allowed_complete_models,
  label = "complete DILS model"
)

code_n = genomic_label_to_code(selected_genomic_Ne, "N")
if (selected_scenario == "SI") {
  expected_name = paste("SI", code_n, sep = "_")
} else {
  code_m = genomic_label_to_code(selected_genomic_migration, "M")
  expected_name = paste(selected_scenario, code_m, code_n, sep = "_")
}
if (!identical(selected_complete_model, expected_name)) {
  stop(
    "TODO-06-2: the name does not match the four allocations. Expected ",
    expected_name, ".",
    call. = FALSE
  )
}

message("Complete DILS model: ", selected_complete_model)

# ---------------------------------------------------------------------------
# Task 3 - TODO-06-3
# Explain the pedagogical construction.
# ---------------------------------------------------------------------------

why_not_a_14_class_forest = student_todo(
  "TODO-06-3",
  "Explain why this name is not mathematically identical to a direct 14-class ABC-RF allocation."
)
if (!is.character(why_not_a_14_class_forest) || nchar(trimws(why_not_a_14_class_forest)) < 40) {
  stop("TODO-06-3: write a short paragraph.", call. = FALSE)
}

# Interpret the output
# You assembled a valid DILS elementary model from four biological
# questions. No new forest was trained.

# What this result does not show
# The name is not guaranteed to match a forest trained on all 14 models
# at once.

# Take-home message
# The teaching sequence keeps each biological question explicit.

# Optional extension
# Write two other valid names from abstract combinations, for example
# AM + heterogeneous migration + homogeneous Ne.

message("Copy this name into selected_complete_model in script 07.")
message("No random forest was trained in this script.")
