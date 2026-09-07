# Small helpers for the guided student scripts.
# These functions never contain empirical workshop answers.

student_todo = function(task_id, what_to_complete = NULL) {
  extra = if (is.null(what_to_complete) || !nzchar(what_to_complete)) {
    "Open the matching teaching page, then replace this call with your choice or code."
  } else {
    what_to_complete
  }
  stop(
    "Incomplete step ", task_id, ".\n",
    extra, "\n",
    "The script parsed, but this scientific choice is still empty. ",
    "Complete the TODO, then rerun this section.",
    call. = FALSE
  )
}

student_checkpoint = function(checkpoint_id, questions = character()) {
  message("CHECKPOINT ", checkpoint_id)
  if (length(questions) > 0) {
    message(paste0("- ", questions, collapse = "\n"))
  }
  message("Write short answers in teaching/student_answer_sheet.md before continuing.")
  invisible(checkpoint_id)
}

check_student_choice = function(task_id, value, allowed, label = "value") {
  if (length(value) != 1) {
    stop(
      task_id, ": please provide a single ", label, ".",
      call. = FALSE
    )
  }
  value_chr = as.character(value)
  if (is.na(value_chr) || !nzchar(value_chr) || identical(value_chr, "TO_COMPLETE")) {
    stop(
      task_id, ": replace TO_COMPLETE or student_todo() with a ", label, ".",
      call. = FALSE
    )
  }
  if (!value_chr %in% allowed) {
    stop(
      task_id, ": ", label, " must be one of: ",
      paste(allowed, collapse = ", "),
      ". Received: ", value_chr, ".",
      call. = FALSE
    )
  }
  invisible(value_chr)
}

check_named_roles = function(task_id, assigned, expected) {
  missing_names = setdiff(names(expected), names(assigned))
  if (length(missing_names) > 0) {
    stop(
      task_id, ": missing names: ",
      paste(missing_names, collapse = ", "),
      call. = FALSE
    )
  }
  mismatches = names(expected)[assigned[names(expected)] != expected]
  if (length(mismatches) > 0) {
    stop(
      task_id, ": at least one column has the wrong role: ",
      paste(mismatches, collapse = ", "),
      ". Use metadata, parameter or summary_statistic.",
      call. = FALSE
    )
  }
  message("OK  ", task_id, ": column roles match the expected scientific categories.")
  invisible(TRUE)
}
