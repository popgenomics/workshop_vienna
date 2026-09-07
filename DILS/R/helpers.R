# Shared constants and small helper functions for the DILS workshop.
# Student scripts should still show the ABC-RF steps explicitly.
# Source this file from the project root:
#   source("R/helpers.R")

metadata_columns = c(
  "migration_status",
  "genomic_migration",
  "genomic_Ne",
  "scenario",
  "model",
  "simulation_id",
  "batch",
  "row_in_batch"
)

classic_statistics = c(
  "bialsites_avg",
  "bialsites_std",
  "sf_avg",
  "sf_std",
  "sxA_avg",
  "sxA_std",
  "sxB_avg",
  "sxB_std",
  "ss_avg",
  "ss_std",
  "piA_avg",
  "piA_std",
  "piB_avg",
  "piB_std",
  "pearson_r_pi",
  "thetaA_avg",
  "thetaA_std",
  "thetaB_avg",
  "thetaB_std",
  "pearson_r_theta",
  "DtajA_avg",
  "DtajA_std",
  "DtajB_avg",
  "DtajB_std",
  "divAB_avg",
  "divAB_std",
  "netdivAB_avg",
  "netdivAB_std",
  "FST_avg",
  "FST_std"
)

jsfs_name_pattern = "^fA[0-9]+_fB[0-9]+$"
jsfs_excluded_cells = c("fA0_fB0", "fA6_fB6")

allowed_complete_models = c(
  "SC_1M_1N", "SC_1M_2N", "SC_2M_1N", "SC_2M_2N",
  "AM_1M_1N", "AM_1M_2N", "AM_2M_1N", "AM_2M_2N",
  "IM_1M_1N", "IM_1M_2N", "IM_2M_1N", "IM_2M_2N",
  "SI_1N", "SI_2N"
)

gof_metadata_columns = c(
  "origin",
  "fitted_model",
  "dataset_id",
  "batch",
  "row_in_batch"
)

gof_fitted_model = "SC_2M_2N"

locus_metadata_columns = c(
  "origin",
  "fitted_model",
  "dataset_id",
  "locus_id",
  "true_class",
  "row_in_source"
)

# Number of trees for ABC-RF. 500 is a reasonable workshop default.
# Lower it (for example to 50) if you only want a quick test.
# Avoid very small values such as 10: with current ranger/abcrf versions,
# some out-of-bag predictions can be missing and abcrf() then fails.
n_trees = 500
if (nzchar(Sys.getenv("DILS_WORKSHOP_NTREES"))) {
  n_trees = as.integer(Sys.getenv("DILS_WORKSHOP_NTREES"))
}
if (
  length(n_trees) != 1 ||
    is.na(n_trees) ||
    !is.finite(n_trees) ||
    n_trees < 50 ||
    n_trees != as.integer(n_trees)
) {
  stop(
    "n_trees must be a single finite integer of at least 50. Received: ",
    paste(n_trees, collapse = ", "),
    call. = FALSE
  )
}
n_trees = as.integer(n_trees)

stop_if_not_project_root = function() {
  required = c(
    "data/reference_summary_statistics.tsv",
    "data/reference_parameters.tsv",
    "data/observed_summary_statistics.tsv",
    "R/helpers.R"
  )
  missing = required[!file.exists(required)]
  if (length(missing) > 0) {
    stop(
      "Please set the working directory to the DILS workshop project root. ",
      "Missing: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
}

ensure_results_dir = function(...) {
  path = file.path("results", ...)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_workshop_tsv = function(path) {
  readr::read_tsv(path, na = c("NA", ""), show_col_types = FALSE)
}

load_reference_summary_statistics = function() {
  read_workshop_tsv("data/reference_summary_statistics.tsv")
}

load_reference_parameters = function() {
  read_workshop_tsv("data/reference_parameters.tsv")
}

load_observed_summary_statistics = function() {
  read_workshop_tsv("data/observed_summary_statistics.tsv")
}

load_goodness_of_fit_data = function() {
  path = "data/goodness_of_fit_data.tsv"
  if (!file.exists(path)) {
    stop(
      "Cannot find data/goodness_of_fit_data.tsv. ",
      "Build it with preparation/02_build_goodness_of_fit_table.R.",
      call. = FALSE
    )
  }
  read_workshop_tsv(path)
}

load_locus_model_choice_data = function() {
  path = "data/locus_model_choice_data.tsv"
  if (!file.exists(path)) {
    stop(
      "Cannot find data/locus_model_choice_data.tsv. ",
      "Build it with preparation/04_build_locus_model_choice_table.R.",
      call. = FALSE
    )
  }
  read_workshop_tsv(path)
}

predictor_names = function(reference_table) {
  setdiff(names(reference_table), metadata_columns)
}

gof_predictor_names = function(gof_table) {
  setdiff(names(gof_table), gof_metadata_columns)
}

locus_predictor_names = function(locus_table) {
  setdiff(names(locus_table), locus_metadata_columns)
}

jsfs_names = function(table) {
  cols = grep(jsfs_name_pattern, names(table), value = TRUE)
  setdiff(cols, jsfs_excluded_cells)
}

count_classes = function(x, label = "class") {
  counts = table(x, useNA = "ifany")
  message("Number of simulations per ", label, ":")
  print(counts)
  invisible(counts)
}

drop_uninformative_predictors = function(training_predictors, observed_predictors) {
  training_predictors = as.data.frame(training_predictors, stringsAsFactors = FALSE)
  observed_predictors = as.data.frame(observed_predictors, stringsAsFactors = FALSE)

  if (!identical(names(training_predictors), names(observed_predictors))) {
    stop(
      "Training and observed predictor names differ. ",
      "They must be identical and in the same order before ABC-RF.",
      call. = FALSE
    )
  }

  keep = vapply(names(training_predictors), function(nm) {
    x = training_predictors[[nm]]
    if (!is.numeric(x)) {
      return(FALSE)
    }
    if (all(is.na(x))) {
      return(FALSE)
    }
    if (any(!is.finite(x))) {
      return(FALSE)
    }
    sdev = stats::sd(x)
    if (is.na(sdev) || sdev < 1e-5) {
      return(FALSE)
    }
    TRUE
  }, logical(1))

  dropped = names(training_predictors)[!keep]
  if (length(dropped) > 0) {
    message(
      "Removed ", length(dropped), " non-informative predictor(s): ",
      paste(dropped, collapse = ", ")
    )
  } else {
    message("No predictor was removed for this comparison.")
  }

  list(
    training = training_predictors[, keep, drop = FALSE],
    observed = observed_predictors[, keep, drop = FALSE],
    dropped = dropped
  )
}

prepare_abcrf_tables = function(reference_subset, observed, response_name) {
  if (!response_name %in% names(reference_subset)) {
    stop("Cannot find the response column '", response_name, "'.", call. = FALSE)
  }

  pred_names = predictor_names(reference_subset)
  if (!identical(pred_names, names(observed))) {
    stop(
      "Observed columns are not identical to the numeric predictors ",
      "of the reference table.",
      call. = FALSE
    )
  }

  cleaned = drop_uninformative_predictors(
    training_predictors = reference_subset[, pred_names, drop = FALSE],
    observed_predictors = observed[, pred_names, drop = FALSE]
  )

  response = reference_subset[[response_name]]
  # abcrf() requires a factor response. A character vector is not enough.
  response = factor(response)

  training_data = data.frame(
    response,
    cleaned$training,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  names(training_data)[1] = response_name

  leaked = intersect(names(training_data)[-1], metadata_columns)
  if (length(leaked) > 0) {
    stop(
      "Metadata columns were accidentally included among predictors: ",
      paste(leaked, collapse = ", "),
      call. = FALSE
    )
  }

  numeric_ok = vapply(training_data[, -1, drop = FALSE], is.numeric, logical(1))
  if (!all(numeric_ok)) {
    stop(
      "Non-numeric predictors found: ",
      paste(names(which(!numeric_ok)), collapse = ", "),
      call. = FALSE
    )
  }

  if (!is.factor(training_data[[response_name]])) {
    stop("The response variable is not a factor.", call. = FALSE)
  }

  if (!identical(names(cleaned$observed), names(cleaned$training))) {
    stop("Observed and training predictors differ after filtering.", call. = FALSE)
  }

  list(
    training_data = training_data,
    observed_predictors = cleaned$observed,
    dropped = cleaned$dropped,
    response_name = response_name
  )
}

print_abcrf_prediction = function(prediction) {
  if (!inherits(prediction, "abcrfpredict")) {
    stop(
      "predict() applied to an abcrf object should return an object ",
      "of class abcrfpredict.",
      call. = FALSE
    )
  }

  cat("\nSelected class (allocation):\n")
  print(prediction$allocation)

  cat("\nVotes of the classification forest (vote):\n")
  print(prediction$vote)

  cat("\nEstimated posterior probability of the selected class (post.prob):\n")
  print(prediction$post.prob)

  cat("\nHow to read these numbers:\n")
  cat("- allocation is the class chosen by the classification forest.\n")
  cat("- vote contains the number of trees supporting each class.\n")
  cat("- post.prob estimates the posterior probability of the selected class.\n")
  cat("- post.prob is not simply the largest vote share. abcrf estimates it\n")
  cat("  with a second, regression random forest trained on out-of-bag errors.\n")
}

save_abcrf_object = function(object, filename) {
  ensure_results_dir()
  path = file.path("results", filename)
  saveRDS(object, path)
  message("Saved ", path)
  path
}

# Variable importance and LDA projection without the interactive pause
# used by plot.abcrf() (which waits for Enter).
plot_abcrf_diagnostics = function(rf_model,
                                   training_data,
                                   observed_predictors,
                                   response_name,
                                   file_prefix) {
  ensure_results_dir()

  importance_file = file.path("results", paste0(file_prefix, "_variable_importance.png"))
  grDevices::png(importance_file, width = 1100, height = 800, res = 120)
  abcrf::variableImpPlot(rf_model)
  grDevices::dev.off()
  message("Saved ", importance_file)

  if (!isTRUE(rf_model$lda) || is.null(rf_model$model.lda)) {
    message("LDA projection is not available for this forest.")
    return(invisible(NULL))
  }

  proj_train = as.data.frame(predict(rf_model$model.lda, training_data)$x)
  proj_obs = as.data.frame(predict(rf_model$model.lda, observed_predictors)$x)
  proj_train$class = training_data[[response_name]]
  n_ld = ncol(proj_obs)

  lda_file = file.path("results", paste0(file_prefix, "_lda.png"))

  if (n_ld == 1) {
    # Two classes: only LD1 is produced. Show class densities plus the observation.
    p = ggplot2::ggplot(proj_train, ggplot2::aes(x = .data$LD1, fill = .data$class)) +
      ggplot2::geom_density(alpha = 0.35, colour = NA) +
      ggplot2::geom_vline(
        xintercept = proj_obs$LD1[1],
        colour = "black",
        linewidth = 0.8
      ) +
      ggplot2::labs(
        title = "LDA projection of the reference table",
        subtitle = "The vertical line is the observed dataset",
        x = "LD1",
        y = "Density",
        fill = response_name
      ) +
      ggplot2::theme_bw()
  } else {
    p = ggplot2::ggplot(
      proj_train,
      ggplot2::aes(x = .data$LD1, y = .data$LD2, colour = .data$class)
    ) +
      ggplot2::geom_point(alpha = 0.15, size = 0.6) +
      ggplot2::geom_point(
        data = proj_obs,
        ggplot2::aes(x = .data$LD1, y = .data$LD2),
        colour = "black",
        shape = 8,
        size = 4,
        inherit.aes = FALSE
      ) +
      ggplot2::labs(
        title = "LDA projection of the reference table",
        subtitle = "The black star is the observed dataset",
        x = "LD1",
        y = "LD2",
        colour = response_name
      ) +
      ggplot2::theme_bw()
  }

  ggplot2::ggsave(lda_file, p, width = 8, height = 5.5, dpi = 120)
  message("Saved ", lda_file)
  print(p)
  invisible(p)
}

genomic_label_to_code = function(value, kind) {
  if (identical(value, "homogeneous")) {
    return(paste0("1", kind))
  }
  if (identical(value, "heterogeneous")) {
    return(paste0("2", kind))
  }
  stop(
    "Expected 'homogeneous' or 'heterogeneous', not: ",
    value,
    call. = FALSE
  )
}
