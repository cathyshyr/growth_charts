# ============================================================
# generate_growth_curves_with_diagnostics.R
#
# Overview:
#   General-purpose SITAR workflow for generating condition-
#   specific growth curve outputs for any condition and sex.
#   This script is written so it can be sourced and called from
#   another script that loops over many conditions.
#
# Main function:
#   generate_growth_curves()
#
# Expected upstream objects:
#   - dx : diagnosis table with person_id, dID, disease_short_title
#   - height_all : prepared height dataset from prep_height_data()
#   - height_no_dx : prepared unaffected dataset
# ============================================================

# ------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

# ------------------------------------------------------------
# Helper: validate core inputs
# ------------------------------------------------------------
validate_growth_curve_inputs <- function(dID,
                                         sex,
                                         dx,
                                         height_all,
                                         height_no_dx) {
  if (missing(dID) || length(dID) != 1) {
    stop("dID must be supplied as a single condition identifier.")
  }

  if (!is.null(sex) && !sex %in% c("M", "F")) {
    stop("sex must be NULL, 'M', or 'F'.")
  }

  required_dx_cols <- c("person_id", "dID", "disease_short_title")
  missing_dx_cols <- setdiff(required_dx_cols, names(dx))
  if (length(missing_dx_cols) > 0) {
    stop("dx is missing required columns: ", paste(missing_dx_cols, collapse = ", "))
  }

  required_height_cols <- c("person_id", "height", "age")
  missing_height_cols <- setdiff(required_height_cols, names(height_all))
  if (length(missing_height_cols) > 0) {
    stop("height_all is missing required columns: ", paste(missing_height_cols, collapse = ", "))
  }

  missing_unaffected_cols <- setdiff(required_height_cols, names(height_no_dx))
  if (length(missing_unaffected_cols) > 0) {
    stop("height_no_dx is missing required columns: ", paste(missing_unaffected_cols, collapse = ", "))
  }

  invisible(TRUE)
}

# ------------------------------------------------------------
# Helper: thin data to 1 data point per patient per age bin
# ------------------------------------------------------------
thin_centile_input_by_age_bin <- function(dat, bin_width = 0.5) {
  dt <- as.data.table(copy(dat))
  
  dt[, age_bin := floor(age / bin_width) * bin_width]
  dt[, age_bin_mid := age_bin + bin_width / 2]
  dt[, dist_to_mid := abs(age - age_bin_mid)]
  
  setorder(dt, person_id, age_bin, dist_to_mid, age)
  
  dt <- dt[, .SD[1], by = .(person_id, age_bin)]
  
  dt[, c("age_bin", "age_bin_mid", "dist_to_mid") := NULL]
  dt[]
}


# ------------------------------------------------------------
# Helper: resolve a readable condition label
# ------------------------------------------------------------
resolve_condition_label <- function(dx, dID, condition_label = NULL) {
  if (!is.null(condition_label) && nzchar(condition_label)) {
    return(condition_label)
  }

  disease_name <- unique(dx$disease_short_title[dx$dID == dID])
  disease_name <- disease_name[!is.na(disease_name)]

  if (length(disease_name) == 0) {
    return(paste0("dID_", dID))
  }

  disease_name[1]
}

# ------------------------------------------------------------
# Helper: make safe file names
# ------------------------------------------------------------
make_safe_slug <- function(x) {
  x |>
    tolower() |>
    gsub("[^a-z0-9]+", "_", x = _) |>
    gsub("^_+|_+$", "", x = _)
}

# ------------------------------------------------------------
# Helper: summarize counts before/after cleaning
# ------------------------------------------------------------
make_cleaning_summary <- function(pipeline_obj) {
  data.frame(
    stage = c("raw", "cleaned"),
    n_ids = c(
      data.table::uniqueN(pipeline_obj$raw_data$person_id),
      data.table::uniqueN(pipeline_obj$cleaned_data$person_id)
    ),
    n_rows = c(
      nrow(pipeline_obj$raw_data),
      nrow(pipeline_obj$cleaned_data)
    ),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Helper: pick age scale used by the best prelim model
# ------------------------------------------------------------
get_best_xvar <- function(pipeline_obj) {
  best_row <- pipeline_obj$prelim_fit_summary_ok[1, ]
  as.character(best_row$xvar[[1]])
}

# ------------------------------------------------------------
# Helper: residual diagnostic plot
# ------------------------------------------------------------
plot_prelim_residuals <- function(prelim_flagged,
                                  title = "Preliminary SITAR residual diagnostics") {
  ggplot(prelim_flagged, aes(x = fitted_height, y = pearson_resid, color = outlier_point)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = c(-4, 4), linetype = 2) +
    theme_minimal() +
    labs(
      title = title,
      x = "Fitted height",
      y = "Pearson residual",
      color = "Outlier"
    )
}


# ------------------------------------------------------------
# Helper: build row-level diagnostic dataset for a fitted model
# ------------------------------------------------------------
make_model_diagnostic_data <- function(fit, dat, stage = c("preliminary", "final")) {
  stage <- match.arg(stage)

  out <- data.table::copy(dat)
  out[, diagnostic_stage := stage]
  out[, fitted_height := as.numeric(fitted(fit, level = 1))]
  out[, residual_response := as.numeric(resid(fit, type = "response"))]
  out[, residual_pearson := as.numeric(resid(fit, type = "pearson"))]
  out[, abs_residual_pearson := abs(residual_pearson)]

  out[]
}

# ------------------------------------------------------------
# Helper: summarize model diagnostics numerically
# ------------------------------------------------------------
make_model_diagnostic_summary <- function(diagnostic_data) {
  if (nrow(diagnostic_data) == 0) {
    return(data.frame())
  }

  rp <- diagnostic_data$residual_pearson
  rr <- diagnostic_data$residual_response

  data.frame(
    diagnostic_stage = unique(diagnostic_data$diagnostic_stage)[1],
    n_rows = nrow(diagnostic_data),
    n_ids = data.table::uniqueN(diagnostic_data$person_id),
    mean_residual_pearson = mean(rp, na.rm = TRUE),
    sd_residual_pearson = stats::sd(rp, na.rm = TRUE),
    mean_abs_residual_pearson = mean(abs(rp), na.rm = TRUE),
    median_abs_residual_pearson = stats::median(abs(rp), na.rm = TRUE),
    q25_residual_pearson = as.numeric(stats::quantile(rp, probs = 0.25, na.rm = TRUE)),
    q75_residual_pearson = as.numeric(stats::quantile(rp, probs = 0.75, na.rm = TRUE)),
    max_abs_residual_pearson = max(abs(rp), na.rm = TRUE),
    n_abs_residual_gt_2 = sum(abs(rp) > 2, na.rm = TRUE),
    n_abs_residual_gt_3 = sum(abs(rp) > 3, na.rm = TRUE),
    n_abs_residual_gt_4 = sum(abs(rp) > 4, na.rm = TRUE),
    cor_abs_residual_vs_fitted = stats::cor(
      abs(rp),
      diagnostic_data$fitted_height,
      use = "complete.obs"
    ),
    mean_residual_response = mean(rr, na.rm = TRUE),
    sd_residual_response = stats::sd(rr, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Helper: extract random effects for diagnostics
# ------------------------------------------------------------
extract_random_effects_diagnostics <- function(fit) {
  re <- tryCatch(nlme::ranef(fit), error = function(e) NULL)

  if (is.null(re)) {
    return(data.frame())
  }

  re_df <- as.data.frame(re)
  re_df$person_id <- rownames(re_df)
  rownames(re_df) <- NULL

  base_cols <- intersect(c("person_id", "a", "b", "c"), names(re_df))
  re_df[, base_cols, drop = FALSE]
}

# ------------------------------------------------------------
# Helper: summarize random effects numerically
# ------------------------------------------------------------
make_random_effects_summary <- function(random_effects_df) {
  param_cols <- intersect(c("a", "b", "c"), names(random_effects_df))

  if (length(param_cols) == 0) {
    return(data.frame())
  }

  dplyr::bind_rows(lapply(param_cols, function(param) {
    x <- random_effects_df[[param]]
    data.frame(
      parameter = param,
      n_ids = sum(!is.na(x)),
      mean = mean(x, na.rm = TRUE),
      sd = stats::sd(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE),
      q25 = as.numeric(stats::quantile(x, probs = 0.25, na.rm = TRUE)),
      q75 = as.numeric(stats::quantile(x, probs = 0.75, na.rm = TRUE)),
      min = min(x, na.rm = TRUE),
      max = max(x, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

# ------------------------------------------------------------
# Helper: QQ plot for model residuals
# ------------------------------------------------------------
plot_model_residual_qq <- function(diagnostic_data,
                                   title = "SITAR residual QQ plot") {
  ggplot(diagnostic_data, aes(sample = residual_pearson)) +
    stat_qq(alpha = 0.4) +
    stat_qq_line() +
    theme_minimal() +
    labs(
      title = title,
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    )
}

# ------------------------------------------------------------
# Helper: residual vs age plot
# ------------------------------------------------------------
plot_model_residuals_vs_age <- function(diagnostic_data,
                                        title = "SITAR residuals vs age") {
  ggplot(diagnostic_data, aes(x = age, y = residual_pearson)) +
    geom_point(alpha = 0.35) +
    geom_hline(yintercept = 0, linetype = 2) +
    theme_minimal() +
    labs(
      title = title,
      x = "Age (years)",
      y = "Pearson residual"
    )
}

# ------------------------------------------------------------
# Helper: observed vs fitted plot
# ------------------------------------------------------------
plot_model_observed_vs_fitted <- function(diagnostic_data,
                                          title = "Observed vs fitted height") {
  ggplot(diagnostic_data, aes(x = fitted_height, y = height)) +
    geom_point(alpha = 0.35) +
    geom_abline(intercept = 0, slope = 1, linetype = 2) +
    theme_minimal() +
    labs(
      title = title,
      x = "Fitted height",
      y = "Observed height"
    )
}

# ------------------------------------------------------------
# Helper: QQ plot for random effects
# ------------------------------------------------------------
plot_random_effects_qq <- function(random_effects_df,
                                   title = "SITAR random-effects QQ plots") {
  param_cols <- intersect(c("a", "b", "c"), names(random_effects_df))

  if (length(param_cols) == 0) {
    return(NULL)
  }

  long_df <- dplyr::bind_rows(lapply(param_cols, function(param) {
    data.frame(
      parameter = param,
      value = random_effects_df[[param]],
      stringsAsFactors = FALSE
    )
  }))

  ggplot(long_df, aes(sample = value)) +
    stat_qq(alpha = 0.4) +
    stat_qq_line() +
    facet_wrap(~ parameter, scales = "free") +
    theme_minimal() +
    labs(
      title = title,
      x = "Theoretical quantiles",
      y = "Sample quantiles"
    )
}

# ------------------------------------------------------------
# Helper: predicted mean curves plot
# ------------------------------------------------------------
plot_affected_unaffected_curves <- function(pred_curves,
                                     condition_label,
                                     sex = NULL) {
  sex_label <- if (is.null(sex)) "all sexes" else if (sex == "M") "males" else "females"

  ggplot(pred_curves, aes(x = age, y = height, color = group)) +
    geom_line(linewidth = 1) +
    theme_minimal() +
    labs(
      title = paste0(condition_label, " vs unaffected (", sex_label, ")"),
      x = "Age (years)",
      y = "Height (cm)",
      color = "Group"
    )
}

# ------------------------------------------------------------
# Main function: generate SITAR growth curve outputs
#
# Overview:
#   Runs the end-to-end SITAR affected-unaffected pipeline for one
#   condition and one optional sex, then returns model objects,
#   summaries, predictions, and plots. Optionally saves outputs.
#
# Inputs:
#   dID : scalar
#     Condition identifier.
#
#   sex : character or NULL, default = NULL
#     Sex filter to apply. Use "M", "F", or NULL.
#
#   dx : data.frame/data.table
#     Diagnosis table.
#
#   height_all : data.table
#     Prepared full height dataset.
#
#   height_no_dx : data.table
#     Prepared unaffected dataset.
#
#   condition_label : character or NULL
#     Optional label used in outputs and file names.
#
#   output_dir : character or NULL
#     If supplied, outputs are written to this folder.
#
#   save_rds : logical, default = TRUE
#     Whether to save the returned object as an .rds file when
#     output_dir is supplied.
#
#   save_plots : logical, default = TRUE
#     Whether to save diagnostic/prediction plots when output_dir
#     is supplied.
#
#   ages_grid : numeric vector
#     Age grid used for curve prediction.
#
#   unaffected_ratio, max_unaffecteds, min_unaffecteds, seed,
#   age_min, age_max, min_obs_per_id, resid_cutoff :
#     Passed to run_condition_pipeline().
#
# Output:
#   A named list containing:
#     - metadata
#     - pipeline
#     - candidate_models
#     - prelim_summary
#     - final_summary
#     - affected_unaffected_effects
#     - sitar_abc
#     - cleaning_summary
#     - outlier_counts
#     - predicted_curves
#     - residual_plot
#     - curve_plot
# ------------------------------------------------------------
generate_growth_curves <- function(dID,
                                   sex = NULL,
                                   dx,
                                   height_all,
                                   height_no_dx,
                                   condition_label = NULL,
                                   variantClass = "None",
                                   subset_ids = NULL,
                                   output_dir = NULL,
                                   save_rds = TRUE,
                                   save_plots = TRUE,
                                   ages_grid = seq(2, 20, by = 0.1),
                                   unaffected_ratio = 10,
                                   max_unaffecteds = 10000,
                                   min_unaffecteds = 1000,
                                   seed = 123,
                                   age_min = 2,
                                   age_max = 20,
                                   min_obs_per_id = 3,
                                   resid_cutoff = 4) {
  
  validate_growth_curve_inputs(
    dID = dID,
    sex = sex,
    dx = dx,
    height_all = height_all,
    height_no_dx = height_no_dx
  )
  
  condition_label <- resolve_condition_label(
    dx = dx,
    dID = dID,
    condition_label = condition_label
  )
  
  sex_label <- if (is.null(sex)) "all" else sex
  
  variant_slug <- if (!is.null(variantClass) && nzchar(variantClass)) {
    make_safe_slug(variantClass)
  } else {
    "none"
  }
  
  slug <- paste(
    make_safe_slug(condition_label),
    tolower(sex_label),
    paste0("did_", dID),
    paste0("vc_", variant_slug),
    sep = "_"
  )
  
  if (condition_label == "Unaffected") {
    slug <- paste(
      "unaffected",
      tolower(sex_label),
      paste0("n_", length(unique(subset_ids))),
      sep = "_"
    )
  }
  
  is_unaffected_target <- !is.null(condition_label) && condition_label == "Unaffected"
  
  if (is_unaffected_target) {
    pipeline <- run_unaffected_pipeline(
      height_no_dx = height_no_dx,
      sex = sex,
      subset_ids = subset_ids,
      age_min = age_min,
      age_max = age_max,
      min_obs_per_id = min_obs_per_id,
      resid_cutoff = resid_cutoff
    )
  } else {
    pipeline <- run_condition_pipeline(
      dID = dID,
      dx = dx,
      height_all = height_all,
      height_no_dx = height_no_dx,
      sex = sex,
      subset_ids = subset_ids,
      condition_label = condition_label,
      unaffected_ratio = unaffected_ratio,
      max_unaffecteds = max_unaffecteds,
      min_unaffecteds = min_unaffecteds,
      seed = seed,
      age_min = age_min,
      age_max = age_max,
      min_obs_per_id = min_obs_per_id,
      resid_cutoff = resid_cutoff
    )
  }
  
  best_xvar <- pipeline$best_xvar
  
  if (is_unaffected_target) {
    predicted_curves <- NULL
    affected_unaffected_effects <- data.frame()
    curve_plot <- NULL
  } else {
    predicted_curves <- predict_affected_unaffected_curves(
      final_fit = pipeline$final_fit,
      ages_grid = ages_grid,
      xvar = best_xvar
    )
    
    affected_unaffected_effects <- extract_affected_unaffected_effects(pipeline$final_fit)
    
    curve_plot <- plot_affected_unaffected_curves(
      pred_curves = predicted_curves,
      condition_label = condition_label,
      sex = sex
    )
  }
  
  sitar_abc <- tryCatch(
    extract_sitar_abc(pipeline$final_fit),
    error = function(e) {
      message("extract_sitar_abc failed: ", conditionMessage(e))
      data.frame(
        parameter = character(),
        group = character(),
        estimate = numeric(),
        se = numeric(),
        lower = numeric(),
        upper = numeric()
      )
    }
  )
  cleaning_summary <- make_cleaning_summary(pipeline)
  outlier_counts <- as.data.frame(table(pipeline$prelim_flagged$outlier_point), stringsAsFactors = FALSE)
  names(outlier_counts) <- c("outlier_point", "n")
  
  prelim_diagnostic_data <- make_model_diagnostic_data(
    fit = pipeline$prelim_fit,
    dat = pipeline$raw_data,
    stage = "preliminary"
  )
  
  final_diagnostic_data <- make_model_diagnostic_data(
    fit = pipeline$final_fit,
    dat = pipeline$cleaned_data,
    stage = "final"
  )
  
  diagnostic_summary <- dplyr::bind_rows(
    make_model_diagnostic_summary(prelim_diagnostic_data),
    make_model_diagnostic_summary(final_diagnostic_data)
  )
  
  random_effects_diagnostics <- extract_random_effects_diagnostics(pipeline$final_fit)
  random_effects_summary <- make_random_effects_summary(random_effects_diagnostics)
  
  residual_plot <- plot_prelim_residuals(
    prelim_flagged = pipeline$prelim_flagged,
    title = paste0(condition_label, " preliminary SITAR residuals")
  )
  
  prelim_residual_qq_plot <- plot_model_residual_qq(
    diagnostic_data = prelim_diagnostic_data,
    title = paste0(condition_label, " preliminary residual QQ plot")
  )
  
  final_residual_qq_plot <- plot_model_residual_qq(
    diagnostic_data = final_diagnostic_data,
    title = paste0(condition_label, " final residual QQ plot")
  )
  
  final_residuals_vs_age_plot <- plot_model_residuals_vs_age(
    diagnostic_data = final_diagnostic_data,
    title = paste0(condition_label, " final residuals vs age")
  )
  
  final_observed_vs_fitted_plot <- plot_model_observed_vs_fitted(
    diagnostic_data = final_diagnostic_data,
    title = paste0(condition_label, " final observed vs fitted")
  )
  
  random_effects_qq_plot <- plot_random_effects_qq(
    random_effects_df = random_effects_diagnostics,
    title = paste0(condition_label, " final random-effects QQ plots")
  )
  
  result <- list(
    metadata = list(
      dID = dID,
      condition_label = condition_label,
      sex = sex,
      variantClass = variantClass,
      subset_ids = subset_ids,
      seed = seed,
      age_min = age_min,
      age_max = age_max,
      min_obs_per_id = min_obs_per_id,
      resid_cutoff = resid_cutoff,
      unaffected_ratio = unaffected_ratio,
      max_unaffecteds = max_unaffecteds,
      min_unaffecteds = min_unaffecteds,
      best_xvar = best_xvar,
      best_df = pipeline$best_df,
      is_unaffected_target = is_unaffected_target
    ),
    pipeline = pipeline,
    candidate_models = pipeline$prelim_fit_summary,
    prelim_summary = summary(pipeline$prelim_fit),
    final_summary = summary(pipeline$final_fit),
    affected_unaffected_effects = affected_unaffected_effects,
    sitar_abc = sitar_abc,
    cleaning_summary = cleaning_summary,
    outlier_counts = outlier_counts,
    diagnostic_summary = diagnostic_summary,
    prelim_diagnostic_data = prelim_diagnostic_data,
    final_diagnostic_data = final_diagnostic_data,
    random_effects_diagnostics = random_effects_diagnostics,
    random_effects_summary = random_effects_summary,
    predicted_curves = predicted_curves,
    residual_plot = residual_plot,
    prelim_residual_qq_plot = prelim_residual_qq_plot,
    final_residual_qq_plot = final_residual_qq_plot,
    final_residuals_vs_age_plot = final_residuals_vs_age_plot,
    final_observed_vs_fitted_plot = final_observed_vs_fitted_plot,
    random_effects_qq_plot = random_effects_qq_plot,
    curve_plot = curve_plot
  )
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    utils::write.csv(
      pipeline$prelim_fit_summary,
      file = file.path(output_dir, paste0(slug, "_prelim_fit_summary.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      pipeline$final_fit_summary,
      file = file.path(output_dir, paste0(slug, "_final_fit_summary.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      sitar_abc,
      file = file.path(output_dir, paste0(slug, "_sitar_abc.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      cleaning_summary,
      file = file.path(output_dir, paste0(slug, "_cleaning_summary.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      outlier_counts,
      file = file.path(output_dir, paste0(slug, "_outlier_counts.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      diagnostic_summary,
      file = file.path(output_dir, paste0(slug, "_diagnostic_summary.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      prelim_diagnostic_data,
      file = file.path(output_dir, paste0(slug, "_prelim_diagnostic_data.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      final_diagnostic_data,
      file = file.path(output_dir, paste0(slug, "_final_diagnostic_data.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      random_effects_diagnostics,
      file = file.path(output_dir, paste0(slug, "_random_effects_diagnostics.csv")),
      row.names = FALSE
    )
    
    utils::write.csv(
      random_effects_summary,
      file = file.path(output_dir, paste0(slug, "_random_effects_summary.csv")),
      row.names = FALSE
    )
    
    if (!is_unaffected_target) {
      utils::write.csv(
        affected_unaffected_effects,
        file = file.path(output_dir, paste0(slug, "_affected_unaffected_effects.csv")),
        row.names = FALSE
      )
      
      utils::write.csv(
        predicted_curves,
        file = file.path(output_dir, paste0(slug, "_predicted_curves.csv")),
        row.names = FALSE
      )
    }
    
    if (isTRUE(save_plots)) {
      ggsave(
        filename = file.path(output_dir, paste0(slug, "_residual_plot.png")),
        plot = residual_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(output_dir, paste0(slug, "_prelim_residual_qq_plot.png")),
        plot = prelim_residual_qq_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(output_dir, paste0(slug, "_final_residual_qq_plot.png")),
        plot = final_residual_qq_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(output_dir, paste0(slug, "_final_residuals_vs_age_plot.png")),
        plot = final_residuals_vs_age_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(output_dir, paste0(slug, "_final_observed_vs_fitted_plot.png")),
        plot = final_observed_vs_fitted_plot,
        width = 8,
        height = 5,
        dpi = 300
      )
      
      if (!is.null(random_effects_qq_plot)) {
        ggsave(
          filename = file.path(output_dir, paste0(slug, "_random_effects_qq_plot.png")),
          plot = random_effects_qq_plot,
          width = 9,
          height = 5,
          dpi = 300
        )
      }
      
      if (!is_unaffected_target && !is.null(curve_plot)) {
        ggsave(
          filename = file.path(output_dir, paste0(slug, "_growth_curves.png")),
          plot = curve_plot,
          width = 8,
          height = 5,
          dpi = 300
        )
      }
    }
    
    if (isTRUE(save_rds)) {
      saveRDS(
        result,
        file = file.path(output_dir, paste0(slug, "_sitar_results.rds"))
      )
    }
  }
  
  result
}
