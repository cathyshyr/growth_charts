# ============================================================
# generate_growth_centiles_with_diagnostics.R
#
# Overview:
#   General-purpose GAMLSS workflow for generating condition-
#   specific height centile curves from the same cleaned data
#   used by the SITAR workflow.
#
# Design:
#   This script is intended to stay separate from the SITAR
#   script, but it can directly reuse the cleaned_data object
#   produced by generate_growth_curves(). That keeps the
#   cleaning logic in one place while letting SITAR and GAMLSS
#   remain modular.
#
# Recommended use:
#   1) Run generate_growth_curves()
#   2) Pass the returned object into generate_growth_centiles()
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(gamlss)
  library(gamlss.dist)
})

# ------------------------------------------------------------
# Helper: make safe file names
# ------------------------------------------------------------
make_safe_slug <- function(x) {
  x <- tolower(as.character(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

# ------------------------------------------------------------
# Helper: human-readable sex labels for outputs
# ------------------------------------------------------------
format_sex_label <- function(sex) {
  if (is.null(sex) || is.na(sex)) return("All")
  if (sex == "M") return("Male")
  if (sex == "F") return("Female")
  as.character(sex)
}

# ------------------------------------------------------------
# Helper: get cleaned data from SITAR results object or file
# ------------------------------------------------------------
resolve_sitar_results <- function(sitar_results = NULL, sitar_rds = NULL) {
  if (!is.null(sitar_results)) {
    return(sitar_results)
  }
  
  if (!is.null(sitar_rds)) {
    return(readRDS(sitar_rds))
  }
  
  stop("Provide either sitar_results or sitar_rds.")
}

# ------------------------------------------------------------
# Helper: prepare cross-sectional GAMLSS input from cleaned SITAR data
# ------------------------------------------------------------
prepare_centile_input <- function(cleaned_data,
                                  cohort = c("affected", "unaffected", "all"),
                                  age_min = 2,
                                  age_max = 20) {
  cohort <- match.arg(cohort)
  dt <- as.data.table(copy(cleaned_data))
  
  required_cols <- c("person_id", "age", "height")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    stop("cleaned_data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if (cohort == "affected") {
    if (!"affected_status" %in% names(dt)) stop("cleaned_data does not contain affected_status.")
    dt <- dt[affected_status == 1L]
  } else if (cohort == "unaffected") {
    if (!"affected_status" %in% names(dt)) stop("cleaned_data does not contain affected_status.")
    dt <- dt[affected_status == 0L]
  }
  
  dt <- dt[
    is.finite(age) &
      is.finite(height) &
      age >= age_min & age <= age_max &
      height > 0
  ]
  
  setorder(dt, age, person_id)
  
  dt[, .(
    person_id = person_id,
    age = as.numeric(age),
    height = as.numeric(height)
  )]
}

# ------------------------------------------------------------
# Candidate transformations used in model search
# ------------------------------------------------------------
y_transforms <- list(
  identity = list(
    forward = function(y) y,
    inverse = function(y) y
  ),
  log = list(
    forward = function(y) log(y),
    inverse = function(y) exp(y)
  ),
  sqrt = list(
    forward = function(y) sqrt(y),
    inverse = function(y) y^2
  )
)

x_transforms <- list(
  identity = function(x) x,
  log = function(x) log(x),
  sqrt = function(x) sqrt(x)
)

# ------------------------------------------------------------
# Helper: enforce non-decreasing centiles with age
# ------------------------------------------------------------
enforce_monotone_centiles <- function(cent_df, centiles) {
  out <- as.data.frame(cent_df)
  centile_cols <- paste0("C", centiles)
  
  for (cc in centile_cols) {
    out[[cc]] <- cummax(out[[cc]])
  }
  
  out
}

# ------------------------------------------------------------
# Helper: convert wide centiles to long plotting format
# ------------------------------------------------------------
make_centiles_long <- function(cent_df, centiles) {
  cent_df %>%
    pivot_longer(
      cols = starts_with("C"),
      names_to = "centile_code",
      values_to = "height_pred"
    ) %>%
    mutate(
      centile = gsub("^C", "", centile_code),
      centile = factor(centile, levels = as.character(centiles))
    )
}

# ------------------------------------------------------------
# Helper: RShiny data format
# ------------------------------------------------------------
make_rshiny_centile_data <- function(cent_df,
                                     centiles,
                                     sex_label,
                                     condition_label,
                                     variantClass = "None",
                                     n_value,
                                     model_label = "best") {
  out <- cent_df %>%
    pivot_longer(
      cols = starts_with("C"),
      names_to = "centile",
      values_to = "height"
    ) %>%
    mutate(
      centile = sub("^C", "c", centile),
      centile = ifelse(grepl("^c\\d$", centile), sub("^c", "c0", centile), centile),
      sex = sex_label,
      condition = condition_label,
      variantClass = variantClass,
      n = as.integer(n_value),
      model = model_label
    ) %>%
    select(age, centile, height, sex, condition, variantClass, n, model)
  
  out
}

# ------------------------------------------------------------
# Helper: label data at maximum age for major centiles
# ------------------------------------------------------------
make_rshiny_label_data <- function(rshiny_data,
                                   label_age = NULL,
                                   label_centiles = c(5, 10, 25, 50, 75, 90, 95)) {
  dat <- as_tibble(rshiny_data)
  
  if (is.null(label_age)) {
    label_age <- max(dat$age, na.rm = TRUE)
  }
  
  label_codes <- paste0("c", sprintf("%02d", label_centiles))
  
  dat %>%
    filter(age == label_age, centile %in% label_codes) %>%
    mutate(centile = sub("^c0?", "", centile)) %>%
    arrange(match(as.numeric(centile), label_centiles))
}

# ------------------------------------------------------------
# Helper: centile plot
# ------------------------------------------------------------
plot_growth_centiles <- function(dat,
                                 cent_long,
                                 model_row,
                                 condition_label,
                                 sex_label,
                                 cohort_label,
                                 point_alpha = 0.35,
                                 title_suffix = "",
                                 subtitle_prefix = "") {
  linewidth_values <- c(
    "5" = 0.5,
    "10" = 0.6,
    "25" = 0.7,
    "50" = 1.1,
    "75" = 0.7,
    "90" = 0.6,
    "95" = 0.5
  )
  
  present_levels <- intersect(names(linewidth_values), as.character(unique(cent_long$centile)))
  linewidth_values <- linewidth_values[present_levels]
  
  subtitle_text <- paste0(
    subtitle_prefix,
    "Cohort: ", cohort_label,
    " | Family: ", model_row$family,
    " | age transform = ", model_row$x_transform,
    " | height transform = ", model_row$y_transform,
    " | BIC = ", round(model_row$bic, 2)
  )
  
  ggplot() +
    geom_point(
      data = dat,
      aes(x = age, y = height),
      color = "grey70",
      alpha = point_alpha,
      size = 1.0
    ) +
    geom_line(
      data = cent_long,
      aes(x = age, y = height_pred, group = centile, linewidth = centile),
      color = "black"
    ) +
    scale_linewidth_manual(values = linewidth_values, guide = "none") +
    labs(
      title = paste0(condition_label, " ", sex_label, " centile curves", title_suffix),
      subtitle = subtitle_text,
      x = "Age (years)",
      y = "Height (cm)"
    ) +
    theme_bw(base_size = 12)
}

# ------------------------------------------------------------
# Helper: diagnostic data from fitted GAMLSS model
# ------------------------------------------------------------
make_gamlss_diagnostic_data <- function(fit_obj, fit_data) {
  data.frame(
    age = fit_data$age,
    height = fit_data$height,
    fitted = fitted(fit_obj),
    rqres = residuals(fit_obj, what = "z-scores")
  )
}

# ------------------------------------------------------------
# Helper: QQ plot of normalized quantile residuals
# ------------------------------------------------------------
plot_gamlss_residual_qq <- function(diagnostic_data,
                                    title = "GAMLSS residual QQ plot") {
  ggplot(diagnostic_data, aes(sample = rqres)) +
    stat_qq(alpha = 0.4) +
    stat_qq_line() +
    theme_bw(base_size = 12) +
    labs(
      title = title,
      x = "Theoretical quantiles",
      y = "Normalized quantile residuals"
    )
}

# ------------------------------------------------------------
# Helper: worm plot
# wp() can use the data already stored in the fitted gamlss object,
# so do not pass data again here.
# ------------------------------------------------------------
# ------------------------------------------------------------
# Helper: worm plot
# Keep the wp() call minimal, because wp() already handles some
# plotting arguments internally.
# ------------------------------------------------------------
plot_gamlss_worm <- function(fit_obj,
                             xvar_name = "age",
                             n_intervals = 6) {
  wp(
    object = fit_obj,
    xvar = as.formula(paste0("~", xvar_name)),
    n.inter = n_intervals
  )
  
  invisible(recordPlot())
}

# ------------------------------------------------------------
# Helper: save recorded base-R plot to file
# ------------------------------------------------------------
save_recorded_plot <- function(recorded_plot,
                               filename,
                               width = 9,
                               height = 6,
                               res = 300) {
  png(filename = filename, width = width, height = height, units = "in", res = res)
  replayPlot(recorded_plot)
  dev.off()
}

# ------------------------------------------------------------
# Helper: fit all candidate GAMLSS models
# ------------------------------------------------------------
fit_gamlss_grid <- function(dat,
                            families = c("NO", "BCCGo", "BCPEo", "BCTo")) {
  model_grid <- expand.grid(
    family = families,
    x_transform = names(x_transforms),
    y_transform = names(y_transforms),
    stringsAsFactors = FALSE
  )
  
  fits_tbl <- purrr::pmap_dfr(
    model_grid,
    ~ fit_one_model(
      data = dat,
      family_name = ..1,
      x_tf_name = ..2,
      y_tf_name = ..3
    )
  )
  
  if (nrow(fits_tbl) == 0) {
    stop("No GAMLSS candidate models converged successfully.")
  }
  
  fits_tbl %>% arrange(bic)
}

# ------------------------------------------------------------
# Helper: extract best NO model, if available
# ------------------------------------------------------------
get_best_no_model_row <- function(fits_tbl) {
  no_tbl <- fits_tbl %>% filter(family == "NO") %>% arrange(bic)
  if (nrow(no_tbl) == 0) {
    return(NULL)
  }
  no_tbl %>% slice(1)
}

# ------------------------------------------------------------
# Helper: build outputs for one chosen model row
# ------------------------------------------------------------
build_model_outputs <- function(dat,
                                model_row,
                                ages_grid,
                                centiles,
                                sex_label,
                                condition_label,
                                cohort_label,
                                variantClass,
                                n_subjects,
                                plot_title_suffix = "",
                                plot_subtitle_prefix = "",
                                rshiny_condition = condition_label,
                                model_label = "best") {
  fit_obj <- model_row$fit[[1]]
  fit_data <- model_row$fit_data[[1]]
  
  cent_df <- predict_centiles(
    fit_obj = fit_obj,
    fit_data = fit_data,
    ages = ages_grid,
    family_name = model_row$family,
    x_tf_name = model_row$x_transform,
    y_tf_name = model_row$y_transform,
    centiles = centiles
  )
  
  cent_df_monotone <- enforce_monotone_centiles(cent_df, centiles)
  cent_long <- make_centiles_long(cent_df_monotone, centiles)
  
  rshiny_data <- make_rshiny_centile_data(
    cent_df = cent_df_monotone,
    centiles = centiles,
    sex_label = sex_label,
    condition_label = rshiny_condition,
    variantClass = variantClass,
    n_value = n_subjects,
    model_label = model_label
  )
  
  rshiny_labels <- make_rshiny_label_data(
    rshiny_data = rshiny_data,
    label_age = max(ages_grid),
    label_centiles = centiles
  )
  
  centile_plot <- plot_growth_centiles(
    dat = dat,
    cent_long = cent_long,
    model_row = model_row,
    condition_label = rshiny_condition,
    sex_label = sex_label,
    cohort_label = cohort_label,
    title_suffix = plot_title_suffix,
    subtitle_prefix = plot_subtitle_prefix
  )
  
  diagnostic_data <- make_gamlss_diagnostic_data(
    fit_obj = fit_obj,
    fit_data = fit_data
  )
  
  residual_qq_plot <- plot_gamlss_residual_qq(
    diagnostic_data = diagnostic_data,
    title = paste0(rshiny_condition, " ", sex_label, plot_title_suffix, " residual QQ plot")
  )
  
  worm_plot <- tryCatch(
    {
      plot_gamlss_worm(
        fit_obj = fit_obj,
        xvar_name = "age",
        n_intervals = 6
      )
    },
    error = function(e) {
      message("Worm plot failed: ", conditionMessage(e))
      NULL
    }
  )
  
  list(
    model_row = model_row,
    fit = fit_obj,
    fit_data = fit_data,
    diagnostic_data = diagnostic_data,
    centiles_wide = cent_df_monotone,
    centiles_long = cent_long,
    RShiny_data = rshiny_data,
    RShiny_data_labels = rshiny_labels,
    centile_plot = centile_plot,
    residual_qq_plot = residual_qq_plot,
    worm_plot = worm_plot
  )
}

# ------------------------------------------------------------
# Main function
# ------------------------------------------------------------
generate_growth_centiles <- function(sitar_results = NULL,
                                     sitar_rds = NULL,
                                     cohort = c("affected", "unaffected", "all"),
                                     condition_label = NULL,
                                     sex = NULL,
                                     variantClass = "None",
                                     ages_grid = seq(2, 20, by = 0.1),
                                     centiles = c(5, 10, 25, 50, 75, 90, 95),
                                     families = c("NO", "BCCGo", "BCPEo", "BCTo"),
                                     age_min = 2,
                                     age_max = 20,
                                     output_dir = NULL,
                                     save_rds = TRUE,
                                     save_plot = TRUE,
                                     save_csv = TRUE,
                                     save_no_model_outputs = TRUE) {
  cohort <- match.arg(cohort)
  
  sitar_obj <- resolve_sitar_results(
    sitar_results = sitar_results,
    sitar_rds = sitar_rds
  )
  
  if (is.null(sitar_obj$pipeline$cleaned_data)) {
    stop("SITAR results object does not contain pipeline$cleaned_data.")
  }
  
  metadata <- sitar_obj$metadata
  if (is.null(condition_label)) condition_label <- metadata$condition_label
  if (is.null(sex)) sex <- metadata$sex
  
  if (is.null(condition_label) || !nzchar(condition_label)) {
    condition_label <- "Unknown condition"
  }
  
  sex_label <- format_sex_label(sex)
  cohort_label <- switch(
    cohort,
    affected = condition_label,
    unaffected = "Unaffected",
    all = "Combined"
  )
  
  dat <- prepare_centile_input(
    cleaned_data = sitar_obj$pipeline$cleaned_data,
    cohort = cohort,
    age_min = age_min,
    age_max = age_max
  )
  
  if (nrow(dat) == 0) {
    stop("No rows available for GAMLSS after cohort and age filtering.")
  }
  
  n_subjects <- uniqueN(dat$person_id)
  n_rows <- nrow(dat)
  
  fits_tbl <- fit_gamlss_grid(dat = dat, families = families)
  
  best_row <- fits_tbl %>% slice(1)
  no_row <- get_best_no_model_row(fits_tbl)
  
  rshiny_condition <- if (cohort == "unaffected") "Unaffected" else condition_label
  
  best_outputs <- build_model_outputs(
    dat = dat,
    model_row = best_row,
    ages_grid = ages_grid,
    centiles = centiles,
    sex_label = sex_label,
    condition_label = condition_label,
    cohort_label = cohort_label,
    variantClass = variantClass,
    n_subjects = n_subjects,
    plot_title_suffix = "",
    plot_subtitle_prefix = "Best GAMLSS model by BIC | ",
    rshiny_condition = rshiny_condition,
    model_label = "best_bic"
  )
  
  no_outputs <- NULL
  if (!is.null(no_row)) {
    no_outputs <- build_model_outputs(
      dat = dat,
      model_row = no_row,
      ages_grid = ages_grid,
      centiles = centiles,
      sex_label = sex_label,
      condition_label = condition_label,
      cohort_label = cohort_label,
      variantClass = variantClass,
      n_subjects = n_subjects,
      plot_title_suffix = " (Normal model)",
      plot_subtitle_prefix = "Best NO family model | ",
      rshiny_condition = rshiny_condition,
      model_label = "best_no"
    )
  }
  
  result <- list(
    metadata = list(
      dID = metadata$dID,
      condition_label = condition_label,
      cohort = cohort,
      cohort_label = cohort_label,
      sex = sex,
      sex_label = sex_label,
      variantClass = variantClass,
      n_subjects = n_subjects,
      n_rows = n_rows,
      age_min = age_min,
      age_max = age_max,
      centiles = centiles,
      ages_grid = ages_grid,
      best_model_family = best_row$family,
      best_no_model_available = !is.null(no_row)
    ),
    gamlss_input = dat,
    candidate_models = fits_tbl,
    best_model = best_outputs$model_row,
    best_fit = best_outputs$fit,
    best_fit_data = best_outputs$fit_data,
    best_diagnostic_data = best_outputs$diagnostic_data,
    centiles_wide = best_outputs$centiles_wide,
    centiles_long = best_outputs$centiles_long,
    RShiny_data = best_outputs$RShiny_data,
    RShiny_data_labels = best_outputs$RShiny_data_labels,
    centile_plot = best_outputs$centile_plot,
    residual_qq_plot = best_outputs$residual_qq_plot,
    worm_plot = best_outputs$worm_plot,
    no_model = if (!is.null(no_outputs)) no_outputs$model_row else NULL,
    no_fit = if (!is.null(no_outputs)) no_outputs$fit else NULL,
    no_fit_data = if (!is.null(no_outputs)) no_outputs$fit_data else NULL,
    no_diagnostic_data = if (!is.null(no_outputs)) no_outputs$diagnostic_data else NULL,
    no_centiles_wide = if (!is.null(no_outputs)) no_outputs$centiles_wide else NULL,
    no_centiles_long = if (!is.null(no_outputs)) no_outputs$centiles_long else NULL,
    no_RShiny_data = if (!is.null(no_outputs)) no_outputs$RShiny_data else NULL,
    no_RShiny_data_labels = if (!is.null(no_outputs)) no_outputs$RShiny_data_labels else NULL,
    no_centile_plot = if (!is.null(no_outputs)) no_outputs$centile_plot else NULL,
    no_residual_qq_plot = if (!is.null(no_outputs)) no_outputs$residual_qq_plot else NULL,
    no_worm_plot = if (!is.null(no_outputs)) no_outputs$worm_plot else NULL
  )
  
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    
    sex_slug <- if (is.null(sex)) "all" else tolower(sex)
    cohort_slug <- make_safe_slug(cohort)
    
    variant_slug <- if (!is.null(variantClass) && !is.na(variantClass) && nzchar(variantClass)) {
      make_safe_slug(variantClass)
    } else {
      "none"
    }
    
    prefix <- paste(
      make_safe_slug(rshiny_condition),
      sex_slug,
      cohort_slug,
      paste0("did_", metadata$dID),
      paste0("vc_", variant_slug),
      sep = "_"
    )
    
    if (isTRUE(save_csv)) {
      write.csv(
        fits_tbl %>% select(-fit, -fit_data),
        file = file.path(output_dir, paste0(prefix, "_gamlss_fit_summary.csv")),
        row.names = FALSE
      )
      
      write.csv(
        best_outputs$centiles_wide,
        file = file.path(output_dir, paste0(prefix, "_centiles_wide.csv")),
        row.names = FALSE
      )
      
      write.csv(
        best_outputs$centiles_long,
        file = file.path(output_dir, paste0(prefix, "_centiles_long.csv")),
        row.names = FALSE
      )
      
      write.csv(
        best_outputs$RShiny_data,
        file = file.path(output_dir, paste0(prefix, "_RShiny_data.csv")),
        row.names = FALSE
      )
      
      write.csv(
        best_outputs$RShiny_data_labels,
        file = file.path(output_dir, paste0(prefix, "_RShiny_data_labels.csv")),
        row.names = FALSE
      )
      
      write.csv(
        best_outputs$diagnostic_data,
        file = file.path(output_dir, paste0(prefix, "_diagnostic_data.csv")),
        row.names = FALSE
      )
      
      write.csv(
        data.frame(metric = c("n_rows", "n_subjects"), value = c(n_rows, n_subjects)),
        file = file.path(output_dir, paste0(prefix, "_sample_sizes.csv")),
        row.names = FALSE
      )
      
      write.csv(
        as.data.frame(best_outputs$model_row) %>% select(-fit, -fit_data),
        file = file.path(output_dir, paste0(prefix, "_best_model_summary.csv")),
        row.names = FALSE
      )
      
      if (isTRUE(save_no_model_outputs) && !is.null(no_outputs)) {
        write.csv(
          no_outputs$centiles_wide,
          file = file.path(output_dir, paste0(prefix, "_no_centiles_wide.csv")),
          row.names = FALSE
        )
        
        write.csv(
          no_outputs$centiles_long,
          file = file.path(output_dir, paste0(prefix, "_no_centiles_long.csv")),
          row.names = FALSE
        )
        
        write.csv(
          no_outputs$RShiny_data,
          file = file.path(output_dir, paste0(prefix, "_no_RShiny_data.csv")),
          row.names = FALSE
        )
        
        write.csv(
          no_outputs$RShiny_data_labels,
          file = file.path(output_dir, paste0(prefix, "_no_RShiny_data_labels.csv")),
          row.names = FALSE
        )
        
        write.csv(
          no_outputs$diagnostic_data,
          file = file.path(output_dir, paste0(prefix, "_no_diagnostic_data.csv")),
          row.names = FALSE
        )
        
        write.csv(
          as.data.frame(no_outputs$model_row) %>% select(-fit, -fit_data),
          file = file.path(output_dir, paste0(prefix, "_no_model_summary.csv")),
          row.names = FALSE
        )
      }
    }
    
    if (isTRUE(save_plot)) {
      ggsave(
        filename = file.path(output_dir, paste0(prefix, "_centile_plot.png")),
        plot = best_outputs$centile_plot,
        width = 9,
        height = 6,
        dpi = 300
      )
      
      ggsave(
        filename = file.path(output_dir, paste0(prefix, "_residual_qq_plot.png")),
        plot = best_outputs$residual_qq_plot,
        width = 7,
        height = 7,
        dpi = 300
      )
      
      if (!is.null(best_outputs$worm_plot)) {
        save_recorded_plot(
          recorded_plot = best_outputs$worm_plot,
          filename = file.path(output_dir, paste0(prefix, "_worm_plot.png")),
          width = 9,
          height = 6,
          res = 300
        )
      }
      
      if (isTRUE(save_no_model_outputs) && !is.null(no_outputs)) {
        ggsave(
          filename = file.path(output_dir, paste0(prefix, "_no_centile_plot.png")),
          plot = no_outputs$centile_plot,
          width = 9,
          height = 6,
          dpi = 300
        )
        
        ggsave(
          filename = file.path(output_dir, paste0(prefix, "_no_residual_qq_plot.png")),
          plot = no_outputs$residual_qq_plot,
          width = 7,
          height = 7,
          dpi = 300
        )
        
        if (!is.null(no_outputs$worm_plot)) {
          save_recorded_plot(
            recorded_plot = no_outputs$worm_plot,
            filename = file.path(output_dir, paste0(prefix, "_no_worm_plot.png")),
            width = 9,
            height = 6,
            res = 300
          )
        }
      }
    }
    
    if (isTRUE(save_rds)) {
      saveRDS(
        result,
        file = file.path(output_dir, paste0(prefix, "_gamlss_results.rds"))
      )
      
      saveRDS(
        best_outputs$fit,
        file = file.path(output_dir, paste0(prefix, "_best_gamlss_fit.rds"))
      )
      
      if (isTRUE(save_no_model_outputs) && !is.null(no_outputs)) {
        saveRDS(
          no_outputs$fit,
          file = file.path(output_dir, paste0(prefix, "_no_gamlss_fit.rds"))
        )
      }
    }
  }
  
  result
}