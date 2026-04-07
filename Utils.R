# Load packages
library(gamlss)
library(gamlss.dist)
library(tidyr)
library(purrr)
library(data.table)
library(dplyr)
library(qs)
library(sitar)
library(nlme)
library(ggplot2)
library(readxl)

############################################################
# 1) Prepare height data for downstream analysis
#
# Overview:
#   Standardize a raw height/vitals dataset so it can be used
#   consistently throughout the pipeline.
#
# Inputs:
#   vitals_height : data.frame or data.table
#     Raw longitudinal height data. Expected columns include:
#     - person_id
#     - height
#     - agedays
#
#   gender : data.frame or data.table
#     Person-level lookup table with at least:
#     - person_id
#     - gender
#
# Output:
#   data.table
#     Cleaned height dataset with standardized columns:
#     - person_id
#     - height
#     - agedays
#     - age
#     - gender
#
# Notes:
#   - age is calculated as agedays / 365.25
#   - rows with missing person_id, age, or height are removed
#   - rows with non-positive heights are removed
############################################################
prep_height_data <- function(vitals_height, gender) {
  dt <- as.data.table(vitals_height)
  dt[, age := agedays / 365.25]
  dt <- merge(dt, gender, by = "person_id", all.x = TRUE)
  
  keep_cols <- intersect(
    c("person_id", "height", "age", "agedays", "gender", "gender.x", "gender.y"),
    names(dt)
  )
  dt <- dt[, ..keep_cols]
  
  if ("gender" %in% names(dt)) {
    # keep as-is
  } else if ("gender.x" %in% names(dt)) {
    setnames(dt, "gender.x", "gender")
  } else if ("gender.y" %in% names(dt)) {
    setnames(dt, "gender.y", "gender")
  }
  
  # If both gender.x and gender.y made it through, keep the first usable one.
  if (sum(names(dt) == "gender") > 1) {
    dt <- dt[, !duplicated(names(dt)), with = FALSE]
  }
  
  dt <- dt[!is.na(person_id) & !is.na(age) & !is.na(height)]
  dt <- dt[height > 0]
  
  dt[]
}

############################################################
# 2) Subset prepared height data to a condition and sex
#
# Overview:
#   Return all height records for individuals with a specified
#   diagnosis. Optionally filter to one sex/gender.
#
# Inputs:
#   dx : data.frame or data.table
#     Diagnosis table with at least:
#     - person_id
#     - dID
#     - disease_short_title
#
#   dID : integer or character scalar
#     Disease identifier to subset.
#
#   height_all : data.table
#     Cleaned height data, typically from prep_height_data().
#
#   sex : character or NULL, default = NULL
#     Optional sex/gender filter (for example "M" or "F").
#     If NULL, no sex-based filtering is applied.
#
# Output:
#   data.table
#     Height records for the selected diagnosis, with added:
#     - dID
#     - dx
############################################################
get_height_dx_data <- function(dx, dID, height_all, sex = NULL, subset_ids = NULL) {
  affected_ids <- unique(dx$person_id[dx$dID == dID])
  
  if (!is.null(subset_ids)) {
    affected_ids <- intersect(affected_ids, subset_ids)
  }
  
  disease_name <- unique(dx$disease_short_title[dx$dID == dID])
  
  dt <- copy(height_all[person_id %in% affected_ids])
  
  if (!is.null(sex)) {
    dt <- dt[gender == sex]
  }
  
  dt[, dID := dID]
  dt[, dx := if (length(disease_name) > 0) disease_name[1] else NA_character_]
  
  dt[]
}

############################################################
# 3) Sample unaffected
#
# Overview:
#   Randomly sample unaffected individuals 
#   using a target unaffected-to-affected ratio, while also respecting
#   a minimum floor and maximum cap. An optional sex/gender
#   filter can be applied before sampling.
#
# Inputs:
#   height_no_dx : data.table
#     Height records for unaffected individuals
#
#   n_affected_ids : integer
#     Number of unique affected individuals.
#
#   ratio : numeric, default = 10
#     Desired unaffected-to-affected ratio.
#
#   max_unaffecteds : integer, default = 10000
#     Maximum number of unique unaffected individuals to sample.
#
#   min_unaffecteds : integer, default = 1000
#     Minimum number of unique unaffected individuals to target,
#     when available.
#
#   sex : character or NULL, default = NULL
#     Optional sex/gender filter (for example "M" or "F").
#     If NULL, unaffecteds are sampled from all available rows.
#
#   seed : integer, default = 123
#     Random seed for reproducibility.
#
# Output:
#   data.table
#     Sampled unaffected records with added columns:
#     - affected_status = 0
#     - group = "Unaffected"
#     - dID = NA
#     - dx = "Unaffected"
############################################################
sample_unaffecteds <- function(height_no_dx,
                            n_affected_ids,
                            ratio = 10,
                            max_unaffecteds = 10000,
                            min_unaffecteds = 1000,
                            sex = NULL,
                            seed = 123) {
  set.seed(seed)
  
  unaffecteds_pool <- copy(height_no_dx)
  
  if (!is.null(sex)) {
    unaffecteds_pool <- unaffecteds_pool[gender == sex]
  }
  
  unaffected_ids_all <- unique(unaffecteds_pool$person_id)
  
  if (length(unaffected_ids_all) == 0) {
    stop("No eligible unaffected IDs found after filtering.")
  }
  
  n_unaffecteds_target <- min(
    max_unaffecteds,
    max(min_unaffecteds, ratio * n_affected_ids),
    length(unaffected_ids_all)
  )
  
  sampled_ids <- sample(unaffected_ids_all, size = n_unaffecteds_target, replace = FALSE)
  unaffecteds <- copy(unaffecteds_pool[person_id %in% sampled_ids])
  
  unaffecteds[, affected_status := 0L]
  unaffecteds[, group := "Unaffected"]
  unaffecteds[, dID := NA_integer_]
  unaffecteds[, dx := "Unaffected"]
  
  unaffecteds[]
}

############################################################
# 4) Combine one condition with sampled unaffecteds
#
# Overview:
#   Build the analysis dataset for a single condition by
#   combining diagnosed affecteds with sampled unaffecteds
#   with optional sex/gender filtering applied to both groups.
#   Basic trajectory filters are then applied before modeling.
#
# Inputs:
#   dID : integer or character scalar
#     Disease identifier to analyze.
#
#   dx : data.frame or data.table
#     Diagnosis lookup table.
#
#   height_all : data.table
#     Prepared full height dataset.
#
#   height_no_dx : data.table
#     Prepared unaffected dataset.
#
#   sex : character or NULL, default = NULL
#     Optional sex/gender filter (for example "M" or "F").
#     If NULL, both affecteds and unaffecteds are drawn from all sexes.
#
#   unaffected_ratio, max_unaffecteds, min_unaffecteds, seed :
#     control sampling parameters
#
#   age_min, age_max : numeric age bounds for inclusion
#   min_obs_per_id : minimum number of observations per person
#
# Output:
#   data.table
#     affected-unaffected analysis dataset with:
#     - affected_status
#     - group
#     - dID
#     - dx
#     - log_age
############################################################
make_affected_unaffected_data <- function(dID,
                                          dx,
                                          height_all,
                                          height_no_dx,
                                          sex = NULL,
                                          subset_ids = NULL,
                                          condition_label = NULL,
                                          unaffected_ratio = 10,
                                          max_unaffecteds = 10000,
                                          min_unaffecteds = 1000,
                                          seed = 123,
                                          age_min = 2,
                                          age_max = 20,
                                          min_obs_per_id = 3) {
  
  is_unaffected_target <- (!is.null(condition_label) && condition_label == "Unaffected") || dID == 0
  
  if (is_unaffected_target) {
    affecteds <- copy(height_no_dx)
    
    if (!is.null(sex)) {
      affecteds <- affecteds[gender == sex]
    }
    
    if (!is.null(subset_ids)) {
      affecteds <- affecteds[person_id %in% subset_ids]
    }
    
    if (nrow(affecteds) == 0) {
      sex_msg <- if (is.null(sex)) "" else paste0(" for sex = '", sex, "'")
      stop(paste0("No unaffected records found", sex_msg, "."))
    }
    
    affecteds[, dID := dID]
    affecteds[, dx := "Unaffected"]
  } else {
    affecteds <- get_height_dx_data(
      dx = dx,
      dID = dID,
      height_all = height_all,
      sex = sex,
      subset_ids = subset_ids
    )
    
    if (nrow(affecteds) == 0) {
      sex_msg <- if (is.null(sex)) "" else paste0(" for sex = '", sex, "'")
      stop(paste0("No affected records found for this dID", sex_msg, "."))
    }
  }
  
  n_affected_ids <- uniqueN(affecteds$person_id)
  
  # Build comparison pool
  unaffecteds_pool <- copy(height_no_dx)
  
  if (!is.null(sex)) {
    unaffecteds_pool <- unaffecteds_pool[gender == sex]
  }
  
  # For the Unaffected target, remove the target IDs from the comparison pool
  if (is_unaffected_target) {
    target_ids <- unique(affecteds$person_id)
    unaffecteds_pool <- unaffecteds_pool[!person_id %in% target_ids]
  }
  
  unaffecteds <- sample_unaffecteds(
    height_no_dx = unaffecteds_pool,
    n_affected_ids = n_affected_ids,
    ratio = unaffected_ratio,
    max_unaffecteds = max_unaffecteds,
    min_unaffecteds = min_unaffecteds,
    sex = NULL,   # already filtered above
    seed = seed
  )
  
  affecteds[, affected_status := 1L]
  affecteds[, group := if (is_unaffected_target) "Unaffected target" else if ("dx" %in% names(affecteds) && !all(is.na(affecteds$dx))) as.character(affecteds$dx[1]) else "affected"]
  
  dat <- rbindlist(list(affecteds, unaffecteds), fill = TRUE, use.names = TRUE)
  
  dat <- dat[age >= age_min & age <= age_max]
  setorder(dat, person_id, age)
  dat <- dat[, .SD[1], by = .(person_id, age)]
  dat <- dat[, if (.N >= min_obs_per_id) .SD, by = person_id]
  dat[, log_age := log(age)]
  
  dat[]
}

make_unaffected_only_data <- function(height_no_dx,
                                      sex = NULL,
                                      subset_ids = NULL,
                                      age_min = 2,
                                      age_max = 20,
                                      min_obs_per_id = 3) {
  dat <- copy(height_no_dx)
  
  if (!is.null(sex)) {
    dat <- dat[gender == sex]
  }
  
  if (!is.null(subset_ids)) {
    dat <- dat[person_id %in% subset_ids]
  }
  
  if (nrow(dat) == 0) {
    sex_msg <- if (is.null(sex)) "" else paste0(" for sex = '", sex, "'")
    stop(paste0("No unaffected records found", sex_msg, "."))
  }
  
  dat[, affected_status := 0L]
  dat[, group := "Unaffected"]
  dat[, dID := NA_integer_]
  dat[, dx := "Unaffected"]
  
  dat <- dat[age >= age_min & age <= age_max]
  setorder(dat, person_id, age)
  dat <- dat[, .SD[1], by = .(person_id, age)]
  dat <- dat[, if (.N >= min_obs_per_id) .SD, by = person_id]
  dat[, log_age := log(age)]
  
  dat[]
}

############################################################
# 5) Fit a grid of candidate SITAR models
#
# Overview:
#   Compare candidate age scales and spline degrees of freedom
#   for a preliminary SITAR fit. All candidate models use
#   a + b + c random effects.
#
# Inputs:
#   dat : data.table
#     Analysis dataset containing height, age/log_age, and person_id.
#
#   age_var_candidates : character vector
#     Candidate age variables to use on the x-axis.
#
#   dfs : integer vector
#     Candidate spline degrees of freedom.
#
# Output:
#   list with:
#     - fits
#     - fit_summary
#     - fit_summary_ok
#     - best_name
#     - best_fit
############################################################
fit_sitar_candidate_grid <- function(dat,
                                     age_var_candidates = c("age", "log_age"),
                                     dfs = 4:9) {
  
  fits <- list()
  fit_rows <- list()
  idx <- 1
  
  for (xvar in age_var_candidates) {
    for (k in dfs) {
      
      fit_name <- paste0("x_", xvar, "_df", k)
      message("Fitting ", fit_name)
      
      dat_model <- copy(dat)
      dat_model[, x_use := get(xvar)]
      
      fit_obj <- tryCatch(
        sitar::sitar(
          x = x_use,
          y = height,
          id = person_id,
          data = dat_model,
          df = k,
          fixed = "a + b + c",
          random = "a + b + c",
          control = nlme::nlmeControl(
            maxIter = 500,
            msMaxIter = 500,
            pnlsMaxIter = 30,
            returnObject = TRUE
          )
        ),
        error = function(e) {
          message("  failed: ", e$message)
          NULL
        }
      )
      
      fits[[fit_name]] <- fit_obj
      
      fit_rows[[idx]] <- data.frame(
        model = fit_name,
        xvar = xvar,
        df = k,
        converged = !is.null(fit_obj),
        AIC = if (!is.null(fit_obj)) AIC(fit_obj) else NA_real_,
        BIC = if (!is.null(fit_obj)) BIC(fit_obj) else NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  fit_summary <- dplyr::bind_rows(fit_rows)
  
  fit_summary_ok <- fit_summary |>
    dplyr::filter(converged, is.finite(BIC)) |>
    dplyr::arrange(BIC)
  
  if (nrow(fit_summary_ok) == 0) {
    stop("No SITAR candidate models converged successfully.")
  }
  
  best_name <- fit_summary_ok$model[1]
  best_fit  <- fits[[best_name]]
  
  list(
    fits = fits,
    fit_summary = fit_summary,
    fit_summary_ok = fit_summary_ok,
    best_name = best_name,
    best_fit = best_fit
  )
}

############################################################
# 6) Flag SITAR outlier points
#
# Overview:
#   Append fitted values and Pearson residual diagnostics to the
#   analysis dataset and flag observations beyond a cutoff.
#
# Inputs:
#   fit : fitted SITAR model
#   dat : data.table used to fit the model
#   resid_cutoff : numeric, default = 4
#
# Output:
#   data.table with residual diagnostics and outlier_point flag
############################################################
flag_sitar_outliers <- function(fit, dat, resid_cutoff = 4) {
  out <- copy(dat)
  
  out[, pearson_resid := as.numeric(resid(fit, type = "pearson"))]
  out[, fitted_height := as.numeric(fitted(fit, level = 1))]
  out[, abs_pearson_resid := abs(pearson_resid)]
  out[, outlier_point := abs_pearson_resid > resid_cutoff]
  
  out[]
}

############################################################
# 7) Clean data using a preliminary SITAR model
#
# Overview:
#   Fit candidate preliminary SITAR models, choose the best one,
#   remove observations with large residuals, and then retain
#   individuals with enough remaining observations.
#
# Inputs:
#   dat : data.table
#   resid_cutoff : numeric, default = 4
#   min_obs_per_id_after : integer, default = 3
#
# Output:
#   list with:
#     - prelim
#     - flagged
#     - cleaned
############################################################
clean_with_prelim_sitar <- function(dat,
                                    resid_cutoff = 4,
                                    min_obs_per_id_after = 3) {
  
  prelim <- fit_sitar_candidate_grid(dat)
  prelim_fit <- prelim$best_fit
  
  flagged <- flag_sitar_outliers(prelim_fit, dat, resid_cutoff = resid_cutoff)
  
  cleaned <- flagged[outlier_point == FALSE]
  cleaned <- cleaned[, if (.N >= min_obs_per_id_after) .SD, by = person_id]
  cleaned[, log_age := log(age)]
  
  list(
    prelim = prelim,
    flagged = flagged,
    cleaned = cleaned
  )
}

############################################################
# 8) Fit final affected-unaffected SITAR models
#
# Overview:
#   Fit final SITAR models in which affected-unaffected status can shift
#   the a, b, and/or c parameters while unaffecteds define the main
#   curve shape. Models are compared by BIC and the best model can
#   optionally be refit with REML.
#
# Inputs:
#   dat_clean : data.table
#   xvar : character, one of the candidate age variables
#   prelim_df : integer, best df from preliminary screening
#   df_window : integer, search window around prelim_df
#   fixed_candidates : character vector of fixed-effect structures
#   refit_best_with_reml : logical
#
# Output:
#   list with fit objects, summaries, and best model information
############################################################
fit_final_affected_unaffected_sitar <- function(dat_clean,
                                         xvar,
                                         prelim_df,
                                         df_window = 1,
                                         fixed_candidates = c("a + b + c"),
                                         refit_best_with_reml = TRUE) {
  
  dat_model <- data.table::copy(dat_clean)
  dat_model[, x_use := get(xvar)]
  
  df_candidates <- sort(unique(pmax(1, prelim_df + (-df_window:df_window))))
  
  fits_ml <- list()
  fit_rows <- list()
  idx <- 1
  
  for (k in df_candidates) {
    for (fixed_str in fixed_candidates) {
      
      fit_name <- paste0(
        "x_", xvar,
        "_df", k,
        "_fixed_", gsub(" ", "", gsub("\\+", "_", fixed_str))
      )
      
      message("Fitting ", fit_name)
      
      fit_obj <- tryCatch({
        
        sitar_args <- list(
          x = quote(x_use),
          y = quote(height),
          id = quote(person_id),
          data = dat_model,
          df = k,
          fixed = fixed_str,
          random = "a + b + c",
          a.formula = ~ affected_status,
          method = "ML",
          control = nlme::nlmeControl(
            maxIter = 500,
            msMaxIter = 500,
            pnlsMaxIter = 30,
            returnObject = TRUE
          )
        )
        
        if (grepl("b", fixed_str, fixed = TRUE)) {
          sitar_args$b.formula <- ~ affected_status
        }
        if (grepl("c", fixed_str, fixed = TRUE)) {
          sitar_args$c.formula <- ~ affected_status
        }
        
        do.call(sitar::sitar, sitar_args)
        
      }, error = function(e) {
        message("  failed: ", e$message)
        NULL
      })
      
      fits_ml[[fit_name]] <- fit_obj
      
      fit_rows[[idx]] <- data.frame(
        model = fit_name,
        xvar = xvar,
        df = k,
        fixed = fixed_str,
        converged = !is.null(fit_obj),
        AIC = if (!is.null(fit_obj)) AIC(fit_obj) else NA_real_,
        BIC = if (!is.null(fit_obj)) BIC(fit_obj) else NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  fit_summary <- dplyr::bind_rows(fit_rows)
  
  fit_summary_ok <- fit_summary |>
    dplyr::filter(converged, is.finite(BIC)) |>
    dplyr::arrange(BIC)
  
  if (nrow(fit_summary_ok) == 0) {
    stop("No final affected-unaffected SITAR models converged successfully.")
  }
  
  best_row <- fit_summary_ok[1, ]
  best_name <- best_row$model
  best_fit_ml <- fits_ml[[best_name]]
  
  best_fit_final <- best_fit_ml
  
  if (isTRUE(refit_best_with_reml)) {
    best_fixed <- best_row$fixed[[1]]
    best_df <- best_row$df[[1]]
    
    best_fit_final <- tryCatch({
      
      sitar_args <- list(
        x = quote(x_use),
        y = quote(height),
        id = quote(person_id),
        data = dat_model,
        df = best_df,
        fixed = best_fixed,
        random = "a + b + c",
        a.formula = ~ affected_status,
        method = "REML",
        control = nlme::nlmeControl(
          maxIter = 500,
          msMaxIter = 500,
          pnlsMaxIter = 30,
          returnObject = TRUE
        )
      )
      
      if (grepl("b", best_fixed, fixed = TRUE)) {
        sitar_args$b.formula <- ~ affected_status
      }else{
        sitar_args$b.formula <- ~1
      }
      if (grepl("c", best_fixed, fixed = TRUE)) {
        sitar_args$c.formula <- ~ affected_status
      }else{
        sitar_args$c.formula <- ~1
      }
      
      do.call(sitar::sitar, sitar_args)
      
    }, error = function(e) {
      message("  REML refit failed, returning ML fit instead: ", e$message)
      best_fit_ml
    })
  }
  
  list(
    fits_ml = fits_ml,
    fit_summary = fit_summary,
    fit_summary_ok = fit_summary_ok,
    best_name = best_name,
    best_row = best_row,
    best_fit_ml = best_fit_ml,
    best_fit = best_fit_final
  )
}

fit_final_unaffected_sitar <- function(dat_clean,
                                       xvar,
                                       prelim_df,
                                       df_window = 1,
                                       fixed_candidates = c("a + b + c"),
                                       refit_best_with_reml = TRUE) {
  dat_model <- data.table::copy(dat_clean)
  dat_model[, x_use := get(xvar)]
  
  df_candidates <- sort(unique(pmax(1, prelim_df + (-df_window:df_window))))
  
  fits_ml <- list()
  fit_rows <- list()
  idx <- 1
  
  for (k in df_candidates) {
    for (fixed_str in fixed_candidates) {
      
      fit_name <- paste0(
        "x_", xvar,
        "_df", k,
        "_fixed_", gsub(" ", "", gsub("\\+", "_", fixed_str))
      )
      
      message("Fitting ", fit_name)
      
      fit_obj <- tryCatch({
        sitar::sitar(
          x = x_use,
          y = height,
          id = person_id,
          data = dat_model,
          df = k,
          fixed = fixed_str,
          random = "a + b + c",
          method = "ML",
          control = nlme::nlmeControl(
            maxIter = 500,
            msMaxIter = 500,
            pnlsMaxIter = 30,
            returnObject = TRUE
          )
        )
      }, error = function(e) {
        message("  failed: ", e$message)
        NULL
      })
      
      fits_ml[[fit_name]] <- fit_obj
      
      fit_rows[[idx]] <- data.frame(
        model = fit_name,
        xvar = xvar,
        df = k,
        fixed = fixed_str,
        converged = !is.null(fit_obj),
        AIC = if (!is.null(fit_obj)) AIC(fit_obj) else NA_real_,
        BIC = if (!is.null(fit_obj)) BIC(fit_obj) else NA_real_,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
  
  fit_summary <- dplyr::bind_rows(fit_rows)
  
  fit_summary_ok <- fit_summary |>
    dplyr::filter(converged, is.finite(BIC)) |>
    dplyr::arrange(BIC)
  
  if (nrow(fit_summary_ok) == 0) {
    stop("No final unaffected-only SITAR models converged successfully.")
  }
  
  best_row <- fit_summary_ok[1, ]
  best_name <- best_row$model
  best_fit_ml <- fits_ml[[best_name]]
  best_fit_final <- best_fit_ml
  
  if (isTRUE(refit_best_with_reml)) {
    best_fixed <- best_row$fixed[[1]]
    best_df <- best_row$df[[1]]
    
    best_fit_final <- tryCatch({
      sitar::sitar(
        x = x_use,
        y = height,
        id = person_id,
        data = dat_model,
        df = best_df,
        fixed = best_fixed,
        random = "a + b + c",
        method = "REML",
        control = nlme::nlmeControl(
          maxIter = 500,
          msMaxIter = 500,
          pnlsMaxIter = 30,
          returnObject = TRUE
        )
      )
    }, error = function(e) {
      message("  REML refit failed, returning ML fit instead: ", e$message)
      best_fit_ml
    })
  }
  
  list(
    fits_ml = fits_ml,
    fit_summary = fit_summary,
    fit_summary_ok = fit_summary_ok,
    best_name = best_name,
    best_row = best_row,
    best_fit_ml = best_fit_ml,
    best_fit = best_fit_final
  )
}
############################################################
# 9) Extract affected-vs-unaffected terms from a final SITAR fit
#
# Overview:
#   Return the subset of coefficient table rows that correspond
#   to affected-unaffected effects.
#
# Inputs:
#   final_fit : fitted SITAR model
#
# Output:
#   data.frame of coefficient rows involving affected_status
############################################################
extract_affected_unaffected_effects <- function(final_fit) {
  tt <- as.data.frame(summary(final_fit)$tTable)
  tt$term <- rownames(tt)
  rownames(tt) <- NULL
  
  tt[dplyr::if_else(grepl("affected_status", tt$term), TRUE, FALSE), ]
}

############################################################
# 10) Predict mean affected and unaffected curves from one fit
#
# Overview:
#   Generate population-level predicted curves for unaffected and
#   affected groups from the same fitted affected-unaffected SITAR model.
#
# Inputs:
#   final_fit : fitted SITAR model
#   ages_grid : numeric vector of ages for prediction
#   xvar : character, either "age" or "log_age"
#
# Output:
#   data.frame with columns: age, height, group
############################################################
predict_affected_unaffected_curves <- function(final_fit,
                                        ages_grid = seq(2, 20, by = 0.1),
                                        xvar = "log_age") {
  
  new0 <- data.frame(age = ages_grid, affected_status = 0L)
  new1 <- data.frame(age = ages_grid, affected_status = 1L)
  
  if (xvar == "log_age") {
    new0$log_age <- log(new0$age)
    new1$log_age <- log(new1$age)
    new0$x_use <- new0$log_age
    new1$x_use <- new1$log_age
  } else if (xvar == "age") {
    new0$x_use <- new0$age
    new1$x_use <- new1$age
  } else {
    stop("xvar must be 'age' or 'log_age'")
  }
  
  pred0 <- data.frame(
    age = ages_grid,
    height = predict(final_fit, newdata = new0, level = 0),
    group = "Unaffected"
  )
  
  pred1 <- data.frame(
    age = ages_grid,
    height = predict(final_fit, newdata = new1, level = 0),
    group = "affected"
  )
  
  rbind(pred0, pred1)
}

############################################################
# 11) Run the full condition pipeline end-to-end
#
# Overview:
#   Build the affected-unaffected dataset for one condition, optionally
#   restrict to one sex/gender, clean it using a preliminary
#   SITAR model, select the best preliminary age scale and df,
#   and then fit the final affected-unaffected SITAR model.
#
# Inputs:
#   dID : integer or character scalar
#   dx : diagnosis lookup table
#   height_all : prepared full height dataset
#   height_no_dx : prepared unaffected unaffected dataset
#   sex : character or NULL, default = NULL
#     Optional sex/gender filter for both affecteds and unaffecteds.
#   unaffected_ratio, max_unaffecteds, min_unaffecteds, seed :
#     unaffected sampling settings
#   age_min, age_max : age restrictions
#   min_obs_per_id : minimum trajectory length
#   resid_cutoff : residual threshold for outlier removal
#
# Output:
#   list containing raw data, cleaned data, preliminary fit,
#   final fit, and associated summaries
############################################################
run_condition_pipeline <- function(dID,
                                   dx,
                                   height_all,
                                   height_no_dx,
                                   sex = NULL,
                                   subset_ids = NULL,
                                   condition_label = NULL,
                                   unaffected_ratio = 10,
                                   max_unaffecteds = 10000,
                                   min_unaffecteds = 1000,
                                   seed = 123,
                                   age_min = 2,
                                   age_max = 20,
                                   min_obs_per_id = 3,
                                   resid_cutoff = 4) {
  
  dat0 <- make_affected_unaffected_data(
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
    min_obs_per_id = min_obs_per_id
  )
  
  cleaned_obj <- clean_with_prelim_sitar(
    dat = dat0,
    resid_cutoff = resid_cutoff,
    min_obs_per_id_after = min_obs_per_id
  )
  
  best_row <- cleaned_obj$prelim$fit_summary_ok[1, ]
  best_xvar <- best_row$xvar
  best_df   <- best_row$df
  
  final_obj <- fit_final_affected_unaffected_sitar(
    dat_clean = cleaned_obj$cleaned,
    xvar = best_xvar,
    prelim_df = best_df,
    df_window = 1,
    fixed_candidates = c("a + b + c"),
    refit_best_with_reml = TRUE
  )
  
  final_fit <- final_obj$best_fit
  final_fit_summary <- final_obj$fit_summary
  final_fit_summary_ok <- final_obj$fit_summary_ok
  
  list(
    raw_data = dat0,
    prelim_flagged = cleaned_obj$flagged,
    cleaned_data = cleaned_obj$cleaned,
    prelim_fit = cleaned_obj$prelim$best_fit,
    prelim_fit_summary = cleaned_obj$prelim$fit_summary,
    prelim_fit_summary_ok = cleaned_obj$prelim$fit_summary_ok,
    best_xvar = best_xvar,
    best_df = best_df,
    final_fit = final_fit,
    final_fit_summary = final_fit_summary,
    final_fit_summary_ok = final_fit_summary_ok
  )
}

run_unaffected_pipeline <- function(height_no_dx,
                                    sex = NULL,
                                    subset_ids = NULL,
                                    age_min = 2,
                                    age_max = 20,
                                    min_obs_per_id = 3,
                                    resid_cutoff = 4) {
  
  dat0 <- make_unaffected_only_data(
    height_no_dx = height_no_dx,
    sex = sex,
    subset_ids = subset_ids,
    age_min = age_min,
    age_max = age_max,
    min_obs_per_id = min_obs_per_id
  )
  
  cleaned_obj <- clean_with_prelim_sitar(
    dat = dat0,
    resid_cutoff = resid_cutoff,
    min_obs_per_id_after = min_obs_per_id
  )
  
  best_row <- cleaned_obj$prelim$fit_summary_ok[1, ]
  best_xvar <- best_row$xvar
  best_df <- best_row$df
  
  final_obj <- fit_final_unaffected_sitar(
    dat_clean = cleaned_obj$cleaned,
    xvar = best_xvar,
    prelim_df = best_df,
    df_window = 1,
    fixed_candidates = c("a + b + c"),
    refit_best_with_reml = TRUE
  )
  
  list(
    raw_data = dat0,
    prelim_flagged = cleaned_obj$flagged,
    cleaned_data = cleaned_obj$cleaned,
    prelim_fit = cleaned_obj$prelim$best_fit,
    prelim_fit_summary = cleaned_obj$prelim$fit_summary,
    prelim_fit_summary_ok = cleaned_obj$prelim$fit_summary_ok,
    best_xvar = best_xvar,
    best_df = best_df,
    final_fit = final_obj$best_fit,
    final_fit_summary = final_obj$fit_summary,
    final_fit_summary_ok = final_obj$fit_summary_ok
  )
}
############################################################
# 12) Extract a, b, and c estimates for affected-unaffected comparison
#
# Overview:
#   Convert the fitted SITAR coefficient table into a tidy summary
#   of unaffected estimates, affected estimates, and affected-unaffected
#   differences for the a, b, and c parameters.
#
# Inputs:
#   fit : fitted SITAR model
#
# Output:
#   data.frame with estimates, SEs, and 95% CIs for:
#   - unaffected
#   - affected
#   - Difference
############################################################

############################################################
# Backward-compatible wrappers
#
# These wrappers preserve the old male-specific entry points
# while delegating to the new generalized functions.
############################################################
get_height_dx_male <- function(dx, dID, height_all) {
  get_height_dx_data(dx = dx, dID = dID, height_all = height_all, sex = "M")
}

sample_male_unaffecteds <- function(height_no_dx_male,
                                 n_affected_ids,
                                 ratio = 10,
                                 max_unaffecteds = 10000,
                                 min_unaffecteds = 1000,
                                 seed = 123) {
  sample_unaffecteds(
    height_no_dx = height_no_dx_male,
    n_affected_ids = n_affected_ids,
    ratio = ratio,
    max_unaffecteds = max_unaffecteds,
    min_unaffecteds = min_unaffecteds,
    sex = "M",
    seed = seed
  )
}

make_affected_unaffected_male_data <- function(dID,
                                        dx,
                                        height_all,
                                        height_no_dx_male,
                                        unaffected_ratio = 10,
                                        max_unaffecteds = 10000,
                                        min_unaffecteds = 1000,
                                        seed = 123,
                                        age_min = 2,
                                        age_max = 20,
                                        min_obs_per_id = 3) {
  make_affected_unaffected_data(
    dID = dID,
    dx = dx,
    height_all = height_all,
    height_no_dx = height_no_dx_male,
    sex = "M",
    unaffected_ratio = unaffected_ratio,
    max_unaffecteds = max_unaffecteds,
    min_unaffecteds = min_unaffecteds,
    seed = seed,
    age_min = age_min,
    age_max = age_max,
    min_obs_per_id = min_obs_per_id
  )
}

run_male_condition_pipeline <- function(dID,
                                        dx,
                                        height_all,
                                        height_no_dx_male,
                                        unaffected_ratio = 10,
                                        max_unaffecteds = 10000,
                                        min_unaffecteds = 1000,
                                        seed = 123,
                                        age_min = 2,
                                        age_max = 20,
                                        min_obs_per_id = 3,
                                        resid_cutoff = 4) {
  run_condition_pipeline(
    dID = dID,
    dx = dx,
    height_all = height_all,
    height_no_dx = height_no_dx_male,
    sex = "M",
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

extract_sitar_abc <- function(fit) {
  
  tt <- as.data.frame(summary(fit)$tTable)
  tt$term <- rownames(tt)
  rownames(tt) <- NULL
  
  vc <- tryCatch(vcov(fit), error = function(e) NULL)
  
  ci <- function(est, se) {
    c(est - 1.96 * se, est + 1.96 * se)
  }
  
  get_param <- function(param) {
    
    intercept_row <- tt[tt$term == param, , drop = FALSE]
    affected_row  <- tt[tt$term == paste0(param, ".affected_status"), , drop = FALSE]
    
    # Case 1: parameter not estimated at all in this model
    if (nrow(intercept_row) == 0) {
      return(data.frame(
        parameter = param,
        group = c("unaffected", "affected", "Difference"),
        estimate = NA_real_,
        se = NA_real_,
        lower = NA_real_,
        upper = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    
    est_ctrl <- intercept_row$Value[1]
    se_ctrl  <- intercept_row$Std.Error[1]
    
    # Case 2: no affected_status term for this parameter
    # then affected == unaffected and difference is NA
    if (nrow(affected_row) == 0) {
      ctrl_ci <- ci(est_ctrl, se_ctrl)
      return(data.frame(
        parameter = param,
        group = c("unaffected", "affected", "Difference"),
        estimate = c(est_ctrl, est_ctrl, NA_real_),
        se = c(se_ctrl, se_ctrl, NA_real_),
        lower = c(ctrl_ci[1], ctrl_ci[1], NA_real_),
        upper = c(ctrl_ci[2], ctrl_ci[2], NA_real_),
        stringsAsFactors = FALSE
      ))
    }
    
    # Case 3: both intercept and affected effect estimated
    est_diff <- affected_row$Value[1]
    se_diff  <- affected_row$Std.Error[1]
    est_affected <- est_ctrl + est_diff
    
    if (!is.null(vc) &&
        param %in% rownames(vc) &&
        paste0(param, ".affected_status") %in% rownames(vc)) {
      
      var_affected <- vc[param, param] +
        vc[paste0(param, ".affected_status"), paste0(param, ".affected_status")] +
        2 * vc[param, paste0(param, ".affected_status")]
      
      se_affected <- sqrt(var_affected)
    } else {
      se_affected <- NA_real_
    }
    
    ctrl_ci <- ci(est_ctrl, se_ctrl)
    diff_ci <- ci(est_diff, se_diff)
    
    if (is.na(se_affected)) {
      aff_ci <- c(NA_real_, NA_real_)
    } else {
      aff_ci <- ci(est_affected, se_affected)
    }
    
    data.frame(
      parameter = param,
      group = c("unaffected", "affected", "Difference"),
      estimate = c(est_ctrl, est_affected, est_diff),
      se = c(se_ctrl, se_affected, se_diff),
      lower = c(ctrl_ci[1], aff_ci[1], diff_ci[1]),
      upper = c(ctrl_ci[2], aff_ci[2], diff_ci[2]),
      stringsAsFactors = FALSE
    )
  }
  
  dplyr::bind_rows(
    get_param("a"),
    get_param("b"),
    get_param("c")
  )
}
# # ------------------------------------------------------------
# # 13. GAMLSS fitting helper
# #
# # Description:
# #   Fits one candidate GAMLSS model for height as a function of age
# #   using specified distribution family and age/height transformations.
# #
# # Inputs:
# #   data : data.frame or tibble
# #       Must contain columns age and height.
# #   family_name : character
# #       One of "NO", "BCCGo", "BCPEo", "BCTo".
# #   x_tf_name : character
# #       Name of age transformation in x_transforms.
# #   y_tf_name : character
# #       Name of height transformation in y_transforms.
# #
# # Output:
# #   tibble with model metadata, BIC, fitted model object, and
# #   transformed fitting data. Returns NULL if fitting fails.
# # ------------------------------------------------------------
# fit_one_model <- function(data, family_name, x_tf_name, y_tf_name) {
#   x_tf <- x_transforms[[x_tf_name]]
#   y_tf <- y_transforms[[y_tf_name]]
#   
#   dd <- data %>%
#     mutate(
#       x = x_tf(age),
#       y = y_tf$forward(height)
#     ) %>%
#     filter(is.finite(x), is.finite(y))
#   
#   fit <- tryCatch({
#     
#     if (family_name == "NO") {
#       gamlss(
#         formula  = y ~ pb(x, control = pb.control(method = "GCV")),
#         sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
#         family   = NO(),
#         data     = dd,
#         trace    = FALSE
#       )
#     } else if (family_name == "BCCGo") {
#       gamlss(
#         formula  = y ~ pb(x, control = pb.control(method = "GCV")),
#         sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
#         nu.fo    = ~ pb(x, control = pb.control(method = "GCV")),
#         family   = BCCGo(),
#         data     = dd,
#         trace    = FALSE
#       )
#     } else if (family_name == "BCPEo") {
#       gamlss(
#         formula  = y ~ pb(x, control = pb.control(method = "GCV")),
#         sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
#         nu.fo    = ~ pb(x, control = pb.control(method = "GCV")),
#         tau.fo   = ~ pb(x, control = pb.control(method = "GCV")),
#         family   = BCPEo(),
#         data     = dd,
#         trace    = FALSE
#       )
#     } else if (family_name == "BCTo") {
#       gamlss(
#         formula  = y ~ pb(x, control = pb.control(method = "GCV")),
#         sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
#         nu.fo    = ~ pb(x, control = pb.control(method = "GCV")),
#         tau.fo   = ~ pb(x, control = pb.control(method = "GCV")),
#         family   = BCTo(),
#         data     = dd,
#         trace    = FALSE
#       )
#     } else {
#       stop("Unknown family")
#     }
#     
#   }, error = function(e) {
#     message(
#       "Model failed: ", family_name,
#       " / x=", x_tf_name,
#       " / y=", y_tf_name,
#       " / error=", e$message
#     )
#     NULL
#   })
#   
#   if (is.null(fit)) return(NULL)
#   
#   tibble(
#     family = family_name,
#     x_transform = x_tf_name,
#     y_transform = y_tf_name,
#     n = nrow(dd),
#     bic = GAIC(fit, k = log(nrow(dd))),
#     fit = list(fit),
#     fit_data = list(dd)
#   )
# }

# ------------------------------------------------------------
# 13. GAMLSS fitting helper
#
# Description:
#   Fits one candidate GAMLSS model for height as a function of age
#   using specified distribution family and age/height transformations.
#
# Inputs:
#   data : data.frame or tibble
#       Must contain columns age and height.
#   family_name : character
#       One of "NO", "BCCGo", "BCPEo", "BCTo".
#   x_tf_name : character
#       Name of age transformation in x_transforms.
#   y_tf_name : character
#       Name of height transformation in y_transforms.
#
# Output:
#   tibble with model metadata, BIC, fitted model object, and
#   transformed fitting data. Returns NULL if fitting fails.
# ------------------------------------------------------------
fit_one_model <- function(data, family_name, x_tf_name, y_tf_name) {
  x_tf <- x_transforms[[x_tf_name]]
  y_tf <- y_transforms[[y_tf_name]]
  
  dd <- data %>%
    mutate(
      x = x_tf(age),
      y = y_tf$forward(height)
    ) %>%
    filter(is.finite(x), is.finite(y))
  
  if (nrow(dd) == 0) {
    message(
      "Model failed: ", family_name,
      " / x=", x_tf_name,
      " / y=", y_tf_name,
      " / error=no finite rows after transformation"
    )
    return(NULL)
  }
  
  fit <- tryCatch({
    
    fit_call <- switch(
      family_name,
      
      "NO" = substitute(
        gamlss(
          formula  = y ~ pb(x, control = pb.control(method = "GCV")),
          sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
          family   = NO(),
          data     = DATA_OBJ,
          trace    = FALSE
        ),
        list(DATA_OBJ = dd)
      ),
      
      "BCCGo" = substitute(
        gamlss(
          formula  = y ~ pb(x, control = pb.control(method = "GCV")),
          sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
          nu.fo    = ~ pb(x, control = pb.control(method = "GCV")),
          family   = BCCGo(),
          data     = DATA_OBJ,
          trace    = FALSE
        ),
        list(DATA_OBJ = dd)
      ),
      
      "BCPEo" = substitute(
        gamlss(
          formula  = y ~ pb(x, control = pb.control(method = "GCV")),
          sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
          nu.fo    = ~ pb(x, control = pb.control(method = "GCV")),
          tau.fo   = ~ pb(x, control = pb.control(method = "GCV")),
          family   = BCPEo(),
          data     = DATA_OBJ,
          trace    = FALSE
        ),
        list(DATA_OBJ = dd)
      ),
      
      "BCTo" = substitute(
        gamlss(
          formula  = y ~ pb(x, control = pb.control(method = "GCV")),
          sigma.fo = ~ pb(x, control = pb.control(method = "GCV")),
          nu.fo    = ~ pb(x, control = pb.control(method = "GCV")),
          tau.fo   = ~ pb(x, control = pb.control(method = "GCV")),
          family   = BCTo(),
          data     = DATA_OBJ,
          trace    = FALSE
        ),
        list(DATA_OBJ = dd)
      ),
      
      stop("Unknown family")
    )
    
    eval(fit_call)
    
  }, error = function(e) {
    message(
      "Model failed: ", family_name,
      " / x=", x_tf_name,
      " / y=", y_tf_name,
      " / error=", conditionMessage(e)
    )
    NULL
  })
  
  if (is.null(fit)) return(NULL)
  
  tibble(
    family = family_name,
    x_transform = x_tf_name,
    y_transform = y_tf_name,
    n = nrow(dd),
    bic = GAIC(fit, k = log(nrow(dd))),
    fit = list(fit),
    fit_data = list(dd)
  )
}

############################################################
# 14) Predict centiles from a fitted GAMLSS model
#
# Overview:
#   Generate predicted centiles on the original height scale for a
#   fitted GAMLSS model across a supplied age grid.
#
# Inputs:
#   fit_obj : fitted GAMLSS model
#   fit_data : transformed fitting data used by predictAll()
#   ages : numeric vector of ages for prediction
#   family_name : character family label
#   x_tf_name : x transformation name
#   y_tf_name : y transformation name
#   centiles : numeric vector of centiles, e.g. c(3, 50, 97)
#
# Output:
#   data.frame with age and one column per requested centile
############################################################
predict_centiles <- function(fit_obj, fit_data, ages, family_name, x_tf_name, y_tf_name, centiles) {
  x_tf <- x_transforms[[x_tf_name]]
  y_tf <- y_transforms[[y_tf_name]]
  
  newdata <- data.frame(x = x_tf(ages))
  
  pars <- predictAll(
    fit_obj,
    newdata = newdata,
    data = fit_data,
    output = "list"
  )
  
  p <- centiles / 100
  
  qmat <- sapply(p, function(pp) {
    if (family_name == "NO") {
      qNO(pp, mu = pars$mu, sigma = pars$sigma)
    } else if (family_name == "BCCGo") {
      qBCCGo(pp, mu = pars$mu, sigma = pars$sigma, nu = pars$nu)
    } else if (family_name == "BCPEo") {
      qBCPEo(pp, mu = pars$mu, sigma = pars$sigma, nu = pars$nu, tau = pars$tau)
    } else if (family_name == "BCTo") {
      qBCTo(pp, mu = pars$mu, sigma = pars$sigma, nu = pars$nu, tau = pars$tau)
    } else {
      stop("Unknown family")
    }
  })
  
  colnames(qmat) <- paste0("C", centiles)
  
  qmat_orig <- apply(qmat, 2, y_tf$inverse)
  
  data.frame(age = ages, qmat_orig, check.names = FALSE)
}
