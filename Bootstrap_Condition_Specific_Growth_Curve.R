#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
i <- as.integer(args[1])  # Get array index from SLURM

library(sitar)
library(data.table)
library(qs)

# Load input data
ages_grid <- seq(2, 20, by = 0.1)
load("height_no_dx_male_mod_df5_abc.RData") # Load final model for unaffected group
height_no_dx_male_mod <- list(height_no_dx_male_mod_df5_abc)
names(height_no_dx_male_mod) <- "mod_fabc_rabc"

demos <- qread("CGDB_demos.qs")
dx <- qread("CGDB_dx.qs")
vitals_height <- qread("peds_study_biometrics_all.qs")
gender <- as.data.table(demos[, c("person_id", "gender")])

# Utility functions
getPeak_test <- function(x, y = NULL, x0 = NULL, peak = TRUE, takeoff = FALSE) {
  xy <- xy.coords(x, y)
  xy <- unique(as.data.frame(xy[1:2])[order(xy$x), ])
  if (!is.null(x0)) {
    xy <- subset(xy, x >= x0)
    if (nrow(xy) < 3) return(c(x = NA, y = NA))
  }
  x <- xy$x
  y <- xy$y
  ddy <- diff(diff(y) > 0)
  tp <- which(ddy == -1) + 1
  if (length(tp) == 0 && (peak || takeoff)) return(c(x = NA, y = NA))
  if (length(tp) > 1) tp <- tp[which.max(y[tp])]
  if (!peak) {
    if (!takeoff) tp <- length(ddy)
    tp <- which(ddy[seq_len(tp)] == 1) + 1
    tp <- tp[which.min(y[tp])]
    if (length(tp) == 0) return(c(x = NA, y = NA))
  }
  n <- 0
  repeat {
    n <- n + 1
    if (tp == n || tp + n > nrow(xy)) break
    curve <- with(xy[(tp - n):(tp + n), ], lm(y ~ poly(x, 2, raw = TRUE)))
    if (curve$rank == 3) break
  }
  x_peak <- -curve$coef[[2]] / (2 * curve$coef[[3]])
  y_peak <- unname(predict(curve, data.frame(x = x_peak)))
  return(c(x = x_peak, y = y_peak))
}

get_Height_Dx_Data <- function(dx, dID){
  person_id_dx <- unique(dx$person_id[dx$dID %in% dID])
  vitals_height_dx <- vitals_height[which(vitals_height$person_id %in% person_id_dx), ]
  vitals_height_dx$age <- vitals_height_dx$agedays/365.25
  vitals_height_dx <- as.data.table(vitals_height_dx)
  vitals_height_dx <- merge(vitals_height_dx, gender, by = "person_id", all.x = TRUE)
  vitals_height_dx[gender.x == "M"]
}

# Prepare data
height_dx <- get_Height_Dx_Data(dx, dID = 190685)
ids_dx <- unique(height_dx$person_id)
id_var <- "person_id"

# Bootstrap iteration
sampled_ids_dx <- sample(ids_dx, size = length(ids_dx), replace = TRUE)
boot_data <- do.call(rbind, lapply(seq_along(sampled_ids_dx), function(j) {
  df <- height_dx[height_dx$person_id == sampled_ids_dx[j], ]
  df[[id_var]] <- paste0("boot", i, "_", j)
  df
}))

# Model with random effects on a and c was chosen as optimal in Generate_Condition_Specific_Growth_Chart.R
boot_model_dx <- try(
  sitar(x = age, y = height, id = person_id, data = boot_data, df = 5,
        fixed = "a + b + c", random = "list(id = pdSymm(a + c ~ 1))",
        control = nlme::nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)),
  silent = TRUE
)

if (!inherits(boot_model_dx, "try-error")) {
  pred_dx <- predict(boot_model_dx, newdata = data.frame(age = ages_grid), level = 0)
  deriv_dx <- predict(boot_model_dx, newdata = data.frame(age = ages_grid), level = 0, deriv = 1)
  peak_dx <- getPeak_test(ages_grid, deriv_dx, x0 = 5)
  
  pred_no_dx <- predict(height_no_dx_male_mod$mod_fabc_rabc, newdata = data.frame(age = ages_grid), level = 0)
  deriv_no_dx <- predict(height_no_dx_male_mod$mod_fabc_rabc, newdata = data.frame(age = ages_grid), level = 0, deriv = 1)
  peak_no_dx <- getPeak_test(ages_grid, deriv_no_dx, x0 = 5)
  
  res <- list(
    sizediff = max(pred_dx) - max(pred_no_dx),
    timediff = peak_dx[[1]] - peak_no_dx[[1]],
    veldiff = peak_dx[[2]] - peak_no_dx[[2]],
    success = TRUE
  )
} else {
  res <- list(success = FALSE)
}

saveRDS(res, file = sprintf("Height_Male_Down_Syndrome/Height_Male_Down_Syndrome_Boot_%04d.rds", i))