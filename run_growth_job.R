#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(qs)
  library(readxl)
})

source("Utils.R")
source("generate_growth_curves_with_diagnostics.R")
source("generate_growth_centiles_with_diagnostics.R")

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript run_growth_job.R <job_index>")
}

job_index <- as.integer(args[1])

if (is.na(job_index) || job_index < 1) {
  stop("job_index must be a positive integer.")
}

message("Starting job_index = ", job_index)

# -----------------------------
# Load data
# -----------------------------
demos <- qread("CGDB_demos.qs")
dx <- qread("CGDB_dx.qs")
vitals_height <- qread("peds_study_biometrics_all.qs")
growth_terms <- qread("GrowthCurves_dx_growth_terms.qs")  # loaded in case downstream code expects it
conditions_list <- read_xlsx("List_of_Conditions.xlsx")
CF_variant_function <- read.csv(
  "12-17-24_CFTR_variant_function_type.csv",
  header = TRUE,
  stringsAsFactors = FALSE
)
MM_id <- unique(CF_variant_function$person_id[CF_variant_function$funcs %in% "M/M"])
MR_id <- unique(CF_variant_function$person_id[CF_variant_function$funcs %in% "M/R"])
gender <- as.data.table(demos[, c("person_id", "gender")])


# -----------------------------
# Validate conditions list
# -----------------------------
required_cols <- c("dID", "Sex", "condition_label", "variantClass")
missing_cols <- setdiff(required_cols, names(conditions_list))
if (length(missing_cols) > 0) {
  stop("conditions_list is missing required columns: ",
       paste(missing_cols, collapse = ", "))
}

n_jobs <- nrow(conditions_list)
if (job_index > n_jobs) {
  stop("job_index = ", job_index, " exceeds nrow(conditions_list) = ", n_jobs)
}

job <- conditions_list[job_index, ]

job_dID <- job$dID[[1]]
job_sex <- toupper(trimws(job$Sex[[1]]))
job_condition <- trimws(job$condition_label[[1]])
job_variantClass <- trimws(job$variantClass[[1]])

if (job_sex %in% c("MALE")) job_sex <- "M"
if (job_sex %in% c("FEMALE")) job_sex <- "F"

if (!job_sex %in% c("M", "F")) {
  stop("Sex must resolve to 'M' or 'F'. Got: ", job_sex)
}

subset_ids <- NULL

if (job_condition == "Cystic fibrosis") {
  if (job_variantClass == "M/M") {
    subset_ids <- MM_id
  } else if (job_variantClass == "M/R") {
    subset_ids <- MR_id
  } else if (job_variantClass == "None") {
    subset_ids <- NULL
  }
}

if (is.na(job_variantClass) || !nzchar(job_variantClass)) {
  job_variantClass <- "None"
}

message("Running row ", job_index, " of ", n_jobs)
message("dID = ", job_dID)
message("Sex = ", job_sex)
message("condition_label = ", job_condition)
message("variantClass = ", job_variantClass)

# -----------------------------
# Prepare height datasets
# -----------------------------
height_all <- prep_height_data(vitals_height, gender)

no_dx_id <- setdiff(unique(height_all$person_id), unique(dx$person_id))
height_no_dx <- height_all[person_id %in% no_dx_id]
set.seed(3)
no_dx_id <- sample(unique(height_no_dx$person_id), 0.1 * length(no_dx_id))

if(job_condition == "Unaffected"){
  subset_ids <- no_dx_id
}

# -----------------------------
# Run SITAR
# -----------------------------
sitar_results <- generate_growth_curves(
  dID = job_dID,
  sex = job_sex,
  dx = dx,
  height_all = height_all,
  height_no_dx = height_no_dx,
  condition_label = job_condition,
  variantClass = job_variantClass,
  subset_ids = subset_ids,
  output_dir = "Results/SITAR",
  seed = 3
)
# -----------------------------
# Run GAMLSS
# -----------------------------
cohort_selector <- ifelse(job_condition == "Unaffected", "unaffected", "affected")

gamlss_results <- generate_growth_centiles(
  sitar_results = sitar_results,
  cohort = cohort_selector,
  condition_label = job_condition,
  sex = job_sex,
  variantClass = job_variantClass,
  output_dir = "Results/GAMLSS"
)

message("Completed job_index = ", job_index)

