#######################################
#           Load Packages             #
#######################################
library("data.table")
library("qs")
library("lme4")
library("ggplot2")
library("dplyr")
library("mgcv")
library("sitar")
library("Hmisc")
library("growthcleanr")
library("doParallel")
library("foreach")
library("fda")
library("pracma")
library("tidyr")
library("boot")
library("tidyr")

# Please change the following path to your data directory
dataDir <- "./path/"

##############################################
#              Load EHR Data                 #
##############################################

######################################################
# `demos` is a data.frame with the following columns:
# person_id - Patient ID (integer)
# gender - Gender (factor - M or F)
######################################################
demos <- qread(file.path(dataDir, 'CGDB_demos.qs'))

######################################################
# `dx` is a data.frame with the following columns:
# person_id - Patient ID (integer)
# dID - Disorder ID (integer)
# disease_short_title - Name of the disorder (char)
######################################################
dx <- qread(file.path(dataDir, 'CGDB_dx.qs'))

######################################################
# `vitals_height` is a data.frame with the following columns:
# person_id - Patient ID (integer)
# agedays - Age in number of days (num)
# measurement_date - Date of measurement (Date)
# gender - Gender (factor - M or F)
# height - Height measurement in centimeters (num)
######################################################
vitals_height <- qread(file.path(dataDir, 'peds_study_biometrics_all.qs'))


##################################################
#              Utility Functions                 #
##################################################
# Function for extracting and formatting male height data for patients with a specific diagnosis ID (dID)
# Replace "M" with "F" for extracting data for females
get_Height_Dx_Data <- function(dx, dID){
  person_id_dx <- unique(dx$person_id[dx$dID %in% dID]); length(person_id_dx)
  vitals_height_dx <- vitals_height[which(vitals_height$person_id %in% person_id_dx), ]; length(unique(vitals_height_dx$person_id))
  vitals_height_dx$age <- vitals_height_dx$agedays/365.25
  vitals_height_dx <- as.data.table(vitals_height_dx)
  vitals_height_dx <- merge(vitals_height_dx, gender, by = "person_id", all.x = TRUE)
  height_dx_male <- vitals_height_dx[which(vitals_height_dx$gender.x %in% "M"), ]
  disease_short_title <- unique(dx$disease_short_title[dx$dID == dID])
  return(list(height_dx_male, disease_short_title = disease_short_title))
}

# Function for fitting SITAR model with spline (df = 3)
model_fit_df3 <- function(data, skip) {
  print("Starting model fit")
  
  # Initialize a named list to store model results
  models <- list()
  
  # Function to conditionally fit a model
  fit_model <- function(model_name, fit_code) {
    if (model_name %in% skip) {
      print(paste(model_name, "skipped"))
      return(NA)
    }
    print(paste("Fitting", model_name))
    tryCatch({
      fit_code
    }, error = function(e) {
      print(paste(model_name, "failed:", e$message))
      return(NA)
    })
  }
  
  # Fit models conditionally
  models$mod_fa_ra <- fit_model("mod_fa_ra", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fa_ra done")
  
  models$mod_fb_rb <- fit_model("mod_fb_rb", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "b", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fb_rb done")
  
  models$mod_fc_rc <- fit_model("mod_fc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fc_rc done")
  
  models$mod_fab_ra <- fit_model("mod_fab_ra", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_ra done")
  
  models$mod_fab_rb <- fit_model("mod_fab_rb", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_rb done")
  
  models$mod_fac_ra <- fit_model("mod_fac_ra", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + c", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_ra done")
  
  models$mod_fac_rc <- fit_model("mod_fac_rc", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_rc done")
  
  models$mod_fbc_rb <- fit_model("mod_fbc_rb", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "b + c", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rb done")
  
  models$mod_fbc_rc <- fit_model("mod_fbc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "b + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rc done")
  
  models$mod_fabc_ra <- fit_model("mod_fabc_ra", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b + c", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_ra done")
  
  models$mod_fabc_rb <- fit_model("mod_fabc_rb", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b + c", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rb done")
  
  models$mod_fabc_rc <- fit_model("mod_fabc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rc done")
  print("Ranef 1 done")
  
  models$mod_fab_rab <- fit_model("mod_fab_rab", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b", random = "list(id = pdSymm(a + b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_rab done")
  
  models$mod_fac_rac <- fit_model("mod_fac_rac", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + c", random = "list(id = pdSymm(a + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_rac done")
  
  models$mod_fbc_rbc <- fit_model("mod_fbc_rbc", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "b + c", random = "list(id = pdSymm(b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rbc done")
  
  models$mod_fabc_rab <- fit_model("mod_fabc_rab", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b + c", random = "list(id = pdSymm(a + b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rab done")
  
  models$mod_fabc_rac <- fit_model("mod_fabc_rac", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b + c", random = "list(id = pdSymm(a + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rac done")
  
  models$mod_fabc_rbc <- fit_model("mod_fabc_rbc", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b + c", random = "list(id = pdSymm(b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rbc done")  
  print("Ranef 2 done")
  
  models$mod_fabc_rabc <- fit_model("mod_fabc_rabc", sitar(x = age, y = height, id = person_id, data = data, df = 3, fixed = "a + b + c", random = "list(id = pdSymm(a + b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rabc done")  
  print("Ranef 3 done")
  
  # Calculate AIC and BIC, handling NAs
  AIC_values <- sapply(models, function(mod) {
    if (all(is.na(mod))) {
      return(NA)
    } else {
      return(AIC(mod))
    }
  })
  
  BIC_values <- sapply(models, function(mod) {
    if (all(is.na(mod))) {
      return(NA)
    } else {
      return(BIC(mod))
    }
  })
  
  # Combine into a data frame
  mod_fit <- data.frame(AIC = AIC_values, BIC = BIC_values)
  
  # Assign row names
  rownames(mod_fit) <- names(models)
  
  return(list(mod_fit = mod_fit, 
              mod_fa_ra = models$mod_fa_ra, 
              mod_fab_ra = models$mod_fab_ra, 
              mod_fab_rb = models$mod_fab_rb, mod_fac_ra = models$mod_fac_ra, mod_fac_rc = models$mod_fac_rc, 
              mod_fabc_ra = models$mod_fabc_ra, mod_fabc_rb = models$mod_fabc_rb, mod_fabc_rc = models$mod_fabc_rc,
              mod_fab_rab = models$mod_fab_rab, mod_fac_rac = models$mod_fac_rac, 
              mod_fabc_rab = models$mod_fabc_rab, mod_fabc_rac = models$mod_fabc_rac, mod_fabc_rbc = models$mod_fabc_rbc,
              mod_fabc_rabc = models$mod_fabc_rabc))
}

# Function for fitting SITAR model with spline (df = 4)
model_fit_df4 <- function(data, skip) {
  print("Starting model fit")
  
  # Initialize a named list to store model results
  models <- list()
  
  # Function to conditionally fit a model
  fit_model <- function(model_name, fit_code) {
    if (model_name %in% skip) {
      print(paste(model_name, "skipped"))
      return(NA)
    }
    print(paste("Fitting", model_name))
    tryCatch({
      fit_code
    }, error = function(e) {
      print(paste(model_name, "failed:", e$message))
      return(NA)
    })
  }
  
  # Fit models conditionally
  models$mod_fa_ra <- fit_model("mod_fa_ra", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fa_ra done")
  
  models$mod_fb_rb <- fit_model("mod_fb_rb", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "b", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fb_rb done")
  
  models$mod_fc_rc <- fit_model("mod_fc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fc_rc done")
  
  models$mod_fab_ra <- fit_model("mod_fab_ra", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_ra done")
  
  models$mod_fab_rb <- fit_model("mod_fab_rb", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_rb done")
  
  models$mod_fac_ra <- fit_model("mod_fac_ra", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + c", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_ra done")
  
  models$mod_fac_rc <- fit_model("mod_fac_rc", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_rc done")
  
  models$mod_fbc_rb <- fit_model("mod_fbc_rb", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "b + c", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rb done")
  
  models$mod_fbc_rc <- fit_model("mod_fbc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "b + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rc done")
  
  models$mod_fabc_ra <- fit_model("mod_fabc_ra", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b + c", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_ra done")
  
  models$mod_fabc_rb <- fit_model("mod_fabc_rb", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b + c", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rb done")
  
  models$mod_fabc_rc <- fit_model("mod_fabc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rc done")
  print("Ranef 1 done")
  
  models$mod_fab_rab <- fit_model("mod_fab_rab", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b", random = "list(id = pdSymm(a + b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_rab done")
  
  models$mod_fac_rac <- fit_model("mod_fac_rac", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + c", random = "list(id = pdSymm(a + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_rac done")
  
  models$mod_fbc_rbc <- fit_model("mod_fbc_rbc", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "b + c", random = "list(id = pdSymm(b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rbc done")
  
  models$mod_fabc_rab <- fit_model("mod_fabc_rab", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b + c", random = "list(id = pdSymm(a + b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rab done")
  
  models$mod_fabc_rac <- fit_model("mod_fabc_rac", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b + c", random = "list(id = pdSymm(a + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rac done")
  
  models$mod_fabc_rbc <- fit_model("mod_fabc_rbc", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b + c", random = "list(id = pdSymm(b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rbc done")  
  print("Ranef 2 done")
  
  models$mod_fabc_rabc <- fit_model("mod_fabc_rabc", sitar(x = age, y = height, id = person_id, data = data, df = 4, fixed = "a + b + c", random = "list(id = pdSymm(a + b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rabc done")  
  print("Ranef 3 done")
  
  # Calculate AIC and BIC, handling NAs
  AIC_values <- sapply(models, function(mod) {
    if (all(is.na(mod))) {
      return(NA)
    } else {
      return(AIC(mod))
    }
  })
  
  BIC_values <- sapply(models, function(mod) {
    if (all(is.na(mod))) {
      return(NA)
    } else {
      return(BIC(mod))
    }
  })
  
  # Combine into a data frame
  mod_fit <- data.frame(AIC = AIC_values, BIC = BIC_values)
  
  # Assign row names
  rownames(mod_fit) <- names(models)
  
  return(list(mod_fit = mod_fit, 
              mod_fa_ra = models$mod_fa_ra, 
              mod_fab_ra = models$mod_fab_ra, 
              mod_fab_rb = models$mod_fab_rb, mod_fac_ra = models$mod_fac_ra, mod_fac_rc = models$mod_fac_rc, 
              mod_fabc_ra = models$mod_fabc_ra, mod_fabc_rb = models$mod_fabc_rb, mod_fabc_rc = models$mod_fabc_rc,
              mod_fab_rab = models$mod_fab_rab, mod_fac_rac = models$mod_fac_rac, 
              mod_fabc_rab = models$mod_fabc_rab, mod_fabc_rac = models$mod_fabc_rac, mod_fabc_rbc = models$mod_fabc_rbc,
              mod_fabc_rabc = models$mod_fabc_rabc))
}

# Function for fitting SITAR model with spline (df = 5)
model_fit_df5 <- function(data, skip) {
  print("Starting model fit")
  
  # Initialize a named list to store model results
  models <- list()
  
  # Function to conditionally fit a model
  fit_model <- function(model_name, fit_code) {
    if (model_name %in% skip) {
      print(paste(model_name, "skipped"))
      return(NA)
    }
    print(paste("Fitting", model_name))
    tryCatch({
      fit_code
    }, error = function(e) {
      print(paste(model_name, "failed:", e$message))
      return(NA)
    })
  }
  
  # Fit models conditionally
  models$mod_fa_ra <- fit_model("mod_fa_ra", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fa_ra done")
  
  models$mod_fb_rb <- fit_model("mod_fb_rb", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "b", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fb_rb done")
  
  models$mod_fc_rc <- fit_model("mod_fc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fc_rc done")
  
  models$mod_fab_ra <- fit_model("mod_fab_ra", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_ra done")
  
  models$mod_fab_rb <- fit_model("mod_fab_rb", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_rb done")
  
  models$mod_fac_ra <- fit_model("mod_fac_ra", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + c", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_ra done")
  
  models$mod_fac_rc <- fit_model("mod_fac_rc", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_rc done")
  
  models$mod_fbc_rb <- fit_model("mod_fbc_rb", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "b + c", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rb done")
  
  models$mod_fbc_rc <- fit_model("mod_fbc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "b + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rc done")
  
  models$mod_fabc_ra <- fit_model("mod_fabc_ra", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b + c", random = "list(id = pdSymm(a ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_ra done")
  
  models$mod_fabc_rb <- fit_model("mod_fabc_rb", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b + c", random = "list(id = pdSymm(b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rb done")
  
  models$mod_fabc_rc <- fit_model("mod_fabc_rc", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b + c", random = "list(id = pdSymm(c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rc done")
  print("Ranef 1 done")
  
  models$mod_fab_rab <- fit_model("mod_fab_rab", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b", random = "list(id = pdSymm(a + b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fab_rab done")
  
  models$mod_fac_rac <- fit_model("mod_fac_rac", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + c", random = "list(id = pdSymm(a + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fac_rac done")
  
  models$mod_fbc_rbc <- fit_model("mod_fbc_rbc", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "b + c", random = "list(id = pdSymm(b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fbc_rbc done")
  
  models$mod_fabc_rab <- fit_model("mod_fabc_rab", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b + c", random = "list(id = pdSymm(a + b ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rab done")
  
  models$mod_fabc_rac <- fit_model("mod_fabc_rac", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b + c", random = "list(id = pdSymm(a + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rac done")
  
  models$mod_fabc_rbc <- fit_model("mod_fabc_rbc", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b + c", random = "list(id = pdSymm(b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rbc done")  
  print("Ranef 2 done")
  
  models$mod_fabc_rabc <- fit_model("mod_fabc_rabc", sitar(x = age, y = height, id = person_id, data = data, df = 5, fixed = "a + b + c", random = "list(id = pdSymm(a + b + c ~ 1))", control = nlmeControl(maxIter = 500, msMaxIter = 500, returnObject = TRUE)))
  # print("Model fabc_rabc done")  
  print("Ranef 3 done")
  
  # Calculate AIC and BIC, handling NAs
  AIC_values <- sapply(models, function(mod) {
    if (all(is.na(mod))) {
      return(NA)
    } else {
      return(AIC(mod))
    }
  })
  
  BIC_values <- sapply(models, function(mod) {
    if (all(is.na(mod))) {
      return(NA)
    } else {
      return(BIC(mod))
    }
  })
  
  # Combine into a data frame
  mod_fit <- data.frame(AIC = AIC_values, BIC = BIC_values)
  
  # Assign row names
  rownames(mod_fit) <- names(models)
  
  return(list(mod_fit = mod_fit, 
              mod_fa_ra = models$mod_fa_ra, 
              mod_fab_ra = models$mod_fab_ra, 
              mod_fab_rb = models$mod_fab_rb, mod_fac_ra = models$mod_fac_ra, mod_fac_rc = models$mod_fac_rc, 
              mod_fabc_ra = models$mod_fabc_ra, mod_fabc_rb = models$mod_fabc_rb, mod_fabc_rc = models$mod_fabc_rc,
              mod_fab_rab = models$mod_fab_rab, mod_fac_rac = models$mod_fac_rac, 
              mod_fabc_rab = models$mod_fabc_rab, mod_fabc_rac = models$mod_fabc_rac, mod_fabc_rbc = models$mod_fabc_rbc,
              mod_fabc_rabc = models$mod_fabc_rabc))
}

# Calculate the centiles for making the growth chart
# model - Final SITAR model
# n_ind - Number of replicates
calculate_centiles_for_plot <- function(model, n_ind) {
  fixed_effects <- fixed.effects(model)  # Population-level (fixed effects)
  residual_sd <- sigma(model)  # Residual standard deviation
  
  # Extract variance-covariance matrix
  vc_matrix <- VarCorr(model)
  
  # Identify random effect names
  re_names <- rownames(vc_matrix)
  re_names <- re_names[re_names != "Residual"]  # Exclude residual variance
  
  num_re <- length(re_names)
  
  if (num_re == 0) {
    stop("The model does not contain any random effects.")
  }
  
  # Extract variances and compute standard deviations
  variances <- as.numeric(vc_matrix[re_names, "Variance"])
  std_devs <- sqrt(variances)
  
  # Initialize covariance matrix
  cov_matrix <- diag(variances, nrow = num_re, ncol = num_re)
  
  # Extract correlations from the last column dynamically
  if (num_re > 1) {
    for (i in 1:(num_re - 1)) {
      for (j in (i + 1):num_re) {
        # Check if correlation exists in the "Corr" column
        corr_value <- suppressWarnings(as.numeric(vc_matrix[j, "Corr"]))
        
        if (!is.na(corr_value)) {
          cov_matrix[i, j] <- corr_value * std_devs[i] * std_devs[j]
          cov_matrix[j, i] <- cov_matrix[i, j]  # Ensure symmetry
        }
      }
    }
  }
  
  # Simulate random effects
  rand_effects <- mvrnorm(n = n_ind, mu = rep(0, num_re), Sigma = cov_matrix)
  
  # Assign names to random effects
  colnames(rand_effects) <- re_names
  
  # Simulate individual growth curves
  individual_curves <- sapply(1:n_ind, function(i) {
    print(i)
    mod <- model
    simulated_random_effects <- rand_effects[i, , drop = FALSE]
    mod$coefficients$random$id[1, ] <- simulated_random_effects
    # Ensure that newdata contains the random effect grouping factor; make predictions for the individual with simulated random effects
    individual_curve <- predict(mod, newdata = data.frame(age = ages_grid, person_id = rownames(mod$coefficients$random$id)[1]), type = "response", level = 1)
    individual_curve <- individual_curve + rnorm(1, mean = 0, sd = residual_sd)
    return(individual_curve)
  })
  
  centiles <- apply(individual_curves, 1, quantile, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95))
  return(centiles)
}

# Function to check monotonicity
check_Monotonicity <- function(df, cols) {
  for (col in cols) {
    for (i in 2:nrow(df)) {
      if (df[i, col] < df[i - 1, col]) {
        df[i, col] <- df[i - 1, col]
      }
    }
  }
  return(df)
}
cols_to_check <- c("c05", "c10", "c25", "c50", "c75", "c90", "c95")

# Age range
ages_grid <- seq(2, 20, by = 0.1)

########################################################################### 
#       Down syndrome dID = 190685
###########################################################################
# dID = 190685
height_dx_Down_syndrome <- get_Height_Dx_Data(dx = dx, dID = 190685)
height_dx_male_Down_syndrome <- height_dx_Down_syndrome[[1]]; length(unique(height_dx_male_Down_syndrome$person_id))
height_dx_male_Down_syndrome$dx <- height_dx_Down_syndrome$disease_short_title

########################################################################### 
#       Model Fitting, Selection, and Sensitivity Analysis
###########################################################################
# Check model fit, perform LRT, check model fit to select optimal model
height_dx_male_Down_syndrome_mod_df3 <- model_fit_df3(height_dx_male_Down_syndrome, skip = NA)
height_dx_male_Down_syndrome_mod_df4 <- model_fit_df4(height_dx_male_Down_syndrome, skip = NA)
height_dx_male_Down_syndrome_mod_df5 <- model_fit_df5(height_dx_male_Down_syndrome, skip = NA)

# Random effect = a
mod3_a <- height_dx_male_Down_syndrome_mod_df3$mod_fabc_ra
mod4_a <- height_dx_male_Down_syndrome_mod_df4$mod_fabc_ra
mod5_a <- height_dx_male_Down_syndrome_mod_df5$mod_fabc_ra
summary(mod3_a) 
summary(mod4_a)
summary(mod5_a)
anova(mod3_a, mod4_a)
anova(mod3_a, mod5_a)
anova(mod4_a, mod5_a) # mod5 has best fit

# Random effect = b
mod3_b <- height_dx_male_Down_syndrome_mod_df3$mod_fabc_rb
mod4_b <- height_dx_male_Down_syndrome_mod_df4$mod_fabc_rb
mod5_b <- height_dx_male_Down_syndrome_mod_df5$mod_fabc_rb
summary(mod3_b) 
summary(mod4_b)
summary(mod5_b)
anova(mod3_b, mod4_b)
anova(mod3_b, mod5_b)
anova(mod4_b, mod5_b) # mod5 has best fit

# Random effect = c
mod3_c <- height_dx_male_Down_syndrome_mod_df3$mod_fabc_rc
mod4_c <- height_dx_male_Down_syndrome_mod_df4$mod_fabc_rc
mod5_c <- height_dx_male_Down_syndrome_mod_df5$mod_fabc_rc
summary(mod3_c) # SE inflated
summary(mod4_c) # Did not converge
summary(mod5_c) # SE inflated
anova(mod3_c, mod4_c)
anova(mod3_c, mod5_c) # mod5 has best fit, but the SE is inflated. So this model will not be considered further.

# Random effect = a, b
mod3_ab <- height_dx_male_Down_syndrome_mod_df3$mod_fabc_rab
mod4_ab <- height_dx_male_Down_syndrome_mod_df4$mod_fabc_rab
mod5_ab <- height_dx_male_Down_syndrome_mod_df5$mod_fabc_rab
summary(mod3_ab) # SE inflated
summary(mod4_ab)
summary(mod5_ab)
anova(mod3_ab, mod4_ab)
anova(mod3_ab, mod5_ab)
anova(mod4_ab, mod5_ab) # mod5 has best fit

# Random effect = a, c
mod3_ac <- height_dx_male_Down_syndrome_mod_df3$mod_fabc_rac
mod4_ac <- height_dx_male_Down_syndrome_mod_df4$mod_fabc_rac
mod5_ac <- height_dx_male_Down_syndrome_mod_df5$mod_fabc_rac
summary(mod3_ac) # SE inflated
summary(mod4_ac)
summary(mod5_ac)
anova(mod3_ac, mod4_ac)
anova(mod3_ac, mod5_ac)
anova(mod4_ac, mod5_ac) # mod5 has best fit

# Random effect = b, c
mod3_bc <- height_dx_male_Down_syndrome_mod_df3$mod_fabc_rbc
mod4_bc <- height_dx_male_Down_syndrome_mod_df4$mod_fabc_rbc
mod5_bc <- height_dx_male_Down_syndrome_mod_df5$mod_fabc_rbc
summary(mod3_bc) 
summary(mod4_bc)
summary(mod5_bc)
anova(mod3_bc, mod4_bc)
anova(mod3_bc, mod5_bc)
anova(mod4_bc, mod5_bc) # mod5 has best fit

# Random effect = a, b, c
mod3_abc <- height_dx_male_Down_syndrome_mod_df3$mod_fabc_rabc
mod4_abc <- height_dx_male_Down_syndrome_mod_df4$mod_fabc_rabc
mod5_abc <- height_dx_male_Down_syndrome_mod_df5$mod_fabc_rabc
summary(mod3_abc) # SE inflated
summary(mod4_abc)
summary(mod5_abc)
anova(mod3_abc, mod4_abc)
anova(mod3_abc, mod5_abc)
anova(mod4_abc, mod5_abc) # mod5 has best fit

# AIC to compare non-nested models
AIC(mod5_a, mod5_b) # mod5_b has better fit
AIC(mod5_b, mod5_c) # mod5_b has better fit
AIC(mod5_b, mod5_ab) # mod5_ab has better fit
AIC(mod5_ab, mod5_ac) # mod5_ac has better fit
AIC(mod5_ac, mod5_bc) # mod5_ac has better fit
AIC(mod5_ac, mod5_abc) # mod5_abc has better fit

# Generate the spaghetti plot to check for adequacy of fit
# against individual trajectories
# Based on visual assessment, mod5_ac has better fit 
predictions <- data.frame(
  age = ages_grid,
  predicted_height = predict(height_dx_male_Down_syndrome_mod_df5$mod_fabc_rac, newdata = data.frame(age = ages_grid), level = 0),
  person_id = 999
)

# Add the predictions as a new layer
ggplot(height_dx_male_Down_syndrome, aes(x = age, y = height, group = person_id, color = as.factor(person_id))) +
  geom_line() +  # Line plot for each person
  geom_line(data = predictions, aes(x = age, y = predicted_height), color = "blue", size = 1) +  # Model predictions
  labs(title = "Down syndrome (Male)",
       x = "Age (years)",
       y = "Height (cm)") +  # Color legend for person_id
  theme_minimal() +
  theme(legend.position = "None")

# Check potential deviations from normality
qqnorm(ranef(mod5_ac)$a, main = "Q-Q plot of size random effect a ")
qqline(ranef(mod5_ac)$a)
qqnorm(ranef(mod5_ac)$c, main = "Q-Q plot of size random effect c")
qqline(ranef(mod5_ac)$c)
qqnorm(residuals(mod5_ac))
qqline(residuals(mod5_ac))

##############################################################
#          Generate Growth Chart at Major Centiles
##############################################################
Down_syndrome_male_centiles <- calculate_centiles_for_plot(model = height_dx_male_Down_syndrome_mod_df5$mod_fabc_rac,
                                                           n_ind = 5000)

# Prepare centile data for plotting
centile_data <- data.frame(
  age = ages_grid,  # Age grid from the predictions
  c05 = Down_syndrome_male_centiles[1, ],  # 5th centile
  c10 = Down_syndrome_male_centiles[2, ],  # 10th centile
  c25 = Down_syndrome_male_centiles[3, ],  # 25th centile
  c50 = Down_syndrome_male_centiles[4, ],  # 50th centile
  c75 = Down_syndrome_male_centiles[5, ],  # 75th centile
  c90 = Down_syndrome_male_centiles[6, ],  # 90th centile
  c95 = Down_syndrome_male_centiles[7, ]   # 95th centile
)
centile_data <- check_Monotonicity(centile_data, cols_to_check)

Down_syndrome_male_centile_data_long <- pivot_longer(centile_data, cols = starts_with("c"), 
                                                     names_to = "centile", values_to = "height")

# Filter to get the last point for each centile
Down_syndrome_male_centile_labels <- Down_syndrome_male_centile_data_long %>%
  group_by(centile) %>%
  filter(age == max(age)) 

Down_syndrome_male_centile_labels <- Down_syndrome_male_centile_labels %>%
  mutate(centile = gsub("^c", "", centile)) %>%
  mutate(centile = as.character(as.numeric(centile))) # Replace "c" at the start with an empty string

Down_syndrome_male_centile_data_long$sex = "Male"
Down_syndrome_male_centile_data_long$condition = "Down Syndrome"
Down_syndrome_male_centile_data_long$variantClass <- NA
Down_syndrome_male_centile_data_long$n <- length(unique(height_dx_male_Down_syndrome$person_id))

Down_syndrome_male_centile_labels$sex = "Male"
Down_syndrome_male_centile_labels$condition = "Down Syndrome"
Down_syndrome_male_centile_labels$variantClass = NA
Down_syndrome_male_centile_labels$n <- length(unique(height_dx_male_Down_syndrome$person_id))

# Generate the plot with blue dashed centile lines
p <- ggplot(height_no_dx_male, aes(x = age, y = height, group = person_id, color = as.factor(person_id))) +
  # Centile lines with a single color (blue) and dashed linetype
  # geom_line() +  
  geom_line(data = Down_syndrome_male_centile_data_long, aes(x = age, y = height, group = centile,
                                                             size = ifelse(centile %in% c("c05", "c50", "c95"), 0.7, 0.4)), 
            color = "blue", linetype = "solid",
  ) +
  # Add labels to the end of the centile lines
  geom_text(data = Down_syndrome_male_centile_labels, 
            aes(x = age, y = height, label = paste0(centile, "th")), 
            color = "blue", hjust = -0.2, size = 3.5, inherit.aes = FALSE) +
  scale_size_identity() +
  scale_x_continuous(
    breaks = seq(0, 20, by = 1),    # Smaller intervals for major grid lines on X
    minor_breaks = seq(0, 20, by = 0.2)  # Even smaller intervals for minor grid lines
  ) +
  scale_y_continuous(
    breaks = seq(55, 200, by = 5),   # Smaller intervals for major grid lines on Y
    minor_breaks = seq(55, 200, by = 1)  # Even smaller intervals for minor grid lines
  ) +
  labs(title = "Unaffected (Male)",
       x = "Age (years)",
       y = "Height (cm)") +
  theme_minimal() +
  theme(
    panel.grid.major = element_line(color = "gray80", size = 0.5),  # Major grid lines
    panel.grid.minor = element_line(color = "gray90", size = 0.2),  # Minor grid lines
    legend.position = "bottom",
    legend.key = element_blank(),  
    panel.background = element_rect(fill = "white") 
  )

p
