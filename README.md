# An EHR-Based Framework for Modeling Growth Curves and Constructing Growth Centile Charts for Genetic Disorders

---

## Scripts

1. **`Utils.R`**  
   Loads packages and helpfer functions

2. **`generate_growth_curves_with_diagnostics.R`**  
   This script runs the full SITAR-based growth-curve workflow for one condition and sex, fits and compares multiple candidate SITAR height-curve models, selects the optimal model via BIC, compares affected versus unaffected cohorts, and generates diagnostic summaries and plots based on the optimal model.

3. **`generate_growth_centiles_with_diagnostics.R`**  
   This script fits and compares multiple candidate GAMLSS height-centile models, selects the optimal model via BIC, and generates centile tables and diagnostic plots based on the optimal model.

4. **`run_growth_job.R`**  
   This script takes one row from a condition list, loads and filters the relevant patient height data for that condition/sex (and CF variant subgroup when applicable), runs a SITAR growth-curve analysis, then runs a GAMLSS centile analysis on the cleaned results, and saves both sets of outputs.

5. **`run_growth_job.sh`**
   This SLURM batch script runs the growth pipeline in parallel by launching an array of jobs (one per condition row), each executing `run_growth_job.R` with a different job index and saving logs for each run.  

---

---

## Data

1. **`List_of_Conditions.xlsx`**  
   List of genetic conditions analyzed in this study, used as input to `run_growth_job.R`

---


