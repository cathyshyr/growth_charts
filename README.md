# Growth Charts for Genetic Disorders in Pediatric Patients Derived from Electronic Health Records

This repository provides code to reproduce the growth curve modeling and bootstrap analysis described in our paper. 

---

## Scripts

1. **`Generate_Condition_Specific_Growth_Chart.R`**  
   Fits, selects, and visualizes sex- and condition-specific height-for-age growth curves (age 2–20) using the SITAR model.  
   - Example provided: **boys with Down syndrome**; adjust `dID` and the sex filter to generate other sex/condition charts.  
   - Produces diagnostic plots (QQ plots, spaghetti overlay) and model-derived centile curves (5th–95th).

2. **`Bootstrap_Condition_Specific_Growth_Curve.R`**  
   Runs bootstrap resampling of the selected SITAR model to estimate:
   - maximum height  
   - age at peak height velocity  
   - peak height velocity  
   Designed for parallel execution on a cluster (array jobs).

3. **`Bootstrap_Condition_Specific_Growth_Curve.sh`**  
   SLURM batch script to submit the bootstrap array (e.g., 5,000 replicates). Please edit SLURM directives as needed for your cluster.

4. **`Calculate_Bootstrap_CI.R`**  
   Aggregates bootstrap outputs, filters successful runs, and computes 95% CIs for height (maximum height), timing (age at peak height velocity), and velocity (peak height velocity) differences. 

---
