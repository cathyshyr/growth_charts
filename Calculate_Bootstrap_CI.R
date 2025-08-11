files <- list.files(pattern = "Height_Male_Down_Syndrome_Boot_\\d+.rds")
results <- lapply(files, readRDS)
results <- Filter(function(x) is.list(x) && x$success, results)

sizediff <- sapply(results, `[[`, "sizediff")
timediff <- sapply(results, `[[`, "timediff")
veldiff <- sapply(results, `[[`, "veldiff")

Down_Syndrome_male_sizediff <- data.frame(sizediff = sizediff, 
                                            dx = "Down Syndrome",
                                            dID = 190685)
Down_Syndrome_male_timediff <- data.frame(timediff = timediff, 
                                            dx = "Down Syndrome",
                                            dID = 190685)
Down_Syndrome_male_veldiff <- data.frame(veldiff = veldiff, 
                                           dx = "Down Syndrome",
                                           dID = 190685)
save(Down_Syndrome_male_sizediff, file = "Down_Syndrome_male_sizediff.RData")
save(Down_Syndrome_male_timediff, file = "Down_Syndrome_male_timediff.RData")
save(Down_Syndrome_male_veldiff, file = "Down_Syndrome_male_veldiff.RData")
quantile(Down_Syndrome_male_sizediff$sizediff, c(0.025, 0.975), na.rm = TRUE)
quantile(Down_Syndrome_male_timediff$timediff, c(0.025, 0.975), na.rm = TRUE)
quantile(Down_Syndrome_male_veldiff$veldiff, c(0.025, 0.975), na.rm = TRUE)