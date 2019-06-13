# Empirical Example with the data from HoeGekisNL

# 0.0 Prepare environment ----
rm(list=ls())
#library(devtools) # if lsttheory of sumplement C of Steyer et al. (2015) has not been installed before.
#install_github("amayer2010/lsttheory", force = TRUE)
#library(lavaan) 
#library(lsttheory)
library(MplusAutomation)
library(MASS)
library(xtable)
library(foreign)
library(psych)
# Easy code to source all the needed files
file.sources <- paste0("R/", list.files(paste0(getwd(), "/R")))
sapply(file.sources,source,.GlobalEnv)
rm(file.sources)
folder <- "HND/" #Folder to store results and all Mplus files

# Read and filter data ----

HND_data <- read.spss(paste0(folder, "190605 HGIN088 Castro.sav"), use.value.label=TRUE, to.data.frame=TRUE)

# Subset the data: Real_period = 1 & time<=90 & meting=1
HND_data_filtered <- subset(HND_data, HND_data$Real_period == 1 & HND_data$meting == 1 & HND_data$time <= 90) # It got reduced to 105362 observations

length(unique(HND_data_filtered$id)) # number of persons that participated at least in the crosectional study. 1396

# listwise deletion in order to check the total number of valid observations.
NA_ind <- apply(HND_data_filtered[, 35:40 ], 1, function(x) length(na.omit(x))) # Indicator for missing values in the PA items

HND_data_filtered_NA <- HND_data_filtered[NA_ind == 6,] # 61417 valid observations.

length(unique(HND_data_filtered_NA$id)) # number of persons with at least 1 observation.1302

# Subset data with person that have enough observations to get personalized feedback 65% ~ 58.5 observations

HND_data_filtered_65 <- HND_data_filtered[HND_data_filtered$Observations_valid90 >= 59, ]

length(unique(HND_data_filtered_65$id)) # number of persons with more than 59 observations. 644

# Select variables of interest

HND_data_filtered_65 <- HND_data_filtered_65[, c("id", "time", "gender", "sex", "age", "country", "mad_diary_5", "mad_diary_7", 
                                                 "mad_diary_9", "mad_diary_11", "mad_diary_13", "mad_diary_15")]

names(HND_data_filtered_65) <- c("id", "time", "gender", "sex", "age", "country",
                                 "PA1", "PA2", "PA3", "PA4", "PA5", "PA6")
### There is one person with 81 observations and two with 87 observations, should I complete these persons with NAs to have every
# person with 90 rows including NAs?

### There are 4 persons from belgium, and another 4 persons from other country, should I deleted them to have only dutch people?

# Fit the MSST ----

file.name <- paste("hnd", "msst",sep = "_")

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(HND_data_filtered_65[,c(1:2, 7:12)], paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(HND_data_filtered_65)[7:12],
                                       cluster = names(HND_data_filtered_65)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000)

ml_syntax <- write.mlmsst.to.Mplus(HND_data_filtered_65[, 7:12])

ml_syntax <- gsub("@0;", "@0.001;", ml_syntax) # replace 0 constraints to 0.001

saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                            ";", "\nOUTPUT: TECH8;")

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax, saveoutput_syntax)

# Run model in Mplus
cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
runModels(paste0(getwd(),"/", folder,file.name,".inp"))
cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0


# Fit the CUTS ----

file.name <- paste("hnd", "cuts",sep = "_")

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(HND_data_filtered_65[,c(1:2, 7:12)], paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(HND_data_filtered_65)[7:12],
                                       cluster = names(HND_data_filtered_65)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000)

ml_syntax <- write.mlcuts.to.Mplus(HND_data_filtered_65[, 7:12])
saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                            ";", "\nOUTPUT: TECH8;")

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax, saveoutput_syntax)

# Run modelin Mplus
cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
runModels(paste0(getwd(),"/", folder,file.name,".inp"))
cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0

# Fit the TSO ----

file.name <- paste("hnd", "tso",sep = "_")

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(HND_data_filtered_65[,c(1:2, 7:12)], paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(HND_data_filtered_65)[7:12],
                                       cluster = names(HND_data_filtered_65)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000)

ml_syntax <- write.mltso.to.Mplus(HND_data_filtered_65[, 7:12])
saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                            ";", "\nOUTPUT: TECH8;")

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax, saveoutput_syntax)

# Run modelin Mplus
cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
runModels(paste0(getwd(),"/", folder,file.name,".inp"))
cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0

# End ----