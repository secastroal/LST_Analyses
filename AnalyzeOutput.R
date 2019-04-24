# This script analyze the results of the simulation analysis
library(plyr)

perf <- read.table(paste(getwd(), "Mplus_Simulation" , "Performance", "performance_1.dat", sep = "/"), header = TRUE)
times <- read.table(paste(getwd(), "Mplus_Simulation" , "Times", "times_1.dat", sep = "/"), header = TRUE)
parameters <- read.table(paste(getwd(), "Mplus_Simulation", "Parameters", "parameters_1.dat", sep = "/"), header = TRUE)
se_psd <- read.table(paste(getwd(), "Mplus_Simulation", "SE_PSD", "se_psd_1.dat", sep = "/"), header = TRUE)
var.coeff <- read.table(paste(getwd(), "Mplus_Simulation", "Var_Coeff", "var_coeff_1.dat", sep = "/"), header = TRUE)
fit.measures <- read.table(paste(getwd(), "Mplus_Simulation", "Fit_Measures", "fit_measures_1.dat", sep = "/"), header = TRUE)
psd.var.coeff <- read.table(paste(getwd(), "Mplus_Simulation", "PSD_Var_Coeff", "psd_var_coeff_1.dat", sep = "/"), header = TRUE)

# Performances ----
files <- paste(getwd(), "Mplus_Simulation" , "Performance", paste0("performance_",  1:15, ".dat"), sep = "/")

performances <- lapply(files, function(x) read.table(file = x, header = TRUE, colClasses = "character"))

performances <- ldply(performances, data.frame)

rm(files)

apply(performances, 2, function(x) length(which(x == "Ok")))

# Running times ----
files <- paste(getwd(), "Mplus_Simulation" , "Times", paste0("times_",  1:15, ".dat"), sep = "/")

running.times <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

running.times <- ldply(running.times, data.frame)

# Parameter estimates ----
files <- paste(getwd(), "Mplus_Simulation" , "Parameters", paste0("parameters_",  1:15, ".dat"), sep = "/")

parameters <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

parameters <- ldply(parameters, data.frame)

# Standard errors and/or posterior sd ----
files <- paste(getwd(), "Mplus_Simulation" , "SE_PSD", paste0("se_psd_",  1:15, ".dat"), sep = "/")

se_psd <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

se_psd <- ldply(se_psd, data.frame)

# Variance coefficients ----
files <- paste(getwd(), "Mplus_Simulation" , "Var_Coeff", paste0("var_coeff_",  1:15, ".dat"), sep = "/")

var_coeff <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

var_coeff <- ldply(var_coeff, data.frame)

# Posterior sd of variance coefficients ----
files <- paste(getwd(), "Mplus_Simulation" , "PSD_Var_Coeff", paste0("psd_var_coeff_",  1:15, ".dat"), sep = "/")

psd_var_coeff <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

psd_var_coeff <- ldply(psd_var_coeff, data.frame)

# Fit measures ----
files <- paste(getwd(), "Mplus_Simulation" , "Fit_Measures", paste0("fit_measures_",  1:15, ".dat"), sep = "/")

fit_measures <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

fit_measures <- ldply(fit_measures, data.frame)

# End ----