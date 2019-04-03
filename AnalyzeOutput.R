
perf <- read.table(paste(getwd(), "Mplus_Simulation" , "Performance", "performance_1.dat", sep = "/"), header = TRUE)
parameters <- read.table(paste(getwd(), "Mplus_Simulation", "Parameters", "parameters_1.dat", sep = "/"), header = TRUE)
se_psd <- read.table(paste(getwd(), "Mplus_Simulation", "SE_PSD", "se_psd_1.dat", sep = "/"), header = TRUE)
var.coeff <- read.table(paste(getwd(), "Mplus_Simulation", "Var_Coeff", "var_coeff_1.dat", sep = "/"), header = TRUE)
fit.measures <- read.table(paste(getwd(), "Mplus_Simulation", "Fit_Measures", "fit_measures_1.dat", sep = "/"), header = TRUE)
psd.var.coeff <- read.table(paste(getwd(), "Mplus_Simulation", "PSD_Var_Coeff", "psd_var_coeff_1.dat", sep = "/"), header = TRUE)

