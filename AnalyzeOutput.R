
perf <- read.table(paste(getwd(), "Mplus_Simulation" , "Performance", "performance_1.dat", sep = "/"), header = TRUE)
times <- read.table(paste(getwd(), "Mplus_Simulation" , "Times", "times_1.dat", sep = "/"), header = TRUE)
parameters <- read.table(paste(getwd(), "Mplus_Simulation", "Parameters", "parameters_1.dat", sep = "/"), header = TRUE)
se_psd <- read.table(paste(getwd(), "Mplus_Simulation", "SE_PSD", "se_psd_1.dat", sep = "/"), header = TRUE)
var.coeff <- read.table(paste(getwd(), "Mplus_Simulation", "Var_Coeff", "var_coeff_1.dat", sep = "/"), header = TRUE)
fit.measures <- read.table(paste(getwd(), "Mplus_Simulation", "Fit_Measures", "fit_measures_1.dat", sep = "/"), header = TRUE)
psd.var.coeff <- read.table(paste(getwd(), "Mplus_Simulation", "PSD_Var_Coeff", "psd_var_coeff_1.dat", sep = "/"), header = TRUE)


files <- paste(getwd(), "Mplus_Simulation" , "Performance", paste0("performance_",  1:6, ".dat"), sep = "/")

performances <- lapply(files, function(x) read.table(file = x, header = TRUE, colClasses = "character"))

files <- paste(getwd(), "Mplus_Simulation" , "Times", paste0("times_",  1:6, ".dat"), sep = "/")

r.times <- lapply(files, function(x) read.table(file = x, header = TRUE))


r.times <- as.data.frame(simplify2array(r.times, higher = FALSE))

r.performances <- as.data.frame(simplify2array(performances, higher = FALSE))


sim.test.results <- cbind(r.performances, r.times)

sim.test.results <- sim.test.results[, c(matrix(1:12, 2, 6, byrow = TRUE))]

sim.test.results <- as.matrix(sim.test.results)

sim.test.results[which(sim.test.results == "Errors and Warnings")] <- "Errors/Warnings"

sim.test.results <- matrix(unlist(sim.test.results), 11, 12)


write.csv(sim.test.results, "Mplus_Simulation/totaltimes.csv")

rm(files, performances, r.times, r.performances, sim.test.results)
rm(times)
