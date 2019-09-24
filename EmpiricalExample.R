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

# NA proportion is 16.72%
sum(c(is.na(HND_data_filtered_65[,7:12])))/length(c(is.na(HND_data_filtered_65[,7:12])))

### There are 4 persons from belgium, and another 4 persons from other country, should I deleted them to have only dutch people?

# Data long to data wide: to sumirize time-invariant variables (sex)

HND_data_filtered_65_w <- reshape(HND_data_filtered_65, v.names = names(HND_data_filtered_65)[7:12],
        timevar = "time", idvar="id", direction="wide")
summary(HND_data_filtered_65_w$gender)
summary(HND_data_filtered_65_w$age)

# plot person's time series

persons <- sample(1:644, 50)
pdf("HND/relaxed.pdf", height = 4)
ts.plot(t(HND_data_filtered_65_w[persons,seq(6,545, by = 6)]), col = gray(0.5), gpars = list(ylab = "Relaxed", las = 1))
lines(1:90, apply(t(HND_data_filtered_65_w[,seq(6,545, by = 6)]), 1, function(x) mean(x, na.rm = TRUE)), 
        col = "black", lty = 1, lwd = 4 )
dev.off()

pdf("HND/energetic.pdf", height = 4)
ts.plot(t(HND_data_filtered_65_w[persons,seq(7,545, by = 6)]), col = gray(0.5), gpars = list(ylab = "Energetic", las = 1))
lines(1:90, apply(t(HND_data_filtered_65_w[,seq(7,545, by = 6)]), 1, function(x) mean(x, na.rm = TRUE)), 
      col = "black", lty = 1, lwd = 4 )
dev.off()

pdf("HND/enthusiastic.pdf", height = 4)
ts.plot(t(HND_data_filtered_65_w[persons,seq(8,545, by = 6)]), col = gray(0.5), gpars = list(ylab = "Enthusiastic", las = 1))
lines(1:90, apply(t(HND_data_filtered_65_w[,seq(8,545, by = 6)]), 1, function(x) mean(x, na.rm = TRUE)), 
      col = "black", lty = 1, lwd = 4 )
dev.off()

pdf("HND/content.pdf", height = 4)
ts.plot(t(HND_data_filtered_65_w[persons,seq(9,545, by = 6)]), col = gray(0.5), gpars = list(ylab = "Content", las = 1))
lines(1:90, apply(t(HND_data_filtered_65_w[,seq(9,545, by = 6)]), 1, function(x) mean(x, na.rm = TRUE)), 
      col = "black", lty = 1, lwd = 4 )
dev.off()

pdf("HND/calm.pdf", height = 4)
ts.plot(t(HND_data_filtered_65_w[persons,seq(10,545, by = 6)]), col = gray(0.5), gpars = list(ylab = "Calm", las = 1))
lines(1:90, apply(t(HND_data_filtered_65_w[,seq(10,545, by = 6)]), 1, function(x) mean(x, na.rm = TRUE)), 
      col = "black", lty = 1, lwd = 4 )
dev.off()

pdf("HND/cheerful.pdf", height = 4)
ts.plot(t(HND_data_filtered_65_w[persons,seq(11,545, by = 6)]), col = gray(0.5), gpars = list(ylab = "Cheerful", las = 1))
lines(1:90, apply(t(HND_data_filtered_65_w[,seq(11,545, by = 6)]), 1, function(x) mean(x, na.rm = TRUE)), 
      col = "black", lty = 1, lwd = 4 )
dev.off()

# plot intraindividual means and intraindividual standard deviations histograms

MEANS <- cbind(apply(HND_data_filtered_65_w[,seq(9,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE)),
               apply(HND_data_filtered_65_w[,seq(10,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE)),
               apply(HND_data_filtered_65_w[,seq(6,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE)),
               apply(HND_data_filtered_65_w[,seq(7,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE)),
               apply(HND_data_filtered_65_w[,seq(8,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE)),
               apply(HND_data_filtered_65_w[,seq(11,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE)))

SDS <- cbind(apply(HND_data_filtered_65_w[,seq(9,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE)),
             apply(HND_data_filtered_65_w[,seq(10,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE)),
             apply(HND_data_filtered_65_w[,seq(6,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE)),
             apply(HND_data_filtered_65_w[,seq(7,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE)),
             apply(HND_data_filtered_65_w[,seq(8,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE)),
             apply(HND_data_filtered_65_w[,seq(11,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE)))

CON <- c(0.41, 0.39, 0.40, 0.35, 0.39, 0.39)
SPE <- c(0.27, 0.28, 0.40, 0.35, 0.45, 0.40)

# plot interindividual means and interindividual standard deviations histograms
pdf(file = "HND/Histograms.pdf", height = 5)
means <- apply(HND_data_filtered_65_w[,seq(9,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE))
hist(means, xlab = "Content Means", xlim = c(0, 100), ylim = c(0, 0.05),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(means),sd=sd(means)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(means) - 2 * sd(means), mean(means) + 2 * sd(means)), y = c(0.001, 0.001), lwd = 3, col = "red")
lines(x = c(mean(means) - 1 * sd(means), mean(means) + 1 * sd(means)), y = c(0.002, 0.002), lwd = 3, col = "blue")
legend("topright", "CON = 0.41")

means <- apply(HND_data_filtered_65_w[,seq(10,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE))
hist(means, xlab = "Calm Means", xlim = c(0, 100), ylim = c(0, 0.05),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(means),sd=sd(means)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(means) - 2 * sd(means), mean(means) + 2 * sd(means)), y = c(0.001, 0.001), lwd = 3, col = "red")
lines(x = c(mean(means) - 1 * sd(means), mean(means) + 1 * sd(means)), y = c(0.002, 0.002), lwd = 3, col = "blue")
legend("topright", "CON = 0.39")

means <- apply(HND_data_filtered_65_w[,seq(6,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE))
hist(means, xlab = "Relaxed Means", xlim = c(0, 100), ylim = c(0, 0.05),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(means),sd=sd(means)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(means) - 2 * sd(means), mean(means) + 2 * sd(means)), y = c(0.001, 0.001), lwd = 3, col = "red")
lines(x = c(mean(means) - 1 * sd(means), mean(means) + 1 * sd(means)), y = c(0.002, 0.002), lwd = 3, col = "blue")
legend("topright", "CON = 0.40")

means <- apply(HND_data_filtered_65_w[,seq(7,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE))
hist(means, xlab = "Energetic Means", xlim = c(0, 100), ylim = c(0, 0.05),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(means),sd=sd(means)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(means) - 2 * sd(means), mean(means) + 2 * sd(means)), y = c(0.001, 0.001), lwd = 3, col = "red")
lines(x = c(mean(means) - 1 * sd(means), mean(means) + 1 * sd(means)), y = c(0.002, 0.002), lwd = 3, col = "blue")
legend("topright", "CON = 0.35")

means <- apply(HND_data_filtered_65_w[,seq(8,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE))
hist(means, xlab = "Enthusiastic Means", xlim = c(0, 100), ylim = c(0, 0.05),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(means),sd=sd(means)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(means) - 2 * sd(means), mean(means) + 2 * sd(means)), y = c(0.001, 0.001), lwd = 3, col = "red")
lines(x = c(mean(means) - 1 * sd(means), mean(means) + 1 * sd(means)), y = c(0.002, 0.002), lwd = 3, col = "blue")
legend("topright", "CON = 0.39")

means <- apply(HND_data_filtered_65_w[,seq(11,545, by = 6)], 1, function(x) mean(x, na.rm = TRUE))
hist(means, xlab = "Cheerful Means", xlim = c(0, 100), ylim = c(0, 0.05),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(means),sd=sd(means)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(means) - 2 * sd(means), mean(means) + 2 * sd(means)), y = c(0.001, 0.001), lwd = 3, col = "red")
lines(x = c(mean(means) - 1 * sd(means), mean(means) + 1 * sd(means)), y = c(0.002, 0.002), lwd = 3, col = "blue")
legend("topright", "CON = 0.39")

sds <- apply(HND_data_filtered_65_w[,seq(9,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE))
hist(sds, xlab = "Content SDs", xlim = c(0, 40), ylim = c(0, 0.1),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(sds),sd=sd(sds)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(sds) - 2 * sd(sds), mean(sds) + 2 * sd(sds)), y = c(0.002, 0.002), lwd = 3, col = "red")
lines(x = c(mean(sds) - 1 * sd(sds), mean(sds) + 1 * sd(sds)), y = c(0.004, 0.004), lwd = 3, col = "blue")
legend("topright", "SPE = 0.27")

sds <- apply(HND_data_filtered_65_w[,seq(10,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE))
hist(sds, xlab = "Calm SDs", xlim = c(0, 40), ylim = c(0, 0.1),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(sds),sd=sd(sds)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(sds) - 2 * sd(sds), mean(sds) + 2 * sd(sds)), y = c(0.002, 0.002), lwd = 3, col = "red")
lines(x = c(mean(sds) - 1 * sd(sds), mean(sds) + 1 * sd(sds)), y = c(0.004, 0.004), lwd = 3, col = "blue")
legend("topright", "SPE = 0.28")

sds <- apply(HND_data_filtered_65_w[,seq(6,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE))
hist(sds, xlab = "Relaxed SDs", xlim = c(0, 40), ylim = c(0, 0.1),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(sds),sd=sd(sds)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(sds) - 2 * sd(sds), mean(sds) + 2 * sd(sds)), y = c(0.002, 0.002), lwd = 3, col = "red")
lines(x = c(mean(sds) - 1 * sd(sds), mean(sds) + 1 * sd(sds)), y = c(0.004, 0.004), lwd = 3, col = "blue")
legend("topright", "SPE = 0.40")

sds <- apply(HND_data_filtered_65_w[,seq(7,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE))
hist(sds, xlab = "Energetic SDs", xlim = c(0, 40), ylim = c(0, 0.1),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(sds),sd=sd(sds)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(sds) - 2 * sd(sds), mean(sds) + 2 * sd(sds)), y = c(0.002, 0.002), lwd = 3, col = "red")
lines(x = c(mean(sds) - 1 * sd(sds), mean(sds) + 1 * sd(sds)), y = c(0.004, 0.004), lwd = 3, col = "blue")
legend("topright", "SPE = 0.35")

sds <- apply(HND_data_filtered_65_w[,seq(8,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE))
hist(sds, xlab = "Enthusiastic SDs", xlim = c(0, 40), ylim = c(0, 0.1),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(sds),sd=sd(sds)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(sds) - 2 * sd(sds), mean(sds) + 2 * sd(sds)), y = c(0.002, 0.002), lwd = 3, col = "red")
lines(x = c(mean(sds) - 1 * sd(sds), mean(sds) + 1 * sd(sds)), y = c(0.004, 0.004), lwd = 3, col = "blue")
legend("topright", "SPE = 0.45")

sds <- apply(HND_data_filtered_65_w[,seq(11,545, by = 6)], 1, function(x) sd(x, na.rm = TRUE))
hist(sds, xlab = "Cheerful SDs", xlim = c(0, 40), ylim = c(0, 0.1),
     breaks = 24, main = NULL, freq = FALSE)
curve(dnorm(x,mean=mean(sds),sd=sd(sds)), add=TRUE,col="black", lwd = 3)
lines(x = c(mean(sds) - 2 * sd(sds), mean(sds) + 2 * sd(sds)), y = c(0.002, 0.002), lwd = 3, col = "red")
lines(x = c(mean(sds) - 1 * sd(sds), mean(sds) + 1 * sd(sds)), y = c(0.004, 0.004), lwd = 3, col = "blue")
legend("topright", "SPE = 0.40")

plot(CON, apply(MEANS, 2, function(x) (mean(x) + sd(x)) - (mean(x) - sd(x))), ylab = "Length Blue Interval", pch = 19, las = 1)
abline(lm(apply(MEANS, 2, function(x) (mean(x) + sd(x)) - (mean(x) - sd(x))) ~ CON))

plot(CON, apply(MEANS, 2, function(x) (mean(x) + 2 * sd(x)) - (mean(x) - 2 * sd(x))), ylab = "Length Red Interval", pch = 19, las = 1)
abline(lm(apply(MEANS, 2, function(x) (mean(x) + 2 * sd(x)) - (mean(x) - 2 * sd(x))) ~ CON))

plot(SPE, apply(SDS, 2, function(x) (mean(x) + sd(x)) - (mean(x) - sd(x))), ylab = "Length Blue Interval", pch = 19, las = 1)
abline(lm(apply(SDS, 2, function(x) (mean(x) + sd(x)) - (mean(x) - sd(x))) ~ SPE))

plot(SPE, apply(SDS, 2, function(x) (mean(x) + 2 * sd(x)) - (mean(x) - 2 * sd(x))), ylab = "Length Red Interval", pch = 19, las = 1)
abline(lm(apply(SDS, 2, function(x) (mean(x) + 2 * sd(x)) - (mean(x) - 2 * sd(x))) ~ SPE))

plot(SPE, apply(SDS, 2, function(x) mean(x)), ylab = "Mean Intraindividual SD", pch = 19, las = 1)
abline(lm(apply(SDS, 2, function(x) mean(x)) ~ SPE))

dev.off()

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

# Run models Final ----
runModels(paste0(getwd(),"/", folder,"hnd_cuts_3itemsa_center_noint_plots.inp"))
runModels(paste0(getwd(),"/", folder,"hnd_cuts_3itemsb_center_noint_plots.inp"))
runModels(paste0(getwd(),"/", folder,"hnd_msst_3itemsa_grandcenter_noint_plots.inp"))
runModels(paste0(getwd(),"/", folder,"hnd_msst_3itemsb_grandcenter_noint_plots.inp"))
runModels(paste0(getwd(),"/", folder,"hnd_tso_3itemsa_grandcenter_noint_plots.inp"))
runModels(paste0(getwd(),"/", folder,"hnd_tso_3itemsb_grandcenter_noint_plots.inp"))

# Compute variance coefficients items deactivation
samples <- read.table(paste0(getwd(), "/", folder, "samples_hnd_tso_3itemsa_grandcenter_noint_plots.dat"))

# Delete burn-in and indicator for chain and draw
samples <- samples[which(samples$V2 >= (max(samples$V2[samples$V1==2])/2 + 1) & samples$V2 <= max(samples$V2[samples$V1==2])),
                   -(1:2)
                   ]

# Create matrix to store variance coefficients draws
pdist.var.coeff <- matrix(NA, dim(samples)[1], 15)

for(i in 1:dim(samples)[1]){
  within.estimates <- list( loadings = c(1, t(samples[i, 1:2])), state.var = samples[i, 7],
                            error.var = c(t(samples[i, 3:5])), ar.effect = samples[i, 6])
  between.estimates <- list( trait.ind.var = c(t(samples[i, c(8, 10, 13)])))
  pdist.var.coeff[i, ] <- tso.var.coeff(I = 3, nT = 90, within.parameters = within.estimates,
                                        between.parameters = between.estimates)[,90]
}

varcoeffD <- data.frame(round(matrix(apply(pdist.var.coeff, 2, median),5, 3, byrow = TRUE),2))

names(varcoeffD) <- c("Relaxed", "Content", "Calm")
row.names(varcoeffD) <- c("\\hspace{1.5cm}Predictability by Trait", "\\hspace{1.5cm}Unpredictability by Trait", 
                          "\\hspace{0.75cm}Consistency", 
                          "\\hspace{0.75cm}Occasion Specificity", "Reliability")

varcoeffD <- varcoeffD[c(5,3,1,2,4),]

print(xtable(varcoeffD, type = "latex", caption = "Variance Coefficients Positive Affect Deactivation",
             label = "tab:varcoeffD", align = c("l", "c", "c", "c")),
      include.colnames=T, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", sanitize.text.function = function(x){x},
      file = "Mplus_Simulation/varcoeffD.txt")


# Compute variance coefficients items activation
samples <- read.table(paste0(getwd(), "/", folder, "samples_hnd_tso_3itemsb_grandcenter_noint_plots.dat"))

# Delete burn-in and indicator for chain and draw
samples <- samples[which(samples$V2 >= (max(samples$V2[samples$V1==2])/2 + 1) & samples$V2 <= max(samples$V2[samples$V1==2])),
                   -(1:2)
                   ]

# Create matrix to store variance coefficients draws
pdist.var.coeff <- matrix(NA, dim(samples)[1], 15)

for(i in 1:dim(samples)[1]){
  within.estimates <- list( loadings = c(1, t(samples[i, 1:2])), state.var = samples[i, 7],
                            error.var = c(t(samples[i, 3:5])), ar.effect = samples[i, 6])
  between.estimates <- list( trait.ind.var = c(t(samples[i, c(8, 10, 13)])))
  pdist.var.coeff[i, ] <- tso.var.coeff(I = 3, nT = 90, within.parameters = within.estimates,
                                        between.parameters = between.estimates)[,90]
}

varcoeffA <- data.frame(round(matrix(apply(pdist.var.coeff, 2, median),5, 3, byrow = TRUE),2))

names(varcoeffA) <- c("Energetic", "Enthusiastic", "Cheerful")
row.names(varcoeffA) <- c("\\hspace{1.5cm}Predictability by Trait", "\\hspace{1.5cm}Unpredictability by Trait", 
                          "\\hspace{0.75cm}Consistency", 
                          "\\hspace{0.75cm}Occasion Specificity", "Reliability")

varcoeffA <- varcoeffA[c(5,3,1,2,4),]

print(xtable(varcoeffA, type = "latex", caption = "Variance Coefficients Positive Affect Activation",
             label = "tab:varcoeffA", align = c("l", "c", "c", "c")),
      include.colnames=T, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", sanitize.text.function = function(x){x},
      file = "Mplus_Simulation/varcoeffA.txt")

# Cross lagged TSO model

WSigma <- matrix(c(108.020, 97.266,97.266,143.410),2, byrow = TRUE)

WR <- solve(diag(sqrt(diag(WSigma))))%*%WSigma%*%solve(diag(sqrt(diag(WSigma))))

BSigma <- matrix(c(125.928, 120.722, 121.035, 98.008, 99.874, 107.330,
                   120.722, 158.245, 119.039, 112.649, 124.776, 132.517,
                   121.035, 119.039, 140.108, 95.981, 95.997, 105.164,
                   98.008, 112.649, 95.981, 146.502, 132.299, 128.355,
                   99.874, 124.776, 95.997, 132.299, 150.499, 139.758,
                   107.330, 132.517, 105.164, 128.355, 139.758, 153.634), 6, byrow = TRUE)
BR <- solve(diag(sqrt(diag(BSigma))))%*%BSigma%*%solve(diag(sqrt(diag(BSigma))))


# End ----