# Contents
# 0.0 Prepare environment
# 1.0 CUTS
# 2.0 MSST
# 3.0 TSO


# 0.0 Prepare environment ----
rm(list=ls())
#library(devtools) # if lsttheory of sumplement C of Steyer et al. (2015) has not been installed before.
#install_github("amayer2010/lsttheory", force = TRUE)
library(lavaan)
library(lsttheory)
library(MplusAutomation)
library(MASS)
source("R/write.Mplus.options.R")
source("R/write.cuts.to.Mplus.R")
source("R/write.mlcuts.to.Mplus.R")
source("R/write.msst.to.Mplus.R")
source("R/write.mlmsst.to.Mplus.R")
source("R/write.tso.to.Mplus.R")
source("R/write.mltso.to.Mplus.R")
source("R/sim.data.cuts.R")
source("R/sim.data.msst.R")
source("R/sim.data.tso.R")
folder <- "Mplus_files_Results/" #Folder to store results

# 1.0 CUTS ----

# 1.1 Set Conditions and true parameters ----

model <- "cuts"
N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
seed <- 123

set.seed(seed)

# Within Parameters

state_loadings <- sample(seq(0.5, 1.2, by = 0.1), size = I, replace = TRUE) # loading parameters for the latent common state
state_loadings[1] <- 1 # fixing first loading at 1
var_CS <- 2 # Variance latent common state
var_US <- sample(seq(1, 3, by = 0.5), size = I, replace = TRUE) # Variance of latent unique states

within.parameters <- list(loadings = state_loadings, CS.var = var_CS, US.var = var_US)

rm(state_loadings, var_CS, var_US)

# Between Paramaters

intercepts <- sample(seq(2, 4, by = 0.5), size = I, replace = TRUE) # intercepts
trait_loadings <- sample(seq(0.5, 1.2, by = 0.1), size = I, replace = TRUE) # loading parametes for the latent common trait
trait_loadings[1] <- 1 # fixing first loading at 1
var_CT <- 1.5 # variance latent common trait
var_UT <- sample(seq(0.5, 2, by = 0.5), size = I, replace = TRUE) # variance latent unique traits

between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, CT.var = var_CT, UT.var = var_UT)

rm(intercepts, trait_loadings, var_CT, var_UT)

# 1.2 Simulate data ----

# Time invariant data
cuts.data <- sim.data.cuts(N, nT, I, within.parameters = within.parameters, 
                           between.parameters = between.parameters, seed = seed)

# Truncate time invariant data long

cuts.data$data.trunc <- trunc(cuts.data$data.long)

# Time variant data

cuts.data.tv <- sim.data.cuts.tv(N, nT, I, within.parameters = within.parameters, time.invariant = FALSE, 
                                 between.parameters = between.parameters, seed = seed)

# 1.3 Model estimation ---- 

# Matrix to store estimated parameters

est.par <- data.frame(matrix(NA, 5*I + 2, 6))
est.par[, 2] <- c(unlist(within.parameters), unlist(between.parameters))
est.par[, 1] <- names(c(unlist(within.parameters), unlist(between.parameters)))
names(est.par) <- c("par", "true", "long", "wide", "trunc", "tv" )

# Matrix to store rmse
rmse <- est.par[, -2] 

# 1.3.4 ML-cuts time-variant data ----

file.name <- paste0(model, "_tv_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(cuts.data.tv$data.long, paste0(folder, file.name, ".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(cuts.data.tv$data.long)[-(1:2)],
                                       cluster = names(cuts.data.tv$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlcuts.to.Mplus(cuts.data.tv$data.long[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

tv.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[, 6] <- tv.fit$unstandardized[,3]

# In this case, the rmse is root mean of the se between the estimated parameter vs the set of 
# time-variant true parameters. 

tv.rmse <- rbind((est.par[1:I, 6] - t(cuts.data.tv$within.parameters$loadings)) ^ 2,
                 (est.par[I + 1, 6] - t(cuts.data.tv$within.parameters$CS.var)) ^ 2,
                 (est.par[(I + 2):(2 * I + 1), 6] - t(cuts.data.tv$within.parameters$US.var)) ^ 2,
                 (est.par[(2 * I + 2):(3 * I + 1), 6] - t(cuts.data.tv$between.parameters$loadings)) ^ 2,
                 (est.par[(3 * I + 2):(4 * I + 1), 6] - t(cuts.data.tv$between.parameters$intercepts)) ^ 2)

rmse[1:(4 * I + 1), 5] <- sqrt(apply(tv.rmse, 1, mean))
rmse[(4*I + 2):(5*I + 2), 5] <- sqrt((est.par[(4*I + 2):(5*I + 2), 6] - est.par[(4*I + 2):(5*I + 2), 2]) ^ 2)

rm(file.name, tv.fit, tv.rmse)

# 1.4 Save output ----

write.table(est.par, paste0(folder, model, "_par_N", N, "_I", I, "_nT", nT, ".dat"))
write.table(rmse, paste0(folder, model, "_rmse_N", N, "_I", I, "_nT", nT, ".dat"))

# 1.5 Clean environment ----

rm(N, nT, I, seed, within.parameters, between.parameters, cuts.data, cuts.data.tv, est.par, model, rmse)


# 2.0 MSST ----

# 2.1 Set Conditions and true parameters ----

model <- "msst" 
N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
seed <- 123

set.seed(seed)

# Within Parameters

state_loadings <- sample(seq(0.5, 1.2, by = 0.1), size = I, replace = TRUE) # loading parameters for the latent state
state_loadings[1] <- 1 # fixing first loading at 1
var_state <- 2 # Variance latent state residual
var_m_error <- sample(seq(1, 3, by = 0.5), size = I, replace = TRUE) # Variance of measurement errors

within.parameters <- list(loadings = state_loadings, state.var = var_state, error.var = var_m_error)

rm(state_loadings, var_state, var_m_error)

# Between Paramaters

intercepts <- sample(seq(0, 1, by = 0.2), size = I, replace = TRUE) # intercepts
intercepts[1] <- 0 # fixing first intercept at 0
trait_loadings <- sample(seq(0.5, 1.2, by = 0.1), size = I, replace = TRUE) # loading parametes for the latent trait
trait_loadings[1] <- 1 # fixing first loading at 1

var_trait <- 2 # variance latent trait variable
mean_trait <- 4 # mean latent trait variable

between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, trait.mean = mean_trait,
                           trait.var = var_trait)

rm(intercepts, trait_loadings, var_trait, mean_trait)

# 2.2 Simulate data ----

#Time invariant data
msst.data <- sim.data.msst(N, nT, I, within.parameters = within.parameters, 
                           between.parameters = between.parameters, seed = seed)

# Truncate time invariant data long

msst.data$data.trunc <- trunc(msst.data$data.long)

# Time variant data

msst.data.tv <- sim.data.msst.tv(N, nT, I,  within.parameters = within.parameters, 
                                 between.parameters = between.parameters, time.invariant = FALSE, seed = seed)

# 2.3 Model estimation ----

# Matrix to store estimated parameters

est.par <- data.frame(matrix(NA, 4*I + 3, 6))
est.par[, 2] <- c(unlist(within.parameters), unlist(between.parameters))
est.par[, 1] <- names(c(unlist(within.parameters), unlist(between.parameters)))
names(est.par) <- c("par", "true", "long", "wide", "trunc", "tv" )

# Matrix to store rmse
rmse <- est.par[, -2]

# 2.3.4 ML-msst time-variant data ----

file.name <- paste0(model, "_tv_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(msst.data.tv$data.long, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(msst.data.tv$data.long)[-(1:2)],
                                       cluster = names(msst.data.tv$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlmsst.to.Mplus(msst.data.tv$data.long[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/", folder,file.name,".inp"))

tv.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[,6] <- tv.fit$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                       (4*I + 2), (6*I +3)),3]

# In this case, the rmse is root mean of the se between the estimated parameter vs the set of 
# time-variant true parameters. 

tv.rmse <- rbind((est.par[1:I, 6] - t(msst.data.tv$within.parameters$loadings)) ^ 2,
                 (est.par[I + 1, 6] - t(msst.data.tv$within.parameters$state.var)) ^ 2,
                 (est.par[(I + 2):(2 * I + 1), 6] - t(msst.data.tv$within.parameters$error.var)) ^ 2,
                 (est.par[(2 * I + 2):(3 * I + 1), 6] - t(msst.data.tv$between.parameters$loadings)) ^ 2,
                 (est.par[(3 * I + 2):(4 * I + 1), 6] - t(msst.data.tv$between.parameters$intercepts)) ^ 2)

rmse[1:(4 * I + 1), 5] <- sqrt(apply(tv.rmse, 1, mean))
rmse[(4*I + 2):(4*I + 3), 5] <- sqrt((est.par[(4*I + 2):(4*I + 3), 6] - est.par[(4*I + 2):(4*I + 3), 2]) ^ 2)

rm(file.name, tv.fit, tv.rmse)

# 2.4 Save output ----

write.table(est.par, paste0(folder, model, "_par_N", N, "_I", I, "_nT", nT, ".dat"))
write.table(rmse, paste0(folder, model, "_rmse_N", N, "_I", I, "_nT", nT, ".dat"))

# 2.5 Clean environment ----

rm(N, nT, I, seed, within.parameters, between.parameters, msst.data, msst.data.tv, est.par, model, rmse)

# 3.0 TSO ----

# 3.1 Set Conditions and true parameters ----

model <- "tso" 
N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
seed <- 123

set.seed(seed)

# Within Parameters

state_loadings <- sample(seq(0.5, 1.2, by = 0.1), size = I, replace = TRUE) # loading parameters for the latent state
state_loadings[1] <- 1 # fixing first loading at 1
var_state <- 2 # Variance latent state residual
var_error <- sample(seq(1, 3, by = 0.5), size = I, replace = TRUE) # Variance of latent measurement errors

ar_effect <- 0.5 # autoregressive effect on the latent state residuals

within.parameters <- list(loadings = state_loadings, ar.effect = ar_effect, error.var = var_error,
                          state.var = var_state)

rm(state_loadings, var_state, var_error, ar_effect)

# Between Paramaters 

intercepts <- sample(seq(2, 4, by = 0.5), size = I, replace = TRUE) # intercepts

var_ind_traits <- sample(seq(1.5, 2.5, by = 0.5), size = I, replace = TRUE) # variance latent indicator trait variables

R <- matrix(sample((6:9)/10, size = I * I, replace = TRUE), I) #correlation matrix trait indicators
R[lower.tri(R)] = t(R)[lower.tri(R)]
diag(R) <- 1

between.parameters <- list(intercepts = intercepts, trait.ind.var = var_ind_traits, 
                           cor.matrix = R)

rm(intercepts, var_ind_traits, R)

# 3.2 Simulate data ----

# Time invariant data
tso.data <- sim.data.tso(N, nT, I, within.parameters = within.parameters, 
                         between.parameters = between.parameters, seed = seed)

# Truncate time invariant data long

tso.data$data.trunc <- trunc(tso.data$data.long)

# Time variant data

tso.data.tv <- sim.data.tso.tv(N, nT, I, within.parameters = within.parameters, time.invariant = FALSE,
                               between.parameters = between.parameters, seed = seed)

# 3.3 Model estimation ----

# Matrix to store estimated parameters

est.par <- data.frame(matrix(NA, 4*I + 2 + (I * (I-1) /2), 6))
# Get true parameters + lower triangle of the true variance-covariance matrix 
est.par[, 2] <- round(c(unlist(within.parameters), unlist(between.parameters)[1:(2*I)],
                        tso.data$between.parameters$Sigma[t(lower.tri(tso.data$between.parameters$Sigma))]),2)
# Get array indices of the variance-covariance matrix
cov.ind <- which(t(lower.tri(tso.data$between.parameters$Sigma))==TRUE, arr.ind = TRUE)
# Define parameter names + covariance names
est.par[, 1] <- c(names(c(unlist(within.parameters), unlist(between.parameters)[1:(2*I)])),
                  paste0("cov", cov.ind[,1], cov.ind[,2]))
names(est.par) <- c("par", "true", "long", "wide", "trunc", "tv" )

# Matrix to store rmse
rmse <- est.par[, -2] 

rm(cov.ind)

# 3.3.4 ML-tso time-variant data ----

file.name <- paste0(model, "_tv_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(tso.data.tv$data.long, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(tso.data.tv$data.long)[-(1:2)],
                                       cluster = names(tso.data.tv$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000,
                                       processors = 2)

ml_syntax <- write.mltso.to.Mplus(tso.data.tv$data.long[, -c(1, 2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

tv.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[, 6] <- tv.fit$unstandardized[c(1:(2*I + 2),
                                        (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2),
                                        (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)), 3]

# In this case, the rmse is root mean of the se between the estimated parameter vs the set of 
# time-variant true parameters. 

tv.rmse <- rbind((est.par[1:I, 6] - t(tso.data.tv$within.parameters$loadings)) ^ 2,
                 (est.par[I + 1, 6] - t(c(tso.data.tv$within.parameters$ar.effect, NA))) ^ 2, # An NA is added to complete a vector of length nT
                 (est.par[(I + 2):(2 * I + 1), 6] - t(tso.data.tv$within.parameters$error.var)) ^ 2,
                 (est.par[(2 * I + 2), 6] - t(tso.data.tv$within.parameters$state.var)) ^ 2,
                 (est.par[(2 * I + 3):(3 * I + 2), 6] - t(tso.data.tv$between.parameters$intercepts)) ^ 2)

rmse[1:(3 * I + 2), 5] <- sqrt(apply(tv.rmse, 1, mean, na.rm = TRUE))
rmse[(3*I + 3):(4*I + 2 + (I * (I-1) /2)), 5] <- sqrt((est.par[(3*I + 3):(4*I + 2 + (I * (I-1) /2)), 6] - 
                                                         est.par[(3*I + 3):(4*I + 2 + (I * (I-1) /2)), 2]) ^ 2)

rm(file.name, tv.fit, tv.rmse)

# 3.4 Save output ----

write.table(est.par, paste0(folder, model, "_par_N", N, "_I", I, "_nT", nT, ".dat"))
write.table(rmse, paste0(folder, model, "_rmse_N", N, "_I", I, "_nT", nT, ".dat"))

# 3.5 Clean environment ----

rm(N, nT, I, seed, within.parameters, between.parameters, tso.data, tso.data.tv, est.par, model, rmse)