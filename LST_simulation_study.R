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

# 1.3.1 ML-cuts data.long ----

file.name <- paste0(model, "_long_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(cuts.data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(cuts.data$data.long)[-(1:2)],
                                       cluster = names(cuts.data$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlcuts.to.Mplus(cuts.data$data.long[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

long.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[, 3] <- long.fit$unstandardized[,3]

rm(file.name, long.fit)

# 1.3.2 cuts data.wide ----

file.name <- paste0(model, "_wide_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(cuts.data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(cuts.data$data.wide)[-1],
                                       analysis_type = "GENERAL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)


mplus_syntax <- write.cuts.to.Mplus(cuts.data$data.wide[-1],  nstate = nT,
                                    method.trait = "om",
                                    scale.invariance = list(int = TRUE, lambda = TRUE),
                                    state.trait.invariance = FALSE,
                                    fixed.method.loadings = TRUE,
                                    homocedasticity.assumption = list(error = TRUE, cs.red = TRUE, ut.red = FALSE))


write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, mplus_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

wide.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[, 4] <- wide.fit$unstandardized[c(1:I,
                                          (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1,
                                          ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I),
                                          (I*nT + 1):(I*nT + I),
                                          ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I),
                                          (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1,
                                          ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)),
                                        3]

rm(file.name, wide.fit)

# 1.3.3 ML-cuts data.trunc ----

file.name <- paste0(model, "_trunc_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(cuts.data$data.trunc, paste0(folder, file.name, ".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(cuts.data$data.trunc)[-(1:2)],
                                       cluster = names(cuts.data$data.trunc)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlcuts.to.Mplus(cuts.data$data.trunc[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

trunc.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[, 5] <- trunc.fit$unstandardized[,3]

rm(file.name, trunc.fit)

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

rm(file.name, tv.fit)

# 1.4 Clean environment ----

rm(N, nT, I, seed, within.parameters, between.parameters, cuts.data, cuts.data.tv, est.par, model)


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

# 2.3.1 ML-msst data.long ----

file.name <- paste0(model, "_long_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(msst.data$data.long, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(msst.data$data.long)[-(1:2)],
                                       cluster = names(msst.data$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlmsst.to.Mplus(msst.data$data.long[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/", folder,file.name,".inp"))

long.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[,3] <- long.fit$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                         (4*I + 2), (6*I +3)),3]

rm(file.name, long.fit)

# 2.3.2 msst data.wide ----

file.name <- paste0(model, "_wide_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(msst.data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(msst.data$data.wide)[-1],
                                       analysis_type = "GENERAL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

mplus_syntax <- write.msst.to.Mplus(msst.data$data.wide[-1], neta = nT, ntheta = 1, 
                                    equiv.assumption = list(tau = "cong", theta = "cong"),
                                    scale.invariance = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                    homocedasticity.assumption = list(error = TRUE, state.red = TRUE),
                                    second.order.trait = FALSE)

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, mplus_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/", folder,file.name,".inp"))

wide.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[,4] <- wide.fit$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                         (4*I + 2), (6*I +3)),3]

rm(file.name, wide.fit)

# 2.3.3 ML-msst data.trunc ----

file.name <- paste0(model, "_trunc_N", N, "_I", I, "_nT", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(msst.data$data.trunc, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(msst.data$data.trunc)[-(1:2)],
                                       cluster = names(msst.data$data.trunc)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlmsst.to.Mplus(msst.data$data.trunc[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/", folder,file.name,".inp"))

trunc.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out"), what = "parameters")$parameters #read Mplus output

est.par[,5] <- trunc.fit$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                         (4*I + 2), (6*I +3)),3]

rm(file.name, trunc.fit)

# 2.3.1 ML-msst time-variant data ----

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

rm(file.name, tv.fit)

# 1.4 Clean environment ----

rm(N, nT, I, seed, within.parameters, between.parameters, msst.data, msst.data.tv, est.par, model)

# 3.0 TSO ----

# Conditions

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

within.parameters <- list(loadings = state_loadings, state.var = var_state, error.var = var_error,
                          ar.effect = ar_effect)

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

#data simulation
tso.data <- sim.data.tso(N, nT, I, within.parameters = within.parameters, 
                         between.parameters = between.parameters, seed = seed)

# model estimation
file.name <- "sim_tso_test1"

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(tso.data$data.long, paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(tso.data$data.long)[-(1:2)],
                                       cluster = names(tso.data$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000,
                                       processors = 2)

ml_syntax <- write.mltso.to.Mplus(tso.data$data.long[, -c(1, 2)])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))
readModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".out"), what = "parameters")$parameters #read Mplus output

rm(N, nT, I, seed, file.name, within.parameters, between.parameters, tso.data)

