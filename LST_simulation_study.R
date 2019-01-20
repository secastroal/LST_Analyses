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

# 1.0 CUTS ----
# Conditions

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

between.parameters <- list(intercepts = intercepts, loadings = trait_loadings, CT.var = var_CT, UT.var = var_UT)

rm(intercepts, trait_loadings, var_CT, var_UT)

#data simulation
cuts.data <- sim.data.cuts(N, nT, I, within.parameters = within.parameters, 
                           between.parameters = between.parameters, seed = seed)

# model estimation
file.name <- "sim_cuts_test1"

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(cuts.data$data.long, paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(cuts.data$data.long)[-(1:2)],
                                       cluster = names(cuts.data$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlcuts.to.Mplus(cuts.data$data.long[, -(1:2)])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))
readModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".out"), what = "parameters")$parameters #read Mplus output

rm(N, nT, I, seed, file.name, within.parameters, between.parameters, cuts.data)

# 2.0 MSST ----

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

between.parameters <- list(intercepts = intercepts, loadings = trait_loadings, trait.mean = mean_trait,
                           trait.var = var_trait)

rm(intercepts, trait_loadings, var_trait, mean_trait)

#data simulation
msst.data <- sim.data.msst(N, nT, I, within.parameters = within.parameters, 
                           between.parameters = between.parameters, seed = seed)

# model estimation
file.name <- "sim_msst_test1"

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(msst.data$data.long, paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(msst.data$data.long)[-(1:2)],
                                       cluster = names(msst.data$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlmsst.to.Mplus(msst.data$data.long[, -(1:2)])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))
readModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".out"), what = "parameters")$parameters #read Mplus output

rm(N, nT, I, seed, file.name, within.parameters, between.parameters, msst.data)

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

