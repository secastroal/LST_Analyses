# Contents
# 0.0 Prepare environment
# 1.0 CUTS
# 2.0 MSST
# 3.0 TSO


# 0.0 Prepare environment ----
rm(list=ls())
#install.packages("lavaan")
#install.packages("MplusAutomation")
#install.packages("devtools")
#library(devtools) # if lsttheory of sumplement C of Steyer et al. (2015) has not been installed before.
#install_github("amayer2010/lsttheory", force = TRUE)
library(lavaan)
library(lsttheory)
library(MplusAutomation)
library(invgamma)
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

lambda_cstate <- round(c(1,runif(I - 1, 0.5, 1.5)), 2) # loading parameters for the latent common state
var_cstate <- rinvgamma(1, shape = 1, scale = 1) # Variance latent common state
var_ustate <- rinvgamma(I, shape = 3, scale = 1) # Variance of latent unique states

within.parameters <- list(loadings = lambda_cstate, CS.var = var_cstate, US.var = var_ustate)

# Between Paramaters

int <- rnorm(I, 3) # intercepts
lambda_ctrait <- c(1,runif(I - 1, 0.5, 1.5)) # loading parametes for the latent common trait

var_ctrait <- rinvgamma(1, shape = 1, scale = 1) # variance latent common trait
var_utrait <- rinvgamma(I, shape = 3, scale = 1) # variance latent unique traits

between.parameters <- list(intercepts = int, loadings = lambda_ctrait, CT.var = var_ctrait, UT.var = var_utrait)

cuts.data <- sim.data.cuts(N, nT, I, within.parameters = within.parameters, between.parameters = between.parameters)

# model estimation
file.name <- "sim_cuts_test"
prepareMplusData(sim_data, paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

analysis_syntax <- "USEVAR = y1 y2 y3 y4;
CLUSTER = subjn;

ANALYSIS:
TYPE = TWOLEVEL;
ESTIMATOR = ML;
ITERATIONS = 500000;" # increase H1 iterations

ml_syntax <- write.mlcuts.to.Mplus(sim_data[, -1])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)


# 2.0 MSST ----

# Conditions

N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
seed <- 123

set.seed(seed)

# Within Parameters ----

lambda_state <- round(c(1,runif(I - 1, 0.5, 1.5)), 2) # loading parameters for the latent state residual
var_state <- rinvgamma(1, shape = 1, scale = 1) # Variance latent state residual
var_m_error <- rinvgamma(I, shape = 3, scale = 1) # Variance of measurement errors

# Between Paramaters ----

int <- c(0, rnorm(I -1, 0, 0.5)) # intercepts
lambda_trait <- c(1,runif(I - 1, 0.5, 1.5)) # loading parametes for the latent trait variable

var_trait <- rinvgamma(1, shape = 1, scale = 1) # variance latent trait variable
mean_trait <- rnorm(1, 3) # mean latent trait variable


# 3.0 TSO ----

# Conditions

N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
seed <- 123

set.seed(seed)

# Within Parameters ----

lambda_state <- round(c(1,runif(I - 1, 0.5, 1.5)), 2) # loading parameters for the latent state residual
var_state <- rinvgamma(1, shape = 1, scale = 1) # Variance latent state residual
var_error <- rinvgamma(I, shape = 3, scale = 1) # Variance of latent measurement errors

ar_effect <- 0.5 # autoregressive effect on the latent state residuals

# Between Paramaters ----

int <- rnorm(I, 3)

var_ind_traits <- rinvgamma(I, shape = 1, scale = 1) # variance latent indicator trait variables

R <- matrix(rbinom(I*I, 10, prob = 0.6)/10, I)
R[lower.tri(R)] = t(R)[lower.tri(R)]
diag(R) <- 1
R
D <- diag(sd_ind_traits)
D
Sigma <- D%*%R%*%D
Sigma



