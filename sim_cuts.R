# Prepare environment---
rm(list=ls())
library(invgamma)

# function to simulate msst data ----

# Conditions

N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
seed <- 123

set.seed(seed)

# Within Parameters ----

lambda_cstate <- round(c(1,runif(I - 1, 0.5, 1.5)), 2) # loading parameters for the latent common state
var_cstate <- rinvgamma(1, shape = 1, scale = 1) # Variance latent common state
var_ustate <- rinvgamma(I, shape = 3, scale = 1) # Variance of latent unique states

sd_cstate <- sqrt(var_cstate) # sd latent common state
sd_ustate <- sqrt(var_ustate) # sd unique states

# Between Paramaters ----

int <- rnorm(I, 3) # intercepts
lambda_ctrait <- c(1,runif(I - 1, 0.5, 1.5)) # loading parametes for the latent common trait

var_ctrait <- rinvgamma(1, shape = 1, scale = 1) # variance latent common trait
var_utrait <- rinvgamma(I, shape = 3, scale = 1) # variance latent unique traits

sd_ctrait <- sqrt(var_ctrait) # sd latent  common trait
sd_utrait <- sqrt(var_utrait) # sd latent  common trait



# Data simulation ----

ctrait_scores <- rnorm(N, 0, sd_ctrait) # common trait scores
ctrait_scores_full <- matrix(rep(ctrait_scores, each = nT), nrow = N * nT, ncol = I, byrow = FALSE) # matrix with common trait scores

# unique trait scores
utrait_scores <- matrix(NA, nrow = N, ncol = I) 
for(i in 1:4){
  utrait_scores[,i] <- rnorm(N, 0, sd_utrait[i])
}
rm(i)
utrait_scores_full <- utrait_scores[rep(1:N, each = nT), ]


cstate_scores <- rnorm(N * nT, 0, sd_cstate) # common state scores
cstate_scores_full <- matrix(cstate_scores, nrow = N * nT, ncol = I , byrow = FALSE) # full matrix with common state scores

# unique state or measurement error
ustate_scores_full <- matrix(NA, N * nT, I)
for(i in 1:I){
  ustate_scores_full[,i] <- rnorm(N * nT, 0, sd_ustate[i])
}
rm(i)

# Complete data
sim_data <- matrix(int, nrow = N *nT, ncol = I, byrow = TRUE) + # intercepts
  ctrait_scores_full * matrix(lambda_ctrait, nrow = N *nT, ncol = I, byrow = TRUE) + # trait scores times trait lambdas
  utrait_scores_full + # unique trait scores 
  cstate_scores_full * matrix(lambda_cstate, nrow = N *nT, ncol = I, byrow = TRUE) + # state scores times state lambdas
  ustate_scores_full # unique states

sim_data <- data.frame(cbind(rep(1:N, each = nT), rep(1:nT, times = N), sim_data), row.names = NULL)

colnames(sim_data) <- c("subjn", "time",  paste0("y", 1:I))



# model estimation ----
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





