# Prepare environment---
library(invgamma)
library(MASS)

# function to simulate msst data ----

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

sd_state <- sqrt(var_state) # sd latent state residual
sd_error <- sqrt(var_error) # sd latent measurement errors

ar_effect <- 0.5 # autoregressive effect on the latent state residuals

# Between Paramaters ----

int <- rnorm(I, 3)

var_ind_traits <- rinvgamma(I, shape = 1, scale = 1) # variance latent indicator trait variables

sd_ind_traits <- sqrt(var_ind_traits) # sd latent indicator trait variables

R <- matrix(rbinom(I*I, 10, prob = 0.6)/10, I)
R[lower.tri(R)] = t(R)[lower.tri(R)]
diag(R) <- 1
R
D <- diag(sqrt(var_ind_traits))
D
Sigma <- D%*%R%*%D
Sigma


# Data simulation ----

trait_scores <- mvrnorm(N, 0, Sigma = Sigma) # factor indicator trait scores
trait_scores_full <- array(trait_scores, dim = c(N,I,nT)) # matrix with factor trait scores

# array with latent state residual in occasion n=1 and latent occasion specific residuals in occasion n>1

state_scores_full <- array(NA, dim = c(N, I, nT)) 
state_scores_full[,,1] <- rnorm(N, 0, sd_state)
for(i in 2:nT){
  state_scores_full[,,i] <- rnorm(N, 0, sd_state) + ar_effect * state_scores_full[,,i-1]
}
rm(i)


# measurement errors
errors <- array(NA,dim = c(N, I , nT))
for(i in 1:I){
  errors[,i,] <- rnorm(N * nT, 0, sd_error[i])
}
rm(i)
errors

# Complete data
sim_data <- array(matrix(int, nrow = N, ncol = I, byrow = TRUE), dim = c(N, I, nT)) + trait_scores_full + # Intercepts and trait scores
  state_scores_full * array(matrix(lambda_state, nrow = N, ncol = I, byrow = TRUE), dim= c(N, I, nT)) + # occasion specific scores times occasion specific lambdas
  errors # errors

sim_data <- aperm(sim_data,c(3,1,2))

sim_data <- matrix(sim_data, N*nT, I)

sim_data <- data.frame(cbind(rep(1:N, each = nT), rep(1:nT, times = N), sim_data), row.names = NULL)

colnames(sim_data) <- c("subjn", "time", paste0("y", 1:I))



# model estimation ----
file.name <- "sim_msst_test"
prepareMplusData(sim_data, paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

analysis_syntax <- "USEVAR = y1 y2 y3 y4;
CLUSTER = subjn;

ANALYSIS:
TYPE = TWOLEVEL;
ESTIMATOR = ML;
ITERATIONS = 500000;" # increase H1 iterations

ml_syntax <- write.mlmsst.to.Mplus(sim_data[, -1])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)





