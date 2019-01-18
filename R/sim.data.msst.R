# Prepare environment---
library(invgamma)

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
var_m_error <- rinvgamma(I, shape = 3, scale = 1) # Variance of measurement errors

sd_state <- sqrt(var_state) # sd latent state residuals
sd_m_error <- sqrt(var_m_error) # sd measurement errors

# Between Paramaters ----

int <- c(0, rnorm(I -1, 0, 0.5)) # intercepts
lambda_trait <- c(1,runif(I - 1, 0.5, 1.5)) # loading parametes for the latent trait variable

var_trait <- rinvgamma(1, shape = 1, scale = 1) # variance latent trait variable
mean_trait <- rnorm(1, 3) # mean latent trait variable

sd_trait <- sqrt(var_trait) # sd latent trait variable

# Data simulation ----

sim.data.msst

trait_scores <- rnorm(N, mean_trait, sd_trait) # factor trait scores
trait_scores_full <- matrix(rep(trait_scores, each = nT), nrow = N * nT, ncol = I, byrow = FALSE) # matrix with factor trait scores

state_scores <- rnorm(N * nT, 0, sd_state) # factor state scores
state_scores_full <- matrix(state_scores, nrow = N * nT, ncol = I , byrow = FALSE) # full matrix with factor state scores

# measurement errors
errors <- matrix(NA, N * nT, I)
for(i in 1:I){
  errors[,i] <- rnorm(N * nT, 0, sd_m_error[i])
}
rm(i)

# Complete data
sim_data <- matrix(int, nrow = N *nT, ncol = I, byrow = TRUE) + # intercepts
  trait_scores_full * matrix(lambda_trait, nrow = N *nT, ncol = I, byrow = TRUE) + # trait scores times trait lambdas
  state_scores_full * matrix(lambda_state, nrow = N *nT, ncol = I, byrow = TRUE) + # state scores times state lambdas
  errors # errors

sim_data <- data.frame(cbind(rep(1:N, each = nT), rep(1:nT, times = N), sim_data), row.names = NULL)

colnames(sim_data) <- c("subjn", "times", paste0("y", 1:I))



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





