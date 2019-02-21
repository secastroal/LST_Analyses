# function to simulate tso data ----

sim.data.tso <- function(N, nT, I, within.parameters, between.parameters, 
                                       na.prop = 0, seed = 123){
  set.seed(seed)
  #Compute standard deviations
  state.sd <- sqrt(within.parameters$state.var) # state residual standard deviation
  error.sd <- sqrt(within.parameters$error.var) # errors sd
  trait.ind.sd <- sqrt(between.parameters$trait.ind.var) # trait indicators sd
  
  #Compute variance covariance matrix
  D <- diag(trait.ind.sd)
  Sigma <- D%*%between.parameters$cor.matrix%*%D
  between.parameters$Sigma <- Sigma
  
  #Latent trait indicators scores
  trait_scores <- mvrnorm(N, rep(0, I), Sigma = Sigma) # factor indicator trait scores
  trait_scores_full <- array(trait_scores, dim = c(N,I,nT)) # array with factor trait scores
  
  # array with latent state residual in occasion n=1 and latent occasion specific residuals in occasion n>1
  
  state_scores_full <- array(NA, dim = c(N, I, nT)) 
  state_scores_full[, 1, 1] <- rnorm(N, 0, state.sd)
  for(i in 2:nT){
    state_scores_full[, 1, i] <- rnorm(N, 0, state.sd) + within.parameters$ar.effect * state_scores_full[, 1, i-1]
    state_scores_full[, , i-1] <- state_scores_full[, 1, i-1] %o% within.parameters$loadings
  }
  rm(i)
  
  state_scores_full[, , nT] <- state_scores_full[, 1, nT] %o% within.parameters$loadings
  
  # measurement errors
  errors_full <- array(NA,dim = c(N, I , nT))
  for(i in 1:I){
    errors_full[,i,] <- rnorm(N * nT, 0, error.sd[i])
  }
  rm(i)
  
  # Complete data
  sim_data <- array(matrix(between.parameters$intercepts, nrow = N, ncol = I, byrow = TRUE), dim = c(N, I, nT)) + # intercepts
    trait_scores_full + # trait scores
    state_scores_full + # latent occasion specific residuals scores
    errors_full # errors
  
  sim_data <- aperm(sim_data,c(3,1,2))
  
  sim_data <- matrix(sim_data, N*nT, I)
  
  sim_data <- data.frame(cbind(rep(1:N, each = nT), rep(1:nT, times = N), sim_data), row.names = NULL)
  
  colnames(sim_data) <- c("subjn", "time", paste0("y", 1:I))
  
  #introduce missing completely at random
  
  if (na.prop != 0) {
    pos.NA   <- matrix(rbinom(N * I * nT, 1, 1 - na.prop), nrow = N * nT)
    pos.NA[pos.NA == 0] <- NA
    sim_data[, 3:(I+2)]     <- sim_data[,3:(I + 2)] * pos.NA
    rm(pos.NA)       # clean
  }
  
  #reshape into wide format
  wide_sim_data <- reshape(sim_data, v.names = paste0("y", 1:I),
                           timevar = "time", idvar="subjn", direction="wide")
  names(wide_sim_data) <- gsub("\\.", "", names(wide_sim_data))
  
  return(list( within.parameters = within.parameters,
               between.parameters = between.parameters,
               data.long = sim_data,
               data.wide = wide_sim_data))
}

# Modified function to allow simulating cuts data with time variant parameters ----
# This function simulates the data in wide format.
# To make it simple, the parameters provided are always time invariant parameters. When time.invariant is false,
#   time variant parameters are generated based on the time invariant parameters provided. Variances are
#   turned time-variant, in such a way that they increase or decrease by 0.1 between the supplied variance and
#   the supplied variance plus one. The autoregressive effect just increase over time between c(0.2, 0.8).
#   Loadings and intercetps remain time invariant.

sim.data.tso.tv <- function(N, nT, I, within.parameters, between.parameters, time.invariant = TRUE, 
                         seed = 123){
  set.seed(seed)
  
  # Put parameters in matrices when time.invariant is TRUE
  if(time.invariant){
    # within parameters
    within.parameters$loadings <- matrix(within.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
    within.parameters$state.var <- rep(within.parameters$state.var, nT)
    within.parameters$error.var <- matrix(within.parameters$error.var, nrow = nT, ncol = I, byrow = TRUE)
    within.parameters$ar.effect <- rep(within.parameters$ar.effect, nT-1)
    
    # between parameters
    between.parameters$intercepts <- matrix(between.parameters$intercepts, nrow = nT, ncol = I, byrow = TRUE)
    between.parameters$loadings <- matrix(1, nrow = nT, ncol = I, byrow = TRUE)
  }
  
  if(!time.invariant){
    # within parameters
    within.parameters$loadings <- matrix(within.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
    within.parameters$ar.effect <- seq(0.2, 0.8, length.out = nT - 1)
    
    var_up <- seq(within.parameters$state.var, within.parameters$state.var+1, by = 0.1)
    var_cycle <- c(var_up, rev(var_up[-c(1, length(var_up))]))
    within.parameters$state.var <- rep(var_cycle, length.out = nT)
    rm(var_up, var_cycle)
    
    old.error.var <- within.parameters$error.var
    
    within.parameters$error.var <- matrix(NA, nrow = nT, ncol = I)
    
    for(i in 1:I){
      var_up <- seq(old.error.var[i], old.error.var[i]+1, by = 0.1)
      var_cycle <- c(var_up, rev(var_up[-c(1, length(var_up))]))
      within.parameters$error.var[,i] <- rep(var_cycle, length.out = nT)
      rm(var_up, var_cycle)
    }
    rm(old.error.var, i)
    
    # between parameters
    between.parameters$intercepts <- matrix(between.parameters$intercepts, nrow = nT, ncol = I, byrow = TRUE)
    between.parameters$loadings <- matrix(1, nrow = nT, ncol = I, byrow = TRUE)
  }
  
  #Compute standard deviations
  state.sd <- sqrt(within.parameters$state.var) # state residual standard deviation
  error.sd <- sqrt(within.parameters$error.var) # errors sd
  trait.ind.sd <- sqrt(between.parameters$trait.ind.var) # trait indicators sd
  
  #Compute variance covariance matrix
  D <- diag(trait.ind.sd)
  Sigma <- D%*%between.parameters$cor.matrix%*%D
  between.parameters$Sigma <- Sigma
  
  
  #Generate simulated data ----
  
  # Latent trait indicators scores
  trait_scores_full <- array(NA, dim = c(N,I,nT)) # array with factor trait scores
  trait_scores <- mvrnorm(N, rep(0, I), Sigma = Sigma) # factor indicator trait scores
  for(i in 1:nT){
    trait_scores_full[ , , i] <- t(between.parameters$loadings[i, ] * t(trait_scores))
  }
  
  # array with latent state residual in occasion n=1 and latent occasion specific residuals in occasion n>1
  
  state_scores_full <- array(NA, dim = c(N, I, nT)) 
  state_scores_full[, 1, 1] <- rnorm(N, 0, state.sd[1])
  for(i in 2:nT){
    state_scores_full[, 1, i] <- rnorm(N, 0, state.sd[i]) + within.parameters$ar.effect[i-1] * state_scores_full[ , 1, i-1]
    state_scores_full[, , i-1] <- state_scores_full[, 1, i-1] %o% within.parameters$loadings[i-1, ]
  }
  rm(i)
  
  state_scores_full[, , nT] <- state_scores_full[, 1, nT] %o% within.parameters$loadings[nT, ]
  
  
  # measurement errors
  errors_full <- array(NA, dim = c(N, I , nT))
  for(i in 1:I){
    for(j in 1:nT){
      errors_full[, i, j] <- rnorm(N, 0, error.sd[j, i])
      }
  }
  rm(i)
  
  # intercepts full array
  intercepts_full <- array(NA, dim = c(N, I, nT))
  for(i in 1:N){
    intercepts_full[i, ,] <- t(between.parameters$intercepts)
  }
  
  # Complete data ----
  sim_data <- intercepts_full + # intercepts
    trait_scores_full + # trait scores
    state_scores_full + # occasion specific scores times occasion specific lambdas
    errors_full # errors
  
  sim_data <- aperm(sim_data,c(3,1,2))
  
  sim_data <- matrix(sim_data, N*nT, I)
  
  sim_data <- data.frame(cbind(rep(1:N, each = nT), rep(1:nT, times = N), sim_data), row.names = NULL)
  
  colnames(sim_data) <- c("subjn", "time", paste0("y", 1:I))
  
  wide_sim_data <- reshape(sim_data, v.names = paste0("y", 1:I),
                           timevar = "time", idvar="subjn", direction="wide")
  names(wide_sim_data) <- gsub("\\.", "", names(wide_sim_data))
  
  return(list( within.parameters = within.parameters,
               between.parameters = between.parameters,
               data.long = sim_data,
               data.wide = wide_sim_data))
}





