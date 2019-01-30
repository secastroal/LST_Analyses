# function to simulate msst data ----

sim.data.msst <- function(N, nT, I, within.parameters, between.parameters, 
                          seed = 123){
  set.seed(seed)
  #Compute standard deviations
  state.sd <- sqrt(within.parameters$state.var) #state sd
  error.sd <- sqrt(within.parameters$error.var) #errors sd
  trait.sd <- sqrt(between.parameters$trait.var) #trait sd
  
  trait_scores <- rnorm(N, between.parameters$trait.mean, trait.sd) # factor trait scores
  trait_scores <- trait_scores %o% between.parameters$loadings # factor trait scores times loadings 
  trait_scores_full <- trait_scores[rep(1:N, each = nT), ] # full matrix with factor trait scores
  
  state_scores <- rnorm(N * nT, 0, state.sd) # factor state scores
  state_scores_full <- state_scores %o% within.parameters$loadings # Factor state scores times state loadings
  
  # measurement errors
  errors_full <- matrix(NA, N * nT, I)
  for(i in 1:I){
    errors_full[,i] <- rnorm(N * nT, 0, error.sd[i])
  }
  rm(i)
  
  # Complete data
  sim_data <- matrix(between.parameters$intercepts, nrow = N *nT, ncol = I, byrow = TRUE) + # intercepts
    trait_scores_full + # trait scores times trait lambdas
    state_scores_full + # state scores times state lambdas
    errors_full # errors
  
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

# Modified function to allow simulating msst data with time variant parameters ----
# This function simulates the data in wide format.
# To make it simple, the parameters provided are always time invariant parameters. When time.invariant is false,
#   time variant parameters are generated based on the time invariant parameters provided. Loadings are always 
#   sampled from the seq(0.5, 1.2, by = 0.1). Intercepts and variances are sampled from an interval centered on
#   the mean of the respective time invariant parameters, which is half the mean long. For example, given the 
#   time-invariant intercetps c(2, 3, 4), the new time-variant intercepts are sampled from ten numbers within
#   the interval c(1.5, 4.5). # Intercepts and variances are sampled in this way because they have a huge 
#   impact on the final simulated data


sim.data.msst.tv <- function(N, nT, I, within.parameters, between.parameters, time.invariant = TRUE,
                          seed = 123){
  set.seed(seed)
  
  # Put parameters in matrices when time.invariant is TRUE
  if(time.invariant){
    # within parameters
    within.parameters$loadings <- matrix(within.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
    within.parameters$state.var <- rep(within.parameters$state.var, nT)
    within.parameters$error.var <- matrix(within.parameters$error.var, nrow = nT, ncol = I, byrow = TRUE)
    
    # between parameters
    between.parameters$intercepts <- matrix(between.parameters$intercepts, nrow = nT, ncol = I, byrow = TRUE)
    between.parameters$loadings <- matrix(between.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
  }
  
  if(!time.invariant){
    #within parameters
    within.parameters$loadings <- matrix(sample(seq(0.5, 1.2, by = 0.1), size = nT * I, replace = TRUE), nrow = nT, 
                                         ncol = I, byrow = TRUE)
    within.parameters$loadings[ , 1] <- 1
    within.parameters$state.var <- sample(seq(within.parameters$state.var * (3/4), within.parameters$state.var * (5/4), 
                                              length.out = 11), size = nT, replace = TRUE)
    within.parameters$error.var <- matrix(sample(seq(mean(within.parameters$error.var) * (3/4), 
                                                     mean(within.parameters$error.var) * (5/4), length.out = 11), 
                                                 size = nT * I, replace = TRUE), nrow = nT, ncol = I, byrow = TRUE)
    
    #between parameters
    between.parameters$intercepts <- matrix(sample(seq(mean(between.parameters$intercepts) *  (3/4), mean(between.parameters$intercepts) *  (5/4), 
                                                       length.out = 11), size = nT * I, replace = TRUE), nrow = nT, 
                                            ncol = I, byrow = TRUE)
    between.parameters$intercepts[1,1] <- 0
    between.parameters$loadings <- matrix(sample(seq(0.5, 1.2, by = 0.1), size = nT * I, replace = TRUE), nrow = nT, 
                                          ncol = I, byrow = TRUE)
    between.parameters$loadings[1,1] <- 1
  }
  
  #Compute standard deviations
  state.sd <- sqrt(within.parameters$state.var) #state sd
  error.sd <- sqrt(within.parameters$error.var) #errors sd
  trait.sd <- sqrt(between.parameters$trait.var) #trait sd
  
  # Generate data in wide format ----
  
  # Latent trait variable
  trait_scores <- rnorm(N, between.parameters$trait.mean, trait.sd) # factor trait scores
  trait_scores_full <- trait_scores %o% c(t(between.parameters$loadings)) # factor trait scores times loadings 
  
  #Latent state residuals
  state_scores_full <- matrix(NA, nrow = N,  ncol = nT * I)
  for(i in 1:nT){ #Generate latent state residuals scores times the corresponding loadings on each time.
    state_scores_full[, (((i - 1) * I) + 1):(i * I)] <- rnorm(N, 0, state.sd[i]) %o% within.parameters$loadings[i, ]
  }
  rm(i)
  
  # measurement errors
  errors_full <- matrix(NA, nrow = N, ncol = nT * I)
  for(j in 1:nT){
    for(i in 1:I){
      errors_full[,((j-1)*I)+i] <- rnorm(N, 0, error.sd[j,i])
    }
  }
  rm(i, j)
  
  # Repeat intercepts in the full matrix
  intercepts_full <- matrix(c(t(between.parameters$intercepts)), nrow = N, ncol = nT * I, byrow = TRUE)
  
  # Complete data ----
  sim_data <- intercepts_full + # intercepts
    trait_scores_full + # trait scores times trait lambdas
    state_scores_full + # state scores times state lambdas
    errors_full # errors
  
  sim_data <- data.frame(cbind(1:N, sim_data))
  
  colnames(sim_data) <- c("subjn",  paste0(paste0("y", 1:I), ".", rep(1:nT, each = I)))
  
  long_sim_data <- reshape(sim_data, varying = colnames(sim_data)[-1],  idvar = "subjn", direction = "long", 
                           sep = ".")
  
  names(sim_data) <- gsub("\\.", "", names(sim_data))
  
  return(list( within.parameters = within.parameters,
               between.parameters = between.parameters,
               data.long = long_sim_data,
               data.wide = sim_data))
}
