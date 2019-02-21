# function to simulate msst data ----

sim.data.msst <- function(N, nT, I, within.parameters, between.parameters, 
                          na.prop = 0, seed = 123){
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
#   time variant parameters are generated based on the time invariant parameters provided. Only variances are
#   turned time-variant, in such a way that they increase or decrease by 0.1 between the supplied variance and
#   the supplied variance plus one.


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
    # within parameters
    within.parameters$loadings <- matrix(within.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
    
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
    between.parameters$loadings <- matrix(between.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
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
