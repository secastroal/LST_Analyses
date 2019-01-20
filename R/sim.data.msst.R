# function to simulate msst data ----

sim.data.msst <- function(N, nT, I, within.parameters, between.parameters, 
                          seed = 123){
  set.seed(seed)
  #Compute standard deviations
  state.sd <- sqrt(within.parameters$state.var) #state sd
  error.sd <- sqrt(within.parameters$error.var) #errors sd
  trait.sd <- sqrt(between.parameters$trait.var) #trait sd
  
  trait_scores <- rnorm(N, between.parameters$trait.mean, trait.sd) # factor trait scores
  trait_scores_full <- matrix(rep(trait_scores, each = nT), nrow = N * nT, ncol = I, byrow = FALSE) # matrix with factor trait scores
  
  state_scores <- rnorm(N * nT, 0, state.sd) # factor state scores
  state_scores_full <- matrix(state_scores, nrow = N * nT, ncol = I , byrow = FALSE) # full matrix with factor state scores
  
  # measurement errors
  errors <- matrix(NA, N * nT, I)
  for(i in 1:I){
    errors[,i] <- rnorm(N * nT, 0, error.sd[i])
  }
  rm(i)
  
  # Complete data
  sim_data <- matrix(between.parameters$intercepts, nrow = N *nT, ncol = I, byrow = TRUE) + # intercepts
    trait_scores_full * matrix(between.parameters$loadings, nrow = N *nT, ncol = I, byrow = TRUE) + # trait scores times trait lambdas
    state_scores_full * matrix(within.parameters$loadings, nrow = N *nT, ncol = I, byrow = TRUE) + # state scores times state lambdas
    errors # errors
  
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


