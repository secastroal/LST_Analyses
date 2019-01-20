# function to simulate tso data ----

sim.data.tso <- function(N, nT, I, within.parameters, between.parameters, 
                                       seed = 123){
  set.seed(seed)
  #Compute standard deviations
  state.sd <- sqrt(within.parameters$state.var) # state residual standard deviation
  error.sd <- sqrt(within.parameters$error.var) # errors sd
  trait.ind.sd <- sqrt(between.parameters$trait.ind.var) # trait indicators sd
  
  #Compute variance covariance matrix
  D <- diag(trait.ind.sd)
  Sigma <- D%*%between.parameters$cor.matrix%*%D
  
  trait_scores <- mvrnorm(N, rep(0, I), Sigma = Sigma) # factor indicator trait scores
  trait_scores_full <- array(trait_scores, dim = c(N,I,nT)) # array with factor trait scores
  
  # array with latent state residual in occasion n=1 and latent occasion specific residuals in occasion n>1
  
  state_scores_full <- array(NA, dim = c(N, I, nT)) 
  state_scores_full[,,1] <- rnorm(N, 0, state.sd)
  for(i in 2:nT){
    state_scores_full[,,i] <- rnorm(N, 0, state.sd) + within.parameters$ar.effect * state_scores_full[,,i-1]
  }
  rm(i)
  
  # measurement errors
  errors <- array(NA,dim = c(N, I , nT))
  for(i in 1:I){
    errors[,i,] <- rnorm(N * nT, 0, error.sd[i])
  }
  rm(i)
  
  # Complete data
  sim_data <- array(matrix(between.parameters$intercepts, nrow = N, ncol = I, byrow = TRUE), dim = c(N, I, nT)) + trait_scores_full + # Intercepts and trait scores
    state_scores_full * array(matrix(within.parameters$loadings, nrow = N, ncol = I, byrow = TRUE), dim= c(N, I, nT)) + # occasion specific scores times occasion specific lambdas
    errors # errors
  
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






