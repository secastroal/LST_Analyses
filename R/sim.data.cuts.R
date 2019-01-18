# function to simulate cuts data ----

sim.data.cuts <- function(N, nT, I, within.parameters, between.parameters, 
                          seed = 123){
  set.seed(seed)
  # Compute standard deviations
  CS.sd <- sqrt(within.parameters$CS.var) #common state sd
  US.sd <- sqrt(within.parameters$US.var) #unique state sd
  CT.sd <- sqrt(between.parameters$CT.var) #common trait sd
  UT.sd <- sqrt(between.parameters$UT.var) #unique trait sd
  
  ctrait_scores <- rnorm(N, 0, CT.sd) # common trait scores
  ctrait_scores_full <- matrix(rep(ctrait_scores, each = nT), nrow = N * nT, ncol = I, byrow = FALSE) # matrix with common trait scores
  
  # unique trait scores
  utrait_scores <- matrix(NA, nrow = N, ncol = I) 
  for(i in 1:4){
    utrait_scores[,i] <- rnorm(N, 0, UT.sd[i])
  }
  rm(i)
  utrait_scores_full <- utrait_scores[rep(1:N, each = nT), ]
  
  
  cstate_scores <- rnorm(N * nT, 0, CS.sd) # common state scores
  cstate_scores_full <- matrix(cstate_scores, nrow = N * nT, ncol = I , byrow = FALSE) # full matrix with common state scores
  
  # unique state or measurement error
  ustate_scores_full <- matrix(NA, N * nT, I)
  for(i in 1:I){
    ustate_scores_full[,i] <- rnorm(N * nT, 0, US.sd[i])
  }
  rm(i)
  
  # Complete data
  sim_data <- matrix(between.parameters$intercepts, nrow = N *nT, ncol = I, byrow = TRUE) + # intercepts
    ctrait_scores_full * matrix(between.parameters$loadings, nrow = N *nT, ncol = I, byrow = TRUE) + # trait scores times trait lambdas
    utrait_scores_full + # unique trait scores 
    cstate_scores_full * matrix(within.parameters$loadings, nrow = N *nT, ncol = I, byrow = TRUE) + # state scores times state lambdas
    ustate_scores_full # unique states
  
  sim_data <- data.frame(cbind(rep(1:N, each = nT), rep(1:nT, times = N), sim_data), row.names = NULL)
  
  colnames(sim_data) <- c("subjn", "time",  paste0("y", 1:I))
  
  wide_sim_data <- reshape(sim_data, v.names = paste0("y", 1:I),
                           timevar = "time", idvar="subjn", direction="wide")
  
  return(list( within.parameters = within.parameters,
               between.parameters = between.parameters,
               data.long = sim_data,
               data.wide = wide_sim_data))
  
}










