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
  ctrait_scores <- ctrait_scores %o% between.parameters$loadings # common trait scores times trait loadings
  ctrait_scores_full <- ctrait_scores[rep(1:N, each = nT), ] # full matrix with common trait scores
  #JL# Why are the loadings across the I items constrained to be equal (= columns are equal)?
  
  # unique trait scores
  utrait_scores <- matrix(NA, nrow = N, ncol = I) 
  for(i in 1:I){ #JL# fixed!
    utrait_scores[,i] <- rnorm(N, 0, UT.sd[i])
  }
  rm(i)
  utrait_scores_full <- utrait_scores[rep(1:N, each = nT), ]
  
  
  cstate_scores <- rnorm(N * nT, 0, CS.sd) # common state scores
  cstate_scores_full <- cstate_scores %o% within.parameters$loadings # full matrix with common state scores times state loadings
  
  # unique state or measurement error
  ustate_scores_full <- matrix(NA, N * nT, I)
  for(i in 1:I){
    ustate_scores_full[,i] <- rnorm(N * nT, 0, US.sd[i])
  }
  rm(i)
  
  # Complete data
  sim_data <- matrix(between.parameters$intercepts, nrow = N *nT, ncol = I, byrow = TRUE) + # intercepts
    ctrait_scores_full + # trait scores times trait lambdas
    utrait_scores_full + # unique trait scores 
    cstate_scores_full + # state scores times state lambdas
    ustate_scores_full # unique states
  #JL# An idea to manipulate: trunc(sim_data)
  
  sim_data <- data.frame(cbind(rep(1:N, each = nT), rep(1:nT, times = N), sim_data), row.names = NULL)
  
  colnames(sim_data) <- c("subjn", "time",  paste0("y", 1:I))
  
  wide_sim_data <- reshape(sim_data, v.names = paste0("y", 1:I),
                           timevar = "time", idvar="subjn", direction="wide")
  names(wide_sim_data) <- gsub("\\.", "", names(wide_sim_data))
  
  return(list( within.parameters = within.parameters,
               between.parameters = between.parameters,
               data.long = sim_data,
               data.wide = wide_sim_data))
  
  #JL#: You could have thought of coding it the other way around: First wide and then long.
  #JL#  The reason is that the matrices would have been much simpler (smaller), with less row and column
  #JL#  repetitions. Easier to debug and code.
}










