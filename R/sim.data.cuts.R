# function to simulate cuts data ----

sim.data.cuts <- function(N, nT, I, within.parameters, between.parameters, 
                          na.prop = 0, seed = 123){
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
  
  #JL#: You could have thought of coding it the other way around: First wide and then long.
  #JL#  The reason is that the matrices would have been much simpler (smaller), with less row and column
  #JL#  repetitions. Easier to debug and code.
}

# Modified function to allow simulating cuts data with time variant parameters ----
# This function simulates the data in wide format.
# To make it simple, the parameters provided are always time invariant parameters. When time.invariant is false,
#   time variant parameters are generated based on the time invariant parameters provided. Loadings are always 
#   sampled from the seq(0.5, 1.2, by = 0.1). Intercepts and variances are sampled from an interval centered on
#   the mean of the respective time invariant parameters, which is half the mean long. For example, given the 
#   time-invariant intercetps c(2, 3, 4), the new time-variant intercepts are sampled from ten numbers within
#   the interval c(1.5, 4.5). # Intercepts and variances are sampled in this way because they have a huge 
#   impact on the final simulated data



sim.data.cuts.tv <- function(N, nT, I, within.parameters, between.parameters, time.invariant = TRUE, 
                          seed = 123){
  set.seed(seed)
  
  # Put parameters in matrices when time.invariant is TRUE
  if(time.invariant){
    # within parameters
    within.parameters$loadings <- matrix(within.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
    within.parameters$CS.var <- rep(within.parameters$CS.var, nT)
    within.parameters$US.var <- matrix(within.parameters$US.var, nrow = nT, ncol = I, byrow = TRUE)
    
    # between parameters
    between.parameters$intercepts <- matrix(between.parameters$intercepts, nrow = nT, ncol = I, byrow = TRUE)
    between.parameters$loadings <- matrix(between.parameters$loadings, nrow = nT, ncol = I, byrow = TRUE)
    between.parameters$loadings.UT <- matrix(1, nrow = nT, ncol = I)
  }
  
  if(!time.invariant){
    #within parameters
    within.parameters$loadings <- matrix(sample(seq(0.5, 1.2, by = 0.1), size = nT * I, replace = TRUE), nrow = nT, 
                                         ncol = I, byrow = TRUE)
    within.parameters$loadings[ , 1] <- 1
    within.parameters$CS.var <- sample(seq(within.parameters$CS.var * (3/4), within.parameters$CS.var * (5/4), 
                                           length.out = 11), size = nT, replace = TRUE)
    within.parameters$US.var <- matrix(sample(seq(mean(within.parameters$US.var) * (3/4), mean(within.parameters$US.var) * (5/4), 
                                                  length.out = 11), size = nT * I, replace = TRUE), nrow = nT, 
                                       ncol = I, byrow = TRUE)
    
    #between parameters
    between.parameters$intercepts <- matrix(sample(seq(mean(between.parameters$intercepts) *  (3/4), mean(between.parameters$intercepts) *  (5/4), 
                                                       length.out = 11), size = nT * I, replace = TRUE), nrow = nT, 
                                            ncol = I, byrow = TRUE)
    between.parameters$loadings <- matrix(sample(seq(0.5, 1.2, by = 0.1), size = nT * I, replace = TRUE), nrow = nT, 
                                          ncol = I, byrow = TRUE)
    between.parameters$loadings[1,1] <- 1
    between.parameters$loadings.UT <- matrix(sample(seq(0.5, 1.2, by = 0.1), size = nT * I, replace = TRUE), nrow = nT, 
                                             ncol = I, byrow = TRUE)
    between.parameters$loadings.UT[1, ] <- 1
  }
  
  # Compute standard deviations
  CS.sd <- sqrt(within.parameters$CS.var) #common state sd
  US.sd <- sqrt(within.parameters$US.var) #unique state sd
  CT.sd <- sqrt(between.parameters$CT.var) #common trait sd
  UT.sd <- sqrt(between.parameters$UT.var) #unique trait sd
  
  #Generate data in wide format----
  # A more optimal way to code this part is to define all the NA matrices before the loops and do only one loop 
  #     for(i in 1:nT){}. Thus, the function would only use three loops instead of 5. However with more loops, 
  #     the code is clearer.
  
  # common trait scores
  ctrait_scores <- rnorm(N, 0, CT.sd) # common trait scores
  ctrait_scores_full <- ctrait_scores %o% c(t(between.parameters$loadings)) # common trait scores times trait loadings
  
  # unique trait scores
  utrait_scores <- matrix(NA, nrow = N, ncol = I) 
  for(i in 1:I){ #JL# fixed!
    utrait_scores[,i] <- rnorm(N, 0, UT.sd[i])
  }
  rm(i)
  
  utrait_scores_full <- matrix(NA, nrow = N, ncol = I * nT)
  for(i in 1:nT){ #unique traits times loadings per each time
    utrait_scores_full[, (((i - 1) * I) + 1):(i * I)] <- t(between.parameters$loadings.UT[i, ] * t(utrait_scores))
  }
  rm(i)
  
  #Common state scores
  cstate_scores_full <- matrix(NA, nrow = N, ncol = I * nT) #matrix to store common state scores
  for(i in 1:nT){ #Generate common state score times the corresponding loadings on each time
    cstate_scores_full[, (((i - 1) * I) + 1):(i * I)] <- rnorm(N, 0, CS.sd[i]) %o% within.parameters$loadings[i, ]
  }
  rm(i)
  
  # unique state or measurement error
  ustate_scores_full <- matrix(NA, nrow = N, ncol = I * nT)
  for(j in 1:nT){
    for(i in 1:I){
      ustate_scores_full[,((j-1)*I)+i] <- rnorm(N, 0, US.sd[j,i])
    }
  }
  rm(i, j)
  
  #intercepts complete matrix
  intercepts_full <- matrix(c(t(between.parameters$intercepts)), nrow = N, ncol = nT * I, byrow = TRUE)
  
  # Complete data ----
  sim_data <- intercepts_full + # intercepts
    ctrait_scores_full + # trait scores times trait lambdas
    utrait_scores_full + # unique trait scores 
    cstate_scores_full + # state scores times state lambdas
    ustate_scores_full # unique states
  #JL# An idea to manipulate: trunc(sim_data)
  
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

