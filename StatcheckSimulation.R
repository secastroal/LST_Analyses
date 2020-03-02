# Stationarity Check for simulated data in first version of the simulation.

# As noted by one of the reviewers, the stationarity assumption should also hold when using 
# state-trait SEMs for intensive longitudinal data. In this document, we check how many subjects
# have at least one of their timeseries non-stationary. 

rm(list=ls())
#library(devtools) # if lsttheory of sumplement C of Steyer et al. (2015) has not been installed before.
#install_github("amayer2010/lsttheory", force = TRUE)
#library(lavaan) 
#library(lsttheory)
library(MplusAutomation)
library(MASS)
library(xtable)
library(tseries)
# Easy code to source all the needed files
file.sources <- paste0("R/", list.files(paste0(getwd(), "/R")))
sapply(file.sources,source,.GlobalEnv)
rm(file.sources)
folder <- "Mplus_Simulation" #Folder to store results and all Mplus files

# Create folders to store output in the working directory ----

if(!dir.exists(paste(getwd(), folder, sep = "/"))){
  dir.create(paste(getwd(), folder, "StationarityCheck", sep = "/"), recursive = TRUE)
}else{
  if(!dir.exists(paste(getwd(), folder, "StationarityCheck", sep = "/"))){
    dir.create(paste(getwd(), folder, "StationarityCheck", sep = "/"))
  }
}

# add / to folder
folder <- paste0(folder, "/")

# 1.0 Set up conditions ----

# Fix conditions
N <- 100 # Sample size
I <- 4 # Number of variables

# Manipulated conditions
timepoints <- c(30, 60, 90) # Number of measurement occasions
missingness <- c(0, 0.1) # Proportion of NAs
dataModel_to_sim <- c("msst", "cuts", "tso") # Model to simulate the data

Cond <- expand.grid(timepoints, missingness, dataModel_to_sim)
names(Cond) <- c("nT", "na.prop", "dataModel")

rm(timepoints, missingness, dataModel_to_sim)

# 2.0 Set up output matrices ----

R <- 100 # Number of replicas

# Create Matrix to count how many timeseries are non-stationary 

statchecksum <- matrix(NA, R, dim(Cond)[1])
colnames(statchecksum) <- paste0("Cond_", 1:18)

# 3.0 Simulation loop ----

time0 <- proc.time()
for(cond in 1:18){
  print(paste("Condition", cond)) # print condition number to track progress
  
  for(r in 1:R){
    print(paste("Replication", r)) #print replication number to track progress
    seed <- 1000 * cond + r
    set.seed(seed)
    nT <- Cond[cond, 1]
    na.prop <- Cond[cond, 2]
    dataModel <- Cond[cond, 3]
    
    # Data simulation base on dataModel ----
    
    if(dataModel == "msst"){
      # Within Parameters
      
      state_loadings <- c(1, 0.5, 1.3, 0.8) # loading parameters for the latent state
      var_state <- 2 # Variance latent state residual
      var_m_error <- c(1, 0.5, 1.5, 0.8) # Variance of measurement errors
      
      within.parameters <- list(loadings = state_loadings, state.var = var_state, error.var = var_m_error)
      
      # Between Paramaters
      
      intercepts <- seq(0, by = 0.2, length.out = I) # intercepts
      trait_loadings <- c(1, 0.8, 1.2, 0.9) # loading parametes for the latent trait
      
      var_trait <- 2 # variance latent trait variable
      mean_trait <- 4 # mean latent trait variable
      
      between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, trait.mean = mean_trait,
                                 trait.var = var_trait)
      
      rm(state_loadings, var_state, var_m_error, intercepts, trait_loadings, var_trait, 
         mean_trait)
      
      # Simulate data
      data <- sim.data.msst(N, nT, I, within.parameters = within.parameters, na.prop = na.prop,
                            between.parameters = between.parameters, seed = seed)
      
      rm(within.parameters, between.parameters)
    }
    
    if(dataModel == "cuts"){
      
      # Within Parameters
      
      state_loadings <- c(1, 0.5, 1.3, 0.8) # loading parameters for the latent common state
      var_CS <- 2 # Variance latent common state
      var_US <- c(1, 0.5, 1.5, 0.8) # Variance of latent unique states
      
      within.parameters <- list(loadings = state_loadings, CS.var = var_CS, US.var = var_US)
      
      # Between Paramaters
      
      intercepts <- seq(2, by = 0.5, length.out = I) # intercepts
      trait_loadings <- c(1, 0.8, 1.2, 0.9) # loading parametes for the latent common trait
      var_CT <- 1.5 # variance latent common trait
      var_UT <- c(0.5, 1, 0.3, 0.8) # variance latent unique traits
      
      between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, CT.var = var_CT, UT.var = var_UT)
      
      rm(state_loadings, var_CS, var_US, intercepts, trait_loadings, var_CT, var_UT)
      
      # Simulate data
      data <- sim.data.cuts(N, nT, I, within.parameters = within.parameters, na.prop = na.prop,
                            between.parameters = between.parameters, seed = seed)
      
      rm(within.parameters, between.parameters)
    }
    
    if(dataModel == "tso"){
      # Within Parameters
      
      state_loadings <- c(1, 0.5, 1.3, 0.8) # loading parameters for the latent state
      var_state <- 2 # Variance latent state residual
      var_error <- c(1, 0.5, 1.5, 0.8) # Variance of latent measurement errors
      
      ar_effect <- 0.5 # autoregressive effect on the latent state residuals
      
      within.parameters <- list(loadings = state_loadings, ar.effect = ar_effect, error.var = var_error,
                                state.var = var_state)
      
      rm(state_loadings, var_state, var_error, ar_effect)
      
      # Between Paramaters 
      
      intercepts <- seq(2, by = 0.5, length.out = I) # intercepts
      
      var_ind_traits <- c(2, 1.5, 2.5, 1.75) # variance latent indicator trait variables
      
      # Create positive definite correlation matrix
      Rcor <- matrix(c(1, .8, .9, .9,
                       .8, 1, .8, .7,
                       .9, .8, 1, .7,
                       .9, .7, .7, 1), I, I, byrow = TRUE)
      
      between.parameters <- list(intercepts = intercepts, trait.ind.var = var_ind_traits, 
                                 cor.matrix = Rcor)
      
      rm(intercepts, var_ind_traits, Rcor)
      
      # Simulate data
      data <- sim.data.tso(N, nT, I, within.parameters = within.parameters, na.prop = na.prop, 
                           between.parameters = between.parameters, seed = seed)
      
      rm(within.parameters, between.parameters)
    }
    
    # Get IDs
    ids <- unique(data$data.long$subjn)
    
    # Create matrix to store the p-value of the KPSS tests.
    statcheck <- matrix(NA, length(ids), I + 1)
    statcheck[, 1] <- ids
    
    for (i in 1:length(ids)){
      tseries_id        <- data$data.long[data$data.long$subjn == ids[i], 3:(I+2)]
      statcheck[i, 2:(I + 1)] <- apply(tseries_id, 2, function(x) kpss.test(na.omit(x), null = "Trend")[[3]])
    }
    rm(i, ids, tseries_id)
    
    # Compute which time series are non-stationary given a p-value lower than 0.05 in the KPSS tests
    statchecksum[r, cond] <- sum(apply(statcheck[, 2:(I+1)] <= 0.05, 1, any))
    rm(statcheck, data, dataModel, na.prop, nT, seed)
  }
}
  
time.total <- proc.time() - time0

write.table(statchecksum, file = paste0(folder, "StationarityCheck/", "statchecksum", ".dat"),
            col.names = TRUE, row.names = FALSE, quote = TRUE)
