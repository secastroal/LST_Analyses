# This script is to perform the simulation study for the paper Using Latent State-Trait Theory 
#   to Analyze Intensive Longitudinal Data.

# CONTENTS
# 0.0 Prepare environment
# 1.0 Set up conditions
# 2.0 Set up output matrices
# 3.0 Simulation loop

# 0.0 Prepare environment ----
rm(list=ls())
#library(devtools) # if lsttheory of sumplement C of Steyer et al. (2015) has not been installed before.
#install_github("amayer2010/lsttheory", force = TRUE)
#library(lavaan) 
#library(lsttheory)
library(MplusAutomation)
library(MASS)
library(xtable)
# Easy code to source all the needed files
file.sources <- paste0("R/", list.files(paste0(getwd(), "/R")))
sapply(file.sources,source,.GlobalEnv)
rm(file.sources)
folder <- "Mplus_Simulation/" #Folder to store results and all Mplus files

# 1.0 Set up conditions ----

# Fix conditions
N <- 100 # Sample size
I <- 4 # Number of variables
timeout <- 3600 # Time limit in seconds to force ending an analysis in Mplus

# Manipulated conditions
times <- c(30, 60, 90) # Number of measurement occasions
missingess <- c(0, 0.1) # Proportion of NAs
dataModel_to_sim <- c("msst", "cuts", "tso") # Model to simulate the data


Cond <- expand.grid(times, missingess, dataModel_to_sim)
names(Cond) <- c("nT", "na.prop", "dataModel")

rm(times, missingness, dataModel_to_sim)

# 2.0 Set up output matrices ----

R <- 50 # Number of replicas

# What do we want to save and how?
# Do we want to save all the parameters?
# Do we want to only save the variance coefficients?
# Are we going to compute bias / absolute bias / RMSE per parameter or per set of parameters (e.g., loadings 
# at the between level)
# Are we saving the standard errors and/or posterior standard deviations?
# what fit measures of what analyses are we saving?


# 3.0 Simulation loop ----

for(cond in 1:18){
  print(cond) # print condition number to track progress
  # Copy table(s) to save results
  
  for(r in 1:R){
    print(r) #print replication number to track progress
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
      
      rm(state_loadings, var_state, var_m_error)
      
      # Between Paramaters
      
      intercepts <- seq(0, by = 0.2, length.out = I) # intercepts
      trait_loadings <- c(1, 0.8, 1.2, 0.9) # loading parametes for the latent trait
      
      var_trait <- 2 # variance latent trait variable
      mean_trait <- 4 # mean latent trait variable
      
      between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, trait.mean = mean_trait,
                                 trait.var = var_trait)
      
      rm(intercepts, trait_loadings, var_trait, mean_trait)
      
      # Simulate data
      
      #Time invariant data
      data <- sim.data.msst(N, nT, I, within.parameters = within.parameters, na.prop = na.prop,
                            between.parameters = between.parameters, seed = seed)
      # Compute true variance coefficients
      
      #var.coeff[(2 * I + 1):(5 * I),2] <- msst.var.coeff(within.parameters = within.parameters, 
       #                                                  between.parameters = between.parameters) 
      
      rm(within.parameters, between.parameters)
    }
    
    if(dataModel == "cuts"){
      
      # Within Parameters
      
      state_loadings <- c(1, 0.5, 1.3, 0.8) # loading parameters for the latent common state
      var_CS <- 2 # Variance latent common state
      var_US <- c(1, 0.5, 1.5, 0.8) # Variance of latent unique states
      
      within.parameters <- list(loadings = state_loadings, CS.var = var_CS, US.var = var_US)
      
      rm(state_loadings, var_CS, var_US)
      
      # Between Paramaters
      
      intercepts <- seq(2, by = 0.5, length.out = I) # intercepts
      trait_loadings <- c(1, 0.8, 1.2, 0.9) # loading parametes for the latent common trait
      var_CT <- 1.5 # variance latent common trait
      var_UT <- c(0.5, 1, 0.3, 0.8) # variance latent unique traits
      
      between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, CT.var = var_CT, UT.var = var_UT)
      
      rm(intercepts, trait_loadings, var_CT, var_UT)
      
      # Simulate data
      
      # Time invariant data
      data <- sim.data.cuts(N, nT, I, within.parameters = within.parameters, na.prop = na.prop,
                            between.parameters = between.parameters, seed = seed)
      
      # Compute true variance coefficients
      
      #var.coeff[,2] <- cuts.var.coeff(within.parameters = within.parameters, 
       #                               between.parameters = between.parameters) 
      
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
      repeat {
        R <- matrix(sample((7:9)/10, size = I * I, replace = TRUE), I) #correlation matrix trait indicators
        R[lower.tri(R)] = t(R)[lower.tri(R)]
        diag(R) <- 1
        print(det(R))
        if (det(R) > 0){
          break
        }
      }
      
      between.parameters <- list(intercepts = intercepts, trait.ind.var = var_ind_traits, 
                                 cor.matrix = R)
      
      rm(intercepts, var_ind_traits, R)
      
      # Simulate data
      
      # Time invariant data
      data <- sim.data.tso(N, nT, I, within.parameters = within.parameters, na.prop = na.prop, 
                           between.parameters = between.parameters, seed = seed)
      # Compute true variance coefficients
      
      #var.coeff[,2] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.parameters, 
       #                              between.parameters = between.parameters)[,nT] 
      
      rm(within.parameters, between.parameters)
    }
    
    # Fitting the models to the simulated data ----
    
    # msst + wide + maximum likelihod ----
    
    file.name <- paste0("msst", "wide", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "ML",
                                           iterations = 50000,
                                           h1iterations = 50000)
    
    mplus_syntax <- write.msst.to.Mplus(data$data.wide[,-1], neta = nT, ntheta = 1, 
                                        equiv.assumption = list(tau = "cong", theta = "cong"),
                                        scale.invariance = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                        homocedasticity.assumption = list(error = TRUE, state.red = TRUE),
                                        second.order.trait = FALSE)
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, mplus_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/", folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      #status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      est.par[,4] <- wide.fit$parameters$unstandardized[c(1:I,
                                                          (I * nT * 3) + ((nT+1) * nT / 2) + 2,
                                                          ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2),
                                                          ((I * nT + 1):(I * nT + I)), 
                                                          ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1),
                                                          ((I * nT * 2) + ((nT+1) * nT / 2) + 1), 
                                                          ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2)),3]
      
      #se[,3] <- wide.fit$parameters$unstandardized[c(1:I,
       #                                              (I * nT * 3) + ((nT+1) * nT / 2) + 2,
        #                                             ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2),
         #                                            ((I * nT + 1):(I * nT + I)), 
          #                                           ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1),
           #                                          ((I * nT * 2) + ((nT+1) * nT / 2) + 1), 
            #                                         ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2)),4]
      
      #bias[, 3] <- (est.par[ , 4] - est.par[ , 2]) 
      
      #rmse[, 3] <- sqrt((est.par[ , 4] - est.par[ , 2]) ^ 2) 
      
      within.estimates <- list( loadings = est.par[1:I, 4], state.var = est.par[I+1, 4],
                                error.var = est.par[(I+2):(2 * I +1), 4])
      between.estimates <- list( loadings = est.par[(2 * I + 2):(3 * I + 1), 4], 
                                 trait.var = est.par[ 4 * I + 3, 4])
      #var.coeff[ , 4] <- msst.var.coeff(within.parameters = within.estimates,
       #                                 between.parameters = between.estimates)
      rm(within.estimates, between.estimates)
      
      
    }else{
      #status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit)
    
    # msst + long + maximun likelihood ----
    
    file.name <- paste0("msst", "long", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "ML",
                                           iterations = 50000,
                                           h1iterations = 50000)
    
    ml_syntax <- write.mlmsst.to.Mplus(data$data.long[, -(1:2)])
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/", folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      status[1] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      est.par[,3] <- fit$parameters$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                                          (4*I + 2), (6*I +3)),3]
      #se[,2] <- fit$parameters$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
       #                                              (4*I + 2), (6*I +3)),4]
      
      #bias[, 2] <- (est.par[ , 3] - est.par[ , 2]) 
      
      #rmse[, 2] <- sqrt((est.par[ , 3] - est.par[ , 2]) ^ 2)
      
      within.estimates <- list( loadings = est.par[1:I, 3], state.var = est.par[I+1, 3],
                                error.var = est.par[(I+2):(2 * I +1), 3])
      between.estimates <- list( loadings = est.par[(2 * I + 2):(3 * I + 1), 3], 
                                 trait.var = est.par[ 4 * I + 3, 3])
      #var.coeff[ , 3] <- msst.var.coeff(within.parameters = within.estimates,
       #                                 between.parameters = between.estimates)
      rm(within.estimates, between.estimates)
      
      
    }else{
      #status[1] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit)
    
    
    # msst + wide + Bayesian ----
    
    file.name <- paste0("msst", "wide", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "BAYES",
                                           iterations = 5000,
                                           processors = 2)
    
    mplus_syntax <- write.msst.to.Mplus(data$data.wide[,-1], neta = nT, ntheta = 1, 
                                        equiv.assumption = list(tau = "cong", theta = "cong"),
                                        scale.invariance = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                        homocedasticity.assumption = list(error = TRUE, state.red = TRUE),
                                        second.order.trait = FALSE)
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, mplus_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/", folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      est.par[,4] <- fit$parameters$unstandardized[c(1:I,
                                                          (I * nT * 3) + ((nT+1) * nT / 2) + 2,
                                                          ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2),
                                                          ((I * nT + 1):(I * nT + I)), 
                                                          ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1),
                                                          ((I * nT * 2) + ((nT+1) * nT / 2) + 1), 
                                                          ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2)),3]
      
      #se[,3] <- fit$parameters$unstandardized[c(1:I,
       #                                              (I * nT * 3) + ((nT+1) * nT / 2) + 2,
        #                                             ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2),
         #                                            ((I * nT + 1):(I * nT + I)), 
          #                                           ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1),
           #                                          ((I * nT * 2) + ((nT+1) * nT / 2) + 1), 
            #                                         ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2)),4]
      
      #bias[, 3] <- (est.par[ , 4] - est.par[ , 2]) 
      
      #rmse[, 3] <- sqrt((est.par[ , 4] - est.par[ , 2]) ^ 2) 
      
      within.estimates <- list( loadings = est.par[1:I, 4], state.var = est.par[I+1, 4],
                                error.var = est.par[(I+2):(2 * I +1), 4])
      between.estimates <- list( loadings = est.par[(2 * I + 2):(3 * I + 1), 4], 
                                 trait.var = est.par[ 4 * I + 3, 4])
      #var.coeff[ , 4] <- msst.var.coeff(within.parameters = within.estimates,
       #                                 between.parameters = between.estimates)
      rm(within.estimates, between.estimates)
      
      
    }else{
      #status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit)
    
    # msst + long + Bayesian ----
    
    file.name <- paste0("msst", "long", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "BAYES",
                                           iterations = 5000,
                                           processors = 2)
    
    ml_syntax <- write.mlmsst.to.Mplus(data$data.long[, -(1:2)])
    
    ml_syntax <- gsub("@0;", "@0.001;", ml_syntax) # replace 0 constraints to 0.001
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/", folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      estimates <- fit$parameters$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                                   (4*I + 2), (6*I +3)),3]
      
      within.estimates <- list( loadings = estimates[1:I], state.var = estimates[I+1],
                                error.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], 
                                 trait.var = estimates[ 4 * I + 3])
      #var.coeff[(2 * I + 1):(5 * I), 4] <- msst.var.coeff(within.parameters = within.estimates,
       #                                                   between.parameters = between.estimates)
      #fit.measures[2, 2:4] <- fit$summaries[, c(11:13)]
      rm(within.estimates, between.estimates)
    }else{
      #stop("Model estimation did not converge or there are warning or error messages in the output")
    }
    
    rm(file.name, fit, estimates)
    
    # cuts + wide + maximum likelihod ----
    
    file.name <- paste0("cuts", "wide", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "ML",
                                           iterations = 50000,
                                           h1iterations = 50000)
    
    
    mplus_syntax <- write.cuts.to.Mplus(data$data.wide[,-1],  nstate = nT,
                                        method.trait = "om",
                                        scale.invariance = list(int = TRUE, lambda = TRUE),
                                        state.trait.invariance = FALSE,
                                        fixed.method.loadings = TRUE,
                                        homocedasticity.assumption = list(error = TRUE, cs.red = TRUE, ut.red = FALSE))
    
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, mplus_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/",folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      est.par[, 4] <- fit$parameters$unstandardized[c(1:I,
                                                           (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1,
                                                           ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I),
                                                           (I*nT + 1):(I*nT + I),
                                                           ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I),
                                                           (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1,
                                                           ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)),
                                                         3]
      
      #se[, 3] <- fit$parameters$unstandardized[c(1:I,
       #                                               (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1,
        #                                              ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I),
         #                                             (I*nT + 1):(I*nT + I),
          #                                            ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I),
           #                                           (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1,
            #                                          ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)),
              #                                      4]
      #bias[, 3] <- (est.par[ , 4] - est.par[ , 2]) 
      
      #rmse[, 3] <- sqrt((est.par[ , 4] - est.par[ , 2]) ^ 2)
      
      within.estimates <- list( loadings = est.par[1:I, 4], CS.var = est.par[I+1, 4],
                                US.var = est.par[(I+2):(2 * I +1), 4])
      between.estimates <- list( loadings = est.par[(2 * I + 2):(3 * I + 1), 4], CT.var = est.par[ 4 * I + 2, 4], 
                                 UT.var = est.par[(4 * I + 3):(5 * I + 2), 4])
      #var.coeff[ , 4] <- cuts.var.coeff(within.parameters = within.estimates,
       #                                 between.parameters = between.estimates)
      rm(within.estimates, between.estimates)
      
      
    }else{
      #status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    
    rm(file.name, fit)
    
    # cuts + long + maximun likelihood ----
    
    file.name <- paste0("cuts", "long", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "ML",
                                           iterations = 50000,
                                           h1iterations = 50000)
    
    ml_syntax <- write.mlcuts.to.Mplus(cuts.data$data.long[, -(1:2)])
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/",folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      status[1] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      est.par[, 3] <- fit$parameters$unstandardized[,3]
      
      #se[, 2] <- fit$parameters$unstandardized[,4]
      
      #bias[, 2] <- (est.par[ , 3] - est.par[ , 2]) 
      
      #rmse[, 2] <- sqrt((est.par[ , 3] - est.par[ , 2]) ^ 2)
      
      within.estimates <- list( loadings = est.par[1:I, 3], CS.var = est.par[I+1, 3],
                                US.var = est.par[(I+2):(2 * I +1), 3])
      between.estimates <- list( loadings = est.par[(2 * I + 2):(3 * I + 1), 3], CT.var = est.par[ 4 * I + 2, 3], 
                                 UT.var = est.par[(4 * I + 3):(5 * I + 2), 3])
      #var.coeff[ , 3] <- cuts.var.coeff(within.parameters = within.estimates,
       #                                 between.parameters = between.estimates)
      rm(within.estimates, between.estimates)
      
    }else{
      #status[1] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    
    rm(file.name, fit)
    
    
    # cuts + wide + Bayesian ----
    
    file.name <- paste0("cuts", "wide", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "BAYES",
                                           iterations = 5000,
                                           processors = 2)
    
    
    mplus_syntax <- write.cuts.to.Mplus(data$data.wide[,-1],  nstate = nT,
                                        method.trait = "om",
                                        scale.invariance = list(int = TRUE, lambda = TRUE),
                                        state.trait.invariance = FALSE,
                                        fixed.method.loadings = TRUE,
                                        homocedasticity.assumption = list(error = TRUE, cs.red = TRUE, ut.red = FALSE))
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, mplus_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/",folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      est.par[, 4] <- fit$parameters$unstandardized[c(1:I,
                                                           (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1,
                                                           ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I),
                                                           (I*nT + 1):(I*nT + I),
                                                           ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I),
                                                           (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1,
                                                           ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)),
                                                         3]
      
      #se[, 3] <- fit$parameters$unstandardized[c(1:I,
       #                                               (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1,
        #                                              ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I),
         #                                             (I*nT + 1):(I*nT + I),
          #                                            ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I),
           #                                           (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1,
            #                                          ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)),
             #                                       4]
      #bias[, 3] <- (est.par[ , 4] - est.par[ , 2]) 
      
      #rmse[, 3] <- sqrt((est.par[ , 4] - est.par[ , 2]) ^ 2)
      
      within.estimates <- list( loadings = est.par[1:I, 4], CS.var = est.par[I+1, 4],
                                US.var = est.par[(I+2):(2 * I +1), 4])
      between.estimates <- list( loadings = est.par[(2 * I + 2):(3 * I + 1), 4], CT.var = est.par[ 4 * I + 2, 4], 
                                 UT.var = est.par[(4 * I + 3):(5 * I + 2), 4])
      #var.coeff[ , 4] <- cuts.var.coeff(within.parameters = within.estimates,
       #                                 between.parameters = between.estimates)
      rm(within.estimates, between.estimates)
      
      
    }else{
      #status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    
    rm(file.name, fit)
    
    # cuts + long + Bayesian ----
    
    file.name <- paste0("cuts", "long", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "BAYES",
                                           iterations = 5000,
                                           processors = 2)
    
    ml_syntax <- write.mlcuts.to.Mplus(data$data.long[, -(1:2)])
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n") 
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      estimates <- fit$parameters$unstandardized[,3]
      
      within.estimates <- list( loadings = estimates[1:I], CS.var = estimates[I+1],
                                US.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], CT.var = estimates[ 4 * I + 2], 
                                 UT.var = estimates[(4 * I + 3):(5 * I + 2)])
      #var.coeff[ , 3] <- cuts.var.coeff(within.parameters = within.estimates,
       #                                 between.parameters = between.estimates)
      #fit.measures[1, 2:4] <- fit$summaries[, c(11:13)]
      rm(within.estimates, between.estimates)
      
    }else{
      # stop("Model estimation did not converge or there are warning or error messages in the output")
    }
    
    
    rm(file.name, fit, estimates)
    
    
    # tso + wide + Bayesian ----
    
    file.name <- paste0("tso", "wide", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "BAYES",
                                           iterations = 5000,
                                           processors = 2)
    
    ml_syntax <- write.tso.to.Mplus(data$data.wide[,-1], nocc = nT, figure = "3b",
                                    equiv.assumption = list(occ = "cong", theta = "equi"),
                                    scale.invariance = list(int = TRUE, lambda = TRUE),
                                    homocedasticity.assumption = list(error = TRUE, occ.red = TRUE),
                                    autoregressive.homogeneity = TRUE)
    
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/",folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      est.par[, 4] <- fit$parameters$unstandardized[c(1:I, #loadings
                                                           (2 * I * nT + 1), # autoregressive effect
                                                           ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                           ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)), # occasion variance
                                                           ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I), # Intercepts
                                                           ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I), # trait indicator variances
                                                           ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)), # trait indicator covariances 
                                                         3]
      
      #post.sd[, 3] <- fit$parameters$unstandardized[c(1:I, #loadings
       #                                                    (2 * I * nT + 1), # autoregressive effect
        #                                                   ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
         #                                                  ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)), # occasion variance
          #                                                 ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I), # Intercepts
           #                                                ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I), # trait indicator variances
            #                                               ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)), # trait indicator covariances 
             #                                            4]
      
      #bias[, 3] <- (est.par[ , 4] - est.par[ , 2])
      
      #rmse[, 3] <- sqrt((est.par[ , 4] - est.par[ , 2]) ^ 2) 
      
      within.estimates <- list( loadings = est.par[1:I, 4], state.var = est.par[2*I+2, 4],
                                error.var = est.par[(I+2):(2 * I +1), 4], ar.effect = est.par[I+1,4])
      between.estimates <- list( trait.ind.var = est.par[(3 * I + 3):(4 * I + 2), 4])
      #var.coeff[ , 4] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
       #                                between.parameters = between.estimates)[,nT]
      rm(within.estimates, between.estimates)
      
    }else{
      #status[2] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit)
    
    # tso + long + Bayesian ----
    
    file.name <- paste0("tso", "long", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "BAYES",
                                           iterations = 5000,
                                           processors = 2)
    
    ml_syntax <- write.mltso.to.Mplus(data$data.long[, -c(1, 2)])
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels(paste0(getwd(),"/",folder,file.name,".inp"))
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      estimates <- fit$parameters$unstandardized[c(1:(2*I + 2),
                                                   (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2),
                                                   (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)), 3]
      
      within.estimates <- list( loadings = estimates[1:I], state.var = estimates[2*I+2],
                                error.var = estimates[(I+2):(2 * I +1)], ar.effect = estimates[I+1])
      between.estimates <- list( trait.ind.var = estimates[(3 * I + 3):(4 * I + 2)])
      #var.coeff[ , 5] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
       #                                between.parameters = between.estimates)[,nT]
      #fit.measures[3, 2:4] <- fit$summaries[, c(11:13)]
      rm(within.estimates, between.estimates)
    }else{
      #stop("Model estimation did not converge or there are warning or error messages in the output")
    }
    
    rm(file.name, fit, estimates)
    
    
    
    
  }
}