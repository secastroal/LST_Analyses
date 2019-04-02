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
missingness <- c(0, 0.1) # Proportion of NAs
dataModel_to_sim <- c("msst", "cuts", "tso") # Model to simulate the data


Cond <- expand.grid(times, missingness, dataModel_to_sim)
names(Cond) <- c("nT", "na.prop", "dataModel")

rm(times, missingness, dataModel_to_sim)

# 2.0 Set up output matrices ----

R <- 50 # Number of replicas

# Create labels to name variables in output matrices 
models <- c("msst", "cuts", "tso")
data_str <- c("wide", "long")
estimator <- c("ml", "bayes")
labels <- expand.grid(data_str, estimator, models)
labels <- labels[-10, c(3,1,2)]
labels <- apply(labels, 1, function(vec) paste(vec, collapse = "+"))
rm(models, data_str, estimator)

# matrix to store check.mplus output (Non-convergence, Warnings and Errors, and Ok)
# parameters, se, variance coefficients, and fit measures are saved if check.mplus output is Ok. 
perf.base <- matrix(NA, R, 11)
colnames(perf.base) <- labels

# matrix to store estimated parameters and standard errors/posterior standard deviations

parameters.base <- matrix(NA, 11 * R + 1,  30)
colnames(parameters.base) <- c(paste0("w_loading", 1:I), "(common)state_var", paste0("unique_state/error_var", 1:I),
                               "ar_effect", paste0("b_loading", 1:I), paste0("intercept", 1:I), "(common)trait_var",
                               "trait_mean", paste0("unique/indicator_trait_var", 1:I), paste0("cov", c(12, 13, 23, 14, 24, 34)))
rownames(parameters.base) <- c("true", paste(rep(labels, times = R), paste0("r", rep(1:R, each = 11)), sep = "+"))

se_psd.base <- parameters.base[-1, ]

# matrix to store variance coefficients
var.coeff.base <- matrix(NA, 11 * R + 1,  I * 5)
colnames(var.coeff.base) <- paste0(rep(c( "ccon/pred_y", "ucon/upred_y", "tcon/con_y", "spe_y", "rel_y"), each = I), 
                              1:I)
rownames(var.coeff.base) <- c("true", paste(rep(labels, times = R), paste0("r", rep(1:R, each = 11)), sep = "+"))

# matrix to store fit measures
fit.measures.base <- matrix(NA, 11 * R,  4)
colnames(fit.measures.base) <- c("AIC", "BIC", "aBIC", "DIC")
rownames(fit.measures.base) <- c(paste(rep(labels, times = R), paste0("r", rep(1:R, each = 11)), sep = "+"))

rm(labels)
# 3.0 Simulation loop ----

for(cond in 1:18){
  print(cond) # print condition number to track progress
  
  # Copy table(s) to save results
  perf <- perf.base
  parameters <- parameters.base
  se_psd <- se_psd.base
  var.coeff <- var.coeff.base
  fit.measures <- fit.measures.base
  
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
      
      # Between Paramaters
      
      intercepts <- seq(0, by = 0.2, length.out = I) # intercepts
      trait_loadings <- c(1, 0.8, 1.2, 0.9) # loading parametes for the latent trait
      
      var_trait <- 2 # variance latent trait variable
      mean_trait <- 4 # mean latent trait variable
      
      between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, trait.mean = mean_trait,
                                 trait.var = var_trait)
      
      # Save True parameters to parameters matrix
      parameters[1, c(1:9, 11:20)] <- c(state_loadings, var_state, var_m_error, trait_loadings,
                                        intercepts, var_trait, mean_trait) 
      
      rm(state_loadings, var_state, var_m_error, intercepts, trait_loadings, var_trait, 
         mean_trait)
      
      # Simulate data
      
      #Time invariant data
      data <- sim.data.msst(N, nT, I, within.parameters = within.parameters, na.prop = na.prop,
                            between.parameters = between.parameters, seed = seed)
      
      # Compute true variance coefficients
      var.coeff[1, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters = within.parameters, 
                                                         between.parameters = between.parameters)) 
      
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
      
      # Save True parameters to parameters matrix
      parameters[1, c(1:9, 11:19, 21:24)] <- c(state_loadings, var_CS, var_US, trait_loadings,
                                        intercepts, var_CT, var_UT) 
      
      rm(state_loadings, var_CS, var_US, intercepts, trait_loadings, var_CT, var_UT)
      
      # Simulate data
      
      # Time invariant data
      data <- sim.data.cuts(N, nT, I, within.parameters = within.parameters, na.prop = na.prop,
                            between.parameters = between.parameters, seed = seed)
      
      # Compute true variance coefficients
      
      var.coeff[1, ] <- t(cuts.var.coeff(within.parameters = within.parameters, 
                                      between.parameters = between.parameters)) 
      
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
      
      # Save True parameters to parameters matrix
      parameters[1, c(1:10, 15:18, 21:30)] <- c(within.parameters$loadings, within.parameters$state.var,
                                                within.parameters$error.var, within.parameters$ar.effect,
                                                between.parameters$intercepts, between.parameters$trait.ind.var,
                                                round(data$between.parameters$Sigma[t(lower.tri(data$between.parameters$Sigma))], 2))
      
      # Compute true variance coefficients
      var.coeff[1, ] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.parameters, 
                                     between.parameters = between.parameters)[,nT] 
      
      rm(within.parameters, between.parameters)
    }
    
    # Fitting the models to the simulated data ----
    
    # msst + wide + maximum likelihod ----
    #give analysis a number
    a <- 1
    
    file.name <- paste("msst", "wide", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "ML",
                                           iterations = 50000)
    
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
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[c(1:I, #loadings
                                                   (I * nT * 3) + ((nT+1) * nT / 2) + 2, # state var 
                                                   ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2), #error variances
                                                   ((I * nT + 1):(I * nT + I)), #trait loadings
                                                   ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1), #intercepts
                                                   ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2), #trait variance
                                                   ((I * nT * 2) + ((nT+1) * nT / 2) + 1)), #trait mean
                                                 3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:20)] <- estimates
      
      # Save standard errors
      stand.errors <- fit$parameters$unstandardized[c(1:I, #loadings
                                                   (I * nT * 3) + ((nT+1) * nT / 2) + 2, # state var 
                                                   ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2), #error variances
                                                   ((I * nT + 1):(I * nT + I)), #trait loadings
                                                   ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1), #intercepts
                                                   ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2), #trait variance
                                                   ((I * nT * 2) + ((nT+1) * nT / 2) + 1)), #trait mean
                                                 4]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:20)] <- stand.errors
      
      # Compute and save variance coefficients
      within.estimates <- list( loadings = estimates[1:I], state.var = estimates[I+1],
                                error.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], 
                                 trait.var = estimates[ 4 * I + 2])
      var.coeff[(11 * (r-1)) + a + 1, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters = within.estimates,
                                                                               between.parameters = between.estimates))
      
      # Save fit measures
      fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                                fit$summaries$aBIC)
      
      rm(within.estimates, between.estimates, estimates, stand.errors)
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit, a)
    
    # msst + long + maximun likelihood ----
    #give analysis a number
    a <- 2
    
    file.name <- paste("msst", "long", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "ML",
                                           iterations = 50000)
    
    ml_syntax <- write.mlmsst.to.Mplus(data$data.long[, -(1:2)])
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[c(1:(2*I+1), #loadings, state var, error variances 
                                                   ((3*I + 2):(4*I + 1)), #loadings 
                                                   ((5*I + 3):(6*I + 2)), #intercepts
                                                   (6*I +3), # trait var
                                                   (4*I + 2)), # trait mean
                                                 3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:20)] <- estimates
      
      # Save standard errors
      stand.errors <- fit$parameters$unstandardized[c(1:(2*I+1), #loadings, state var, error variances 
                                                      ((3*I + 2):(4*I + 1)), #loadings 
                                                      ((5*I + 3):(6*I + 2)), #intercepts
                                                      (6*I +3), # trait var
                                                      (4*I + 2)), # trait mean
                                                    4]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:20)] <- stand.errors
      
      # Compute and save variance coefficients
      within.estimates <- list( loadings = estimates[1:I], state.var = estimates[I+1],
                                error.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], 
                                 trait.var = estimates[ 4 * I + 2])
      var.coeff[(11 * (r-1)) + a + 1, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters = within.estimates,
                                                                               between.parameters = between.estimates))
      
      # Save fit measures
      fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                               fit$summaries$aBIC)
      rm(within.estimates, between.estimates, estimates, stand.errors)
      
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit, a)
    
    
    # msst + wide + Bayesian ----
    #give analysis a number
    a <- 3
    
    file.name <- paste("msst", "wide", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "BAYES",
                                           iterations = 5000)
    
    mplus_syntax <- write.msst.to.Mplus(data$data.wide[,-1], neta = nT, ntheta = 1, 
                                        equiv.assumption = list(tau = "cong", theta = "cong"),
                                        scale.invariance = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                        homocedasticity.assumption = list(error = TRUE, state.red = TRUE),
                                        second.order.trait = FALSE)
    saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                                ";", "\nOUTPUT: TECH8;")
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)
    write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, mplus_syntax, saveoutput_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[c(1:I, #loadings
                                                   (I * nT * 3) + ((nT+1) * nT / 2) + 2, # state var 
                                                   ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2), #error variances
                                                   ((I * nT + 1):(I * nT + I)), #trait loadings
                                                   ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1), #intercepts
                                                   ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2), #trait variance
                                                   ((I * nT * 2) + ((nT+1) * nT / 2) + 1)), #trait mean
                                                 3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:20)] <- estimates
      
      # Save posterior standard deviations
      post.sd <- fit$parameters$unstandardized[c(1:I, #loadings
                                                      (I * nT * 3) + ((nT+1) * nT / 2) + 2, # state var 
                                                      ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT+1) * nT / 2) + nT + I + 2), #error variances
                                                      ((I * nT + 1):(I * nT + I)), #trait loadings
                                                      ((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 2) + ((nT+1) * nT / 2) + I + 1), #intercepts
                                                      ((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2), #trait variance
                                                      ((I * nT * 2) + ((nT+1) * nT / 2) + 1)), #trait mean
                                                    4]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:20)] <- post.sd
      
      # Compute and save variance coefficients
      within.estimates <- list( loadings = estimates[1:I], state.var = estimates[I+1],
                                error.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], 
                                 trait.var = estimates[ 4 * I + 2])
      var.coeff[(11 * (r-1)) + a + 1, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters = within.estimates,
                                                                               between.parameters = between.estimates))
      
      # Save fit measures
      if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
      
      rm(within.estimates, between.estimates, estimates, post.sd)
    
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit, a)
    
    # msst + long + Bayesian ----
    #give analysis a number
    a <- 4
    
    file.name <- paste("msst", "long", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "BAYES",
                                           iterations = 5000)
    
    ml_syntax <- write.mlmsst.to.Mplus(data$data.long[, -(1:2)])
    
    ml_syntax <- gsub("@0;", "@0.001;", ml_syntax) # replace 0 constraints to 0.001
    
    saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                                ";", "\nOUTPUT: TECH8;")
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[c(1:(2*I+1), #loadings, state var, error variances 
                                                   ((3*I + 2):(4*I + 1)), #loadings 
                                                   ((5*I + 3):(6*I + 2)), #intercepts
                                                   (6*I +3), # trait var
                                                   (4*I + 2)), # trait mean
                                                 3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:20)] <- estimates
      
      # Save posterior standard deviations
      post.sd <- fit$parameters$unstandardized[c(1:(2*I+1), #loadings, state var, error variances 
                                                      ((3*I + 2):(4*I + 1)), #loadings 
                                                      ((5*I + 3):(6*I + 2)), #intercepts
                                                      (6*I +3), # trait var
                                                      (4*I + 2)), # trait mean
                                                    4]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:20)] <- post.sd
      
      # Compute and save variance coefficients
      within.estimates <- list( loadings = estimates[1:I], state.var = estimates[I+1],
                                error.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], 
                                 trait.var = estimates[ 4 * I + 2])
      var.coeff[(11 * (r-1)) + a + 1, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters = within.estimates,
                                                                               between.parameters = between.estimates))
      
      # Save fit measures
      if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
      
      rm(within.estimates, between.estimates, estimates, post.sd)
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit, a)
    
    # cuts + wide + maximum likelihod ----
    #give analysis a number
    a <- 5
    
    file.name <- paste("cuts", "wide", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "ML",
                                           iterations = 50000)
    
    
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
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[c(1:I, # CS loadings
                                                   (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1, #CS var
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I), #UT variances
                                                   (I*nT + 1):(I*nT + I), # CT loadings
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I), #intercepts
                                                   (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1, #CT var
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)), # UT var
                                                 3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- estimates
      
      # Save standard errors
      stand.errors <- fit$parameters$unstandardized[c(1:I, # CS loadings
                                                   (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1, #CS var
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I), #UT variances
                                                   (I*nT + 1):(I*nT + I), # CT loadings
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I), #intercepts
                                                   (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1, #CT var
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)), # UT var
                                                 4]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:19, 21:24)] <- stand.errors
      
      # Compute variance coefficients
      within.estimates <- list( loadings = estimates[1:I], CS.var = estimates[I+1],
                                US.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], CT.var = estimates[ 4 * I + 2], 
                                 UT.var = estimates[(4 * I + 3):(5 * I + 2)])
      var.coeff[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                                                between.parameters = between.estimates))
      
      # Save fit measures
      fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                               fit$summaries$aBIC)
      rm(within.estimates, between.estimates)
      
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    
    rm(file.name, fit, a)
    
    # cuts + long + maximun likelihood ----
    #give analysis a number
    a <- 6
    
    file.name <- paste("cuts", "long", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "ML",
                                           iterations = 50000)
    
    ml_syntax <- write.mlcuts.to.Mplus(cuts.data$data.long[, -(1:2)])
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[,3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- estimates
      
      # Save standard errors
      stand.errors <- fit$parameters$unstandardized[,3]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:19, 21:24)] <- stand.errors
      
      # Compute variance coefficients
      within.estimates <- list( loadings = estimates[1:I], CS.var = estimates[I+1],
                                US.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], CT.var = estimates[ 4 * I + 2], 
                                 UT.var = estimates[(4 * I + 3):(5 * I + 2)])
      var.coeff[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                                                between.parameters = between.estimates))
      
      # Save fit measures
      fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                               fit$summaries$aBIC)
      rm(within.estimates, between.estimates)
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    
    rm(file.name, fit, a)
    
    # cuts + wide + Bayesian ----
    #give analysis a number
    a <- 7
    
    file.name <- paste("cuts", "wide", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "BAYES",
                                           iterations = 5000)
    
    
    mplus_syntax <- write.cuts.to.Mplus(data$data.wide[,-1],  nstate = nT,
                                        method.trait = "om",
                                        scale.invariance = list(int = TRUE, lambda = TRUE),
                                        state.trait.invariance = FALSE,
                                        fixed.method.loadings = TRUE,
                                        homocedasticity.assumption = list(error = TRUE, cs.red = TRUE, ut.red = FALSE))
    saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                                ";", "\nOUTPUT: TECH8;")
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)
    write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, mplus_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[c(1:I, # CS loadings
                                                   (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1, #CS var
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I), #UT variances
                                                   (I*nT + 1):(I*nT + I), # CT loadings
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I), #intercepts
                                                   (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1, #CT var
                                                   ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)), # UT var
                                                 3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- estimates
      
      # Save posterior standard deviations
      post.sd <- fit$parameters$unstandardized[c(1:I, # CS loadings
                                                 (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1, #CS var
                                                 ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + 2*I), #UT variances
                                                 (I*nT + 1):(I*nT + I), # CT loadings
                                                 ((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + I), #intercepts
                                                 (((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1, #CT var
                                                 ((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)), # UT var
                                               4]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:19, 21:24)] <- post.sd
      
      # Compute variance coefficients
      within.estimates <- list( loadings = estimates[1:I], CS.var = estimates[I+1],
                                US.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], CT.var = estimates[ 4 * I + 2], 
                                 UT.var = estimates[(4 * I + 3):(5 * I + 2)])
      var.coeff[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                                                between.parameters = between.estimates))
      
      # Save fit measures
      if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
      
      rm(within.estimates, between.estimates, estimates, post.sd)
      
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    
    rm(file.name, fit, a)
    
    # cuts + long + Bayesian ----
    #give analysis a number
    a <- 8
    
    file.name <- paste("cuts", "long", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "BAYES",
                                           iterations = 5000)
    
    ml_syntax <- write.mlcuts.to.Mplus(data$data.long[, -(1:2)])
    saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                                ";", "\nOUTPUT: TECH8;")
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n") 
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      # Save estimates
      estimates <- fit$parameters$unstandardized[,3]
      
      parameters[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- estimates
      
      # Save posterior standard deviations
      post.sd <- fit$parameters$unstandardized[,3]
      
      se_psd[(11 * (r-1)) + a, c(1:9, 11:19, 21:24)] <- post.sd
      
      # Compute variance coefficients
      within.estimates <- list( loadings = estimates[1:I], CS.var = estimates[I+1],
                                US.var = estimates[(I+2):(2 * I +1)])
      between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], CT.var = estimates[ 4 * I + 2], 
                                 UT.var = estimates[(4 * I + 3):(5 * I + 2)])
      var.coeff[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                                                between.parameters = between.estimates))
      
      # Save fit measures
      if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
      
      rm(within.estimates, between.estimates, estimates, post.sd)
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    
    rm(file.name, fit, a)
    
    # tso + wide + maximum likelihood ----
    #give analysis a number
    a <- 9
    
    file.name <- paste("tso", "wide", "ml", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "ML",
                                           iterations = 50000)
    
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
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
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
      
      # Save fit measures
      fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                               fit$summaries$aBIC)
      
      rm(within.estimates, between.estimates)
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit, a)
    
    # tso + wide + Bayesian ----
    #give analysis a number
    a <- 10
    
    file.name <- paste("tso", "wide", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.wide)[-1],
                                           analysis_type = "GENERAL",
                                           estimator = "BAYES",
                                           iterations = 5000)
    
    ml_syntax <- write.tso.to.Mplus(data$data.wide[,-1], nocc = nT, figure = "3b",
                                    equiv.assumption = list(occ = "cong", theta = "equi"),
                                    scale.invariance = list(int = TRUE, lambda = TRUE),
                                    homocedasticity.assumption = list(error = TRUE, occ.red = TRUE),
                                    autoregressive.homogeneity = TRUE)
    saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                                ";", "\nOUTPUT: TECH8;")
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
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
      # Save fit measures
      if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
      
      rm(within.estimates, between.estimates)
      
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit, a)
    
    # tso + long + Bayesian ----
    #give analysis a number
    a <- 11
    
    file.name <- paste("tso", "long", "bayes", cond, r,sep = "_")
    
    # Prepare data: Write data in Mplus format and write input file template
    prepareMplusData(data$data.long, paste0(folder,file.name,".dat"), inpfile = T)
    
    # Complete Mplus syntax
    analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                           cluster = names(data$data.long)[1],
                                           analysis_type = "TWOLEVEL",
                                           estimator = "BAYES",
                                           iterations = 5000)
    
    ml_syntax <- write.mltso.to.Mplus(data$data.long[, -c(1, 2)])
    saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", paste0("samples_", file.name, ".dat"),
                                ";", "\nOUTPUT: TECH8;")
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    write(saveoutput_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n")
    runModels_2(paste0(getwd(),"/",folder,file.name,".inp"), timeout = timeout) 
    cat("\n"); print(Sys.time()); cat("\n")
    
    fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
    
    if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
      
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
      
      estimates <- fit$parameters$unstandardized[c(1:(2*I + 2),
                                                   (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2),
                                                   (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)), 3]
      
      within.estimates <- list( loadings = estimates[1:I], state.var = estimates[2*I+2],
                                error.var = estimates[(I+2):(2 * I +1)], ar.effect = estimates[I+1])
      between.estimates <- list( trait.ind.var = estimates[(3 * I + 3):(4 * I + 2)])
      #var.coeff[ , 5] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
       #                                between.parameters = between.estimates)[,nT]
      #fit.measures[3, 2:4] <- fit$summaries[, c(11:13)]
      # Save fit measures
      if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
      
      rm(within.estimates, between.estimates)
    }else{
      perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
    }
    
    rm(file.name, fit, a)
    
    
    
    
  }
}


#read all analysis to check
#fit <- readModels(paste0(getwd(),"/","Mplus_files_Results/","cuts_long_n100_i4_nt5_na0",".out"))
