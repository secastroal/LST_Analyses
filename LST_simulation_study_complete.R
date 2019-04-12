# This script is to perform the simulation study for the paper Using Latent State-Trait Theory 
#   to Analyze Intensive Longitudinal Data.

# CONTENTS
# 0.0 Prepare environment
# 1.0 Set up conditions
# 2.0 Set up output matrices
# 3.0 Simulation loop
# 4.0 Save output

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
folder <- "Mplus_Simulation" #Folder to store results and all Mplus files

# Create folders to store output in the working directory ----

if(!dir.exists(paste(getwd(), folder, sep = "/"))){
  dir.create(paste(getwd(), folder, "Performance", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Times", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Parameters", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "SE_PSD", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Var_Coeff", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Fit_Measures", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "PSD_Var_Coeff", sep = "/"), recursive = TRUE)
}else{
  if(!dir.exists(paste(getwd(), folder, "Performance", sep = "/"))){
    dir.create(paste(getwd(), folder, "Performance", sep = "/"))
  }
  if(!dir.exists(paste(getwd(), folder, "Times", sep = "/"))){
    dir.create(paste(getwd(), folder, "Times", sep = "/"))
  }
  if(!dir.exists(paste(getwd(), folder, "Parameters", sep = "/"))){
    dir.create(paste(getwd(), folder, "Parameters", sep = "/"))
  }
  if(!dir.exists(paste(getwd(), folder, "SE_PSD", sep = "/"))){
    dir.create(paste(getwd(), folder, "SE_PSD", sep = "/"))
  }
  if(!dir.exists(paste(getwd(), folder, "Var_Coeff", sep = "/"))){
    dir.create(paste(getwd(), folder, "Var_Coeff", sep = "/"))
  }
  if(!dir.exists(paste(getwd(), folder, "Fit_Measures", sep = "/"))){
    dir.create(paste(getwd(), folder, "Fit_Measures", sep = "/"))
  }
  if(!dir.exists(paste(getwd(), folder, "PSD_Var_Coeff", sep = "/"))){
    dir.create(paste(getwd(), folder, "PSD_Var_Coeff", sep = "/"))
  }
}

# Remove files in output folder
unlink(paste(getwd(), folder, "*", sep = "/"))

# add / to folder
folder <- paste0(folder, "/")

# 1.0 Set up conditions ----

# Fix conditions
N <- 100 # Sample size
I <- 4 # Number of variables
timeout <- 120 # Time limit in seconds to force ending an analysis in Mplus

# Manipulated conditions
timepoints <- c(30, 60, 90) # Number of measurement occasions
missingness <- c(0, 0.1) # Proportion of NAs
dataModel_to_sim <- c("msst", "cuts", "tso") # Model to simulate the data

Cond <- expand.grid(timepoints, missingness, dataModel_to_sim)
names(Cond) <- c("nT", "na.prop", "dataModel")

rm(timepoints, missingness, dataModel_to_sim)

# 2.0 Set up output matrices ----

R <- 1 # Number of replicas

# Create labels to name variables in output matrices 
models <- c("msst", "cuts", "tso")
data_str <- c("wide", "long")
estimator <- c("ml", "bayes")
labels <- expand.grid(data_str, estimator, models)
labels <- labels[-10, c(3,1,2)]
labels <- apply(labels, 1, function(vec) paste(vec, collapse = "_"))
rm(models, data_str, estimator)

# matrix to store check.mplus output (Non-convergence, Warnings and Errors, and Ok)
# parameters, se, variance coefficients, and fit measures are saved if check.mplus output is Ok. 
perf.base <- matrix(NA, R, 11)
colnames(perf.base) <- labels
times.base <- perf.base

# matrix to store estimated parameters and standard errors/posterior standard deviations

parameters.base <- matrix(NA, 11 * R + 1,  30)
colnames(parameters.base) <- c(paste0("w_loading", 1:I), "common_state_var", paste0("unique_state_error_var", 1:I),
                               "ar_effect", paste0("b_loading", 1:I), paste0("intercept", 1:I), "common_trait_var",
                               "trait_mean", paste0("unique_indicator_trait_var", 1:I), paste0("cov", c(12, 13, 23, 14, 24, 34)))
rownames(parameters.base) <- c("true", paste(rep(labels, times = R), paste0("r", rep(1:R, each = 11)), sep = "_"))

se_psd.base <- parameters.base[-1, ]

# matrix to store variance coefficients
var.coeff.base <- matrix(NA, 11 * R + 1,  I * 5)
colnames(var.coeff.base) <- paste0(rep(c( "ccon_pred_y", "ucon_upred_y", "tcon_con_y", "spe_y", "rel_y"), each = I), 
                              1:I)
rownames(var.coeff.base) <- c("true", paste(rep(labels, times = R), paste0("r", rep(1:R, each = 11)), sep = "_"))

# matrix to store posterior standard deviations of variance coefficients
psd.var.coeff.base <- matrix(NA, 6 * R,  I * 5)
colnames(psd.var.coeff.base) <- paste0(rep(c( "ccon_pred_y", "ucon_upred_y", "tcon_con_y", "spe_y", "rel_y"), each = I), 
                                   1:I)
rownames(psd.var.coeff.base) <- c(paste(rep(labels[grep(pattern = "bayes", labels)], times = R),
                                        paste0("r", rep(1:R, each = 6)), sep = "_"))

# matrix to store fit measures
fit.measures.base <- matrix(NA, 11 * R,  4)
colnames(fit.measures.base) <- c("AIC", "BIC", "aBIC", "DIC")
rownames(fit.measures.base) <- c(paste(rep(labels, times = R), paste0("r", rep(1:R, each = 11)), sep = "_"))

rm(labels)

# 3.0 Simulation loop ----

time0 <- proc.time()
for(cond in c(9, 11, 12)){
  print(paste("Condition", cond)) # print condition number to track progress
  
  # Copy table(s) to save results
  perf <- perf.base
  times <- times.base
  parameters <- parameters.base
  se_psd <- se_psd.base
  var.coeff <- var.coeff.base
  fit.measures <- fit.measures.base
  psd.var.coeff <- psd.var.coeff.base
  
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
      
      # Save True parameters to parameters matrix
      parameters[1, c(1:9, 11:20)] <- c(state_loadings, var_state, var_m_error, trait_loadings,
                                        intercepts, var_trait, mean_trait) 
      
      rm(state_loadings, var_state, var_m_error, intercepts, trait_loadings, var_trait, 
         mean_trait)
      
      # Simulate data
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
    
    if(cond %in% c(5, 6, 11, 12, 17, 18)){
      perf[r, a] <- "notrun"
      rm(a)
    }else{
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
      cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
      runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
      cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
      
      if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
        system("taskkill /im Mplus.exe /f")
        
        perf[r, a] <- "timeout"
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        
        rm(file.name, a, t0, tf)
      }else{
        fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
        
        if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
          
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          
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
          times[r, a] <- round(tf[3], 2)
          parameters[(11 * (r-1)) + a + 1, ] <- -999
          se_psd[(11 * (r-1)) + a, ] <- -999 
          var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
          fit.measures[(11 * (r-1)) + a, ] <- -999 
        }
        
        rm(file.name, fit, a, t0, tf)
      }
    }
    
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
    cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
    
    if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
      system("taskkill /im Mplus.exe /f")
      
      perf[r, a] <- "timeout"
      times[r, a] <- round(tf[3], 2)
      parameters[(11 * (r-1)) + a + 1, ] <- -999
      se_psd[(11 * (r-1)) + a, ] <- -999 
      var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
      fit.measures[(11 * (r-1)) + a, ] <- -999
      
      rm(file.name, a, t0, tf)
    }else{
      fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
      
      if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
        
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        
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
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999 
      }
      
      rm(file.name, fit, a, t0, tf)
    }
    
    # msst + wide + Bayesian ----
    #give analysis a number
    a <- 3
    
    if(cond %in% c(3, 6, 9, 12, 15, 18)){
      perf[r, a] <- "notrun"
      rm(a)
    }else{
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
      cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
      runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
      cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
      
      if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
        system("taskkill /im Mplus.exe /f")
        
        perf[r, a] <- "timeout"
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
        
        rm(file.name, a, t0, tf)
      }else{
        fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
        
        if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
          
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          
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
          
          # Read MCMC draws
          samples <- read.table(paste0(getwd(), "/", folder, "samples_", file.name, ".dat"))
          
          # Delete burn-in and indicator for chain and draw
          samples <- samples[which(samples$V2 >= (max(samples$V2)/2 + 1)), -(1:2)]
          
          # Create matrix to store variance coefficients draws
          pdist.var.coeff <- matrix(NA, dim(samples)[1], 12)
          
          for(i in 1:dim(samples)[1]){
            within.estimates <- list( loadings = c(1, t(samples[i, c(4, 6, 8)])), state.var = samples[i, 15],
                                      error.var = c(t(samples[i, 10:13])))
            between.estimates <- list( loadings = c(1, t(samples[i, c(5, 7, 9)])), 
                                       trait.var = samples[i, 16])
            pdist.var.coeff[i, ] <- t(msst.var.coeff(within.parameters = within.estimates,
                                                     between.parameters = between.estimates))
          }
          
          var.coeff[(11 * (r-1)) + a + 1, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, median)
          psd.var.coeff[(6 * (r-1)) + a - 2, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, sd)
          
          # Save fit measures
          if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
          
          rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
          
        }else{
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          parameters[(11 * (r-1)) + a + 1, ] <- -999
          se_psd[(11 * (r-1)) + a, ] <- -999 
          var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
          fit.measures[(11 * (r-1)) + a, ] <- -999 
          psd.var.coeff[(6 * (r-1)) + a - 2, ] <- -999
        }
        rm(file.name, fit, a, t0, tf)
      }
    }
    
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
    
    rm(analysis_syntax, ml_syntax, saveoutput_syntax)
    
    # Run model in Mplus
    cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
    
    if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
      system("taskkill /im Mplus.exe /f")
      
      perf[r, a] <- "timeout"
      times[r, a] <- round(tf[3], 2)
      parameters[(11 * (r-1)) + a + 1, ] <- -999
      se_psd[(11 * (r-1)) + a, ] <- -999 
      var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
      fit.measures[(11 * (r-1)) + a, ] <- -999
      psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
      
      rm(file.name, a, t0, tf)
    }else{
      fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
      
      if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
        
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        
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
        
        # Read MCMC draws
        samples <- read.table(paste0(getwd(), "/", folder, "samples_", file.name, ".dat"))
        
        # Delete burn-in and indicator for chain and draw
        samples <- samples[which(samples$V2 >= (max(samples$V2)/2 + 1)), -(1:2)]
        
        # Create matrix to store variance coefficients draws
        pdist.var.coeff <- matrix(NA, dim(samples)[1], 12)
        
        for(i in 1:dim(samples)[1]){
          within.estimates <- list( loadings = c(1, t(samples[i, 1:3])), state.var = samples[i, 8],
                                    error.var = c(t(samples[i, 4:7])))
          between.estimates <- list( loadings = c(1, t(samples[i, c(13:15)])), 
                                     trait.var = samples[i, 16])
          pdist.var.coeff[i, ] <- t(msst.var.coeff(within.parameters = within.estimates,
                                                   between.parameters = between.estimates))
        }
        
        var.coeff[(11 * (r-1)) + a + 1, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, median)
        psd.var.coeff[(6 * (r-1)) + a - 2, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, sd)
        
        # Save fit measures
        if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
        
        rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
      }else{
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        psd.var.coeff[(6 * (r-1)) + a - 2, ] <- -999
      }
      
      rm(file.name, fit, a, t0, tf)
    }
    
    # cuts + wide + maximum likelihod ----
    #give analysis a number
    a <- 5
    
    if(cond %in% c(5, 6, 11, 12, 17, 18)){
      perf[r, a] <- "notrun"
      rm(a)
    }else{
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
      cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
      runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
      cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
      
      if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
        system("taskkill /im Mplus.exe /f")
        
        perf[r, a] <- "timeout"
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        
        rm(file.name, a, t0, tf)
      }else{
        fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
        
        if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
          
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          
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
          var.coeff[(11 * (r-1)) + a + 1, ] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                                between.parameters = between.estimates))
          
          # Save fit measures
          fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                                   fit$summaries$aBIC)
          rm(within.estimates, between.estimates, estimates, stand.errors)
          
          
        }else{
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          parameters[(11 * (r-1)) + a + 1, ] <- -999
          se_psd[(11 * (r-1)) + a, ] <- -999 
          var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
          fit.measures[(11 * (r-1)) + a, ] <- -999 
        }
        rm(file.name, fit, a, t0, tf)
      }
    }
    
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
    
    ml_syntax <- write.mlcuts.to.Mplus(data$data.long[, -(1:2)])
    
    write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
    write(ml_syntax, paste0(folder,file.name,".inp"), append = T)
    
    rm(analysis_syntax, ml_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
    
    if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
      system("taskkill /im Mplus.exe /f")
      
      perf[r, a] <- "timeout"
      times[r, a] <- round(tf[3], 2)
      parameters[(11 * (r-1)) + a + 1, ] <- -999
      se_psd[(11 * (r-1)) + a, ] <- -999 
      var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
      fit.measures[(11 * (r-1)) + a, ] <- -999
      
      rm(file.name, a, t0, tf)
    }else{
      fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
      
      if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
        
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        
        # Save estimates
        estimates <- fit$parameters$unstandardized[,3]
        
        parameters[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- estimates
        
        # Save standard errors
        stand.errors <- fit$parameters$unstandardized[,4]
        
        se_psd[(11 * (r-1)) + a, c(1:9, 11:19, 21:24)] <- stand.errors
        
        # Compute variance coefficients
        within.estimates <- list( loadings = estimates[1:I], CS.var = estimates[I+1],
                                  US.var = estimates[(I+2):(2 * I +1)])
        between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], CT.var = estimates[ 4 * I + 2], 
                                   UT.var = estimates[(4 * I + 3):(5 * I + 2)])
        var.coeff[(11 * (r-1)) + a + 1, ] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                              between.parameters = between.estimates))
        
        # Save fit measures
        fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                                 fit$summaries$aBIC)
        rm(within.estimates, between.estimates, estimates, stand.errors)
        
      }else{
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999 
      }
      rm(file.name, fit, a, t0, tf)
    }
    
    # cuts + wide + Bayesian ----
    #give analysis a number
    a <- 7
    
    if(cond %in% c(3, 6, 9, 12, 15, 18)){
      perf[r, a] <- "notrun"
      rm(a)
    }else{
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
      
      rm(analysis_syntax, mplus_syntax, saveoutput_syntax)
      
      # Run model in Mplus
      cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
      runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
      cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
      
      if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
        system("taskkill /im Mplus.exe /f")
        
        perf[r, a] <- "timeout"
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
        
        rm(file.name, a, t0, tf)
      }else{
        fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
        
        if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
          
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          
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
          
          # Read MCMC draws
          samples <- read.table(paste0(getwd(), "/", folder, "samples_", file.name, ".dat"))
          
          # Delete burn-in and indicator for chain and draw
          samples <- samples[which(samples$V2 >= (max(samples$V2)/2 + 1)), -(1:2)]
          
          # Create matrix to store variance coefficients draws
          pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
          
          for(i in 1:dim(samples)[1]){
            within.estimates <- list( loadings = c(1, t(samples[i, c(5, 7, 9)])), CS.var = samples[i, 15],
                                      US.var = c(t(samples[i, 11:14])))
            between.estimates <- list( loadings = c(1, t(samples[i, c(6, 8, 10)])), CT.var = samples[i, 16], 
                                       UT.var = c(t(samples[i, 17:20])))
            pdist.var.coeff[i, ] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                     between.parameters = between.estimates))
          }
          
          var.coeff[(11 * (r-1)) + a + 1, ] <- apply(pdist.var.coeff, 2, median)
          psd.var.coeff[(6 * (r-1)) + a - 4, ] <- apply(pdist.var.coeff, 2, sd)
          
          # Save fit measures
          if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
          
          rm(within.estimates, between.estimates, estimates, post.sd, pdist.var.coeff, samples,i)
          
          
        }else{
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          parameters[(11 * (r-1)) + a + 1, ] <- -999
          se_psd[(11 * (r-1)) + a, ] <- -999 
          var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
          fit.measures[(11 * (r-1)) + a, ] <- -999
          psd.var.coeff[(6 * (r-1)) + a - 4, ] <- -999
        }
        rm(file.name, fit, a, t0, tf)
      }
    }
    
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
    
    rm(analysis_syntax, ml_syntax, saveoutput_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
    
    if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
      system("taskkill /im Mplus.exe /f")
      
      perf[r, a] <- "timeout"
      times[r, a] <- round(tf[3], 2)
      parameters[(11 * (r-1)) + a + 1, ] <- -999
      se_psd[(11 * (r-1)) + a, ] <- -999 
      var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
      fit.measures[(11 * (r-1)) + a, ] <- -999
      psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
      
      rm(file.name, a, t0, tf)
    }else{
      fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
      
      if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
        
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        
        # Save estimates
        estimates <- fit$parameters$unstandardized[,3]
        
        parameters[(11 * (r-1)) + a + 1, c(1:9, 11:19, 21:24)] <- estimates
        
        # Save posterior standard deviations
        post.sd <- fit$parameters$unstandardized[,4]
        
        se_psd[(11 * (r-1)) + a, c(1:9, 11:19, 21:24)] <- post.sd
        
        # Compute variance coefficients
        
        # Read MCMC draws
        samples <- read.table(paste0(getwd(), "/", folder, "samples_", file.name, ".dat"))
        
        # Delete burn-in and indicator for chain and draw
        samples <- samples[which(samples$V2 >= (max(samples$V2)/2 + 1)), -(1:2)]
        
        # Create matrix to store variance coefficients draws
        pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
        
        for(i in 1:dim(samples)[1]){
          within.estimates <- list( loadings = c(1, t(samples[i, 1:3])), CS.var = samples[i, 8],
                                    US.var = c(t(samples[i, 4:7])))
          between.estimates <- list( loadings = c(1, t(samples[i, 13:15])), CT.var = samples[i, 20], 
                                     UT.var = c(t(samples[i, 16:19])))
          pdist.var.coeff[i, ] <- t(cuts.var.coeff(within.parameters = within.estimates,
                                                   between.parameters = between.estimates))
        }
        
        var.coeff[(11 * (r-1)) + a + 1, ] <- apply(pdist.var.coeff, 2, median)
        psd.var.coeff[(6 * (r-1)) + a - 4, ] <- apply(pdist.var.coeff, 2, sd)
        
        # Save fit measures
        if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
        
        rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
        
      }else{
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        psd.var.coeff[(6 * (r-1)) + a - 4, ] <- -999
      }
      
      rm(file.name, fit, a, t0, tf)
    }
    
    # tso + wide + maximum likelihood ----
    #give analysis a number
    a <- 9
    
    if(cond %in% c(5, 6, 11, 12, 17, 18)){
      perf[r, a] <- "notrun"
      rm(a)
    }else{
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
      cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
      runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
      cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
      
      if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
        system("taskkill /im Mplus.exe /f")
        
        perf[r, a] <- "timeout"
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        
        rm(file.name, a, t0, tf)
      }else{
        fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
        
        if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
          
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          
          # Save estimates
          estimates <- fit$parameters$unstandardized[c(1:I, #loadings
                                                       ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)), # occasion variance
                                                       ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                       (2 * I * nT + 1), # autoregressive effect
                                                       ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I), # Intercepts
                                                       ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I), # trait indicator variances
                                                       ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)), # trait indicator covariances 
                                                     3]
          
          parameters[(11 * (r-1)) + a + 1, c(1:10, 15:18, 21:30)] <- estimates
          
          # Save standard errors
          stand.errors <- fit$parameters$unstandardized[c(1:I, #loadings
                                                          ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)), # occasion variance
                                                          ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                          (2 * I * nT + 1), # autoregressive effect
                                                          ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I), # Intercepts
                                                          ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I), # trait indicator variances
                                                          ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)), # trait indicator covariances 
                                                        4]
          
          se_psd[(11 * (r-1)) + a, c(1:10, 15:18, 21:30)] <- stand.errors
          
          # Compute and save variance coefficients
          within.estimates <- list( loadings = estimates[1:I], state.var = estimates[I+1],
                                    error.var = estimates[(I+2):(2 * I +1)], ar.effect = estimates[(2 * I + 2)])
          between.estimates <- list( trait.ind.var = estimates[(3 * I + 3):(4 * I + 2)])
          var.coeff[(11 * (r-1)) + a + 1, ] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
                                                             between.parameters = between.estimates)[,nT]
          
          # Save fit measures
          fit.measures[(11 * (r-1)) + a, 1:3] <- c(fit$summaries$AIC, fit$summaries$BIC, 
                                                   fit$summaries$aBIC)
          
          rm(within.estimates, between.estimates, estimates, stand.errors)
          
        }else{
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          parameters[(11 * (r-1)) + a + 1, ] <- -999
          se_psd[(11 * (r-1)) + a, ] <- -999 
          var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
          fit.measures[(11 * (r-1)) + a, ] <- -999 
        }
        rm(file.name, fit, a, t0, tf)
      }
    }
    
    # tso + wide + Bayesian ----
    #give analysis a number
    a <- 10
    
    if(cond %in% c(3, 6, 9, 12, 15, 18)){
      perf[r, a] <- "notrun"
      rm(a)
    }else{
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
      
      rm(analysis_syntax, ml_syntax, saveoutput_syntax)
      
      # Run modelin Mplus
      cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
      runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
      cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
      
      if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
        system("taskkill /im Mplus.exe /f")
        
        perf[r, a] <- "timeout"
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
        
        rm(file.name, a, t0, tf)
      }else{
        fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
        
        if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
          
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          
          # Save estimates
          estimates <- fit$parameters$unstandardized[c(1:I, #loadings
                                                       ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)), # occasion variance
                                                       ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                       (2 * I * nT + 1), # autoregressive effect
                                                       ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I), # Intercepts
                                                       ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I), # trait indicator variances
                                                       ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)), # trait indicator covariances 
                                                     3]
          
          parameters[(11 * (r-1)) + a + 1, c(1:10, 15:18, 21:30)] <- estimates
          
          # Save posterior standard deviations
          post.sd <- fit$parameters$unstandardized[c(1:I, #loadings
                                                     ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)), # occasion variance
                                                     ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                     (2 * I * nT + 1), # autoregressive effect
                                                     ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I), # Intercepts
                                                     ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I), # trait indicator variances
                                                     ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)), # trait indicator covariances 
                                                   4]
          
          se_psd[(11 * (r-1)) + a, c(1:10, 15:18, 21:30)] <- post.sd
          
          # Compute and save variance coefficients
          
          # Read MCMC draws
          samples <- read.table(paste0(getwd(), "/", folder, "samples_", file.name, ".dat"))
          
          # Delete burn-in and indicator for chain and draw
          samples <- samples[which(samples$V2 >= (max(samples$V2)/2 + 1)), -(1:2)]
          
          # Create matrix to store variance coefficients draws
          pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
          
          for(i in 1:dim(samples)[1]){
            within.estimates <- list( loadings = c(1, t(samples[i, 5:7])), state.var = samples[i, 13],
                                      error.var = c(t(samples[i, 8:11])), ar.effect = samples[i, 12])
            between.estimates <- list( trait.ind.var = c(t(samples[i, c(14, 16, 19, 23)])))
            pdist.var.coeff[i, ] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
                                                  between.parameters = between.estimates)[,nT]
          }
          
          var.coeff[(11 * (r-1)) + a + 1, ] <- apply(pdist.var.coeff, 2, median)
          psd.var.coeff[(6 * (r-1)) + a - 5, ] <- apply(pdist.var.coeff, 2, sd)
          
          # Save fit measures
          if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
          
          rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
          
        }else{
          perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
          times[r, a] <- round(tf[3], 2)
          parameters[(11 * (r-1)) + a + 1, ] <- -999
          se_psd[(11 * (r-1)) + a, ] <- -999 
          var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
          fit.measures[(11 * (r-1)) + a, ] <- -999
          psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
        }
        rm(file.name, fit, a, t0, tf)
      }
    }
    
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
    
    rm(analysis_syntax, ml_syntax, saveoutput_syntax)
    
    # Run modelin Mplus
    cat("\n"); print(Sys.time()); cat("\n"); t0 <- proc.time()
    runModels_2(paste0(getwd(),"/", folder,file.name,".inp"), timeout = timeout)
    cat("\n"); print(Sys.time()); cat("\n"); tf <- proc.time() - t0
    
    if(length(grep("^Mplus.exe", system2( 'tasklist' , stdout = TRUE ))) != 0){
      system("taskkill /im Mplus.exe /f")
      
      perf[r, a] <- "timeout"
      times[r, a] <- round(tf[3], 2)
      parameters[(11 * (r-1)) + a + 1, ] <- -999
      se_psd[(11 * (r-1)) + a, ] <- -999 
      var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
      fit.measures[(11 * (r-1)) + a, ] <- -999
      psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
      
      rm(file.name, a, data, dataModel, nT, na.prop, seed, t0, tf)
    }else{
      fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output
      
      if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
        
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        
        # Save estimates
        estimates <- fit$parameters$unstandardized[c(1:I, #loadings
                                                     (2 * I + 2), # state var
                                                     ((I + 2):(2 * I + 1)), # error variances
                                                     (I + 1), # autoregressive effect
                                                     (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2), # intercepts and indicator trait variances
                                                     (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)), #Covariances 
                                                   3]
        
        parameters[(11 * (r-1)) + a + 1, c(1:10, 15:18, 21:30)] <- estimates
        
        # Save posterior standard deviations
        post.sd <- fit$parameters$unstandardized[c(1:I, #loadings
                                                   (2 * I + 2), # state var
                                                   ((I + 2):(2 * I + 1)), # error variances
                                                   (I + 1), # autoregressive effect
                                                   (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2), # intercepts and indicator trait variances
                                                   (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)), #Covariances 
                                                 4]
        
        se_psd[(11 * (r-1)) + a, c(1:10, 15:18, 21:30)] <- post.sd
        
        # Compute and save variance coefficients
        
        # Read MCMC draws
        samples <- read.table(paste0(getwd(), "/", folder, "samples_", file.name, ".dat"))
        
        # Delete burn-in and indicator for chain and draw
        samples <- samples[which(samples$V2 >= (max(samples$V2)/2 + 1)), -(1:2)]
        
        # Create matrix to store variance coefficients draws
        pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
        
        for(i in 1:dim(samples)[1]){
          within.estimates <- list( loadings = c(1, t(samples[i, 1:3])), state.var = samples[i, 9],
                                    error.var = c(t(samples[i, 4:7])), ar.effect = samples[i, 8])
          between.estimates <- list( trait.ind.var = c(t(samples[i, c(14, 16, 19, 23)])))
          pdist.var.coeff[i, ] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
                                                between.parameters = between.estimates)[,nT]
        }
        
        var.coeff[(11 * (r-1)) + a + 1, ] <- apply(pdist.var.coeff, 2, median)
        psd.var.coeff[(6 * (r-1)) + a - 5, ] <- apply(pdist.var.coeff, 2, sd)
        
        # Save fit measures
        if(!is.null(fit$summaries$DIC)){fit.measures[(11 * (r-1)) + a, 4] <- c(fit$summaries$DIC)}
        
        rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
      }else{
        perf[r, a] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
        times[r, a] <- round(tf[3], 2)
        parameters[(11 * (r-1)) + a + 1, ] <- -999
        se_psd[(11 * (r-1)) + a, ] <- -999 
        var.coeff[(11 * (r-1)) + a + 1, ] <- -999 
        fit.measures[(11 * (r-1)) + a, ] <- -999
        psd.var.coeff[(6 * (r-1)) + a - 5, ] <- -999
      }
      
      rm(file.name, fit, a, data, dataModel, nT, na.prop, seed, t0, tf)
    }
    
    # clean directory
    unlink(paste0(getwd(), "/", folder, "*"))
    }
  
  # 4.0 Save output tables ----
  
  write.table(perf, file = paste0(folder, "Performance/", paste("performance", cond, sep = "_"), ".dat"),
              col.names = TRUE, row.names = FALSE, quote = TRUE)
  write.table(times, file = paste0(folder, "Times/", paste("times", cond, sep = "_"), ".dat"),
              col.names = TRUE, row.names = FALSE, quote = TRUE)
  write.table(parameters, file = paste0(folder, "Parameters/", paste("parameters", cond, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(se_psd, file = paste0(folder, "SE_PSD/", paste("se_psd", cond, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(var.coeff, file = paste0(folder, "Var_Coeff/", paste("var_coeff", cond, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(psd.var.coeff, file = paste0(folder, "PSD_Var_Coeff/", paste("psd_var_coeff", cond, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  write.table(fit.measures, file = paste0(folder, "Fit_Measures/", paste("fit_measures", cond, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  rm(perf, parameters, se_psd, var.coeff, psd.var.coeff, fit.measures, times)
}
time <- proc.time() - time0

# 5.0 Clean enviroment ----
rm(Cond, fit.measures.base, parameters.base, perf.base, psd.var.coeff.base, var.coeff.base, 
   se_psd.base, times.base, folder, I, N, timeout, R, cond, r)
