# This script is to perform the simulation study for the paper Using Latent State-Trait Theory 
#   to Analyze Intensive Longitudinal Data.

# CONTENTS
# 0.0 Prepare environment
# 1.0 Set up conditions
# 2.0 Set up output matrices
# 3.0 Simulation foreach loops
# 4.0 Save output tables
# 5.0 Stop cluster and clean enviroment

# 0.0 Prepare environment ----
rm(list=ls())
library(doParallel)
library(foreach)
library(MplusAutomation)
library(MASS)
library(bayesplot)
library(coda)
library(mcmcr)
# Easy code to source all the needed files
file.sources <- paste0("R/", list.files(paste0(getwd(), "/R")))
sapply(file.sources,source,.GlobalEnv)
rm(file.sources)
folder <- "Mplus_Simulation" #Folder to store results and all Mplus files

# Linux does stop mplus when timeout is reached
#  then to identify when an analysis timeout we should check and clean the warnings in r
# For example:
# runModels_2("ex3.1Bayes.inp", Mplus_command = "/home/p280329/mplusdemo/mpdemo")
# if(length(warnings())!=0){}else{}
# assign("last.warning", NULL, envir = baseenv())


# Create folders to store output in the working directory ----

if (!dir.exists(paste(getwd(), folder, sep = "/"))) {
  dir.create(paste(getwd(), folder, "Performance", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Times", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Parameters", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "SE_PSD", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Var_Coeff", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Fit_Measures", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "PSD_Var_Coeff", sep = "/"), recursive = TRUE)
  dir.create(paste(getwd(), folder, "Plots", sep = "/"), recursive = TRUE)
} else {
  if (!dir.exists(paste(getwd(), folder, "Performance", sep = "/"))) {
    dir.create(paste(getwd(), folder, "Performance", sep = "/"))
  }
  if (!dir.exists(paste(getwd(), folder, "Times", sep = "/"))) {
    dir.create(paste(getwd(), folder, "Times", sep = "/"))
  }
  if (!dir.exists(paste(getwd(), folder, "Parameters", sep = "/"))) {
    dir.create(paste(getwd(), folder, "Parameters", sep = "/"))
  }
  if (!dir.exists(paste(getwd(), folder, "SE_PSD", sep = "/"))) {
    dir.create(paste(getwd(), folder, "SE_PSD", sep = "/"))
  }
  if (!dir.exists(paste(getwd(), folder, "Var_Coeff", sep = "/"))) {
    dir.create(paste(getwd(), folder, "Var_Coeff", sep = "/"))
  }
  if (!dir.exists(paste(getwd(), folder, "Fit_Measures", sep = "/"))) {
    dir.create(paste(getwd(), folder, "Fit_Measures", sep = "/"))
  }
  if (!dir.exists(paste(getwd(), folder, "PSD_Var_Coeff", sep = "/"))) {
    dir.create(paste(getwd(), folder, "PSD_Var_Coeff", sep = "/"))
  }
  if (!dir.exists(paste(getwd(), folder, "Plots", sep = "/"))) {
    dir.create(paste(getwd(), folder, "Plots", sep = "/"))
  }
}

# Create function needed to combine foreach output ----
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}


# Remove files in output folder
#unlink(paste(getwd(), folder, "*", sep = "/"))

# add / to folder
folder <- paste0(folder, "/")

# 1.0 Set up conditions ----

# Fix conditions
N       <- 100   # Sample size
I       <- 4     # Number of variables
timeout <- 3600  # Time limit in seconds to force ending an analysis in Mplus
dplots  <- TRUE  # Save diagnosis plots for bayesian analyses.

# Manipulated conditions
timepoints       <- c(10, 15, 20, 30, 40, 50, 60, 90) # Number of measurement occasions
missingness      <- c(0.0, 0.1, 0.2)                  # Proportion of NAs
dataModel_to_sim <- c("msst", "cuts", "tso")          # Model to simulate the data
var_ratio        <- c("1:3", "1:1", "3:1")            # Ratio between consistency and occasion-specificity

Cond        <- expand.grid(timepoints, missingness, dataModel_to_sim, var_ratio)
names(Cond) <- c("nT", "na.prop", "dataModel", "ratio")

rm(timepoints, missingness, dataModel_to_sim, var_ratio)

# 2.0 Set up output matrices and vectors ----

# Create labels to name variables in output matrices 
models    <- c("msst", "cuts", "tso")
data_str  <- c("wide", "long")
estimator <- c("ml", "bayes")
labels    <- expand.grid(data_str, estimator, models)
labels    <- labels[-10, c(3,1,2)]
labels    <- apply(labels, 1, function(vec) paste(vec, collapse = "_"))
rm(models, data_str, estimator)

# Vector to store check.mplus output (Non-convergence, Warnings and Errors, Ok, and timeouts)
# Parameters, se, variance coefficients, and fit measures are saved if check.mplus output is Ok.
# The running time of each analisis is also saved.
perf.base           <- matrix(NA, 1, 11)
colnames(perf.base) <- labels
times.base          <- perf.base

# Matrix to store estimated parameters and standard errors/posterior standard deviations
parameters.base <- matrix(NA, 11,  30)
colnames(parameters.base) <- c(paste0("w_loading", 1:I), "common_state_var", paste0("unique_state_error_var", 1:I),
                               "ar_effect", paste0("b_loading", 1:I), paste0("intercept", 1:I), "common_trait_var",
                               "trait_mean", paste0("unique_indicator_trait_var", 1:I), paste0("cov", c(12, 13, 23, 14, 24, 34)))
se_psd.base <- parameters.base

# Matrix to store variance coefficients
var.coeff.base <- matrix(NA, 11,  I * 5)
colnames(var.coeff.base) <- paste0(rep(c( "ccon_pred_y", "ucon_upred_y", "tcon_con_y", "spe_y", "rel_y"), each = I), 
                              1:I)

# Matrix to store posterior standard deviations of variance coefficients
psd.var.coeff.base <- matrix(NA, 6,  I * 5)
colnames(psd.var.coeff.base) <- paste0(rep(c( "ccon_pred_y", "ucon_upred_y", "tcon_con_y", "spe_y", "rel_y"), each = I), 
                                   1:I)

# Matrix to store fit measures. When estimation is MLE, number of parameters, df, CFI, TLI, RMSEA, AIC,
# BIC, and adjusted BIC are saved if available. When estimation is Bayesian, number of parameters, DIC, 
# estimated number of parameters(pD), and the posterior predictive p-value (ppp) are saved if available.
fit.measures.base <- matrix(NA, 11,  13)
colnames(fit.measures.base) <- c("\\# Parameters", "Chitest", "df", "Chipval", 
                                 "CFI", "TLI", "RMSEA", "AIC", "BIC", "aBIC", 
                                 "DIC", "pD", "ppp")

# 3.0 Simulation foreach loops ----

# Setup parallel backend to use 6 parallel tasks between 6 to 24(when bayes) processors (cores):
cl <- makeCluster(8)
registerDoParallel(cl, cores = 8)

# Get conditions and replications from the batch file
args <- commandArgs(trailingOnly = TRUE)

time0 <- proc.time()
outcome.simulation <- foreach(cond = args[1]:args[2], .combine = 'list', .multicombine = TRUE) %:%
  foreach(r = args[3]:args[4], .combine = 'comb', .multicombine = TRUE, 
          .packages = c("MplusAutomation", "MASS", "bayesplot", "coda", "mcmcr")) %dopar% {
            
            # Copy table(s) to save results
            perf          <- perf.base
            times         <- times.base
            parameters    <- parameters.base
            se_psd        <- se_psd.base
            var.coeff     <- var.coeff.base
            fit.measures  <- fit.measures.base
            psd.var.coeff <- psd.var.coeff.base
            
            # Define rownames of output tables
            rownames(parameters)    <- c(paste(labels, paste0("r", r), sep = "_"))
            rownames(se_psd)        <- c(paste(labels, paste0("r", r), sep = "_"))
            rownames(var.coeff)     <- c(paste(labels, paste0("r", r), sep = "_"))
            rownames(fit.measures)  <- c(paste(labels, paste0("r", r), sep = "_"))
            rownames(psd.var.coeff) <- c(paste(labels[grep(pattern = "bayes", labels)],
                                                    paste0("r", r), sep = "_"))
            
            # Define seed and conditions to be use
            nT        <- Cond[cond, 1]
            na.prop   <- Cond[cond, 2]
            dataModel <- Cond[cond, 3]
            ratio     <- Cond[cond, 4]
            seed      <- 1000 * cond + r
            set.seed(seed)
            
            # Data simulation given dataModel and ratio ----
            
            if (dataModel == "msst") {
              
              if (ratio == "1:3") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent state
                var_state      <- 1.8                    # Variance latent state residual
                var_m_error    <- c(0.6, 0.25, 0.7, 0.5) # Variance of measurement errors
                
                within.parameters <- list(loadings  = state_loadings, 
                                          state.var = var_state, 
                                          error.var = var_m_error)
                
                # Between Paramaters
                intercepts     <- seq(0, by = 0.2, length.out = I) # intercepts
                trait_loadings <- c(1, 0.5, 1.3, 0.8)              # loading parametes for the latent trait
                var_trait      <- 0.6                                # variance latent trait variable
                mean_trait     <- 4                                # mean latent trait variable
                
                between.parameters <- list(loadings   = trait_loadings, 
                                           intercepts = intercepts, 
                                           trait.mean = mean_trait,
                                           trait.var  = var_trait)
                
                rm(state_loadings, var_state, var_m_error, intercepts, trait_loadings, var_trait, 
                   mean_trait)
                
                # Simulate data
                data <- sim.data.msst(N, nT, I, 
                                      within.parameters  = within.parameters, 
                                      between.parameters = between.parameters, 
                                      na.prop            = na.prop,
                                      seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
              
              if (ratio == "1:1") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent state
                var_state      <- 1.2                    # Variance latent state residual
                var_m_error    <- c(0.6, 0.25, 0.7, 0.5) # Variance of measurement errors
                
                within.parameters <- list(loadings  = state_loadings, 
                                          state.var = var_state, 
                                          error.var = var_m_error)
                
                # Between Paramaters
                intercepts     <- seq(0, by = 0.2, length.out = I) # intercepts
                trait_loadings <- c(1, 0.5, 1.3, 0.8)              # loading parametes for the latent trait
                var_trait      <- 1.2                              # variance latent trait variable
                mean_trait     <- 4                                # mean latent trait variable
                
                between.parameters <- list(loadings   = trait_loadings, 
                                           intercepts = intercepts, 
                                           trait.mean = mean_trait,
                                           trait.var  = var_trait)
                
                rm(state_loadings, var_state, var_m_error, intercepts, trait_loadings, var_trait, 
                   mean_trait)
                
                # Simulate data
                data <- sim.data.msst(N, nT, I, 
                                      within.parameters  = within.parameters, 
                                      between.parameters = between.parameters, 
                                      na.prop            = na.prop,
                                      seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
              
              if (ratio == "3:1") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent state
                var_state      <- 0.6                    # Variance latent state residual
                var_m_error    <- c(0.6, 0.25, 0.7, 0.5) # Variance of measurement errors
                
                within.parameters <- list(loadings  = state_loadings, 
                                          state.var = var_state, 
                                          error.var = var_m_error)
                
                # Between Paramaters
                intercepts     <- seq(0, by = 0.2, length.out = I) # intercepts
                trait_loadings <- c(1, 0.5, 1.3, 0.8)              # loading parametes for the latent trait
                var_trait      <- 1.8                              # variance latent trait variable
                mean_trait     <- 4                                # mean latent trait variable
                
                between.parameters <- list(loadings   = trait_loadings, 
                                           intercepts = intercepts, 
                                           trait.mean = mean_trait,
                                           trait.var  = var_trait)
                
                rm(state_loadings, var_state, var_m_error, intercepts, trait_loadings, var_trait, 
                   mean_trait)
                
                # Simulate data
                data <- sim.data.msst(N, nT, I, 
                                      within.parameters  = within.parameters, 
                                      between.parameters = between.parameters, 
                                      na.prop            = na.prop,
                                      seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
            }
            
            if (dataModel == "cuts") {
              
              if (ratio == "1:3") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent common state
                var_CS         <- 1.8                    # Variance latent common state
                var_US         <- c(0.6, 0.25, 0.8, 0.5) # Variance of latent unique states
                
                within.parameters <- list(loadings = state_loadings, 
                                          CS.var   = var_CS, 
                                          US.var   = var_US)
                
                # Between Paramaters
                intercepts     <- seq(2, by = 0.5, length.out = I) # intercepts
                trait_loadings <- c(1, 0.5, 1.3, 0.8)              # loading parametes for the latent common trait
                var_CT         <- 0.4                              # variance latent common trait
                var_UT         <- c(0.2, 0.1, 0.25, 0.15)          # variance latent unique traits
                
                between.parameters <- list(loadings   = trait_loadings, 
                                           intercepts = intercepts, 
                                           CT.var     = var_CT, 
                                           UT.var     = var_UT)
                
                rm(state_loadings, var_CS, var_US, intercepts, trait_loadings, var_CT, var_UT)
                
                # Simulate data
                data <- sim.data.cuts(N, nT, I, 
                                      within.parameters  = within.parameters, 
                                      between.parameters = between.parameters, 
                                      na.prop            = na.prop,
                                      seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
              
              if (ratio == "1:1") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent common state
                var_CS         <- 1.2                    # Variance latent common state
                var_US         <- c(0.6, 0.25, 0.8, 0.5) # Variance of latent unique states
                
                within.parameters <- list(loadings = state_loadings, 
                                          CS.var   = var_CS, 
                                          US.var   = var_US)
                
                # Between Paramaters
                intercepts     <- seq(2, by = 0.5, length.out = I) # intercepts
                trait_loadings <- c(1, 0.5, 1.3, 0.8)              # loading parametes for the latent common trait
                var_CT         <- 1.2                              # variance latent common trait
                var_UT         <- c(0.2, 0.1, 0.25, 0.15)          # variance latent unique traits
                
                between.parameters <- list(loadings   = trait_loadings, 
                                           intercepts = intercepts, 
                                           CT.var     = var_CT, 
                                           UT.var     = var_UT)
                
                rm(state_loadings, var_CS, var_US, intercepts, trait_loadings, var_CT, var_UT)
                
                # Simulate data
                data <- sim.data.cuts(N, nT, I, 
                                      within.parameters  = within.parameters, 
                                      between.parameters = between.parameters, 
                                      na.prop            = na.prop,
                                      seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
              
              if (ratio == "3:1") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent common state
                var_CS         <- 0.6                    # Variance latent common state
                var_US         <- c(0.6, 0.25, 0.8, 0.5) # Variance of latent unique states
                
                within.parameters <- list(loadings = state_loadings, 
                                          CS.var   = var_CS, 
                                          US.var   = var_US)
                
                # Between Paramaters
                intercepts     <- seq(2, by = 0.5, length.out = I) # intercepts
                trait_loadings <- c(1, 0.5, 1.3, 0.8)              # loading parametes for the latent common trait
                var_CT         <- 1.8                              # variance latent common trait
                var_UT         <- c(0.2, 0.1, 0.25, 0.15)          # variance latent unique traits
                
                between.parameters <- list(loadings   = trait_loadings, 
                                           intercepts = intercepts, 
                                           CT.var     = var_CT, 
                                           UT.var     = var_UT)
                
                rm(state_loadings, var_CS, var_US, intercepts, trait_loadings, var_CT, var_UT)
                
                # Simulate data
                data <- sim.data.cuts(N, nT, I, 
                                      within.parameters  = within.parameters, 
                                      between.parameters = between.parameters, 
                                      na.prop            = na.prop,
                                      seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
            }
            
            if (dataModel == "tso") {
              
              if (ratio == "1:3") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent state
                var_state      <- 1.8                    # Variance latent state residual
                var_error      <- c(0.6, 0.25, 0.7, 0.5) # Variance of latent measurement errors
                ar_effect      <- 0.5                    # autoregressive effect on the latent state residuals
                
                within.parameters <- list(loadings  = state_loadings, 
                                          ar.effect = ar_effect, 
                                          error.var = var_error,
                                          state.var = var_state)
                
                rm(state_loadings, var_state, var_error, ar_effect)
                
                # Between Paramaters 
                intercepts     <- seq(2, by = 0.5, length.out = I) # intercepts
                var_ind_traits <- c(0.6, 0.2, 0.8, 0.5)            # variance latent indicator trait variables
                
                # Correlation matrix of the trait indicators
                # Correlation are high between 0.7 and 0.9. This correlation matrix was
                # generated with the seed 13001 during the first run of the simulation.
                Rcor <- matrix(c(1.0, 0.8, 0.9, 0.9,
                                 0.8, 1.0, 0.8, 0.7,
                                 0.9, 0.8, 1.0, 0.7,
                                 0.9, 0.7, 0.7, 1.0), 
                               nrow = I, ncol = I, byrow = TRUE)
                
                between.parameters <- list(intercepts    = intercepts, 
                                           trait.ind.var = var_ind_traits, 
                                           cor.matrix    = Rcor)
                
                rm(intercepts, var_ind_traits, Rcor)
                
                # Simulate data
                data <- sim.data.tso(N, nT, I, 
                                     within.parameters  = within.parameters, 
                                     between.parameters = between.parameters, 
                                     na.prop            = na.prop, 
                                     seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
              
              if (ratio == "1:1") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent state
                var_state      <- 1.2                    # Variance latent state residual
                var_error      <- c(0.6, 0.25, 0.7, 0.5) # Variance of latent measurement errors
                ar_effect      <- 0.5                    # autoregressive effect on the latent state residuals
                
                within.parameters <- list(loadings  = state_loadings, 
                                          ar.effect = ar_effect, 
                                          error.var = var_error,
                                          state.var = var_state)
                
                rm(state_loadings, var_state, var_error, ar_effect)
                
                # Between Paramaters 
                intercepts     <- seq(2, by = 0.5, length.out = I) # intercepts
                var_ind_traits <- c(1.2, 0.35, 1.8, 0.8)           # variance latent indicator trait variables
                
                # Correlation matrix of the trait indicators
                # Correlation are high between 0.7 and 0.9. This correlation matrix was
                # generated with the seed 13001 during the first run of the simulation.
                Rcor <- matrix(c(1.0, 0.8, 0.9, 0.9,
                                 0.8, 1.0, 0.8, 0.7,
                                 0.9, 0.8, 1.0, 0.7,
                                 0.9, 0.7, 0.7, 1.0), 
                               nrow = I, ncol = I, byrow = TRUE)
                
                between.parameters <- list(intercepts    = intercepts, 
                                           trait.ind.var = var_ind_traits, 
                                           cor.matrix    = Rcor)
                
                rm(intercepts, var_ind_traits, Rcor)
                
                # Simulate data
                data <- sim.data.tso(N, nT, I, 
                                     within.parameters  = within.parameters, 
                                     between.parameters = between.parameters, 
                                     na.prop            = na.prop, 
                                     seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
              
              if (ratio == "3:1") {
                # Within Parameters
                state_loadings <- c(1, 0.5, 1.3, 0.8)    # loading parameters for the latent state
                var_state      <- 0.6                    # Variance latent state residual
                var_error      <- c(0.6, 0.25, 0.7, 0.5) # Variance of latent measurement errors
                ar_effect      <- 0.5                    # autoregressive effect on the latent state residuals
                
                within.parameters <- list(loadings  = state_loadings, 
                                          ar.effect = ar_effect, 
                                          error.var = var_error,
                                          state.var = var_state)
                
                rm(state_loadings, var_state, var_error, ar_effect)
                
                # Between Paramaters 
                intercepts     <- seq(2, by = 0.5, length.out = I) # intercepts
                var_ind_traits <- c(1.8, 0.5, 2.5, 1.1)            # variance latent indicator trait variables
                
                # Correlation matrix of the trait indicators
                # Correlation are high between 0.7 and 0.9. This correlation matrix was
                # generated with the seed 13001 during the first run of the simulation.
                Rcor <- matrix(c(1.0, 0.8, 0.9, 0.9,
                                 0.8, 1.0, 0.8, 0.7,
                                 0.9, 0.8, 1.0, 0.7,
                                 0.9, 0.7, 0.7, 1.0), 
                               nrow = I, ncol = I, byrow = TRUE)
                
                between.parameters <- list(intercepts    = intercepts, 
                                           trait.ind.var = var_ind_traits, 
                                           cor.matrix    = Rcor)
                
                rm(intercepts, var_ind_traits, Rcor)
                
                # Simulate data
                data <- sim.data.tso(N, nT, I, 
                                     within.parameters  = within.parameters, 
                                     between.parameters = between.parameters, 
                                     na.prop            = na.prop, 
                                     seed               = seed)
                
                rm(within.parameters, between.parameters)
              }
            }
            
            # Fitting the models to the simulated data ----
            # Next, the different models are fitted to the data in both wide and long 
            # format and in both MLE and bayesian estimation. In total, the simulated
            # dataset is analyzed 11 times because the TSO in long format with MLE
            # is not possible.
            
            # msst + wide + maximum likelihod ----
            # Give analysis a number
            a <- 1
            
            file.name <- paste("msst", "wide", "ml", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.wide, paste0(folder,file.name,".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.wide)[-1],
                                                   analysis_type = "GENERAL",
                                                   estimator     = "ML",
                                                   iterations    = 50000)
            
            mplus_syntax <- write.msst.to.Mplus(data$data.wide[,-1], 
                                                neta               = nT, 
                                                ntheta             = 1, 
                                                equiv.assumption   = list(tau = "cong", theta = "cong"),
                                                scale.invariance   = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                                homocedasticity.assumption = list(error = TRUE, state.red = TRUE),
                                                second.order.trait = FALSE)
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T) 
            write(mplus_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, mplus_syntax)
            
            # Run model in Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]        <- "timeout"
              times[1, a]       <- round(tf[3], 2)
              parameters[a, ]   <- -999
              se_psd[a, ]       <- -999 
              var.coeff[a, ]    <- -999 
              fit.measures[a, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:I, #loadings
                                                             (I * nT * 3) + ((nT + 1) * nT / 2) + 2,                                                          # state var 
                                                             ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT + 1) * nT / 2) + nT + I + 2), # error variances
                                                             ((I * nT + 1):(I * nT + I)),                                                                     # trait loadings
                                                             ((I * nT * 2) + ((nT + 1) * nT / 2) + 2):((I * nT * 2) + ((nT + 1) * nT / 2) + I + 1),           # intercepts
                                                             ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 2),                                                   # trait variance
                                                             ((I * nT * 2) + ((nT + 1) * nT / 2) + 1)),                                                       # trait mean
                                                           3]
                
                parameters[a, c(1:9, 11:20)] <- estimates
                
                # Save standard errors
                stand.errors <- fit$parameters$unstandardized[c(1:I,                                                                                             # loadings
                                                                (I * nT * 3) + ((nT + 1) * nT / 2) + 2,                                                          # state var 
                                                                ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT + 1) * nT / 2) + nT + I + 2), # error variances
                                                                ((I * nT + 1):(I * nT + I)),                                                                     # trait loadings
                                                                ((I * nT * 2) + ((nT + 1) * nT / 2) + 2):((I * nT * 2) + ((nT + 1) * nT / 2) + I + 1),           # intercepts
                                                                ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 2),                                                   # trait variance
                                                                ((I * nT * 2) + ((nT + 1) * nT / 2) + 1)),                                                       # trait mean
                                                              4]
                
                se_psd[a, c(1:9, 11:20)] <- stand.errors
                
                # Compute and save variance coefficients
                
                within.estimates  <- list(loadings  = estimates[1:I], 
                                          state.var = estimates[I + 1],
                                          error.var = estimates[(I + 2):(2 * I + 1)])
                between.estimates <- list(loadings  = estimates[(2 * I + 2):(3 * I + 1)], 
                                          trait.var = estimates[ 4 * I + 2])
                var.coeff[a, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters  = within.estimates,
                                                                      between.parameters = between.estimates))
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$ChiSqM_Value)) {
                  fit.measures[a, 2] <- c(fit$summaries$ChiSqM_Value)
                }
                if (!is.null(fit$summaries$ChiSqM_DF)) {
                  fit.measures[a, 3] <- c(fit$summaries$ChiSqM_DF)
                }
                if (!is.null(fit$summaries$ChiSqM_PValue)) {
                  fit.measures[a, 4] <- c(fit$summaries$ChiSqM_PValue)
                }
                if (!is.null(fit$summaries$CFI)) {
                  fit.measures[a, 5] <- c(fit$summaries$CFI)
                  }
                if (!is.null(fit$summaries$TLI)) {
                  fit.measures[a, 6] <- c(fit$summaries$TLI)
                  }
                if (!is.null(fit$summaries$RMSEA_Estimate)) {
                  fit.measures[a, 7] <- c(fit$summaries$RMSEA_Estimate)
                  }
                if (!is.null(fit$summaries$AIC)) {
                  fit.measures[a, 8] <- c(fit$summaries$AIC)
                  }
                if (!is.null(fit$summaries$BIC)) {
                  fit.measures[a, 9] <- c(fit$summaries$BIC)
                  }
                if (!is.null(fit$summaries$aBIC)) {
                  fit.measures[a, 10] <- c(fit$summaries$aBIC)
                  }
                
                rm(within.estimates, between.estimates, estimates, stand.errors)
                
              } else {
                perf[1, a]        <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]       <- round(tf[3], 2)
                parameters[a, ]   <- -999
                se_psd[a, ]       <- -999 
                var.coeff[a, ]    <- -999 
                fit.measures[a, ] <- -999 
              }
              
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # msst + long + maximun likelihood ----
            # Give analysis a number
            a <- 2
            
            file.name <- paste("msst", "long", "ml", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.long)[-(1:2)],
                                                   cluster       = names(data$data.long)[1],
                                                   analysis_type = "TWOLEVEL",
                                                   estimator     = "ML",
                                                   iterations    = 50000)
            
            ml_syntax <- write.mlmsst.to.Mplus(data$data.long[, -(1:2)])
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(ml_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, ml_syntax)
            
            # Run model in Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]        <- "timeout"
              times[1, a]       <- round(tf[3], 2)
              parameters[a, ]   <- -999
              se_psd[a, ]       <- -999 
              var.coeff[a, ]    <- -999 
              fit.measures[a, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:(2 * I + 1),             # loadings, state var, error variances 
                                                             ((3 * I + 2):(4 * I + 1)), # loadings 
                                                             ((5 * I + 3):(6 * I + 2)), # intercepts
                                                             (6 * I + 3),               # trait var
                                                             (4 * I + 2)),              # trait mean
                                                           3]
                
                parameters[a, c(1:9, 11:20)] <- estimates
                
                # Save standard errors
                stand.errors <- fit$parameters$unstandardized[c(1:(2 * I + 1),             # loadings, state var, error variances 
                                                                ((3 * I + 2):(4 * I + 1)), # loadings 
                                                                ((5 * I + 3):(6 * I + 2)), # intercepts
                                                                (6 * I + 3),               # trait var
                                                                (4 * I + 2)),              # trait mean
                                                              4]
                
                se_psd[a, c(1:9, 11:20)] <- stand.errors
                
                # Compute and save variance coefficients
                within.estimates  <- list(loadings  = estimates[1:I], 
                                          state.var = estimates[I + 1],
                                          error.var = estimates[(I + 2):(2 * I + 1)])
                between.estimates <- list(loadings  = estimates[(2 * I + 2):(3 * I + 1)], 
                                          trait.var = estimates[ 4 * I + 2])
                var.coeff[a, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters  = within.estimates,
                                                                      between.parameters = between.estimates))
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$ChiSqM_Value)) {
                  fit.measures[a, 2] <- c(fit$summaries$ChiSqM_Value)
                }
                if (!is.null(fit$summaries$ChiSqM_DF)) {
                  fit.measures[a, 3] <- c(fit$summaries$ChiSqM_DF)
                }
                if (!is.null(fit$summaries$ChiSqM_PValue)) {
                  fit.measures[a, 4] <- c(fit$summaries$ChiSqM_PValue)
                }
                if (!is.null(fit$summaries$CFI)) {
                  fit.measures[a, 5] <- c(fit$summaries$CFI)
                }
                if (!is.null(fit$summaries$TLI)) {
                  fit.measures[a, 6] <- c(fit$summaries$TLI)
                }
                if (!is.null(fit$summaries$RMSEA_Estimate)) {
                  fit.measures[a, 7] <- c(fit$summaries$RMSEA_Estimate)
                }
                if (!is.null(fit$summaries$AIC)) {
                  fit.measures[a, 8] <- c(fit$summaries$AIC)
                }
                if (!is.null(fit$summaries$BIC)) {
                  fit.measures[a, 9] <- c(fit$summaries$BIC)
                }
                if (!is.null(fit$summaries$aBIC)) {
                  fit.measures[a, 10] <- c(fit$summaries$aBIC)
                }
                
                rm(within.estimates, between.estimates, estimates, stand.errors)
                
              } else {
                perf[1, a]        <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]       <- round(tf[3], 2)
                parameters[a, ]   <- -999
                se_psd[a, ]       <- -999 
                var.coeff[a, ]    <- -999 
                fit.measures[a, ] <- -999 
              }
              
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # msst + wide + Bayesian ----
            # Give analysis a number
            a <- 3
            
            file.name <- paste("msst", "wide", "bayes", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.wide)[-1],
                                                   analysis_type = "GENERAL",
                                                   estimator     = "BAYES",
                                                   iterations    = 5000)
            
            mplus_syntax <- write.msst.to.Mplus(data$data.wide[,-1], 
                                                neta               = nT, 
                                                ntheta             = 1, 
                                                equiv.assumption   = list(tau = "cong", theta = "cong"),
                                                scale.invariance   = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                                homocedasticity.assumption = list(error = TRUE, state.red = TRUE),
                                                second.order.trait = FALSE)
            saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", 
                                        "samples_", 
                                        file.name, 
                                        ".dat;", 
                                        "\nOUTPUT: TECH8;",
                                        "\nPLOT: TYPE = PLOT2;")
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(mplus_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(saveoutput_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, mplus_syntax, saveoutput_syntax)
            
            # Run model in Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]             <- "timeout"
              times[1, a]            <- round(tf[3], 2)
              parameters[a, ]        <- -999
              se_psd[a, ]            <- -999 
              var.coeff[a, ]         <- -999 
              fit.measures[a, ]      <- -999
              psd.var.coeff[a - 2, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:I,                                                                                             # loadings
                                                             (I * nT * 3) + ((nT + 1) * nT / 2) + 2,                                                          # state var 
                                                             ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT + 1) * nT / 2) + nT + I + 2), # error variances
                                                             ((I * nT + 1):(I * nT + I)),                                                                     # trait loadings
                                                             ((I * nT * 2) + ((nT + 1) * nT / 2) + 2):((I * nT * 2) + ((nT + 1) * nT / 2) + I + 1),           # intercepts
                                                             ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 2),                                                   # trait variance
                                                             ((I * nT * 2) + ((nT + 1) * nT / 2) + 1)),                                                       # trait mean
                                                           3]
                
                parameters[a, c(1:9, 11:20)] <- estimates
                
                # Save posterior standard deviations
                post.sd <- fit$parameters$unstandardized[c(1:I,                                                                                             # loadings
                                                           (I * nT * 3) + ((nT + 1) * nT / 2) + 2,                                                          # state var 
                                                           ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 3):((I * nT * 3) + ((nT + 1) * nT / 2) + nT + I + 2), # error variances
                                                           ((I * nT + 1):(I * nT + I)),                                                                     # trait loadings
                                                           ((I * nT * 2) + ((nT + 1) * nT / 2) + 2):((I * nT * 2) + ((nT + 1) * nT / 2) + I + 1),           # intercepts
                                                           ((I * nT * 3) + ((nT + 1) * nT / 2) + nT + 2),                                                   # trait variance
                                                           ((I * nT * 2) + ((nT + 1) * nT / 2) + 1)),                                                       # trait mean
                                                         4]
                
                se_psd[a, c(1:9, 11:20)] <- post.sd
                
                # Compute and save variance coefficients
                
                # Put valid samples in an array 
                fit_samples <- as.array(fit$bparameters$valid_draw)
                fit_samples <- aperm(fit_samples, perm = c(1, 3, 2))
                fit_samples <- fit_samples[, , -(1:2)] # Take out chain and iteration number.
                
                # Rearrange samples in a matrix
                samples <- matrix(fit_samples, prod(dim(fit_samples)[1:2]), dim(fit_samples)[3])
                
                # Create matrix to store variance coefficients draws
                pdist.var.coeff <- matrix(NA, dim(samples)[1], 12)
                
                for (i in 1:dim(samples)[1]) {
                  within.estimates  <- list(loadings  = c(1, t(samples[i, c(4, 6, 8)])), 
                                            state.var = samples[i, 15],
                                            error.var = c(t(samples[i, 10:13])))
                  between.estimates <- list(loadings  = c(1, t(samples[i, c(5, 7, 9)])), 
                                            trait.var = samples[i, 16])
                  pdist.var.coeff[i, ] <- t(msst.var.coeff(within.parameters  = within.estimates,
                                                           between.parameters = between.estimates))
                }
                
                var.coeff[a, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, median)
                psd.var.coeff[a - 2, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, sd)
                
                # Save mcmc diagnosis plots
                if (dplots) {
                  pdf(file = paste0(folder, "Plots/", file.name, ".pdf"))
                  print(mcmc_rhat(apply(fit_samples, 3, function (x) mcmcr::rhat(coda::as.mcmc(x)))))
                  # Diagnosis within loadings
                  print(mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(4, 6, 8)]))
                  print(mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(4, 6, 8)]))
                  # Diagnosis between loadings 
                  print(mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(5, 7, 9)]))
                  print(mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(5, 7, 9)]))
                  # Diagnosis error variances and state variance
                  print(mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(10:13, 15)]))
                  print(mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(10:13, 15)]))
                  # Diagnosis trait mean and trait variance
                  print(mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(14, 16)]))
                  print(mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(14, 16)]))
                  # Diagnosis intercepts
                  print(mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[1:3]))
                  print(mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[1:3]))
                  dev.off()
                }
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$DIC)) {
                  fit.measures[a, 11] <- c(fit$summaries$DIC)
                }
                if (!is.null(fit$summaries$pD)) {
                  fit.measures[a, 12] <- c(fit$summaries$pD)
                }
                if (!is.null(fit$summaries$PostPred_PValue)) {
                  fit.measures[a, 13] <- c(fit$summaries$PostPred_PValue)
                }
                
                rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
                
              } else {
                perf[1, a]             <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]            <- round(tf[3], 2)
                parameters[a, ]        <- -999
                se_psd[a, ]            <- -999 
                var.coeff[a, ]         <- -999 
                fit.measures[a, ]      <- -999 
                psd.var.coeff[a - 2, ] <- -999
              }
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, file.name, ".gh5"))) {unlink(paste0(folder, file.name, ".gh5"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # msst + long + Bayesian ----
            # Give analysis a number
            a <- 4
            
            file.name <- paste("msst", "long", "bayes", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.long)[-(1:2)],
                                                   cluster       = names(data$data.long)[1],
                                                   analysis_type = "TWOLEVEL",
                                                   estimator     = "BAYES",
                                                   iterations    = 5000)
            
            ml_syntax <- write.mlmsst.to.Mplus(data$data.long[, -(1:2)])
            
            ml_syntax <- gsub("@0;", "@0.001;", ml_syntax) # replace 0 constraints to 0.001
            
            saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", 
                                        "samples_", 
                                        file.name, 
                                        ".dat;", 
                                        "\nOUTPUT: TECH8;",
                                        "\nPLOT: TYPE = PLOT2;")
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(ml_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(saveoutput_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, ml_syntax, saveoutput_syntax)
            
            # Run model in Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]             <- "timeout"
              times[1, a]            <- round(tf[3], 2)
              parameters[a, ]        <- -999
              se_psd[a, ]            <- -999 
              var.coeff[a, ]         <- -999 
              fit.measures[a, ]      <- -999
              psd.var.coeff[a - 2, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:(2 * I + 1),             # loadings, state var, error variances 
                                                             ((3 * I + 2):(4 * I + 1)), # loadings 
                                                             ((5 * I + 3):(6 * I + 2)), # intercepts
                                                             (6 * I + 3),               # trait var
                                                             (4 * I + 2)),              # trait mean
                                                           3]
                
                parameters[a, c(1:9, 11:20)] <- estimates
                
                # Save posterior standard deviations
                post.sd <- fit$parameters$unstandardized[c(1:(2 * I + 1),             # loadings, state var, error variances 
                                                           ((3 * I + 2):(4 * I + 1)), # loadings 
                                                           ((5 * I + 3):(6 * I + 2)), # intercepts
                                                           (6 * I + 3),               # trait var
                                                           (4 * I + 2)),              # trait mean
                                                         4]
                
                se_psd[a, c(1:9, 11:20)] <- post.sd
                
                # Compute and save variance coefficients
                
                # Put valid samples in an array 
                fit_samples <- as.array(fit$bparameters$valid_draw)
                fit_samples <- aperm(fit_samples, perm = c(1, 3, 2))
                fit_samples <- fit_samples[, , -(1:2)] # Take out chain and iteration number.
                
                # Rearrange samples in a matrix
                samples <- matrix(fit_samples, prod(dim(fit_samples)[1:2]), dim(fit_samples)[3])
                
                # Create matrix to store variance coefficients draws
                pdist.var.coeff <- matrix(NA, dim(samples)[1], 12)
                
                for (i in 1:dim(samples)[1]) {
                  within.estimates  <- list(loadings  = c(1, t(samples[i, 1:3])), 
                                            state.var = samples[i, 8],
                                            error.var = c(t(samples[i, 4:7])))
                  between.estimates <- list(loadings  = c(1, t(samples[i, c(13:15)])), 
                                            trait.var = samples[i, 16])
                  pdist.var.coeff[i, ] <- t(msst.var.coeff(within.parameters  = within.estimates,
                                                           between.parameters = between.estimates))
                }
                
                var.coeff[a, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, median)
                psd.var.coeff[a - 2, (2 * I + 1):(5 * I)] <- apply(pdist.var.coeff, 2, sd)
                
                # Save mcmc diagnosis plots
                pdf(file = paste0(folder, "Plots/", file.name, ".pdf"))
                if (dplots) {
                  mcmc_rhat(apply(fit_samples, 3, function (x) mcmcr::rhat(coda::as.mcmc(x))))
                  # Diagnosis within loadings
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(1:3)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(1:3)])
                  # Diagnosis between loadings 
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[13:15])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[13:15])
                  # Diagnosis error variances and state variance
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[4:8])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[4:8])
                  # Diagnosis trait mean and trait variance
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(12, 16)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(12, 16)])
                  # Diagnosis intercepts
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[9:11])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[9:11])
                }
                dev.off()
                if (!dplots) {unlink(paste0(folder, "Plots/", file.name, ".pdf"))}
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$DIC)) {
                  fit.measures[a, 11] <- c(fit$summaries$DIC)
                }
                if (!is.null(fit$summaries$pD)) {
                  fit.measures[a, 12] <- c(fit$summaries$pD)
                }
                if (!is.null(fit$summaries$PostPred_PValue)) {
                  fit.measures[a, 13] <- c(fit$summaries$PostPred_PValue)
                }
                
                rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
                
              } else {
                perf[1, a]             <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]            <- round(tf[3], 2)
                parameters[a, ]        <- -999
                se_psd[a, ]            <- -999 
                var.coeff[a, ]         <- -999 
                fit.measures[a, ]      <- -999
                psd.var.coeff[a - 2, ] <- -999
              }
              
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, file.name, ".gh5"))) {unlink(paste0(folder, file.name, ".gh5"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # cuts + wide + maximum likelihod ----
            # Give analysis a number
            a <- 5
            
            file.name <- paste("cuts", "wide", "ml", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.wide)[-1],
                                                   analysis_type = "GENERAL",
                                                   estimator     = "ML",
                                                   iterations    = 50000)
            
            
            mplus_syntax <- write.cuts.to.Mplus(data$data.wide[,-1],  
                                                nstate                 = nT,
                                                method.trait           = "om",
                                                scale.invariance       = list(int = TRUE, lambda = TRUE),
                                                state.trait.invariance = FALSE,
                                                fixed.method.loadings  = TRUE,
                                                homocedasticity.assumption = list(error = TRUE, cs.red = TRUE, ut.red = FALSE))
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(mplus_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, mplus_syntax)
            
            # Run model in Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]        <- "timeout"
              times[1, a]       <- round(tf[3], 2)
              parameters[a, ]   <- -999
              se_psd[a, ]       <- -999 
              var.coeff[a, ]    <- -999 
              fit.measures[a, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # CS loadings
                                                             (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + 1,                                                                              # CS var
                                                             ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2 + I):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + 2 * I), # UT variances
                                                             (I * nT + 1):(I * nT + I),                                                                                                       # CT loadings
                                                             ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + 1):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + I),                       # intercepts
                                                             (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1,                                                                         # CT var
                                                             ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + I)),        # UT var
                                                           3]
                
                parameters[a, c(1:9, 11:19, 21:24)] <- estimates
                
                # Save standard errors
                stand.errors <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # CS loadings
                                                                (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + 1,                                                                              # CS var
                                                                ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2 + I):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + 2 * I), # UT variances
                                                                (I * nT + 1):(I * nT + I),                                                                                                       # CT loadings
                                                                ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + 1):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + I),                       # intercepts
                                                                (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1,                                                                         # CT var
                                                                ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + I)),        # UT var
                                                              4]
                
                se_psd[a, c(1:9, 11:19, 21:24)] <- stand.errors
                
                # Compute variance coefficients
                within.estimates  <- list(loadings = estimates[1:I], 
                                          CS.var   = estimates[I + 1],
                                          US.var   = estimates[(I + 2):(2 * I + 1)])
                between.estimates <- list(loadings = estimates[(2 * I + 2):(3 * I + 1)], 
                                          CT.var   = estimates[4 * I + 2], 
                                          UT.var   = estimates[(4 * I + 3):(5 * I + 2)])
                var.coeff[a, ] <- t(cuts.var.coeff(within.parameters  = within.estimates,
                                                   between.parameters = between.estimates))
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$ChiSqM_Value)) {
                  fit.measures[a, 2] <- c(fit$summaries$ChiSqM_Value)
                }
                if (!is.null(fit$summaries$ChiSqM_DF)) {
                  fit.measures[a, 3] <- c(fit$summaries$ChiSqM_DF)
                }
                if (!is.null(fit$summaries$ChiSqM_PValue)) {
                  fit.measures[a, 4] <- c(fit$summaries$ChiSqM_PValue)
                }
                if (!is.null(fit$summaries$CFI)) {
                  fit.measures[a, 5] <- c(fit$summaries$CFI)
                }
                if (!is.null(fit$summaries$TLI)) {
                  fit.measures[a, 6] <- c(fit$summaries$TLI)
                }
                if (!is.null(fit$summaries$RMSEA_Estimate)) {
                  fit.measures[a, 7] <- c(fit$summaries$RMSEA_Estimate)
                }
                if (!is.null(fit$summaries$AIC)) {
                  fit.measures[a, 8] <- c(fit$summaries$AIC)
                }
                if (!is.null(fit$summaries$BIC)) {
                  fit.measures[a, 9] <- c(fit$summaries$BIC)
                }
                if (!is.null(fit$summaries$aBIC)) {
                  fit.measures[a, 10] <- c(fit$summaries$aBIC)
                }
                
                rm(within.estimates, between.estimates, estimates, stand.errors)
                
                
              } else {
                perf[1, a]        <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
                times[1, a]       <- round(tf[3], 2)
                parameters[a, ]   <- -999
                se_psd[a, ]       <- -999 
                var.coeff[a, ]    <- -999 
                fit.measures[a, ] <- -999 
              }
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # cuts + long + maximun likelihood ----
            # Give analysis a number
            a <- 6
            
            file.name <- paste("cuts", "long", "ml", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.long)[-(1:2)],
                                                   cluster       = names(data$data.long)[1],
                                                   analysis_type = "TWOLEVEL",
                                                   estimator     = "ML",
                                                   iterations    = 50000)
            
            ml_syntax <- write.mlcuts.to.Mplus(data$data.long[, -(1:2)])
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(ml_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, ml_syntax)
            
            # Run modelin Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]        <- "timeout"
              times[1, a]       <- round(tf[3], 2)
              parameters[a, ]   <- -999
              se_psd[a, ]       <- -999 
              var.coeff[a, ]    <- -999 
              fit.measures[a, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name,".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[,3]
                
                parameters[a, c(1:9, 11:19, 21:24)] <- estimates
                
                # Save standard errors
                stand.errors <- fit$parameters$unstandardized[,4]
                
                se_psd[a, c(1:9, 11:19, 21:24)] <- stand.errors
                
                # Compute variance coefficients
                within.estimates  <- list(loadings = estimates[1:I], 
                                          CS.var   = estimates[I + 1],
                                          US.var   = estimates[(I + 2):(2 * I + 1)])
                between.estimates <- list(loadings = estimates[(2 * I + 2):(3 * I + 1)], 
                                          CT.var   = estimates[4 * I + 2], 
                                          UT.var   = estimates[(4 * I + 3):(5 * I + 2)])
                var.coeff[a, ] <- t(cuts.var.coeff(within.parameters  = within.estimates,
                                                   between.parameters = between.estimates))
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$ChiSqM_Value)) {
                  fit.measures[a, 2] <- c(fit$summaries$ChiSqM_Value)
                }
                if (!is.null(fit$summaries$ChiSqM_DF)) {
                  fit.measures[a, 3] <- c(fit$summaries$ChiSqM_DF)
                }
                if (!is.null(fit$summaries$ChiSqM_PValue)) {
                  fit.measures[a, 4] <- c(fit$summaries$ChiSqM_PValue)
                }
                if (!is.null(fit$summaries$CFI)) {
                  fit.measures[a, 5] <- c(fit$summaries$CFI)
                }
                if (!is.null(fit$summaries$TLI)) {
                  fit.measures[a, 6] <- c(fit$summaries$TLI)
                }
                if (!is.null(fit$summaries$RMSEA_Estimate)) {
                  fit.measures[a, 7] <- c(fit$summaries$RMSEA_Estimate)
                }
                if (!is.null(fit$summaries$AIC)) {
                  fit.measures[a, 8] <- c(fit$summaries$AIC)
                }
                if (!is.null(fit$summaries$BIC)) {
                  fit.measures[a, 9] <- c(fit$summaries$BIC)
                }
                if (!is.null(fit$summaries$aBIC)) {
                  fit.measures[a, 10] <- c(fit$summaries$aBIC)
                }
                
                rm(within.estimates, between.estimates, estimates, stand.errors)
                
              } else {
                perf[1, a]        <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]       <- round(tf[3], 2)
                parameters[a, ]   <- -999
                se_psd[a, ]       <- -999 
                var.coeff[a, ]    <- -999 
                fit.measures[a, ] <- -999 
              }
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # cuts + wide + Bayesian ----
            # Give analysis a number
            a <- 7
            
            file.name <- paste("cuts", "wide", "bayes", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.wide)[-1],
                                                   analysis_type = "GENERAL",
                                                   estimator     = "BAYES",
                                                   iterations    = 5000)
            
            
            mplus_syntax <- write.cuts.to.Mplus(data$data.wide[,-1],  
                                                nstate                 = nT,
                                                method.trait           = "om",
                                                scale.invariance       = list(int = TRUE, lambda = TRUE),
                                                state.trait.invariance = FALSE,
                                                fixed.method.loadings  = TRUE,
                                                homocedasticity.assumption = list(error = TRUE, cs.red = TRUE, ut.red = FALSE))
            saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", 
                                        "samples_", 
                                        file.name, 
                                        ".dat;", 
                                        "\nOUTPUT: TECH8;",
                                        "\nPLOT: TYPE = PLOT2;")
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(mplus_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(saveoutput_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, mplus_syntax, saveoutput_syntax)
            
            # Run model in Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]             <- "timeout"
              times[1, a]            <- round(tf[3], 2)
              parameters[a, ]        <- -999
              se_psd[a, ]            <- -999 
              var.coeff[a, ]         <- -999 
              fit.measures[a, ]      <- -999
              psd.var.coeff[a - 4, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # CS loadings
                                                             (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + 1,                                                                              # CS var
                                                             ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2 + I):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + 2 * I), # UT variances
                                                             (I * nT + 1):(I * nT + I),                                                                                                       # CT loadings
                                                             ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + 1):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + I),                       # intercepts
                                                             (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1,                                                                         # CT var
                                                             ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + I)),        # UT var
                                                           3]
                
                parameters[a, c(1:9, 11:19, 21:24)] <- estimates
                
                # Save posterior standard deviations
                post.sd <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # CS loadings
                                                           (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + 1,                                                                              # CS var
                                                           ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2 + I):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + 2 * I), # UT variances
                                                           (I * nT + 1):(I * nT + I),                                                                                                       # CT loadings
                                                           ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + 1):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 3) + I),                       # intercepts
                                                           (((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1,                                                                         # CT var
                                                           ((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 2):((((nT + I + 1) * (nT + I)) / 2) + (nT * I * 4) + nT + 1 + I)),        # UT var
                                                         4]
                
                se_psd[a, c(1:9, 11:19, 21:24)] <- post.sd
                
                # Compute variance coefficients
                
                # Put valid samples in an array 
                fit_samples <- as.array(fit$bparameters$valid_draw)
                fit_samples <- aperm(fit_samples, perm = c(1, 3, 2))
                fit_samples <- fit_samples[, , -(1:2)] # Take out chain and iteration number.
                
                # Rearrange samples in a matrix
                samples <- matrix(fit_samples, prod(dim(fit_samples)[1:2]), dim(fit_samples)[3])
                
                # Create matrix to store variance coefficients draws
                pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
                
                for(i in 1:dim(samples)[1]){
                  within.estimates  <- list(loadings = c(1, t(samples[i, c(5, 7, 9)])), 
                                            CS.var   = samples[i, 15],
                                            US.var   = c(t(samples[i, 11:14])))
                  between.estimates <- list(loadings = c(1, t(samples[i, c(6, 8, 10)])), 
                                            CT.var   = samples[i, 16], 
                                            UT.var   = c(t(samples[i, 17:20])))
                  pdist.var.coeff[i, ] <- t(cuts.var.coeff(within.parameters  = within.estimates,
                                                           between.parameters = between.estimates))
                }
                
                var.coeff[a, ] <- apply(pdist.var.coeff, 2, median)
                psd.var.coeff[a - 4, ] <- apply(pdist.var.coeff, 2, sd)
                
                # Save mcmc diagnosis plots
                pdf(file = paste0(folder, "Plots/", file.name, ".pdf"))
                if (dplots) {
                  mcmc_rhat(apply(fit_samples, 3, function (x) mcmcr::rhat(coda::as.mcmc(x))))
                  # Diagnosis within loadings
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(5, 7, 9)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(5, 7, 9)])
                  # Diagnosis between loadings 
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(6, 8, 10)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(6, 8, 10)])
                  # Diagnosis unique and common state variances
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[11:15])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[11:15])
                  # Diagnosis unique and common trait variances
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[16:20])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[16:20])
                  # Diagnosis intercepts
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[1:4])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[1:4])
                }
                dev.off()
                if (!dplots) {unlink(paste0(folder, "Plots/", file.name, ".pdf"))}
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$DIC)) {
                  fit.measures[a, 11] <- c(fit$summaries$DIC)
                }
                if (!is.null(fit$summaries$pD)) {
                  fit.measures[a, 12] <- c(fit$summaries$pD)
                }
                if (!is.null(fit$summaries$PostPred_PValue)) {
                  fit.measures[a, 13] <- c(fit$summaries$PostPred_PValue)
                }
                
                rm(within.estimates, between.estimates, estimates, post.sd, pdist.var.coeff, samples,i)
                
              } else {
                perf[1, a]             <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]            <- round(tf[3], 2)
                parameters[a, ]        <- -999
                se_psd[a, ]            <- -999 
                var.coeff[a, ]         <- -999 
                fit.measures[a, ]      <- -999
                psd.var.coeff[a - 4, ] <- -999
              }
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, file.name, ".gh5"))) {unlink(paste0(folder, file.name, ".gh5"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # cuts + long + Bayesian ----
            # Give analysis a number
            a <- 8
            
            file.name <- paste("cuts", "long", "bayes", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.long)[-(1:2)],
                                                   cluster       = names(data$data.long)[1],
                                                   analysis_type = "TWOLEVEL",
                                                   estimator     = "BAYES",
                                                   iterations    = 5000)
            
            ml_syntax <- write.mlcuts.to.Mplus(data$data.long[, -(1:2)])
            saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", 
                                        "samples_", 
                                        file.name, 
                                        ".dat;", 
                                        "\nOUTPUT: TECH8;",
                                        "\nPLOT: TYPE = PLOT2;")
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(ml_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(saveoutput_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, ml_syntax, saveoutput_syntax)
            
            # Run modelin Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]             <- "timeout"
              times[1, a]            <- round(tf[3], 2)
              parameters[a, ]        <- -999
              se_psd[a, ]            <- -999 
              var.coeff[a, ]         <- -999 
              fit.measures[a, ]      <- -999
              psd.var.coeff[a - 4, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[,3]
                
                parameters[a, c(1:9, 11:19, 21:24)] <- estimates
                
                # Save posterior standard deviations
                post.sd <- fit$parameters$unstandardized[,4]
                
                se_psd[a, c(1:9, 11:19, 21:24)] <- post.sd
                
                # Compute variance coefficients
                
                # Put valid samples in an array 
                fit_samples <- as.array(fit$bparameters$valid_draw)
                fit_samples <- aperm(fit_samples, perm = c(1, 3, 2))
                fit_samples <- fit_samples[, , -(1:2)] # Take out chain and iteration number.
                
                # Rearrange samples in a matrix
                samples <- matrix(fit_samples, prod(dim(fit_samples)[1:2]), dim(fit_samples)[3])
                
                # Create matrix to store variance coefficients draws
                pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
                
                for(i in 1:dim(samples)[1]){
                  within.estimates  <- list(loadings = c(1, t(samples[i, 1:3])), 
                                            CS.var   = samples[i, 8],
                                            US.var   = c(t(samples[i, 4:7])))
                  between.estimates <- list(loadings = c(1, t(samples[i, 13:15])), 
                                            CT.var   = samples[i, 20], 
                                            UT.var   = c(t(samples[i, 16:19])))
                  pdist.var.coeff[i, ] <- t(cuts.var.coeff(within.parameters  = within.estimates,
                                                           between.parameters = between.estimates))
                }
                
                var.coeff[a, ] <- apply(pdist.var.coeff, 2, median)
                psd.var.coeff[a - 4, ] <- apply(pdist.var.coeff, 2, sd)
                
                # Save mcmc diagnosis plots
                pdf(file = paste0(folder, "Plots/", file.name, ".pdf"))
                if (dplots) {
                  mcmc_rhat(apply(fit_samples, 3, function (x) mcmcr::rhat(coda::as.mcmc(x))))
                  # Diagnosis within loadings
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[1:3])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[1:3])
                  # Diagnosis between loadings 
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[13:15])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[13:14])
                  # Diagnosis unique and common state variances
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[4:8])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[4:8])
                  # Diagnosis unique and common trait variances
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[16:20])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[16:20])
                  # Diagnosis intercepts
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[9:12])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[9:12])
                }
                dev.off()
                if (!dplots) {unlink(paste0(folder, "Plots/", file.name, ".pdf"))}
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$DIC)) {
                  fit.measures[a, 11] <- c(fit$summaries$DIC)
                }
                if (!is.null(fit$summaries$pD)) {
                  fit.measures[a, 12] <- c(fit$summaries$pD)
                }
                if (!is.null(fit$summaries$PostPred_PValue)) {
                  fit.measures[a, 13] <- c(fit$summaries$PostPred_PValue)
                }
                
                rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
                
              } else {
                perf[1, a]             <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]            <- round(tf[3], 2)
                parameters[a, ]        <- -999
                se_psd[a, ]            <- -999 
                var.coeff[a, ]         <- -999 
                fit.measures[a, ]      <- -999
                psd.var.coeff[a - 4, ] <- -999
              }
              
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, file.name, ".gh5"))) {unlink(paste0(folder, file.name, ".gh5"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # tso + wide + maximum likelihood ----
            # Give analysis a number
            a <- 9
            
            file.name <- paste("tso", "wide", "ml", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.wide)[-1],
                                                   analysis_type = "GENERAL",
                                                   estimator     = "ML",
                                                   iterations    = 50000)
            
            ml_syntax <- write.tso.to.Mplus(data$data.wide[,-1], 
                                            nocc             = nT, 
                                            figure           = "3b",
                                            equiv.assumption = list(occ = "cong", theta = "equi"),
                                            scale.invariance = list(int = TRUE, lambda = TRUE),
                                            homocedasticity.assumption = list(error = TRUE, occ.red = TRUE),
                                            autoregressive.homogeneity = TRUE)
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(ml_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, ml_syntax)
            
            # Run modelin Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]        <- "timeout"
              times[1, a]       <- round(tf[3], 2)
              parameters[a, ]   <- -999
              se_psd[a, ]       <- -999 
              var.coeff[a, ]    <- -999 
              fit.measures[a, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # loadings
                                                             ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)),                                                                         # occasion variance
                                                             ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                             (2 * I * nT + 1),                                                                                                                # autoregressive effect
                                                             ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I),                 # Intercepts
                                                             ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I),         # trait indicator variances
                                                             ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)),    # trait indicator covariances 
                                                           3]
                
                parameters[a, c(1:10, 15:18, 21:30)] <- estimates
                
                # Save standard errors
                stand.errors <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # loadings
                                                                ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)),                                                                         # occasion variance
                                                                ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                                (2 * I * nT + 1),                                                                                                                # autoregressive effect
                                                                ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I),                 # Intercepts
                                                                ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I),         # trait indicator variances
                                                                ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)),    # trait indicator covariances 
                                                              4]
                
                se_psd[a, c(1:10, 15:18, 21:30)] <- stand.errors
                
                # Compute and save variance coefficients
                within.estimates  <- list(loadings      = estimates[1:I], 
                                          state.var     = estimates[I+1],
                                          error.var     = estimates[(I+2):(2 * I +1)], 
                                          ar.effect     = estimates[(2 * I + 2)])
                between.estimates <- list(trait.ind.var = estimates[(3 * I + 3):(4 * I + 2)])
                var.coeff[a, ] <- tso.var.coeff(I = I, nT = nT, 
                                                within.parameters  = within.estimates,
                                                between.parameters = between.estimates)[,nT]
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$ChiSqM_Value)) {
                  fit.measures[a, 2] <- c(fit$summaries$ChiSqM_Value)
                }
                if (!is.null(fit$summaries$ChiSqM_DF)) {
                  fit.measures[a, 3] <- c(fit$summaries$ChiSqM_DF)
                }
                if (!is.null(fit$summaries$ChiSqM_PValue)) {
                  fit.measures[a, 4] <- c(fit$summaries$ChiSqM_PValue)
                }
                if (!is.null(fit$summaries$CFI)) {
                  fit.measures[a, 5] <- c(fit$summaries$CFI)
                }
                if (!is.null(fit$summaries$TLI)) {
                  fit.measures[a, 6] <- c(fit$summaries$TLI)
                }
                if (!is.null(fit$summaries$RMSEA_Estimate)) {
                  fit.measures[a, 7] <- c(fit$summaries$RMSEA_Estimate)
                }
                if (!is.null(fit$summaries$AIC)) {
                  fit.measures[a, 8] <- c(fit$summaries$AIC)
                }
                if (!is.null(fit$summaries$BIC)) {
                  fit.measures[a, 9] <- c(fit$summaries$BIC)
                }
                if (!is.null(fit$summaries$aBIC)) {
                  fit.measures[a, 10] <- c(fit$summaries$aBIC)
                }
                
                rm(within.estimates, between.estimates, estimates, stand.errors)
                
              } else {
                perf[1, a]        <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]       <- round(tf[3], 2)
                parameters[a, ]   <- -999
                se_psd[a, ]       <- -999 
                var.coeff[a, ]    <- -999 
                fit.measures[a, ] <- -999 
              }
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # tso + wide + Bayesian ----
            # Give analysis a number
            a <- 10
            
            file.name <- paste("tso", "wide", "bayes", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.wide)[-1],
                                                   analysis_type = "GENERAL",
                                                   estimator     = "BAYES",
                                                   iterations    = 5000)
            
            ml_syntax <- write.tso.to.Mplus(data$data.wide[,-1], 
                                            nocc             = nT, 
                                            figure           = "3b",
                                            equiv.assumption = list(occ = "cong", theta = "equi"),
                                            scale.invariance = list(int = TRUE, lambda = TRUE),
                                            homocedasticity.assumption = list(error = TRUE, occ.red = TRUE),
                                            autoregressive.homogeneity = TRUE)
            saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", 
                                        "samples_", 
                                        file.name, 
                                        ".dat;", 
                                        "\nOUTPUT: TECH8;",
                                        "\nPLOT: TYPE = PLOT2;")
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(ml_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(saveoutput_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, ml_syntax, saveoutput_syntax)
            
            # Run modelin Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]             <- "timeout"
              times[1, a]            <- round(tf[3], 2)
              parameters[a, ]        <- -999
              se_psd[a, ]            <- -999 
              var.coeff[a, ]         <- -999 
              fit.measures[a, ]      <- -999
              psd.var.coeff[a - 5, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # loadings
                                                             ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)),                                                                         # occasion variance
                                                             ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                             (2 * I * nT + 1),                                                                                                                # autoregressive effect
                                                             ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I),                 # Intercepts
                                                             ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I),         # trait indicator variances
                                                             ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)),    # trait indicator covariances 
                                                           3]
                
                parameters[a, c(1:10, 15:18, 21:30)] <- estimates
                
                # Save posterior standard deviations
                post.sd <- fit$parameters$unstandardized[c(1:I,                                                                                                                             # loadings
                                                           ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)),                                                                         # occasion variance
                                                           ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 2 * I), # Error variances
                                                           (2 * I * nT + 1),                                                                                                                # autoregressive effect
                                                           ((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + I),                 # Intercepts
                                                           ((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I),         # trait indicator variances
                                                           ((2 * I * nT + nT) + (((nT + I) * (nT + I - 1) - I * (I - 1)) / 2)):((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) - 1)),    # trait indicator covariances 
                                                         4]
                
                se_psd[a, c(1:10, 15:18, 21:30)] <- post.sd
                
                # Compute and save variance coefficients
                
                # Put valid samples in an array 
                fit_samples <- as.array(fit$bparameters$valid_draw)
                fit_samples <- aperm(fit_samples, perm = c(1, 3, 2))
                fit_samples <- fit_samples[, , -(1:2)] # Take out chain and iteration number.
                
                # Rearrange samples in a matrix
                samples <- matrix(fit_samples, prod(dim(fit_samples)[1:2]), dim(fit_samples)[3])
                
                # Create matrix to store variance coefficients draws
                pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
                
                for(i in 1:dim(samples)[1]){
                  within.estimates  <- list(loadings  = c(1, t(samples[i, 5:7])), 
                                            state.var = samples[i, 13],
                                            error.var = c(t(samples[i, 8:11])), 
                                            ar.effect = samples[i, 12])
                  between.estimates <- list(trait.ind.var = c(t(samples[i, c(14, 16, 19, 23)])))
                  pdist.var.coeff[i, ] <- tso.var.coeff(I = I, nT = nT, 
                                                        within.parameters  = within.estimates,
                                                        between.parameters = between.estimates)[,nT]
                }
                
                var.coeff[a, ] <- apply(pdist.var.coeff, 2, median)
                psd.var.coeff[a - 5, ] <- apply(pdist.var.coeff, 2, sd)
                
                # Save mcmc diagnosis plots
                pdf(file = paste0(folder, "Plots/", file.name, ".pdf"))
                if (dplots) {
                  mcmc_rhat(apply(fit_samples, 3, function (x) mcmcr::rhat(coda::as.mcmc(x))))
                  # Diagnosis within loadings
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[5:7])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[5:7])
                  # Diagnosis error variances and state variance
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(9:11, 13)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(9:11, 13)])
                  # Diagnosis autoregressive effect
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[12])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[12])
                  # Diagnosis trait indicator variances
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(14, 16, 19, 23)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(14, 16, 19, 23)])
                  # Diagnosis intercepts
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[1:4])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[1:4])
                }
                dev.off()
                if (!dplots) {unlink(paste0(folder, "Plots/", file.name, ".pdf"))}
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$DIC)) {
                  fit.measures[a, 11] <- c(fit$summaries$DIC)
                }
                if (!is.null(fit$summaries$pD)) {
                  fit.measures[a, 12] <- c(fit$summaries$pD)
                }
                if (!is.null(fit$summaries$PostPred_PValue)) {
                  fit.measures[a, 13] <- c(fit$summaries$PostPred_PValue)
                }
                
                rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
                
              } else {
                perf[1, a]             <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]            <- round(tf[3], 2)
                parameters[a, ]        <- -999
                se_psd[a, ]            <- -999 
                var.coeff[a, ]         <- -999 
                fit.measures[a, ]      <- -999
                psd.var.coeff[a - 5, ] <- -999
              }
              rm(fit, a, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, file.name, ".gh5"))) {unlink(paste0(folder, file.name, ".gh5"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # tso + long + Bayesian ----
            # Give analysis a number
            a <- 11
            
            file.name <- paste("tso", "long", "bayes", cond, r, sep = "_")
            
            # Prepare data: Write data in Mplus format and write input file template
            prepareMplusData(data$data.long, paste0(folder, file.name, ".dat"), inpfile = T)
            
            # Complete Mplus syntax
            analysis_syntax <- write.Mplus.options(usevariables  = names(data$data.long)[-(1:2)],
                                                   cluster       = names(data$data.long)[1],
                                                   analysis_type = "TWOLEVEL",
                                                   estimator     = "BAYES",
                                                   iterations    = 5000)
            
            ml_syntax <- write.mltso.to.Mplus(data$data.long[, -c(1, 2)])
            saveoutput_syntax <- paste0("\nSAVEDATA: BPARAMETERS = ", 
                                        "samples_", 
                                        file.name, 
                                        ".dat;", 
                                        "\nOUTPUT: TECH8;",
                                        "\nPLOT: TYPE = PLOT2;")
            
            # Write analysis specifications in input file
            write(analysis_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(ml_syntax, paste0(folder, file.name, ".inp"), append = T)
            write(saveoutput_syntax, paste0(folder, file.name, ".inp"), append = T)
            
            rm(analysis_syntax, ml_syntax, saveoutput_syntax)
            
            # Run modelin Mplus
            assign("last.warning", NULL, envir = baseenv())
            t0  <- proc.time()
            run <- runModels_2(paste0(getwd(), "/", folder, file.name, ".inp"), timeout = timeout)
            tf  <- proc.time() - t0
            
            if (run == "timeout") {
              
              perf[1, a]             <- "timeout"
              times[1, a]            <- round(tf[3], 2)
              parameters[a, ]        <- -999
              se_psd[a, ]            <- -999 
              var.coeff[a, ]         <- -999 
              fit.measures[a, ]      <- -999
              psd.var.coeff[a - 5, ] <- -999
              assign("last.warning", NULL, envir = baseenv())
              
              rm(a, data, dataModel, nT, na.prop, seed, t0, tf)
            } else {
              fit <- readModels(paste0(getwd(), "/", folder, file.name, ".out")) #read Mplus output
              
              if (check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out")) == "Ok") {
                
                perf[1, a]  <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a] <- round(tf[3], 2)
                
                # Save estimates
                estimates <- fit$parameters$unstandardized[c(1:I,                   # loadings
                                                             (2 * I + 2),           # state var
                                                             ((I + 2):(2 * I + 1)), # error variances
                                                             (I + 1),               # autoregressive effect
                                                             (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2), # intercepts and indicator trait variances
                                                             (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)),                    # Covariances 
                                                           3]
                
                parameters[a, c(1:10, 15:18, 21:30)] <- estimates
                
                # Save posterior standard deviations
                post.sd <- fit$parameters$unstandardized[c(1:I,                   # loadings
                                                           (2 * I + 2),           # state var
                                                           ((I + 2):(2 * I + 1)), # error variances
                                                           (I + 1),               # autoregressive effect
                                                           (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2), # intercepts and indicator trait variances
                                                           (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)),                    # Covariances 
                                                         4]
                
                se_psd[a, c(1:10, 15:18, 21:30)] <- post.sd
                
                # Compute and save variance coefficients
                
                # Put valid samples in an array 
                fit_samples <- as.array(fit$bparameters$valid_draw)
                fit_samples <- aperm(fit_samples, perm = c(1, 3, 2))
                fit_samples <- fit_samples[, , -(1:2)] # Take out chain and iteration number.
                
                # Rearrange samples in a matrix
                samples <- matrix(fit_samples, prod(dim(fit_samples)[1:2]), dim(fit_samples)[3])
                
                # Create matrix to store variance coefficients draws
                pdist.var.coeff <- matrix(NA, dim(samples)[1], 20)
                
                for(i in 1:dim(samples)[1]){
                  within.estimates  <- list(loadings  = c(1, t(samples[i, 1:3])), 
                                            state.var = samples[i, 9],
                                            error.var = c(t(samples[i, 4:7])), 
                                            ar.effect = samples[i, 8])
                  between.estimates <- list(trait.ind.var = c(t(samples[i, c(14, 16, 19, 23)])))
                  pdist.var.coeff[i, ] <- tso.var.coeff(I = I, nT = nT, 
                                                        within.parameters  = within.estimates,
                                                        between.parameters = between.estimates)[,nT]
                }
                
                var.coeff[a, ] <- apply(pdist.var.coeff, 2, median)
                psd.var.coeff[a - 5, ] <- apply(pdist.var.coeff, 2, sd)
                
                # Save mcmc diagnosis plots
                pdf(file = paste0(folder, "Plots/", file.name, ".pdf"))
                if (dplots) {
                  mcmc_rhat(apply(fit_samples, 3, function (x) mcmcr::rhat(coda::as.mcmc(x))))
                  # Diagnosis within loadings
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[1:3])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[1:3])
                  # Diagnosis error variances and state variance
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(4:7, 9)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(4:7, 9)])
                  # Diagnosis autoregressive effect
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[8])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[8])
                  # Diagnosis trait indicator variances
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[c(14, 16, 19, 23)])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[c(14, 16, 19, 23)])
                  # Diagnosis intercepts
                  mcmc_trace(fit_samples, pars = dimnames(fit_samples)$var[10:13])
                  mcmc_acf(fit_samples, pars = dimnames(fit_samples)$var[10:13])
                }
                dev.off()
                if (!dplots) {unlink(paste0(folder, "Plots/", file.name, ".pdf"))}
                
                # Save fit measures
                if (!is.null(fit$summaries$Parameters)) {
                  fit.measures[a, 1] <- c(fit$summaries$Parameters)
                }
                if (!is.null(fit$summaries$DIC)) {
                  fit.measures[a, 11] <- c(fit$summaries$DIC)
                }
                if (!is.null(fit$summaries$pD)) {
                  fit.measures[a, 12] <- c(fit$summaries$pD)
                }
                if (!is.null(fit$summaries$PostPred_PValue)) {
                  fit.measures[a, 13] <- c(fit$summaries$PostPred_PValue)
                }
                
                rm(within.estimates, between.estimates, estimates, post.sd, samples, pdist.var.coeff, i)
                
              } else {
                perf[1, a]             <- check.mplus(fit, paste0(getwd(), "/", folder, file.name, ".out"))
                times[1, a]            <- round(tf[3], 2)
                parameters[a, ]        <- -999
                se_psd[a, ]            <- -999 
                var.coeff[a, ]         <- -999 
                fit.measures[a, ]      <- -999
                psd.var.coeff[a - 5, ] <- -999
              }
              
              rm(fit, a, data, dataModel, nT, na.prop, seed, t0, tf)
            }
            
            # Remove Mplus files of previous analysis 
            if (file.exists(paste0(folder, file.name, ".dat"))) {unlink(paste0(folder, file.name, ".dat"))}
            if (file.exists(paste0(folder, file.name, ".inp"))) {unlink(paste0(folder, file.name, ".inp"))}
            if (file.exists(paste0(folder, file.name, ".out"))) {unlink(paste0(folder, file.name, ".out"))}
            if (file.exists(paste0(folder, file.name, ".gh5"))) {unlink(paste0(folder, file.name, ".gh5"))}
            if (file.exists(paste0(folder, "samples_", file.name, ".dat"))) {unlink(paste0(folder, "samples_", file.name, ".dat"))}
            rm(file.name, run)
            
            # Export output to outer foreach ----
            result <- list(perf, times, parameters, se_psd, var.coeff, psd.var.coeff, fit.measures)
            rm(perf, times, parameters, se_psd, var.coeff, psd.var.coeff, fit.measures)
            result
            
  }
time <- proc.time() - time0

# 4.0 Save output tables ----

for(i in args[1]:args[2]){
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[1]]
  write.table(tmp, file = paste0(folder, "Performance/", paste("performance", i, sep = "_"), ".dat"),
              col.names = TRUE, row.names = FALSE, quote = TRUE)
  rm(tmp)
  
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[2]]
  write.table(tmp, file = paste0(folder, "Times/", paste("times", i, sep = "_"), ".dat"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  rm(tmp)
  
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[3]]
  write.table(tmp, file = paste0(folder, "Parameters/", paste("parameters", i, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  rm(tmp)
  
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[4]]
  write.table(tmp, file = paste0(folder, "SE_PSD/", paste("se_psd", i, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  rm(tmp)
  
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[5]]
  write.table(tmp, file = paste0(folder, "Var_Coeff/", paste("var_coeff", i, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  rm(tmp)
  
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[6]]
  write.table(tmp, file = paste0(folder, "PSD_Var_Coeff/", paste("psd_var_coeff", i, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  rm(tmp)
  
  tmp           <- outcome.simulation[[i - as.numeric(args[1]) + 1]][[7]]
  write.table(tmp, file = paste0(folder, "Fit_Measures/", paste("fit_measures", i, sep = "_"), ".dat"),
              col.names = TRUE, row.names = TRUE, quote = FALSE)
  rm(tmp)
  
}

# 5.0 Stop cluster and clean enviroment ----
stopCluster(cl)

rm(list = ls())
