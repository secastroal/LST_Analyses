# Contents
# 0.0 Prepare environment
# 1.0 CUTS
# 2.0 MSST
# 3.0 TSO


# 0.0 Prepare environment ----
rm(list=ls())
#library(devtools) # if lsttheory of sumplement C of Steyer et al. (2015) has not been installed before.
#install_github("amayer2010/lsttheory", force = TRUE)
library(lavaan)
library(lsttheory)
library(MplusAutomation)
library(MASS)
library(xtable)
source("R/write.Mplus.options.R")
source("R/write.cuts.to.Mplus.R")
source("R/write.mlcuts.to.Mplus.R")
source("R/write.msst.to.Mplus.R")
source("R/write.mlmsst.to.Mplus.R")
source("R/write.tso.to.Mplus.R")
source("R/write.mltso.to.Mplus.R")
source("R/sim.data.cuts.R")
source("R/sim.data.msst.R")
source("R/sim.data.tso.R")
source("R/check.Mplus.R")
source("R/var.coeff.R")
source("R/var.coeff.tv.R")
source("R/runModels_2.R")
folder <- "Mplus_files_Cross/" #Folder to store results

# 1.0 Set Conditions ----

model <- "tso" # model to simulate data
N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
na.prop <- 0 #proportion of missingness
seed <- 123

set.seed(seed)

# 1.1 Matrices to store output ----
# Matrix of variance coefficients

var.coeff <- data.frame(matrix(NA, I * 5, 5))
var.coeff[ , 1] <- paste0(rep(c( "ccon/pred_y", "ucon/upred_y", "tcon/con_y", "spe_y", "rel_y"), each = I), 1:I)
names(var.coeff) <- c("coefficient", "true", "CUTS", "MSST", "TSO")


# Matrix of fit measures
fit.measures <- data.frame(matrix(NA, 3, 4))
names(fit.measures) <- c("Data","Parameters","DIC","pD")
fit.measures[,1] <- c("CUTS", "MSST", "TSO")

# 2.0 Set true paramaters and simulate data given the model ----
# 2.1 simulate CUTS data ----
if(model == "cuts"){
  
  # Within Parameters
  
  state_loadings <- c(1, 0.5, 1.3, 0.8, 1.1)[1:I] # loading parameters for the latent common state
  var_CS <- 2 # Variance latent common state
  var_US <- c(1, 0.5, 1.5, 0.8, 1.2)[1:I] # Variance of latent unique states
  
  within.parameters <- list(loadings = state_loadings, CS.var = var_CS, US.var = var_US)
  
  rm(state_loadings, var_CS, var_US)
  
  # Between Paramaters
  
  intercepts <- seq(2, by = 0.5, length.out = I) # intercepts
  trait_loadings <- c(1, 0.8, 1.2, 0.9, 1.1)[1:I] # loading parametes for the latent common trait
  var_CT <- 1.5 # variance latent common trait
  var_UT <- c(0.5, 1, 0.3, 0.8, 0.5)[1:I] # variance latent unique traits
  
  between.parameters <- list(loadings = trait_loadings, intercepts = intercepts, CT.var = var_CT, UT.var = var_UT)
  
  rm(intercepts, trait_loadings, var_CT, var_UT)
  
  # Simulate data
  
  # Time invariant data
  data <- sim.data.cuts(N, nT, I, within.parameters = within.parameters, na.prop = na.prop,
                             between.parameters = between.parameters, seed = seed)
  
  # Compute true variance coefficients
  
  var.coeff[,2] <- cuts.var.coeff(within.parameters = within.parameters, 
                                  between.parameters = between.parameters) 
  
  rm(within.parameters, between.parameters)
}
# 2.2 simulate MSST data ----
if(model == "msst"){
  # Within Parameters
  
  state_loadings <- c(1, 0.5, 1.3, 0.8, 1.1)[1:I] # loading parameters for the latent state
  var_state <- 2 # Variance latent state residual
  var_m_error <- c(1, 0.5, 1.5, 0.8, 1.2)[1:I] # Variance of measurement errors
  
  within.parameters <- list(loadings = state_loadings, state.var = var_state, error.var = var_m_error)
  
  rm(state_loadings, var_state, var_m_error)
  
  # Between Paramaters
  
  intercepts <- seq(0, by = 0.2, length.out = I) # intercepts
  trait_loadings <- c(1, 0.8, 1.2, 0.9, 1.1)[1:I] # loading parametes for the latent trait
  
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
  
  var.coeff[(2 * I + 1):(5 * I),2] <- msst.var.coeff(within.parameters = within.parameters, 
                                  between.parameters = between.parameters) 
  
  rm(within.parameters, between.parameters)
}
# 2.3 simulate TSO data ----
if(model == "tso"){
  # Within Parameters
  
  state_loadings <- c(1, 0.5, 1.3, 0.8, 1.1)[1:I] # loading parameters for the latent state
  var_state <- 2 # Variance latent state residual
  var_error <- c(1, 0.5, 1.5, 0.8, 1.2)[1:I] # Variance of latent measurement errors
  
  ar_effect <- 0.5 # autoregressive effect on the latent state residuals
  
  within.parameters <- list(loadings = state_loadings, ar.effect = ar_effect, error.var = var_error,
                            state.var = var_state)
  
  rm(state_loadings, var_state, var_error, ar_effect)
  
  # Between Paramaters 
  
  intercepts <- seq(2, by = 0.5, length.out = I) # intercepts
  
  var_ind_traits <- c(2, 1.5, 2.5, 1.75, 2.25)[1:I] # variance latent indicator trait variables
  
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
  
  var.coeff[,2] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.parameters, 
                                  between.parameters = between.parameters)[,nT] 
  
  rm(within.parameters, between.parameters)
}


# 3.0 Model estimation ----

# 3.1 cuts ----
file.name <- paste0(model, "dataXcuts_n", N, "_i", I, "_nt", nT, "_na", na.prop)

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
runModels(paste0(getwd(),"/",folder,file.name,".inp"))
cat("\n"); print(Sys.time()); cat("\n")

fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  
  estimates <- fit$parameters$unstandardized[,3]
  
  within.estimates <- list( loadings = estimates[1:I], CS.var = estimates[I+1],
                            US.var = estimates[(I+2):(2 * I +1)])
  between.estimates <- list( loadings = estimates[(2 * I + 2):(3 * I + 1)], CT.var = estimates[ 4 * I + 2], 
                             UT.var = estimates[(4 * I + 3):(5 * I + 2)])
  var.coeff[ , 3] <- cuts.var.coeff(within.parameters = within.estimates,
                                    between.parameters = between.estimates)
  fit.measures[1, 2:4] <- fit$summaries[, c(11:13)]
  rm(within.estimates, between.estimates)
  
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}


rm(file.name, fit, estimates)

# 3.2 msst ----

file.name <- paste0(model, "dataXmsst_n", N, "_i", I, "_nt", nT, "_na", na.prop)

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

ml_syntax <- gsub("@0;", "@0.001;", ml_syntax)

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
  var.coeff[(2 * I + 1):(5 * I), 4] <- msst.var.coeff(within.parameters = within.estimates,
                                    between.parameters = between.estimates)
  fit.measures[2, 2:4] <- fit$summaries[, c(11:13)]
  rm(within.estimates, between.estimates)
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

rm(file.name, fit, estimates)

# 3.3 tso ----
file.name <- paste0(model, "dataXtso_n", N, "_i", I, "_nt", nT, "_na", na.prop)

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
  var.coeff[ , 5] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
                                   between.parameters = between.estimates)[,nT]
  fit.measures[3, 2:4] <- fit$summaries[, c(11:13)]
  rm(within.estimates, between.estimates)
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

rm(file.name, fit, estimates)

# 4.0 Plot variance coefficients and save fit measures ----

jpeg(paste0(getwd(), "/", folder, "OutputPlots/varcoeff_", model, "data_n", N, "_i", I, "_nt", nT, ".jpg" ))
plot(c(seq(1,1.16, by =0.04), seq(2,2.16, by =0.04), seq(3,3.16, by =0.04)), unlist(var.coeff[seq(3,20, by = I), 3:5]), ylim = c(0,1),
     col = c("green", "orange", "red", "blue", "black"), xlab = "Model", ylab = "Explained Variance", pch = 19, xlim = c(0.75, 3.75), 
     main = paste0(model, " data: Variance coefficients."),xaxt = 'n')
abline(h = var.coeff[seq(3,20, by = I),2], col = c("green", "orange", "red", "blue", "black"), lty = 5:1)
abline(h = 0)
axis(side = 1, at = 1.1:3.1, labels = c("CUTS", "MSST", "TSO"))
legend("bottomright", legend = c("CCon/Pred", "UCon/Upred", "TCon/Con", "Spe", "Rel"), col = c("green", "orange", "red", "blue", "black"),
       pch = 19, cex = 0.8)
dev.off()


print(xtable(fit.measures, type = "latex", caption = paste0("Fit measures: ", model, "data with N =", N, ", I = ", I, ", and nT =", nT, ".")), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/FitMeasures_", model, "data_n", N, "_i", I, "_nt", nT, ".txt"))



# 5.0 Clean environment ----

rm(N, nT, I, seed, data, model, fit.measures, var.coeff, na.prop)

