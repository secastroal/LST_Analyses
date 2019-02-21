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

model <- "cuts" # model to simulate data
N <- 100 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
na.prop <- 0 #proportion of missingness
seed <- 123

set.seed(seed)

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
}


# 3.0 Matrices to store output ----

# 4.0 Model estimation ----

# 4.1 cuts ----
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
  status[1] <- check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out"))
}


rm(file.name, long.fit)

# 4.2 msst ----

file.name <- paste0(model, "_long_n", N, "_i", I, "_nt", nT, "_na", na.prop)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(msst.data$data.long, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                       cluster = names(data$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000,
                                       processors = 2)

ml_syntax <- write.mlmsst.to.Mplus(msst.data$data.long[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
cat("\n"); print(Sys.time()); cat("\n")
runModels(paste0(getwd(),"/", folder,file.name,".inp"))
cat("\n"); print(Sys.time()); cat("\n")

long.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(long.fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  
  status[1] <- check.mplus(long.fit, paste0(getwd(),"/",folder,file.name,".out"))
  
  est.par[,3] <- long.fit$parameters$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                                      (4*I + 2), (6*I +3)),3]
  se[,2] <- long.fit$parameters$unstandardized[c(1:(2*I+1), ((3*I + 2):(4*I + 1)), ((5*I + 3):(6*I + 2)),
                                                 (4*I + 2), (6*I +3)),4]
  
  bias[, 2] <- (est.par[ , 3] - est.par[ , 2]) 
  
  rmse[, 2] <- sqrt((est.par[ , 3] - est.par[ , 2]) ^ 2)
  
  within.estimates <- list( loadings = est.par[1:I, 3], state.var = est.par[I+1, 3],
                            error.var = est.par[(I+2):(2 * I +1), 3])
  between.estimates <- list( loadings = est.par[(2 * I + 2):(3 * I + 1), 3], 
                             trait.var = est.par[ 4 * I + 3, 3])
  var.coeff[ , 3] <- msst.var.coeff(within.parameters = within.estimates,
                                    between.parameters = between.estimates)
  rm(within.estimates, between.estimates)
  
  
}else{
  status[1] <- check.mplus(long.fit, paste0(getwd(),"/",folder,file.name,".out"))
}

rm(file.name, long.fit)

# 4.3 tso ----
file.name <- paste0(model, "_long_n", N, "_i", I, "_nt", nT, "_na", na.prop)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(tso.data$data.long, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data$data.long)[-(1:2)],
                                       cluster = names(data$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000,
                                       processors = 2)

ml_syntax <- write.mltso.to.Mplus(tso.data$data.long[, -c(1, 2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
cat("\n"); print(Sys.time()); cat("\n")
runModels(paste0(getwd(),"/",folder,file.name,".inp"))
cat("\n"); print(Sys.time()); cat("\n")

long.fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(long.fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  
  status[1] <- check.mplus(long.fit, paste0(getwd(),"/",folder,file.name,".out"))
  
  est.par[, 3] <- long.fit$parameters$unstandardized[c(1:(2*I + 2),
                                                       (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2),
                                                       (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)), 3]
  
  post.sd[, 2] <- long.fit$parameters$unstandardized[c(1:(2*I + 2),
                                                       (4 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2),
                                                       (3 * I + 3):(3 * I + (I * (I - 1) / 2) + 2)), 4]
  
  bias[, 2] <- (est.par[ , 3] - est.par[ , 2]) 
  
  rmse[, 2] <- sqrt((est.par[ , 3] - est.par[ , 2]) ^ 2)
  
  within.estimates <- list( loadings = est.par[1:I, 3], state.var = est.par[2*I+2, 3],
                            error.var = est.par[(I+2):(2 * I +1), 3], ar.effect = est.par[I+1,3])
  between.estimates <- list( trait.ind.var = est.par[(3 * I + 3):(4 * I + 2), 3])
  var.coeff[ , 3] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.estimates,
                                   between.parameters = between.estimates)[,nT]
  rm(within.estimates, between.estimates)
  
  
}else{
  status[1] <- check.mplus(long.fit, paste0(getwd(),"/",folder,file.name,".out"))
}

rm(file.name, long.fit)




# 1.3 Model estimation ---- 

# Matrix to store estimated parameters

est.par <- data.frame(matrix(NA, 5*I + 2, 6))
est.par[, 2] <- c(unlist(within.parameters), unlist(between.parameters))
est.par[, 1] <- names(c(unlist(within.parameters), unlist(between.parameters)))
names(est.par) <- c("par", "true", "long", "wide", "Ltrunc", "Wtrunc" )

# Matrix to store rmse and mse
se <- rmse <- bias <- est.par[, -2] 

# Matrix to store variance coefficients

var.coeff <- data.frame(matrix(NA, I * 5, 6))
var.coeff[ , 1] <- paste0(rep(c( "ccon_y", "ucon_y", "tcon_y", "spe_y", "rel_y"), each = I), 1:I)
var.coeff[ , 2] <- cuts.var.coeff(within.parameters = within.parameters, 
                                between.parameters = between.parameters)
names(var.coeff) <- c("coefficient", "true", "long", "wide", "Ltrunc", "Wtrunc" )

# Save convergence check

status <- rep(NA, 4)

# 1.3.1 ML-cuts data.long ----


# 1.4 Save output ----

est.par[dim(est.par)[1] + 1, ] <- c("Status", "NA", status)

est.par <- est.par[c(dim(est.par)[1], 1:(dim(est.par)[1]-1)), ]

file.name <- paste0(model, "_n", N, "_i", I, "_nt", nT, "_na", na.prop, ".txt")
caption <- paste("of model", model, "with N =", N,", I =", I, ", nT =", nT,
                 ", and missingness proportion =", na.prop, ".", sep = " ")

print(xtable(est.par, type = "latex", caption = paste0("Parameter estimates ", caption)), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/param_", file.name))
print(xtable(bias, type = "latex", caption = paste0("Bias ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/bias_", file.name ))
print(xtable(rmse, type = "latex", caption = paste0("RMSE ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/rmse_", file.name))
print(xtable(se, type = "latex", caption = paste0("Standard errors ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/se_", file.name))
print(xtable(var.coeff, type = "latex", caption = paste0("LST variance coefficients ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/varcoeff_", file.name))


# 1.5 Clean environment ----

rm(within.parameters, between.parameters, cuts.data, est.par, model, rmse, 
   se, status, var.coeff, file.name, caption, bias)



# 2.3 Model estimation ----

# Matrix to store estimated parameters

est.par <- data.frame(matrix(NA, 4*I + 3, 6))
est.par[, 2] <- c(unlist(within.parameters), unlist(between.parameters))
est.par[, 1] <- names(c(unlist(within.parameters), unlist(between.parameters)))
names(est.par) <- c("par", "true", "long", "wide", "Ltrunc", "Wtrunc")

# Matrix to store rmse and mse
se <- rmse <- bias <- est.par[, -2]

# Matrix to store variance coefficients

var.coeff <- data.frame(matrix(NA, I * 3, 6))
var.coeff[ , 1] <- paste0(rep(c("con_y", "spe_y", "rel_y"), each = I), 1:I)
var.coeff[ , 2] <- msst.var.coeff(within.parameters = within.parameters, 
                                  between.parameters = between.parameters)
names(var.coeff) <- c("coefficient", "true", "long", "wide", "Ltrunc", "Wtrunc" )


# Save convergence check

status <- rep(NA, 4)

# 2.3.1 ML-msst data.long ----


# 2.4 Save output ----

est.par[dim(est.par)[1] + 1, ] <- c("Status", "NA", status)

est.par <- est.par[c(dim(est.par)[1], 1:(dim(est.par)[1]-1)), ]

file.name <- paste0(model, "_n", N, "_i", I, "_nt", nT, "_na", na.prop, ".txt")
caption <- paste("of model", model, "with N =", N,", I =", I, ", nT =", nT,
                 ", and missingness proportion =", na.prop, ".", sep = " ")

print(xtable(est.par, type = "latex", caption = paste0("Parameter estimates ", caption)), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/param_", file.name))
print(xtable(bias, type = "latex", caption = paste0("Bias ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/bias_", file.name ))
print(xtable(rmse, type = "latex", caption = paste0("RMSE ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/rmse_", file.name))
print(xtable(se, type = "latex", caption = paste0("Standard errors ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/se_", file.name))
print(xtable(var.coeff, type = "latex", caption = paste0("LST variance coefficients ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/varcoeff_", file.name))

# 2.5 Clean environment ----

rm(within.parameters, between.parameters, msst.data, est.par, model, rmse, 
   se, status, var.coeff, file.name, caption, bias)



# 3.3 Model estimation ----

# Matrix to store estimated parameters

est.par <- data.frame(matrix(NA, 4*I + 2 + (I * (I-1) /2), 6))
# Get true parameters + lower triangle of the true variance-covariance matrix 
est.par[, 2] <- round(c(unlist(within.parameters), unlist(between.parameters)[1:(2*I)],
                  tso.data$between.parameters$Sigma[t(lower.tri(tso.data$between.parameters$Sigma))]),2)
# Get array indices of the variance-covariance matrix
cov.ind <- which(t(lower.tri(tso.data$between.parameters$Sigma))==TRUE, arr.ind = TRUE)
# Define parameter names + covariance names
est.par[, 1] <- c(names(c(unlist(within.parameters), unlist(between.parameters)[1:(2*I)])),
                  paste0("cov", cov.ind[,1], cov.ind[,2]))
names(est.par) <- c("par", "true", "long", "wide", "Ltrunc", "Wtrunc")

# Matrix to store rmse
post.sd <- rmse <- bias <- est.par[, -2] 

rm(cov.ind)

# Matrix to store MEAN variance coefficients

var.coeff <- data.frame(matrix(NA, I * 5, 6))
var.coeff[ , 1] <- paste0(rep(c( "pred_y", "upred_y", "con_y", "spe_y", "rel_y"), each = I), 1:I)
var.coeff[ , 2] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.parameters, 
                                  between.parameters = between.parameters)[, nT]
names(var.coeff) <- c("coefficient", "true", "long", "wide", "Ltrunc", "Wtrunc" )

# Save convergence check

status <- rep(NA, 4)

# 3.3.1 ML-tso data.long ----


# 3.4 Save output ----

est.par[dim(est.par)[1] + 1, ] <- c("Status", "NA", status)

est.par <- est.par[c(dim(est.par)[1], 1:(dim(est.par)[1]-1)), ]

file.name <- paste0(model, "_n", N, "_i", I, "_nt", nT, "_na", na.prop, ".txt")
caption <- paste("of model", model, "with N =", N,", I =", I, ", nT =", nT,
                 ", and missingness proportion =", na.prop, ".", sep = " ")

print(xtable(est.par, type = "latex", caption = paste0("Parameter estimates ", caption)), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/param_", file.name))
print(xtable(bias, type = "latex", caption = paste0("Bias ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/bias_", file.name ))
print(xtable(rmse, type = "latex", caption = paste0("RMSE ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/rmse_", file.name))
print(xtable(post.sd, type = "latex", caption = paste0("Posterior standard deviations ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/postsd_", file.name))
print(xtable(var.coeff, type = "latex", caption = paste0("LST variance coefficients ", caption), digits = 3), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/varcoeff_", file.name))
# 3.5 Clean environment ----

rm(N, nT, I, seed, within.parameters, between.parameters, tso.data, est.par, model, rmse,  
   bias, post.sd, var.coeff, caption, file.name, na.prop, status)

