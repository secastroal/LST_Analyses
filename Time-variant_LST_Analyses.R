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
folder <- "Mplus_files_TV/" #Folder to store results

# 1.0 CUTS ----

# 1.1 Set Conditions and true parameters ----

model <- "cuts"
N <- 200 # number of persons
nT <- 30 # number of times // measurement occasions
I <- 4 # number of variables // items
seed <- 123

set.seed(seed)

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

# 1.2 Simulate data ----

# Time variant data
data.tv <- sim.data.cuts.tv(N, nT, I, within.parameters = within.parameters, time.invariant = FALSE, 
                                 between.parameters = between.parameters, seed = seed)

rm(within.parameters, between.parameters)

# 1.3 Model estimation ----

# Compute true variance coefficients

true.var.coeff <- cuts.var.coeff.tv(I, nT, within.parameters = data.tv$within.parameters, 
                                    between.parameters = data.tv$between.parameters)
# Matrix to save fit measures

fit.measures <- data.frame(matrix(NA, 2, 12))
names(fit.measures) <- c("Data","LL","Parameters","ChiSqM_Value","ChiSqM_DF","ChiSqM_PValue",      
                         "CFI","TLI","AIC",                 
                         "BIC","aBIC","RMSEA_Estimate")
fit.measures[,1] <- c("Long", "Wide")

# 1.3.1 ML-cuts with time-variant data ----

file.name <- paste0(model, "long_tv_n", N, "_i", I, "_nt", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(data.tv$data.long, paste0(folder, file.name, ".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data.tv$data.long)[-(1:2)],
                                       cluster = names(data.tv$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlcuts.to.Mplus(data.tv$data.long[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  estimates <- fit$parameters$unstandardized[,3]
  
  est.within <- list(loadings = estimates[1:I], #within loadings
                     CS.var = estimates[I + 1], #CS.var
                     US.var = estimates[(I + 2):(2 * I + 1)]) # US.var
  
  est.between <- list(loadings = estimates[(2 * I + 2):(3 * I + 1)], #loadings
                      intercepts = estimates[(3 * I + 2):(4 * I + 1)], # intercepts
                      CT.var = estimates[(4 * I + 2)], # CT.var
                      UT.var = estimates[(4 * I + 3):(5 * I + 2)]) # UT.var
  
  est.var.coeff.long <- cuts.var.coeff(within.parameters = est.within,
                                       between.parameters = est.between)
  
  fit.measures[1, 2:12] <- fit$summaries[, c(18, 11:14, 20:25)]
  
  rm(file.name, fit, estimates, est.between, est.within)
  
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

# 1.3.2 cuts data.wide ("Perfect fit") ----

file.name <- paste0(model, "_wide_tv_n", N, "_i", I, "_nt", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(data.tv$data.wide, paste0(folder, file.name, ".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data.tv$data.wide)[-1],
                                       analysis_type = "GENERAL",
                                       estimator = "ML",
                                       iterations = 50000,
                                       h1iterations = 50000)


mplus_syntax <- write.cuts.to.Mplus(data.tv$data.wide[,-1],  nstate = nT,
                                    method.trait = "om",
                                    scale.invariance = list(int = TRUE, lambda = TRUE),
                                    state.trait.invariance = FALSE,
                                    fixed.method.loadings = TRUE,
                                    homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE))


write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, mplus_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  estimates <- fit$parameters$unstandardized[ ,3]
  
  wloadings <- estimates[1:(I*nT)] # within loadings
  CS.var <- estimates[((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT)] #CS variances
  US.var <- estimates[((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2 + I):((((nT + I + 1)*(nT + I))/2) + (nT*I*5) + nT + 1 + I)] #US variances
  bloadings <- estimates[(I*nT + 1):(I*nT*2)] #between loadings CT
  loadings.UT <- estimates[(I*nT*2 + 1):(I*nT*3)] #between loadings UT
  intercepts <- estimates[((((nT + I + 1)*(nT + I))/2) + (nT*I*3) + 1):((((nT + I + 1)*(nT + I))/2) + (nT*I*4))] #intercepts
  CT.var <- estimates[(((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1] # CT variance
  UT.var <- estimates[((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 2):((((nT + I + 1)*(nT + I))/2) + (nT*I*4) + nT + 1 + I)] #UT variances
  
  est.within <- list(loadings = matrix(wloadings, nrow = nT, ncol = I, byrow = TRUE),
                     CS.var = CS.var,
                     US.var = matrix(US.var,  nrow = nT, ncol = I, byrow = TRUE))
  est.between <- list(loadings = matrix(bloadings, nrow = nT, ncol = I, byrow = TRUE),
                      intercepts = matrix(intercepts, nrow = nT, ncol = I, byrow = TRUE),
                      loadings.UT = matrix(loadings.UT, nrow = nT, ncol = I),
                      CT.var = CT.var,
                      UT.var = UT.var)
  
  est.var.coeff.wide <- cuts.var.coeff.tv(I=I, nT= nT, within.parameters = est.within,
                                          between.parameters = est.between)
  
  fit.measures[2, 2:12] <- fit$summaries[, c(18, 11:14, 20:25)]
  
  rm(file.name, fit, estimates, est.between, est.within, wloadings, CS.var, US.var, 
     bloadings, loadings.UT, intercepts, CT.var, UT.var)
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

# 1.4 Plot variance coefficients and save fit measures ---- 

jpeg(paste0(getwd(), "/", folder, "OutputPlots/varcoeff_", model, "_n", N, "_i", I, "_nt", nT, ".jpg" ))
matplot(cbind(t(true.var.coeff[seq(1,20, by = I),]), t(est.var.coeff.wide[seq(1,20, by = I),])), type = c("l"), ylim = c(0, 1),
        col = c("green", "orange", "red", "blue", "black"), xlab = "Time", ylab = "Explained Variance", lty = rep(1:2, each = 5), lwd = 2,
        main = paste0(model, ": Variance coefficients."))
abline(h = est.var.coeff.long[seq(1,20, by = I),], col = c("green", "orange", "red", "blue", "black"), lty = 3)
legend("bottomright", legend = c("CCon", "UCon", "TCon", "Spe", "Rel"), col = c("green", "orange", "red", "blue", "black"),
       lty = 1, lwd = 2, pch = c(NA,NA,16,16,16), cex = 0.8)
dev.off()

print(xtable(fit.measures, type = "latex", caption = paste0("Fit measures: multi-level and single-level ", model, " N =", N, ", I = ", I, ", and nT =", nT, ".")), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/FitMeasures_", model, "_n", N, "_i", I, "_nt", nT, ".txt"))


# 1.5 Clean environment ----

rm(data.tv, est.var.coeff.long, est.var.coeff.wide, true.var.coeff, fit.measures, model)


# 2.0 MSST ----

# 2.1 Set Conditions and true parameters ----

model <- "msst" 
set.seed(seed)

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

# 2.2 Simulate data ----

# Time variant data

data.tv <- sim.data.msst.tv(N, nT, I,  within.parameters = within.parameters, 
                                 between.parameters = between.parameters, time.invariant = FALSE, seed = seed)
rm(within.parameters, between.parameters)

# 2.3 Model estimation ----

# Compute true variance coefficients

true.var.coeff <- msst.var.coeff.tv(I, nT, within.parameters = data.tv$within.parameters, 
                                    between.parameters = data.tv$between.parameters)
# Matrix to save fit measures

fit.measures <- data.frame(matrix(NA, 2, 12))
names(fit.measures) <- c("Data","LL","Parameters","ChiSqM_Value","ChiSqM_DF","ChiSqM_PValue",      
                         "CFI","TLI","AIC",                 
                         "BIC","aBIC","RMSEA_Estimate")
fit.measures[,1] <- c("Long", "Wide")

# 2.3.1 ML-msst time-variant data ----

file.name <- paste0(model, "long_tv_n", N, "_i", I, "_nt", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(data.tv$data.long, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data.tv$data.long)[-(1:2)],
                                       cluster = names(data.tv$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "ML",
                                       iterations = 500000,
                                       h1iterations = 500000)

ml_syntax <- write.mlmsst.to.Mplus(data.tv$data.long[, -(1:2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/", folder,file.name,".inp"))

fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  estimates <- fit$parameters$unstandardized[,3]
  
  est.within <- list(loadings = estimates[1:I], #within loadings
                     state.var = estimates[I + 1], #state.var
                     error.var = estimates[(I + 2):(2 * I + 1)]) # error.var
  
  est.between <- list(loadings = estimates[((3*I + 2):(4*I + 1))], #loadings
                      intercepts = estimates[((5*I + 3):(6*I + 2))], # intercepts
                      trait.mean = estimates[(4 * I + 2)], # trait.mean
                      trait.var = estimates[(6*I +3)]) # trait.var
  
  est.var.coeff.long <- msst.var.coeff(within.parameters = est.within,
                                       between.parameters = est.between)
  
  fit.measures[1, 2:12] <- fit$summaries[, c(18, 11:14, 20:25)]
  
  rm(file.name, fit, est.between, est.within, estimates)
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

# 2.3.2 msst data.wide ("Perfect fit") ----

file.name <- paste0(model, "_wide_tv_n", N, "_i", I, "_nt", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(data.tv$data.wide, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data.tv$data.wide)[-1],
                                       analysis_type = "GENERAL",
                                       estimator = "ML",
                                       iterations = 50000,
                                       h1iterations = 50000)

mplus_syntax <- write.msst.to.Mplus(data.tv$data.wide[,-1], neta = nT, ntheta = 1, 
                                    equiv.assumption = list(tau = "cong", theta = "cong"),
                                    scale.invariance = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                    homocedasticity.assumption = list(error = FALSE, state.red = FALSE),
                                    second.order.trait = FALSE)

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(mplus_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, mplus_syntax)

# Run model in Mplus
runModels(paste0(getwd(),"/", folder,file.name,".inp"))

fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  estimates <- fit$parameters$unstandardized[ ,3]
  
  wloadings <- estimates[1:(I*nT)] # within loadings
  state.var <- estimates[((I * nT * 3) + ((nT+1) * nT / 2) + 2):((I * nT * 3) + ((nT+1) * nT / 2) + nT + 1)] #state variances
  error.var <- estimates[((I * nT * 3) + ((nT+1) * nT / 2) + nT + 3):((I * nT * 4) + ((nT+1) * nT / 2) + nT + 2)] #error variances
  
  bloadings <- estimates[((I * nT + 1):(I * nT * 2))] #between loadings
  intercepts <- estimates[((I * nT * 2) + ((nT+1) * nT / 2) + 2):((I * nT * 3) + ((nT+1) * nT / 2) + 1)] #intercepts
  trait.mean <- estimates[((I * nT * 2) + ((nT+1) * nT / 2) + 1)] #trait mean
  trait.var <- estimates[((I * nT * 3) + ((nT+1) * nT / 2) + nT + 2)] # trait variance
  
  est.within <- list(loadings = matrix(wloadings, nrow = nT, ncol = I, byrow = TRUE),
                     state.var = state.var,
                     error.var = matrix(error.var,  nrow = nT, ncol = I, byrow = TRUE))
  est.between <- list(loadings = matrix(bloadings, nrow = nT, ncol = I, byrow = TRUE),
                      intercepts = matrix(intercepts, nrow = nT, ncol = I, byrow = TRUE),
                      trait.mean = trait.mean,
                      trait.var = trait.var)
  
  est.var.coeff.wide <- msst.var.coeff.tv(I=I, nT= nT, within.parameters = est.within,
                                          between.parameters = est.between)
  
  fit.measures[2, 2:12] <- fit$summaries[, c(18, 11:14, 20:25)]
  
  rm(file.name, fit, estimates, est.within, est.between, bloadings, error.var, intercepts, state.var,
     trait.mean, trait.var, wloadings)
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

# 2.4 Plot variance coefficients and save fit measures ---- 

jpeg(paste0(getwd(), "/", folder, "OutputPlots/varcoeff_", model, "_n", N, "_i", I, "_nt", nT, ".jpg" ))
matplot(cbind(t(true.var.coeff[seq(4,12, by = I),]), t(est.var.coeff.wide[seq(4,12, by = I),])), type = c("l"), ylim = c(0, 1),
        col = c("red", "blue", "black"), xlab = "Time", ylab = "Explained Variance", lty = rep(1:2, each = 3), lwd = 2,
        main = paste0(model, ": Variance coefficients."))
abline(h = est.var.coeff.long[seq(4,12, by = I),], col = c("red", "blue", "black"), lty = 3)
legend("bottomright", legend = c("Con", "Spe", "Rel"), col = c("red", "blue", "black"),
       lty = 1, lwd = 2, pch = c(16,16,16), cex = 0.8)
dev.off()

print(xtable(fit.measures, type = "latex", caption = paste0("Fit measures: multi-level and single-level ", model, " N =", N, ", I = ", I, ", and nT =", nT, ".")), 
      include.rownames = FALSE, file = paste0(folder, "OutputTables/FitMeasures_", model, "_n", N, "_i", I, "_nt", nT, ".txt"))

# 2.5 Clean environment ----

rm(data.tv, est.var.coeff.long, est.var.coeff.wide, true.var.coeff, fit.measures, model)

# 3.0 TSO ----

# 3.1 Set Conditions and true parameters ----

model <- "tso" 
set.seed(seed)

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

# 3.2 Simulate data ----

# Time variant data

data.tv <- sim.data.tso.tv(N, nT, I, within.parameters = within.parameters, time.invariant = FALSE,
                               between.parameters = between.parameters, seed = seed)
rm(within.parameters, between.parameters)

# 3.3 Model estimation ----

# Compute true variance coefficients

true.var.coeff <- tso.var.coeff.tv(I, nT, within.parameters = data.tv$within.parameters, 
                                    between.parameters = data.tv$between.parameters)
# Matrix to save fit measures
#change to bayesian  fit measures
fit.measures <- data.frame(matrix(NA, 2, 4))
names(fit.measures) <- c("Data","Parameters","DIC","pD")
fit.measures[,1] <- c("Long", "Wide")

# 3.3.1 ML-tso time-variant data ----

file.name <- paste0(model, "long_tv_n", N, "_i", I, "_nt", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(data.tv$data.long, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data.tv$data.long)[-(1:2)],
                                       cluster = names(data.tv$data.long)[1],
                                       analysis_type = "TWOLEVEL",
                                       estimator = "BAYES",
                                       iterations = 5000,
                                       processors = 2)

ml_syntax <- write.mltso.to.Mplus(data.tv$data.long[, -c(1, 2)])

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  estimates <- fit$parameters$unstandardized[,3]
  
  est.within <- list(loadings = estimates[1:I], #within loadings
                     ar.effect = estimates[I + 1], # autoregresive effect
                     error.var = estimates[(I + 2):(2 * I + 1)], #error.var
                     state.var = estimates[(2 * I + 2)]) # state.var
  
  est.between <- list(intercepts = estimates[(4 * I + (I * (I - 1) / 2) + 3):(5 * I + (I * (I - 1) / 2) + 2)], # intercepts
                      trait.ind.var = estimates[(5 * I + (I * (I - 1) / 2) + 3):(6 * I + (I * (I - 1) / 2) + 2)]) # trait indicator variances
  
  
  est.var.coeff.long <- tso.var.coeff(I = I, nT = nT, within.parameters = est.within,
                                      between.parameters = est.between)
  
  fit.measures[1, 2:4] <- fit$summaries[, c(11:13)]
  
  rm(file.name, fit, est.between, est.within, estimates)
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

# 3.3.2 tso data.wide ("Perfect fit") ----

file.name <- paste0(model, "_wide_tv_n", N, "_i", I, "_nt", nT)

# Prepare data: Write data in Mplus format and write input file template
prepareMplusData(data.tv$data.wide, paste0(folder,file.name,".dat"), inpfile = T)

# Complete Mplus syntax
analysis_syntax <- write.Mplus.options(usevariables = names(data.tv$data.wide)[-1],
                                       analysis_type = "GENERAL",
                                       estimator = "BAYES",
                                       iterations = 30000,
                                       processors = 2)

ml_syntax <- write.tso.to.Mplus(data.tv$data.wide[,-1], nocc = nT, figure = "3b",
                                equiv.assumption = list(occ = "cong", theta = "equi"),
                                scale.invariance = list(int = TRUE, lambda = TRUE),
                                homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                autoregressive.homogeneity = FALSE)

write(analysis_syntax, paste0(folder,file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0(folder,file.name,".inp"), append = T)

rm(analysis_syntax, ml_syntax)

# Run modelin Mplus
runModels(paste0(getwd(),"/",folder,file.name,".inp"))

fit <- readModels(paste0(getwd(),"/",folder,file.name,".out")) #read Mplus output

if(check.mplus(fit, paste0(getwd(),"/",folder,file.name,".out")) == "Ok"){
  estimates <- fit$parameters$unstandardized[ ,3]
  
  wloadings <- estimates[1:(I*nT)] # within loadings
  ar.effects <- estimates[(2 * I * nT + 1):(2 * I * nT + nT - 1)] #autoregressive effects
  error.var <- estimates[((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((4 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I)] #error variances
  state.var <- estimates[c(((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2)), 
                           ((4 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + 1):((4 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I + nT - 1))] #state variances
  
  bloadings <- estimates[(I*nT + 1):(I*nT*2)] #between loadings
  intercepts <- estimates[((2 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + nT) + ((nT + I) * (nT + I - 1) / 2))] #intercepts
  trait.ind.var <- estimates[((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + 1):((3 * I * nT + 2 * nT) + ((nT + I) * (nT + I - 1) / 2) + I)] # trait indicators variances
  
  est.within <- list(loadings = matrix(wloadings, nrow = nT, ncol = I, byrow = TRUE),
                     ar.effects = ar.effects,
                     state.var = state.var,
                     error.var = matrix(error.var,  nrow = nT, ncol = I, byrow = TRUE))
  est.between <- list(loadings = matrix(bloadings, nrow = nT, ncol = I),
                      intercepts = matrix(intercepts, nrow = nT, ncol = I, byrow = TRUE),
                      trait.ind.var = trait.ind.var)
  
  est.var.coeff.wide <- tso.var.coeff.tv(I=I, nT= nT, within.parameters = est.within,
                                         between.parameters = est.between)
  
  #fit.measures[2, 2:4] <- fit$summaries[, c(11:13)]
  
  
  rm(file.name, fit, est.between, est.within, ar.effects, bloadings,
     error.var, estimates, intercepts, state.var, trait.ind.var, wloadings)
}else{
  stop("Model estimation did not converge or there are warning or error messages in the output")
}

# 3.4 Plot variance coefficients and save fit measures ---- 

jpeg(paste0(getwd(), "/", folder, "OutputPlots/varcoeff_", model, "_n", N, "_i", I, "_nt", nT, ".jpg" ))
matplot(cbind(t(true.var.coeff[seq(4,20, by = I),]), t(est.var.coeff.wide[seq(4,20, by = I),])), type = c("l"), ylim = c(0, 1),
        col = c("green", "orange", "red", "blue", "black"), xlab = "Time", ylab = "Explained Variance", lty = rep(1:2, each = 5), lwd = 2,
        main = paste0(model, ": Variance coefficients."))
abline(h = est.var.coeff.long[seq(4,20, by = I),nT], col = c("green", "orange", "red", "blue", "black"), lty = 3)
legend("bottomright", legend = c("Pred", "UPred", "Con", "Spe", "Rel"), col = c("green", "orange", "red", "blue", "black"),
       lty = 1, lwd = 2, pch = c(NA,NA,16,16,16), cex = 0.8)
dev.off()

#print(xtable(fit.measures, type = "latex", caption = paste0("Fit measures: multi-level and single-level ", model, " N =", N, ", I = ", I, ", and nT =", nT, ".")), 
 #     include.rownames = FALSE, file = paste0(folder, "OutputTables/FitMeasures_", model, "_n", N, "_i", I, "_nt", nT, ".txt"))


# 3.5 Clean environment ----

rm(N, nT, I, seed, data.tv, model, est.var.coeff.long, est.var.coeff.wide, fit.measures,
   true.var.coeff)
