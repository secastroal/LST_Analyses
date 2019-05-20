# This script analyze the results of the simulation analysis
#Prepare environment ----
library(plyr)
library(xtable)
library(ggplot2)
file.sources <- paste0("R/", list.files(paste0(getwd(), "/R")))
sapply(file.sources,source,.GlobalEnv)
rm(file.sources)

# Create true parameters ----
I <- 4
nT <- 30
R <- 100

# matrix to store True parameters
 
parameters.true <- matrix(NA, 3, 30)
colnames(parameters.true) <- c(paste0("w_loading", 1:I), "common_state_var", paste0("unique_state_error_var", 1:I),
                               "ar_effect", paste0("b_loading", 1:I), paste0("intercept", 1:I), "common_trait_var",
                               "trait_mean", paste0("unique_indicator_trait_var", 1:I), paste0("cov", c(12, 13, 23, 14, 24, 34)))
row.names(parameters.true) <- c("MSST", "CUTS", "TSO")

var_coeff.true <- matrix(NA, 3, 20)
colnames(var_coeff.true) <- paste0(rep(c( "ccon_pred_y", "ucon_upred_y", "tcon_con_y", "spe_y", "rel_y"), each = I), 
                                   1:I)
row.names(var_coeff.true) <- c("MSST", "CUTS", "TSO")

# MSST true parameters

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
parameters.true[1, c(1:9, 11:20)] <- c(state_loadings, var_state, var_m_error, trait_loadings,
                                  intercepts, var_trait, mean_trait) 

rm(state_loadings, var_state, var_m_error, intercepts, trait_loadings, var_trait, 
   mean_trait)

# Compute true variance coefficients
var_coeff.true[1, (2 * I + 1):(5 * I)] <- t(msst.var.coeff(within.parameters = within.parameters, 
                                                      between.parameters = between.parameters)) 

rm(within.parameters, between.parameters)

# CUTS true parameters

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
parameters.true[2, c(1:9, 11:19, 21:24)] <- c(state_loadings, var_CS, var_US, trait_loadings,
                                         intercepts, var_CT, var_UT) 

rm(state_loadings, var_CS, var_US, intercepts, trait_loadings, var_CT, var_UT)

# Compute true variance coefficients

var_coeff.true[2, ] <- t(cuts.var.coeff(within.parameters = within.parameters, 
                                   between.parameters = between.parameters)) 

rm(within.parameters, between.parameters)

# TSO true parameters

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

set.seed(13001)
Rcor <- matrix(sample((7:9)/10, size = I * I, replace = TRUE), I) #correlation matrix trait indicators
Rcor[lower.tri(Rcor)] = t(Rcor)[lower.tri(Rcor)]
diag(Rcor) <- 1

D <- diag(sqrt(var_ind_traits))
Sigma <- D%*%Rcor%*%D #Variance-Covariance matrix of the trait indicators

between.parameters <- list(intercepts = intercepts, trait.ind.var = var_ind_traits, 
                           cor.matrix = Rcor, Sigma = Sigma)

rm(intercepts, var_ind_traits, Rcor, D, Sigma)

# Save True parameters to parameters matrix
parameters.true[3, c(1:10, 15:18, 21:30)] <- c(within.parameters$loadings, within.parameters$state.var,
                                          within.parameters$error.var, within.parameters$ar.effect,
                                          between.parameters$intercepts, between.parameters$trait.ind.var,
                                          round(between.parameters$Sigma[t(lower.tri(between.parameters$Sigma))], 2))

# Compute true variance coefficients
var_coeff.true[3, ] <- tso.var.coeff(I = I, nT = nT, within.parameters = within.parameters, 
                                between.parameters = between.parameters)[,nT] 

rm(within.parameters, between.parameters)

#Export true parameters to LaTeX

parameters.true.ex <- parameters.true
colnames(parameters.true.ex) <- c(paste0("Within loading $\\lambda_{S_{", 1:I,"}}$"), 
                                  "(Common) State variance $var(\\zeta)$", 
                               paste0("Unique state/ Error variance $var(\\varepsilon_", 1:I, ")$"), 
                               "Autoregressive effect $\\beta$", 
                               paste0("Between loading $\\lambda_{T_{", 1:I, "}}$"), 
                               paste0("Intercept $\\alpha_", 1:I, "$"), 
                               "(Common) Trait variance $var(\\xi)$",
                               "Trait mean $\\hat{\\xi}$", 
                               paste0("Unique trait/ Indicator trait variance $var(\\vartheta/\\xi_", 1:I, ")$"), 
                               paste0("$Cov(", c("\\xi_1, \\xi_2)$", "\\xi_1, \\xi_3)$", 
                                                "\\xi_2, \\xi_3)$", "\\xi_1, \\xi_4)$", 
                                                "\\xi_2, \\xi_4)$", "\\xi_3, \\xi_4)$")))

print(xtable(t(parameters.true.ex), type = "latex", caption = "True Parameters Used per Model",
             label = "tab:Truepar", align = c("l", "c", "c", "c")),
      include.colnames=T, sanitize.rownames.function = identity,
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", file = "Mplus_Simulation/trueparameters.txt")


var_coeff.true.ex <- var_coeff.true
colnames(var_coeff.true.ex) <- paste0(rep(c( "Common consistency/Trait predictability Y", 
                                          "Unique consistency/Trait unpredictability Y", 
                                          "(Total) Consistency Y", 
                                          "Occasion Specificity Y", 
                                          "Reliability Y"), each = I), 
                                   1:I)

print(xtable(t(var_coeff.true.ex), type = "latex", caption = "True Variance Coefficient Components",
             label = "tab:Truevar", align = c("l", "c", "c", "c"), digits = c(0,2,2,2)), 
      include.rownames = TRUE, NA.string = "-", caption.placement = "top", file = "Mplus_Simulation/truevarcoeff.txt")






# Create factor variables for anova analyses ----
b.model <- factor(rep(1:3, each = R * 6 * 11), levels = c("1", "2", "3"), labels = c("b.msst", "b.cuts", "b.tso"))
NA.prop  <- factor(rep(rep(1:2, 3), each = R * 3 * 11), levels = c("1", "2"), labels = c("0%", "10%"))
times   <- factor(rep(rep(1:3, 6), each = R * 11), levels = c("1", "2", "3"), labels = c("30", "60", "90"))
model   <- factor(rep(rep(1:3, each = 4)[1:11], R * 3 * 6), levels = c("1", "2", "3"), labels = c("msst", "cuts", "tso"))
est.method <- factor(rep(c(1,1,2,2,1,1,2,2,1,2,2), R * 3 * 6), levels = c("1", "2"), labels = c("ml", "bayes"))
format  <- factor(rep(c(1,2,1,2,1,2,1,2,1,1,2), R * 3 * 6), levels = c("1", "2"), labels = c("wide", "long"))

factor.var <- data.frame(b.model, NA.prop, times, model, est.method, format)
rm(b.model, NA.prop, times, model, est.method, format)

# Performances ----
files <- paste(getwd(), "Mplus_Simulation" , "Performance", paste0("performance_",  1:18, ".dat"), sep = "/")

performances <- lapply(files, function(x) read.table(file = x, header = TRUE, colClasses = "character"))

performances <- ldply(performances, data.frame)

rm(files)

perf.table <- rbind(apply(performances, 2, function(x) length(which(x == "Ok"))),
      apply(performances, 2, function(x) length(which(x == "Errors/Warnings"))),
      apply(performances, 2, function(x) length(which(x == "Non-convergence"))),
      apply(performances, 2, function(x) length(which(x == "timeout"))))

perf.prop <- t(round(prop.table(as.table(perf.table), margin = 2)*100, 2))
dimnames(perf.prop)[[2]] <- c("Ok", "Warnings/Errors", "Non-convergence", "timeout")

print(xtable(perf.prop, type = "latex", caption = "Performance percentanges across all conditions"), 
include.rownames = TRUE, file = "Mplus_Simulation/performances.txt")

# Running times ----
files <- paste(getwd(), "Mplus_Simulation" , "Times", paste0("times_",  1:18, ".dat"), sep = "/")

running.times <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

running.times <- ldply(running.times, data.frame)

# Parameter estimates ----
files <- paste(getwd(), "Mplus_Simulation" , "Parameters", paste0("parameters_",  1:18, ".dat"), sep = "/")

parameters <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

parameters <- ldply(parameters, data.frame)
parameters[which(parameters == -999, arr.ind = TRUE)] <- NA

parameters.bias <- parameters - parameters.true

fit <- lm(parameters.bias$common_state_var ~ factor.var$b.model + factor.var$NA.prop + factor.var$times + 
            factor.var$model + factor.var$est.method + factor.var$format, na.action =  "na.exclude")
anova(fit)

# Standard errors and/or posterior sd ----
files <- paste(getwd(), "Mplus_Simulation" , "SE_PSD", paste0("se_psd_",  1:18, ".dat"), sep = "/")

se_psd <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

se_psd <- ldply(se_psd, data.frame)

# Variance coefficients ----
files <- paste(getwd(), "Mplus_Simulation" , "Var_Coeff", paste0("var_coeff_",  1:18, ".dat"), sep = "/")

var_coeff <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

var_coeff <- ldply(var_coeff, data.frame)
var_coeff[which(var_coeff == -999, arr.ind = TRUE)] <- NA

var_coeff.bias <- var_coeff - var_coeff.true

fit <- lm(var_coeff.est$rel_y2 ~ factor.var$b.model + factor.var$NA.prop + factor.var$times + 
            factor.var$model + factor.var$est.method + factor.var$format, na.action =  "na.exclude")
anova(fit)

jpeg("Mplus_Simulation/rel1_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), rel_y1)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/rel2_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), rel_y2)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/rel3_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), rel_y3)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/rel4_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), rel_y4)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/con1_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), tcon_con_y1)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/con2_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), tcon_con_y2)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/con3_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), tcon_con_y3)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/con4_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), tcon_con_y4)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/spe1_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), spe_y1)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/spe2_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), spe_y2)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/spe3_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), spe_y3)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()

jpeg("Mplus_Simulation/spe4_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(b.model, model), spe_y4)) + geom_violin()+geom_boxplot(width = 0.1)
dev.off()



# Posterior sd of variance coefficients ----
files <- paste(getwd(), "Mplus_Simulation" , "PSD_Var_Coeff", paste0("psd_var_coeff_",  1:18, ".dat"), sep = "/")

psd_var_coeff <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

psd_var_coeff <- ldply(psd_var_coeff, data.frame)

# Fit measures ----
files <- paste(getwd(), "Mplus_Simulation" , "Fit_Measures", paste0("fit_measures_",  1:18, ".dat"), sep = "/")

fit_measures <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

fit_measures <- ldply(fit_measures, data.frame)

# End ----