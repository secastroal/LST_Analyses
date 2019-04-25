# This script analyze the results of the simulation analysis
library(plyr)
library(xtable)
library(ggplot2)

# Create factor variables for anova analyses
b.model <- factor(rep(1:3, each = 50 * 6 * 11)[1:8250], levels = c("1", "2", "3"), labels = c("b.msst", "b.cuts", "b.tso"))
NA.prop  <- factor(rep(rep(1:2, 3), each = 50 * 3 * 11)[1:8250], levels = c("1", "2"), labels = c("0%", "10%"))
times   <- factor(rep(rep(1:3, 6), each = 50 * 11)[1:8250], levels = c("1", "2", "3"), labels = c("30", "60", "90"))
model   <- factor(rep(rep(1:3, each = 4)[1:11], 50 * 3 * 6)[1:8250], levels = c("1", "2", "3"), labels = c("msst", "cuts", "tso"))
est.method <- factor(rep(c(1,1,2,2,1,1,2,2,1,2,2), 50 * 3 * 6)[1:8250], levels = c("1", "2"), labels = c("ml", "bayes"))
format  <- factor(rep(c(1,2,1,2,1,2,1,2,1,1,2), 50 * 3 * 6)[1:8250], levels = c("1", "2"), labels = c("wide", "long"))

factor.var <- data.frame(b.model, NA.prop, times, model, est.method, format)
rm(b.model, NA.prop, times, model, est.method, format)

# Performances ----
files <- paste(getwd(), "Mplus_Simulation" , "Performance", paste0("performance_",  1:15, ".dat"), sep = "/")

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
include.rownames = FALSE, file = "Mplus_Simulation/performances.txt")

# Running times ----
files <- paste(getwd(), "Mplus_Simulation" , "Times", paste0("times_",  1:15, ".dat"), sep = "/")

running.times <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

running.times <- ldply(running.times, data.frame)

# Parameter estimates ----
files <- paste(getwd(), "Mplus_Simulation" , "Parameters", paste0("parameters_",  1:15, ".dat"), sep = "/")

parameters <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

parameters <- ldply(parameters, data.frame)
parameters[which(parameters == -999, arr.ind = TRUE)] <- NA

parameters.true <- parameters[seq(1, 8265, by = 11 * 50 + 1), ]
parameters.true <- parameters.true[rep(1:15, each = 11 * 50), ]

parameters.est <- parameters[-seq(1, 8265, by = 11 * 50 + 1),]

parameters.bias <- parameters.est - parameters.true

fit <- lm(parameters.bias$common_state_var ~ factor.var$b.model + factor.var$NA.prop + factor.var$times + 
            factor.var$model + factor.var$est.method + factor.var$format, na.action =  "na.exclude")
anova(fit)

# Standard errors and/or posterior sd ----
files <- paste(getwd(), "Mplus_Simulation" , "SE_PSD", paste0("se_psd_",  1:15, ".dat"), sep = "/")

se_psd <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

se_psd <- ldply(se_psd, data.frame)

# Variance coefficients ----
files <- paste(getwd(), "Mplus_Simulation" , "Var_Coeff", paste0("var_coeff_",  1:15, ".dat"), sep = "/")

var_coeff <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

var_coeff <- ldply(var_coeff, data.frame)
var_coeff[which(var_coeff == -999, arr.ind = TRUE)] <- NA

var_coeff.true <- var_coeff[seq(1, 8265, by = 11 * 50 + 1), ]
var_coeff.true <- var_coeff.true[rep(1:15, each = 11 * 50), ]

var_coeff.est <- var_coeff[-seq(1, 8265, by = 11 * 50 + 1),]

var_coeff.bias <- var_coeff.est - var_coeff.true

fit <- lm(var_coeff.est$rel_y2 ~ factor.var$b.model + factor.var$NA.prop + factor.var$times + 
            factor.var$model + factor.var$est.method + factor.var$format, na.action =  "na.exclude")
anova(fit)

jpeg("Mplus_Simulation/rel1_bias.jpg")
ggplot(data.frame(factor.var, var_coeff.bias), 
       aes(interaction(est.method, times), rel_y1)) + geom_violin()+geom_boxplot(width = 0.1)
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
files <- paste(getwd(), "Mplus_Simulation" , "PSD_Var_Coeff", paste0("psd_var_coeff_",  1:15, ".dat"), sep = "/")

psd_var_coeff <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

psd_var_coeff <- ldply(psd_var_coeff, data.frame)

# Fit measures ----
files <- paste(getwd(), "Mplus_Simulation" , "Fit_Measures", paste0("fit_measures_",  1:15, ".dat"), sep = "/")

fit_measures <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

fit_measures <- ldply(fit_measures, data.frame)

# End ----