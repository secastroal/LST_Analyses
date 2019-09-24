# This script analyze the results of the simulation analysis
#Prepare environment ----
rm(list=ls())
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

# Create labels for the 11 types of analyses 
models <- c("msst", "cuts", "tso")
data_str <- c("wide", "long")
estimator <- c("ml", "bayes")
labels <- expand.grid(data_str, estimator, models)
labels <- labels[-10, c(3,1,2)]
labels <- apply(labels, 1, function(vec) paste(vec, collapse = "_"))
rm(models, data_str, estimator)

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

# Table: True Parameters ----

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


var_coeff.true.ex <- var_coeff.true[, c(17:20, 9:12, 1:8, 1:8, 13:16)]
var_coeff.true.ex[2, 17:24] <- NA
var_coeff.true.ex[3, 9:16] <- NA
colnames(var_coeff.true.ex) <- paste0(rep(c( "Reliability $Y_",
                                             "\\hspace{0.5cm}(Total) Consistency $Y_",
                                             "\\hspace{1cm}Common Consistency $Y_", 
                                             "\\hspace{1cm}Unique Consistency $Y_", 
                                             "\\hspace{1cm}Predictability by Trait $Y_", 
                                             "\\hspace{1cm}Unpredictability by Trait $Y_", 
                                             "\\hspace{0.5cm}Occasion Specificity $Y_" 
                                          ), each = I), 
                                   1:I, "$")

print(xtable(t(var_coeff.true.ex), type = "latex", caption = "True Variance Coefficient Components",
             label = "tab:Truevar", align = c("l", "c", "c", "c"), digits = c(0,2,2,2)), 
      include.rownames = TRUE, NA.string = "-", sanitize.text.function = function(x){x}, 
      caption.placement = "top", file = "Mplus_Simulation/truevarcoeff.txt")
rm(var_coeff.true.ex, parameters.true.ex)

# Repeat true parameters in a matrix with 19800 rows ----
parameters.true <- parameters.true[rep(1:3, each = 11 * R * 6),]

var_coeff.true <- var_coeff.true[rep(1:3, each = 11 * R * 6),]

# Create factor variables for anova analyses ----
cond <- factor(rep(1:18, each = R * 11))
b.model <- factor(rep(1:3, each = R * 6 * 11), levels = c("1", "2", "3"), labels = c("b.msst", "b.cuts", "b.tso"))
NA.prop  <- factor(rep(rep(1:2, 3), each = R * 3 * 11), levels = c("1", "2"), labels = c("0%", "10%"))
times   <- factor(rep(rep(1:3, 6), each = R * 11), levels = c("1", "2", "3"), labels = c("30", "60", "90"))
model   <- factor(rep(rep(1:3, each = 4)[1:11], R * 3 * 6), levels = c("1", "2", "3"), labels = c("msst", "cuts", "tso"))
est.method <- factor(rep(c(1,1,2,2,1,1,2,2,1,2,2), R * 3 * 6), levels = c("1", "2"), labels = c("ml", "bayes"))
format  <- factor(rep(c(1,2,1,2,1,2,1,2,1,1,2), R * 3 * 6), levels = c("1", "2"), labels = c("wide", "long"))

factor.var <- data.frame(cond, b.model, NA.prop, times, model, est.method, format)
rm(cond, b.model, NA.prop, times, model, est.method, format)

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

perf.table <- cbind(c("MSST", "ML-MSST", "MSST", "ML-MSST", "CUTS", "ML-CUTS", "CUTS", "ML-CUTS",
                      "TSO", "TSO", "ML-TSO"), 
                    c("MLE", "MLE", "Bayes", "Bayes", "MLE", "MLE", "Bayes", "Bayes", "MLE", "Bayes", "Bayes"), 
                    t(perf.table))
dimnames(perf.table)[[2]] <- c("Model", "Est.Method", "Successful", "Warnings/Errors", "Non-convergence", "Timeout")

perf.prop <- cbind(c("MSST", "ML-MSST", "MSST", "ML-MSST", "CUTS", "ML-CUTS", "CUTS", "ML-CUTS",
                      "TSO", "TSO", "ML-TSO"), 
                    c("MLE", "MLE", "Bayes", "Bayes", "MLE", "MLE", "Bayes", "Bayes", "MLE", "Bayes", "Bayes"), 
                    perf.prop)

dimnames(perf.prop)[[2]] <- c("Model", "Est.Method", "Successful", "Warnings/Errors", "Non-convergence", "Timeout")

perf.cond <- matrix(NA, 18, 11*4)
colnames(perf.cond) <- c(paste(rep(labels, 4),rep(c("Ok", "WE", "NC", "TOut"), each = 11), 
                               sep = " "))

for(i in 1:18){
  perf.subset <- performances[((100 *(i-1))+1):(100*i),]
  perf.cond[i, ] <- c(apply(perf.subset, 2, function(x) length(which(x == "Ok"))),
                      apply(perf.subset, 2, function(x) length(which(x == "Errors/Warnings"))),
                      apply(perf.subset, 2, function(x) length(which(x == "Non-convergence"))),
                      apply(perf.subset, 2, function(x) length(which(x == "timeout"))))
}
rm(i, perf.subset)

# Table: Performance ----
print(xtable(perf.table, type = "latex", align = c("l", "l", "l", "c", "c", "c", "c"), label = "tab:perf", 
             caption = "Number of Successes and Failures per Type of Analysis on 1800 Analyses Each"), 
include.rownames = FALSE, caption.placement = "top", file = "Mplus_Simulation/performancesnumbers.txt")

print(xtable(perf.prop, type = "latex", align = c("l", "l", "l", "c", "c", "c", "c"), label = "tab:perf", 
             caption = "Percentages of Successes and Failures per Type of Analysis on 1800 Analyses Each"), 
      include.rownames = FALSE, caption.placement = "top", file = "Mplus_Simulation/performancespercentages.txt")


# Plot: Successful analysis ----
pdf("Mplus_Simulation/Okplot.pdf")
par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
matplot(1:6,perf.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(-2,-40,c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
       lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
mtext("Number of Successful Analyses", 2, outer=TRUE, line=2.5)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
dev.off()

# Plot: Analyses with warnings and error messages ----
pdf("Mplus_Simulation/warningplot.pdf")
par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
matplot(1:6,perf.cond[1:6,12:15],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,16:19],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,20:22],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[7:12,12:15],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,16:19],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,20:22],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[13:18,12:15],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,16:19],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,20:22],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(-2,-40,c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
       lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
mtext("Number of Analyses with Warnings or Errors", 2, outer=TRUE, line=2.5)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
dev.off()

# Plot: Analyses that did not converge ----
pdf("Mplus_Simulation/nonconvergenceplot.pdf")
par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
matplot(1:6,perf.cond[1:6,23:26],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,27:30],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,31:33],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[7:12,23:26],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,27:30],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,31:33],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[13:18,23:26],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,27:30],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,31:33],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(-2,-40,c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
       lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
mtext("Number of Analyses that did not Converge", 2, outer=TRUE, line=2.5)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
dev.off()

# Plot: Analyses that time out ----
pdf("Mplus_Simulation/timeoutplot.pdf")
par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
matplot(1:6,perf.cond[1:6,34:37],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,38:41],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[1:6,42:44],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[7:12,34:37],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,38:41],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
matplot(1:6,perf.cond[7:12,42:44],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
matplot(1:6,perf.cond[13:18,34:37],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,38:41],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,perf.cond[13:18,42:44],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
        xlim=c(0.7,6.2),ylim=c(0,110),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]-6,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(-2,-40,c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
       lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
mtext("Number of Analyses that Timeout", 2, outer=TRUE, line=2.5)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
dev.off()

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

# Plots: Bias Parameters----
ylabel <- c(expression(paste("Bias ", lambda["S"[2]])),
            expression(paste("Bias ", lambda["S"[3]])),
            expression(paste("Bias ", lambda["S"[4]])),
            expression(paste("Bias var(", zeta, ")")),
            expression(paste("Bias var(", epsilon[1], ")")),
            expression(paste("Bias var(", epsilon[2], ")")),
            expression(paste("Bias var(", epsilon[3], ")")),
            expression(paste("Bias var(", epsilon[4], ")")))
ylimits.down <- c(-0.11, -0.165, -0.05, -0.1, -0.25, -0.25, -0.1, -0.1 )
ylimits.up <- c(0.025, 0.025, 0.05, 1.75, 0.65, 0.25, 0.75, 1.1 )
xlabel.factor <- c(0.075,0.075,0.1,1,0.25,0.15,0.5,0.85)
legend.factor <- c(1.45,1.4,1.7,7.5,2.25,1.7,4,5.35)

for(j in 2:9){
  bias.cond <- matrix(NA, 18, 11)
  colnames(bias.cond) <- labels
  for(i in 1:18){
    bias.subset <- subset(parameters.bias, factor.var$cond == i)[,j]
    bias.subset <- matrix(bias.subset, 100, 11, byrow = TRUE)
    bias.cond[i,] <- apply(bias.subset, 2, function(x) mean(x, na.rm = TRUE))
  }
  rm(i, bias.subset)
  bias.cond[which(perf.cond[,1:11] < 10, arr.ind = TRUE)] <- NA
  
  pdf(paste0("Mplus_Simulation/biasparameter", j, "plot.pdf"))
  par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
  matplot(1:6,bias.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-1]*xlabel.factor[j-1],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-1]*xlabel.factor[j-1],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-1]*xlabel.factor[j-1],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  
  legend(-2,ylimits.down[j-1]*legend.factor[j-1],c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
         lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
  mtext(ylabel[j-1], 2, outer=TRUE, line=3.5)
  mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
  mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
  mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
  mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  dev.off()
}
rm(j, bias.cond, legend.factor, xlabel.factor, ylabel, ylimits.down, ylimits.up)

# Plots: AbBias Parameters----
ylabel <- c(expression(paste("AbBias ", lambda["S"[2]])),
            expression(paste("AbBias ", lambda["S"[3]])),
            expression(paste("AbBias ", lambda["S"[4]])),
            expression(paste("AbBias var(", zeta, ")")),
            expression(paste("AbBias var(", epsilon[1], ")")),
            expression(paste("AbBias var(", epsilon[2], ")")),
            expression(paste("AbBias var(", epsilon[3], ")")),
            expression(paste("AbBias var(", epsilon[4], ")")))
ylimits <- c(0.16, 0.16, 0.065, 1.75, 0.6, 0.26, 0.75, 1.1 )

for(j in 2:9){
  bias.cond <- matrix(NA, 18, 11)
  colnames(bias.cond) <- labels
  for(i in 1:18){
    bias.subset <- subset(parameters.bias, factor.var$cond == i)[,j]
    bias.subset <- abs(bias.subset)
    bias.subset <- matrix(bias.subset, 100, 11, byrow = TRUE)
    bias.cond[i,] <- apply(bias.subset, 2, function(x) mean(x, na.rm = TRUE))
  }
  rm(i, bias.subset)
  bias.cond[which(perf.cond[,1:11] < 10, arr.ind = TRUE)] <- NA
  
  pdf(paste0("Mplus_Simulation/abbiasparameter", j, "plot.pdf"))
  par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
  matplot(1:6,bias.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]-ylimits[j-1]*0.05,
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]-ylimits[j-1]*0.05,
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]-ylimits[j-1]*0.05,
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  
  legend(-2,-ylimits[j-1]*0.40,c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
         lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
  mtext(ylabel[j-1], 2, outer=TRUE, line=3.5)
  mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
  mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
  mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
  mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  dev.off()
}
rm(bias.cond, j, ylabel, ylimits)

# Plots: RMSE Parameters----
ylabel <- c(expression(paste("RMSE ", lambda["S"[2]])),
            expression(paste("RMSE ", lambda["S"[3]])),
            expression(paste("RMSE ", lambda["S"[4]])),
            expression(paste("RMSE var(", zeta, ")")),
            expression(paste("RMSE var(", epsilon[1], ")")),
            expression(paste("RMSE var(", epsilon[2], ")")),
            expression(paste("RMSE var(", epsilon[3], ")")),
            expression(paste("RMSE var(", epsilon[4], ")")))
ylimits <- c(0.16, 0.16, 0.065, 1.75, 0.6, 0.26, 0.75, 1.1 )

for(j in 2:9){
  bias.cond <- matrix(NA, 18, 11)
  colnames(bias.cond) <- labels
  for(i in 1:18){
    bias.subset <- subset(parameters.bias, factor.var$cond == i)[,j]
    bias.subset <- bias.subset*bias.subset
    bias.subset <- matrix(bias.subset, 100, 11, byrow = TRUE)
    bias.cond[i,] <- apply(bias.subset, 2, function(x) sqrt(mean(x, na.rm = TRUE)))
  }
  rm(i, bias.subset)
  bias.cond[which(perf.cond[,1:11] < 10, arr.ind = TRUE)] <- NA
  
  pdf(paste0("Mplus_Simulation/RMSEparameter", j, "plot.pdf"))
  par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
  matplot(1:6,bias.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]-ylimits[j-1]*0.05,
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]-ylimits[j-1]*0.05,
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(-ylimits[j-1]*0.05,ylimits[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]-ylimits[j-1]*0.05,
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  
  legend(-2,-ylimits[j-1]*0.40,c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
         lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
  mtext(ylabel[j-1], 2, outer=TRUE, line=3.5)
  mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
  mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
  mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
  mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  dev.off()
}
rm(bias.cond, j, ylabel, ylimits)

# Standard errors and/or posterior sd ----
files <- paste(getwd(), "Mplus_Simulation" , "SE_PSD", paste0("se_psd_",  1:18, ".dat"), sep = "/")

se_psd <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

se_psd <- ldply(se_psd, data.frame)

se_psd[which(se_psd == -999, arr.ind = TRUE)] <- NA


# Plots: SE/PSD----
ylabel <- c(expression(paste("SE/PSD ", lambda["S"[2]])),
            expression(paste("SE/PSD ", lambda["S"[3]])),
            expression(paste("SE/PSD ", lambda["S"[4]])),
            expression(paste("SE/PSD var(", zeta, ")")),
            expression(paste("SE/PSD var(", epsilon[1], ")")),
            expression(paste("SE/PSD var(", epsilon[2], ")")),
            expression(paste("SE/PSD var(", epsilon[3], ")")),
            expression(paste("SE/PSD var(", epsilon[4], ")")))
ylimits.down <- c(-0.005, -0.005, -0.005, -0.01, -0.005, -0.005, -0.01, -0.005 )
ylimits.up <- c(0.065, 0.065, 0.065, 0.165, 0.065, 0.065, 0.165, 0.065 )
xlabel.factor <- c(1,1,1,1,1,1,1,1)
legend.factor <- c(6,6,6,7.5,6,6,7.5,6)

for(j in 2:9){
  se.psd.cond <- matrix(NA, 18, 11)
  colnames(se.psd.cond) <- labels
  for(i in 1:18){
    se.psd.subset <- subset(se_psd, factor.var$cond == i)[,j]
    se.psd.subset <- matrix(se.psd.subset, 100, 11, byrow = TRUE)
    se.psd.cond[i,] <- apply(se.psd.subset, 2, function(x) mean(x, na.rm = TRUE))
  }
  rm(i, se.psd.subset)
  se.psd.cond[which(perf.cond[,1:11] < 10, arr.ind = TRUE)] <- NA
  
  pdf(paste0("Mplus_Simulation/mean.se.psd.par", j, "plot.pdf"))
  par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
  matplot(1:6,se.psd.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,se.psd.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,se.psd.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,se.psd.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,se.psd.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,se.psd.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,se.psd.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-1]*xlabel.factor[j-1],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,se.psd.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-1]*xlabel.factor[j-1],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,se.psd.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-1],ylimits.up[j-1]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-1]*xlabel.factor[j-1],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  
  legend(-2,ylimits.down[j-1]*legend.factor[j-1],c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
         lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
  mtext(ylabel[j-1], 2, outer=TRUE, line=3.5)
  mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
  mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
  mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
  mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  dev.off()
}
rm(j, se.psd.cond, legend.factor, xlabel.factor, ylabel, ylimits.down, ylimits.up)

# Variance coefficients ----
files <- paste(getwd(), "Mplus_Simulation" , "Var_Coeff", paste0("var_coeff_",  1:18, ".dat"), sep = "/")

var_coeff <- lapply(files, function(x) read.table(file = x, header = TRUE))

rm(files)

var_coeff <- ldply(var_coeff, data.frame)
var_coeff[which(var_coeff == -999, arr.ind = TRUE)] <- NA

var_coeff.bias <- var_coeff - var_coeff.true

# Plot: Bias Variance Coefficients ----
ylabel <- paste0(rep("Bias ", 3 * I),rep(c("Consistency Y", "Specificity Y", "Reliability Y"), each = I), rep(1:I, 3))
ylimits.down <- c(-0.65, -0.65, -0.65, -0.65, -0.1, -0.1, -0.1, -0.1, -0.25, -0.25, -0.25, -0.25)
ylimits.up <- c(0.1, 0.1, 0.1, 0.1, 0.65, 0.65, 0.65, 0.65, 0.075, 0.075, 0.075, 0.075)
xlabel.factor <- c(0.065, 0.065, 0.065, 0.065, 0.45, 0.45, 0.45, 0.45, 0.065, 0.065, 0.065, 0.065)
legend.factor <- c(1.4, 1.4, 1.4, 1.4, 3.65, 3.65, 3.65, 3.65, 1.45, 1.45, 1.45, 1.45)

for(j in 9:20){
  bias.cond <- matrix(NA, 18, 11)
  colnames(bias.cond) <- labels
  for(i in 1:18){
    bias.subset <- subset(var_coeff.bias, factor.var$cond == i)[,j]
    bias.subset <- matrix(bias.subset, 100, 11, byrow = TRUE)
    bias.cond[i,] <- apply(bias.subset, 2, function(x) mean(x, na.rm = TRUE))
  }
  rm(i, bias.subset)
  bias.cond[which(perf.cond[,1:11] < 10, arr.ind = TRUE)] <- NA
  
  pdf(paste0("Mplus_Simulation/", gsub("[[:space:]]", "", ylabel[j-8]), "plot.pdf"))
  par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
  matplot(1:6,bias.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  
  legend(-2,ylimits.down[j-8]*legend.factor[j-8],c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
         lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
  mtext(ylabel[j-8], 2, outer=TRUE, line=3.5)
  mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
  mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
  mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
  mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  dev.off()
}
rm(j, bias.cond, legend.factor, xlabel.factor, ylabel, ylimits.down, ylimits.up)

# Plot: Absolute Bias Variance Coefficients ----
ylabel <- paste0(rep("AbBias ", 3 * I),rep(c("Consistency Y", "Specificity Y", "Reliability Y"), each = I), rep(1:I, 3))
ylimits.down <- c(-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.075, -0.075, -0.075, -0.075)
ylimits.up <- c(0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.25, 0.25, 0.25, 0.25)
xlabel.factor <- c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.25, 0.25, 0.25, 0.25)
legend.factor <- c(3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 2.6, 2.6, 2.6, 2.6)

for(j in 9:20){
  bias.cond <- matrix(NA, 18, 11)
  colnames(bias.cond) <- labels
  for(i in 1:18){
    bias.subset <- subset(var_coeff.bias, factor.var$cond == i)[,j]
    bias.subset <- abs(bias.subset)
    bias.subset <- matrix(bias.subset, 100, 11, byrow = TRUE)
    bias.cond[i,] <- apply(bias.subset, 2, function(x) mean(x, na.rm = TRUE))
  }
  rm(i, bias.subset)
  bias.cond[which(perf.cond[,1:11] < 10, arr.ind = TRUE)] <- NA
  
  pdf(paste0("Mplus_Simulation/", gsub("[[:space:]]", "", ylabel[j-8]), "plot.pdf"))
  par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
  matplot(1:6,bias.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  
  legend(-2,ylimits.down[j-8]*legend.factor[j-8],c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
         lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
  mtext(ylabel[j-8], 2, outer=TRUE, line=3.5)
  mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
  mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
  mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
  mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  dev.off()
}
rm(j, bias.cond, legend.factor, xlabel.factor, ylabel, ylimits.down, ylimits.up)

# Plot: RMSEs Variance Coefficients ----
ylabel <- paste0(rep("RMSE ", 3 * I),rep(c("Consistency Y", "Specificity Y", "Reliability Y"), each = I), rep(1:I, 3))
ylimits.down <- c(-0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.075, -0.075, -0.075, -0.075)
ylimits.up <- c(0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.65, 0.25, 0.25, 0.25, 0.25)
xlabel.factor <- c(0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45, 0.25, 0.25, 0.25, 0.25)
legend.factor <- c(3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 3.65, 2.6, 2.6, 2.6, 2.6)

for(j in 9:20){
  bias.cond <- matrix(NA, 18, 11)
  colnames(bias.cond) <- labels
  for(i in 1:18){
    bias.subset <- subset(var_coeff.bias, factor.var$cond == i)[,j]
    bias.subset <- bias.subset*bias.subset
    bias.subset <- matrix(bias.subset, 100, 11, byrow = TRUE)
    bias.cond[i,] <- apply(bias.subset, 2, function(x) sqrt(mean(x, na.rm = TRUE)))
  }
  rm(i, bias.subset)
  bias.cond[which(perf.cond[,1:11] < 10, arr.ind = TRUE)] <- NA
  
  pdf(paste0("Mplus_Simulation/", gsub("[[:space:]]", "", ylabel[j-8]), "plot.pdf"))
  par(mfrow=c(3,3),mar=c(0,0,0,0),oma=c(8,6,4,6),xpd=NA)
  matplot(1:6,bias.cond[1:6,1:4],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[1:6,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8), 
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,1:4],type="l",xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[7:12,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  matplot(1:6,bias.cond[13:18,1:4],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,5:8],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  matplot(1:6,bias.cond[13:18,9:11],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray(c(1,3,4)/8),
          xlim=c(0.7,6.2),ylim=c(ylimits.down[j-8],ylimits.up[j-8]),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
  abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
  axis(1, at=1:6, labels=FALSE)
  text(x=1:6, y=par()$usr[3]+ylimits.down[j-8]*xlabel.factor[j-8],
       labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
  
  legend(-2,ylimits.down[j-8]*legend.factor[j-8],c("Wide-ML","Long-ML", "Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
         lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1.1)
  mtext(ylabel[j-8], 2, outer=TRUE, line=3.5)
  mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
  mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
  mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
  mtext("Model", 3, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 3, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 3, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 3, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("Base Model", 4, at=3/6,cex=1, outer=TRUE, line=2)
  mtext("MSST", 4, at=5/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("CUTS", 4, at=3/6,cex=0.75, outer=TRUE, line=0.5)
  mtext("TSO", 4, at=1/6,cex=0.75, outer=TRUE, line=0.5)
  dev.off()
}
rm(j, bias.cond, legend.factor, xlabel.factor, ylabel, ylimits.down, ylimits.up)

# Plot: Unique Consistency Base Model MSST----

Ucon.cond <- matrix(NA, 6 * I, 3)
colnames(Ucon.cond) <- labels[6:8]
for(i in 1:6){
  for(j in 5:8){
    Ucon.subset <- subset(var_coeff, factor.var$cond == i)[,j]
    Ucon.subset <- matrix(Ucon.subset, 100, 11, byrow = TRUE)
    Ucon.cond[((j-5)*6)+i,] <- apply(Ucon.subset, 2, function(x) mean(x, na.rm = TRUE))[6:8]
  }
}
rm(i,j, Ucon.subset)
Ucon.cond[1:6,][which(perf.cond[1:6,6:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[7:12,][which(perf.cond[1:6,6:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[13:18,][which(perf.cond[1:6,6:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[19:24,][which(perf.cond[1:6,6:8] < 10, arr.ind = TRUE)] <- NA

ylabel <- "Unique Consistency Estimate"
ylimits.down <- -0.0125
ylimits.up <- 0.0125
xlabel.factor <- 0.1
legend.factor <- 1.54

pdf(paste0("Mplus_Simulation/Ucon_msst_plot.pdf"))
par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(8,6,2,2),xpd=NA)
matplot(1:6,Ucon.cond[1:6,],type="l",lty = c(5,2,6),xaxt="n",xlab="",ylab="", col = gray((2:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up), cex = 1.5, lwd = c(4.5,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y1", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Ucon.cond[7:12,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((2:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(5,2,6), cex = 1.5, lwd = c(4.5,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y2", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Ucon.cond[13:18,],type="l",xaxt="n",xlab="",ylab="",col = gray((2:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(5,2,6), cex = 1.5, lwd = c(4.5,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y3", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,Ucon.cond[19:24,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((2:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(5,2,6), cex = 1.5, lwd = c(4.5,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y4", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(1.5,ylimits.down*legend.factor,c("Long-ml", "Wide-Bayes","Long-Bayes"), col = gray((2:4)/8),
       lty = c(5,2,6),lwd = c(4.5,3,4),ncol=1, seg.len = 5, cex = 1)
mtext(ylabel, 2, outer=TRUE, line=4)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
dev.off()

rm(Ucon.cond,legend.factor,xlabel.factor,ylabel, ylimits.down, ylimits.up)

# Plot: Unique Consistency Base Model TSO----

Ucon.cond <- matrix(NA, 6 * I, 2)
colnames(Ucon.cond) <- labels[7:8]
for(i in 13:18){
  for(j in 5:8){
    Ucon.subset <- subset(var_coeff, factor.var$cond == i)[,j]
    Ucon.subset <- matrix(Ucon.subset, 100, 11, byrow = TRUE)
    Ucon.cond[((j-5)*6)+(i-12),] <- apply(Ucon.subset, 2, function(x) mean(x, na.rm = TRUE))[7:8]
  }
}
rm(i,j, Ucon.subset)
Ucon.cond[1:6,][which(perf.cond[13:18,7:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[7:12,][which(perf.cond[13:18,7:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[13:18,][which(perf.cond[13:18,7:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[19:24,][which(perf.cond[13:18,7:8] < 10, arr.ind = TRUE)] <- NA

ylabel <- "Unique Consistency Estimate"
ylimits.down <- -0.01
ylimits.up <- 0.225
xlabel.factor <- 0.9
legend.factor <- 8

pdf(paste0("Mplus_Simulation/Ucon_tso_plot.pdf"))
par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(8,6,2,2),xpd=NA)
matplot(1:6,Ucon.cond[1:6,],type="l",lty = c(2,6),xaxt="n",xlab="",ylab="", col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y1", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Ucon.cond[7:12,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(2,6), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y2", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Ucon.cond[13:18,],type="l",xaxt="n",xlab="",ylab="",col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(2,6), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y3", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,Ucon.cond[19:24,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(2,6), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y4", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(1.5,ylimits.down*legend.factor,c("Wide-Bayes","Long-Bayes"), col = gray((3:4)/8),
       lty = c(2,6),lwd = c(3,4),ncol=1, seg.len = 5, cex = 1)
mtext(ylabel, 2, outer=TRUE, line=4)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
dev.off()

rm(Ucon.cond,legend.factor,xlabel.factor,ylabel, ylimits.down, ylimits.up)

# Plot: Unpredictability Base Model MSST----

Upred.cond <- matrix(NA, 6 * I, 2)
colnames(Upred.cond) <- labels[10:11]
for(i in 1:6){
  for(j in 5:8){
    Upred.subset <- subset(var_coeff, factor.var$cond == i)[,j]
    Upred.subset <- matrix(Upred.subset, 100, 11, byrow = TRUE)
    Upred.cond[((j-5)*6)+i,] <- apply(Upred.subset, 2, function(x) mean(x, na.rm = TRUE))[10:11]
  }
}
rm(i,j, Upred.subset)
Upred.cond[1:6,][which(perf.cond[1:6,10:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[7:12,][which(perf.cond[1:6,10:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[13:18,][which(perf.cond[1:6,10:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[19:24,][which(perf.cond[1:6,10:11] < 10, arr.ind = TRUE)] <- NA

ylabel <- "Trait Unpredictability Estimate"
ylimits.down <- -0.000125
ylimits.up <- 0.00125
xlabel.factor <- 0.55
legend.factor <- 4.5

pdf(paste0("Mplus_Simulation/Upred_msst_plot.pdf"))
par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(8,6,2,2),xpd=NA)
matplot(1:6,Upred.cond[1:6,],type="l",lty = c(2,6),xaxt="n",xlab="",ylab="", col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y1", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Upred.cond[7:12,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(2,6), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y2", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Upred.cond[13:18,],type="l",xaxt="n",xlab="",ylab="",col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(2,6), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y3", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,Upred.cond[19:24,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((3:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(2,6), cex = 1.5, lwd = c(3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y4", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(1.5,ylimits.down*legend.factor,c("Wide-Bayes","Long-Bayes"), col = gray((3:4)/8),
       lty = c(2,6),lwd = c(3,4),ncol=1, seg.len = 5, cex = 1)
mtext(ylabel, 2, outer=TRUE, line=4)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
dev.off()

rm(Upred.cond,legend.factor,xlabel.factor,ylabel, ylimits.down, ylimits.up)

# Plot: Unpredictability Base Model CUTS----

Upred.cond <- matrix(NA, 6 * I, 3)
colnames(Upred.cond) <- labels[9:11]
for(i in 7:12){
  for(j in 5:8){
    Upred.subset <- subset(var_coeff, factor.var$cond == i)[,j]
    Upred.subset <- matrix(Upred.subset, 100, 11, byrow = TRUE)
    Upred.cond[((j-5)*6)+(i-6),] <- apply(Upred.subset, 2, function(x) mean(x, na.rm = TRUE))[9:11]
  }
}
rm(i,j, Upred.subset)
Upred.cond[1:6,][which(perf.cond[7:12,9:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[7:12,][which(perf.cond[7:12,9:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[13:18,][which(perf.cond[7:12,9:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[19:24,][which(perf.cond[7:12,9:11] < 10, arr.ind = TRUE)] <- NA

ylabel <- "Trait Unpredictability Estimate"
ylimits.down <- -0.000125
ylimits.up <- 0.00125
xlabel.factor <- 0.55
legend.factor <- 3.9

pdf(paste0("Mplus_Simulation/Upred_cuts_plot.pdf"))
par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(8,6,2,2),xpd=NA)
matplot(1:6,Upred.cond[1:6,],type="l",lty = c(1,2,6),xaxt="n",xlab="",ylab="", col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y1", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Upred.cond[7:12,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y2", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Upred.cond[13:18,],type="l",xaxt="n",xlab="",ylab="",col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y3", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,Upred.cond[19:24,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y4", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(1.5,ylimits.down*legend.factor,c("Wide-ml", "Wide-Bayes","Long-Bayes"), col = gray((c(1,3,4))/8),
       lty = c(1,2,6),lwd = c(3,3,4),ncol=1, seg.len = 5, cex = 1)
mtext(ylabel, 2, outer=TRUE, line=4)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
dev.off()

rm(Upred.cond,legend.factor,xlabel.factor,ylabel, ylimits.down, ylimits.up)

# Plot: Unique Consistency Bias----

Ucon.cond <- matrix(NA, 6 * I, 4)
colnames(Ucon.cond) <- labels[5:8]
for(i in 7:12){
  for(j in 5:8){
    Ucon.subset <- subset(var_coeff.bias, factor.var$cond == i)[,j]
    Ucon.subset <- matrix(Ucon.subset, 100, 11, byrow = TRUE)
    Ucon.cond[((j-5)*6)+(i-6),] <- apply(Ucon.subset, 2, function(x) mean(x, na.rm = TRUE))[5:8]
  }
}
rm(i,j, Ucon.subset)
Ucon.cond[1:6,][which(perf.cond[7:12,5:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[7:12,][which(perf.cond[7:12,5:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[13:18,][which(perf.cond[7:12,5:8] < 10, arr.ind = TRUE)] <- NA
Ucon.cond[19:24,][which(perf.cond[7:12,5:8] < 10, arr.ind = TRUE)] <- NA

ylabel <- "Bias Unique Consistency"
ylimits.down <- -0.0165
ylimits.up <- 0.0165
xlabel.factor <- 0.1
legend.factor <- 1.625

pdf(paste0("Mplus_Simulation/Ucon_bias_plot.pdf"))
par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(8,6,2,2),xpd=NA)
matplot(1:6,Ucon.cond[1:6,],type="l",lty = c(1,5,2,6),xaxt="n",xlab="",ylab="", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up), cex = 1.5, lwd = c(3,4.5,3,4), las =1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y1", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Ucon.cond[7:12,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las =1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y2", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Ucon.cond[13:18,],type="l",xaxt="n",xlab="",ylab="",col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las =1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y3", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,Ucon.cond[19:24,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((1:4)/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,5,2,6), cex = 1.5, lwd = c(3,4.5,3,4), las =1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y4", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(-0.3,ylimits.down*legend.factor,c("Wide-ml","Long-ml","Wide-Bayes","Long-Bayes"), col = gray((1:4)/8),
       lty = c(1,5,2,6),lwd = c(3,4.5,3,4),ncol=2, seg.len = 5, cex = 1)
mtext(ylabel, 2, outer=TRUE, line=4)
mtext("Number of Measurement Times", 1, at=1/5, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/5, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/5, outer=TRUE, line=6, cex=0.8)
dev.off()

rm(Ucon.cond,legend.factor,xlabel.factor,ylabel, ylimits.down, ylimits.up)

# Plot: Unpredictability Bias----

Upred.cond <- matrix(NA, 6 * I, 3)
colnames(Upred.cond) <- labels[9:11]
for(i in 13:18){
  for(j in 5:8){
    Upred.subset <- subset(var_coeff.bias, factor.var$cond == i)[,j]
    Upred.subset <- matrix(Upred.subset, 100, 11, byrow = TRUE)
    Upred.cond[((j-5)*6)+(i-12),] <- apply(Upred.subset, 2, function(x) mean(x, na.rm = TRUE))[9:11]
  }
}
rm(i,j, Upred.subset)
Upred.cond[1:6,][which(perf.cond[13:18,9:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[7:12,][which(perf.cond[13:18,9:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[13:18,][which(perf.cond[13:18,9:11] < 10, arr.ind = TRUE)] <- NA
Upred.cond[19:24,][which(perf.cond[13:18,9:11] < 10, arr.ind = TRUE)] <- NA

ylabel <- "Bias Trait Unpredictability"
ylimits.down <- -0.0165
ylimits.up <- 0.0165
xlabel.factor <- 0.1
legend.factor <- 1.525

pdf(paste0("Mplus_Simulation/Upred_bias_plot.pdf"))
par(mfrow=c(2,2),mar=c(0,0,0,0),oma=c(8,6,2,2),xpd=NA)
matplot(1:6,Upred.cond[1:6,],type="l",lty = c(1,2,6),xaxt="n",xlab="",ylab="", col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y1", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Upred.cond[7:12,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y2", 3, outer = FALSE, line = -1.5, cex = 1)
matplot(1:6,Upred.cond[13:18,],type="l",xaxt="n",xlab="",ylab="",col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y3", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)
matplot(1:6,Upred.cond[19:24,],type="l",xaxt="n",xlab="",ylab="",yaxt="n", col = gray((c(1,3,4))/8),
        xlim=c(0.7,6.2),ylim=c(ylimits.down,ylimits.up),lty = c(1,2,6), cex = 1.5, lwd = c(3,3,4), las = 1)
abline(h = 0, xpd = FALSE, col = rgb(.211, .211, .211, .25))
mtext("Y4", 3, outer = FALSE, line = -1.5, cex = 1)
axis(1, at=1:6, labels=FALSE)
text(x=1:6, y=par()$usr[3]+ylimits.down*xlabel.factor,
     labels=c("30-0%","60-0%","90-0%", "30-10%","60-10%","90-10%"), srt=45, adj=1)

legend(1.5,ylimits.down*legend.factor,c("Wide-ml", "Wide-Bayes","Long-Bayes"), col = gray((c(1,3,4))/8),
       lty = c(1,2,6),lwd = c(3,3,4),ncol=1, seg.len = 5, cex = 1.1)
mtext(ylabel, 2, outer=TRUE, line=4)
mtext("Number of Measurement Times", 1, at=1/4, outer=TRUE, line=4, cex=0.8)
mtext(expression(""%*%""), 1, at=1/4, outer=TRUE, line=5, cex=0.8)
mtext("Percentage of Missingness", 1, at=1/4, outer=TRUE, line=6, cex=0.8)
dev.off()

rm(Upred.cond,legend.factor,xlabel.factor,ylabel, ylimits.down, ylimits.up)


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
fit_measures[which(fit_measures == -999, arr.ind = TRUE)] <- NA

# Plot: DIC selection ----

fit.AIC <- matrix(fit_measures[,1], R * 18, 11, byrow = TRUE) 
fit.AIC <- fit.AIC[,c(1,2,5,6,9)]
AIC.select <- unlist(apply(fit.AIC, 1, function(x) which(x == min(x, na.rm = TRUE))[1]))
AIC.select <- factor(AIC.select, levels = 1:5, labels = c("MSST", "ML-MSST", "CUTS", "ML-CUTS", "TSO"))

fit.BIC <- matrix(fit_measures[,2], R * 18, 11, byrow = TRUE) 
fit.BIC <- fit.BIC[,c(1,2,5,6,9)]
BIC.select <- unlist(apply(fit.BIC, 1, function(x) which(x == min(x, na.rm = TRUE))[1]))
BIC.select <- factor(BIC.select, levels = 1:5, labels = c("MSST", "ML-MSST", "CUTS", "ML-CUTS", "TSO"))

fit.aBIC <- matrix(fit_measures[,2], R * 18, 11, byrow = TRUE) 
fit.aBIC <- fit.aBIC[,c(1,2,5,6,9)]
aBIC.select <- unlist(apply(fit.aBIC, 1, function(x) which(x == min(x, na.rm = TRUE))[1]))
aBIC.select <- factor(aBIC.select, levels = 1:5, labels = c("MSST", "ML-MSST", "CUTS", "ML-CUTS", "TSO"))

fit.DIC <- matrix(fit_measures[,4], R * 18, 11, byrow = TRUE) 
fit.DIC <- fit.DIC[,c(4,8,11)]
DIC.select <- unlist(apply(fit.DIC, 1, function(x) which(x == min(x, na.rm = TRUE))))
DIC.select <- factor(DIC.select, levels = 1:3, labels = c("ML-MSST", "ML-CUTS", "ML-TSO"))

mle.models <- sort(c("MSST", "ML-MSST", "CUTS", "ML-CUTS", "TSO"))

AIC.select <- matrix(AIC.select, 600, 3, byrow = FALSE)
X <- apply(AIC.select, 2, function(x) summary(as.factor(x)))

AIC.cond <- matrix(0, 3, 5)
colnames(AIC.cond) <- mle.models

for(i in 1:3){
  AIC.cond[i, which(mle.models%in%names(X[[i]]))]<-X[[i]]
}
rm(X)

BIC.select <- matrix(BIC.select, 600, 3, byrow = FALSE)
X <- apply(BIC.select, 2, function(x) summary(as.factor(x)))

BIC.cond <- matrix(0, 3, 5)
colnames(BIC.cond) <- mle.models

for(i in 1:3){
  BIC.cond[i, which(mle.models%in%names(X[[i]]))]<-X[[i]]
}
rm(X)

aBIC.select <- matrix(aBIC.select, 600, 3, byrow = FALSE)
X <- apply(aBIC.select, 2, function(x) summary(as.factor(x)))

aBIC.cond <- matrix(0, 3, 5)
colnames(aBIC.cond) <- mle.models

for(i in 1:3){
  aBIC.cond[i, which(mle.models%in%names(X[[i]]))]<-X[[i]]
}
rm(X)

mcmc.models <- sort(c("ML-MSST", "ML-CUTS", "ML-TSO"))

DIC.select <- matrix(DIC.select, 600, 3, byrow = FALSE)
X <- apply(DIC.select, 2, function(x) summary(as.factor(x)))

DIC.cond <- matrix(0, 3, 3)
colnames(DIC.cond) <- mcmc.models

DIC.cond[,3] <- X
rm(X)

AIC.cond <- AIC.cond[, c(4,3,1,2,5)]
BIC.cond <- BIC.cond[, c(4,3,1,2,5)]
aBIC.cond <- aBIC.cond[, c(4,3,1,2,5)]

DIC.cond <- DIC.cond[, c(2,1,3)]

mle.ic <- rbind(AIC.cond, BIC.cond, aBIC.cond)
mle.ic <- round(prop.table(mle.ic, 1)*100, 1) 
mle.ic <- mle.ic[c(1,4,7,2,5,8,3,6,9),]
mle.ic <- cbind(rep(c("MSST", "CUTS", "TSO"), each = 3),
                rep(c("AIC", "BIC", "aBIC"), times = 3),
                mle.ic)
colnames(mle.ic)[1:2]<- c("Base Model", "Information Criterion")

print(xtable(mle.ic, type = "latex", align = c("l", "l", "l", "c", "c", "c", "c", "c"), label = "tab:fitmle", 
             caption = "Percentage of Times a Model was Selected as Best Model According to the Information Criteria Indices Available with Maximum Likelihood Estimation"), 
      include.rownames = FALSE, caption.placement = "top", file = "Mplus_Simulation/fitmle.txt")


AIC.delta <- matrix(apply(fit.AIC, 1, function(x) max(x, na.rm = TRUE)[1] - min(x, na.rm = TRUE)[1]), 100, 18, byrow = FALSE)

apply(AIC.delta, 2, mean)

# End ----