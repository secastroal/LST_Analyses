# Contents
# 0.0 Prepare environment
# 1.0 Read data
# 1.1 Selection of measurements and variables, turn data into wide format
# 2.0 Fit multistate-singletrait with the lsttheory package
# 2.1 Plotting variance components and factor scores
# 3.0 Fit multistate-singletrait with Mplus
# 4.0 Fit trait-state-occasion model with Mplus
# 5.0 Fit common and unique trait state model with Mplus
# 6.0 mplus to lavaan


# 0.0 Prepare environment ----
rm(list=ls())
#library(devtools) # if lsttheory of sumplement C of Steyer et al. (2015) has not been installed before.
#install_github("amayer2010/lsttheory", force = TRUE)
library(lavaan)
library(lsttheory)
library(MplusAutomation)
source("msst_to_Mplus.R")
source("tso_to_Mplus.R")
source("cuts_to_Mplus.R")


# 1.0 Read data  ----
# This data is from an ESM study (Geschwind N, Peeters F, Drukker M, van Os J, Wichers M, 2011).
# The study followed 129 participants during 6 days in the baseline and during 6 days
# after intervention. Respondents had to fill the ESM self-assessment form 10 times a day.
# The variables assessed in the ESM form were cheerful, relaxed, worry, sad, fearful, and
# pleasantness of the event.

PlosOne <- read.table("Data S1.txt", header = TRUE, sep = ",", row.names = NULL)
PlosOne <- PlosOne[, 2:13]

# rename variables
names(PlosOne) <- c("subjn", "dayn", "beepn", "inform", "period",
                    "cheer", "eventp", "worry", "fear", "sad", "relax", "neur")

# 1.1 Selection of measurements and variables, turn data into wide format ----

# All 11 rows were discarded because most of them were NA
PlosOne_C <- PlosOne[-seq(from = 11, to = 28600, by = 11),]
# For these analyses, we are only interested on the baseline measurements.
PlosOne_C <- subset(PlosOne_C, PlosOne_C$period == 0)
# Measurements after day 6 are discarded.
PlosOne_C <- subset(PlosOne_C, PlosOne_C$dayn <= 6)
# A variable time is added, in order to be able to turm the data into wide format
time <- rep(1:60, 130)
PlosOne_C <- cbind(PlosOne_C, time)
rm(time)

# Only the variables related to negative affect are selected to try LST models.
PlosOne_W <- reshape(PlosOne_C[, c(1, 8:10, 13)], v.names=c("worry", "fear", "sad"),
                     timevar = "time", idvar="subjn", direction="wide")
# Person 74 is dicarded due to complete non-response, also the variable subject number is discarded.
PlosOne_W <- PlosOne_W[-74, -1]

# Change variable names to delete "." given that this is not allowed in Mplus syntax.
names(PlosOne_W) <- gsub("\\.", "", names(PlosOne_W))

# 2.0 Fit multistate-singletrait with the lsttheory package ----

# data subset to allow for several analyses with diferent number of measurement occasions.

m <- 30 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

# multistate-singletrait analysis
time0 <- proc.time()
msst <- lsttheory(m, 1, data = PlosOne_Wr, 
                  equiv.assumption=list(tau="cong", theta="equi"), 
                  scale.invariance = list(lait0 = T, lait1 = T, lat0 = T, lat1 = T),
                  missing = "ml")
time.msst.lavaan <- proc.time() - time0
rm(time0)

# print lsttheory object returns the variance components: reliability, consistency, and 
# occasion specifity   
msst 

# returns the parameter estimates matrix with se, significance test, and lower 
# and upper ci bounds.
parameterEstimates(msst@lavaanres) 

# write lavaan syntax (optional)
writeLines(msst@lavaansyntax, "lavaan_files/PlosOne_msst1.txt")

# Predict factor scores, it returns a matrix of dimensions nxm. With n equal to 
# the sample size and m equal to the number of states (measurement occasions) plus 1.
# These are the factor score for the latent state variables and for the latent
# trait variable.
lavPredict(msst@lavaanres)

# 2.1 Plotting variance components and factor scores ----

# count number of missings per person
count_na <- as.vector(apply(PlosOne_Wr, 1, function(x) sum(is.na(x))))
cbind(1:129, count_na)


i <- 12 # select person to plot
# plot worry observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m], 
     type = "l", ylim = c(0,7.5), col = "blue")# traceplot of the latent state scores of one person
points(1:m, PlosOne_W[i, seq(1,t.m, by = 3)], type = "b", col = "red")
abline(h =lavPredict(msst@lavaanres)[i,(m+1)]) #plot latent trait factor score
abline(h = mean(as.numeric(PlosOne_W[i, seq(1,t.m, by = 3)]), na.rm = T), col="green") # plot mean of worry of person i
abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange") # plot combined mean of all observed variables of person i

# plot fear observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m] * parameterEstimates(msst@lavaanres)[seq(2, 90, by = 3), 5], 
     type = "l", ylim = c(0,7.5), col = "blue")# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(2,t.m, by = 3)], type = "b", col = "red")
abline(h =lavPredict(msst@lavaanres)[i,(m+1)])
abline(h = mean(as.numeric(PlosOne_W[i, seq(2,t.m, by = 3)]), na.rm = T), col="green")
abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")


# plot sad observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m] * parameterEstimates(msst@lavaanres)[seq(3, 90, by = 3), 5], 
     type = "l", ylim = c(0,7.5), col = "blue")# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(3,t.m, by = 3)], type = "b", col = "red")
abline(h =lavPredict(msst@lavaanres)[i,(m+1)])
abline(h = mean(as.numeric(PlosOne_W[i, seq(3,t.m, by = 3)]), na.rm = T), col="green")
abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")


# plot variance components by variable.
var.comp <- as.data.frame(msst@lsttheory)
matplot(var.comp[seq(1, t.m, by= 3),], type = "l") # worry variance components
matplot(var.comp[seq(2, t.m, by= 3),], type = "l") # fear variance components
matplot(var.comp[seq(3, t.m, by= 3),], type = "l") # sad variance components


# 3.0 Fit multistate-singletrait with Mplus ----
# To work on this section, the complete version of Mplus should be installed.

#  Prepare data and syntax

m <- 30 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

msst.file.name <- paste0("msst", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                  equiv.assumption = list(tau = "cong", theta = "equi"))

write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus <- proc.time() - time0
rm(time0)

# Read Mplus output into R
fit.msst <- readModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".out"))
fit.msst$parameters$unstandardized # parameter estimates

# 4.0 Fit trait-state-occasion model with Mplus ----
# To work on this section, the complete version of Mplus should be installed.

#  Prepare data and syntax

m <- 30 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

tso.file.name <- paste0("tso", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",tso.file.name,".dat"), inpfile = T)

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, trait.indicators = TRUE,
                                 equiv.assumption = list(occ = "cong", theta = "equi"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)

write(tso_syntax, paste0("Mplus_files/",tso.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"))
time.tso.mplus <- proc.time() - time0
rm(time0)

# Read Mplus output into R
fit.tso <- readModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".out"))
fit.tso$parameters$unstandardized # parameter estimates

# 5.0 Fit common and unique trait state model with Mplus ----
# To work on this section, the complete version of Mplus should be installed.

#  Prepare data and syntax

m <- 30 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

cuts.file.name <- paste0("cuts", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
                                   scale.invariance = list(int = FALSE, lambda = FALSE),
                                   state.trait.invariance = FALSE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )
  
write(cuts_syntax, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",cuts.file.name,".inp"))
time.cuts.mplus <- proc.time() - time0
rm(time0)

# Read Mplus output into R
fit.cuts <- readModels(paste0(getwd(),"/Mplus_files/",cuts.file.name,".out"))
fit.cuts$parameters$unstandardized # parameter estimates


# 7.0 mplus to lavaan ----
# In this section, we use the function mplus2lavaan (lavaan package) to verify 
# that both programs lead to the same solutions for the TSO and the CUTS models.
# This section is buggy, mplus2lavaan function is failing in reading the data 
# from the mplus input file, and also it seems it is not converting mplus syntax
# successfully into lavaan syntax.

# TSO

fit.tso.lavaan <- mplus2lavaan(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"),
                               run = TRUE)

# CUTS

fit.cuts.lavaan <- mplus2lavaan(paste0(getwd(),"/Mplus_files/",cuts.file.name,".inp"),
                               run = TRUE)






