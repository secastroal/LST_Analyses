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
# 6.1 fitting the TSO with lavaan
# 6.2 fitting the CUTS with lavaan
# 7.0 Computing and plotting variance components of TSO and CUTS
# 7.1 Selected TSO model 
# 7.2 Selected CUTS model with om 




# 0.0 Prepare environment ----
rm(list=ls())
#install.packages("lavaan")
#install.packages("MplusAutomation")
#install.packages("devtools")
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
# h1: congenericity in states and congenericity in traits
# h2: congenericity in states and essential equivalence in traits
# h3: congenericity in states and equivalence in traits
# h4: essential equivalence in states and congenericity in traits
# h5: essential equivalence in states and essential equivalence in traits
# h6: essential equivalence in states and equivalence in traits
# h7: equivalence in states and congenericity in traits
# h8: equivalence in states and essential equivalence in traits
# h9: equivalence in states and equivalence in traits


time0 <- proc.time()
msst.h1 <- lsttheory(m, 1, data = PlosOne_Wr, 
                  equiv.assumption=list(tau="cong", theta="cong"), 
                  scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                  missing = "ml")
time.msst.h1 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h2 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="cong", theta="ess"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h2 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h3 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="cong", theta="equi"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h3 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h4 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="ess", theta="cong"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h4 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h5 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="ess", theta="ess"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h5 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h6 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="ess", theta="equi"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h6 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h7 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="equi", theta="cong"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h7 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h8 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="equi", theta="ess"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h8 <- proc.time() - time0
rm(time0)

time0 <- proc.time()
msst.h9 <- lsttheory(m, 1, data = PlosOne_Wr, 
                     equiv.assumption=list(tau="equi", theta="equi"), 
                     scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F),
                     missing = "ml")
time.msst.h9 <- proc.time() - time0
rm(time0)


# When 60 measurements are included, the analysis gives the following warning message:
# In lav_data_full(data = data, group = group, cluster = cluster,  :
#                   lavaan WARNING: small number of observations (nobs < nvar)
#                 nobs = 129 nvar = 180
# However, the model actually converge and give estimates.

### Optionally the analysis can be saved to don't run it again

save(msst.h9, file = paste0("lavaan_files/msst_h9_m",m, ".R"))

### load all 9 msst analyses
for(i in 1:9){
  load( file = paste0("lavaan_files/msst_h", i, "_m",m, ".R"))
}
rm(i)

# print lsttheory object returns the variance components: reliability, consistency, and 
# occasion specifity   
summary(msst.h9@lavaanres, fit.measures = TRUE, estimates = FALSE)

anova(msst.h1@lavaanres,msst.h9@lavaanres)

# The best msst model is the msst h2: congenericity in states and 
# essential equivalence in traits. Therefore, the factor scores and
# variance coeffcients are computed based on this model.
msst <- msst.h2

# returns the parameter estimates matrix with se, significance test, and lower 
# and upper ci bounds.
parameterEstimates(msst@lavaanres) 

# Predict factor scores, it returns a matrix of dimensions nxm. With n equal to 
# the sample size and m equal to the number of states (measurement occasions) plus 1.
# These are the factor score for the latent state variables and for the latent
# trait variable.
lavPredict(msst@lavaanres)

# 2.1 Plotting variance components and factor scores ----

# count number of missings per person
count_na <- as.vector(apply(PlosOne_Wr, 1, function(x) sum(is.na(x))))
cbind(1:129, count_na)
(sum(is.na(PlosOne_Wr)))/(dim(PlosOne_Wr)[1]*dim(PlosOne_Wr)[2]) # % of missingess

sample(1:129, 1)
i <- 53 # select person to plot
# plot worry observed vs latent state variables person i
par(mfrow= c(3,1),mar= c(0,5,0,2), oma= c(4,0,1,8), xpd = NA)

plot(lavPredict(msst@lavaanres)[i,1:m], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xaxt="n",xlab="", cex.lab = 1.5, ylab = "Worry", las = 1)# traceplot of the latent state scores of one person
points(1:m, PlosOne_W[i, seq(1,t.m, by = 3)], type = "b", col = "red", pch = 16)
abline(h =lavPredict(msst@lavaanres)[i,(m+1)], xpd = FALSE) #plot latent trait factor score
#abline(h = mean(as.numeric(PlosOne_W[i, seq(1,t.m, by = 3)]), na.rm = T), 
 #      col="green", xpd = TRUE) # plot mean of worry of person i
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange", 
 #      xpd = TRUE) # plot combined mean of all observed variables of person i

# plot fear observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m] * parameterEstimates(msst@lavaanres)[seq(2, t.m, by = 3), 5], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xaxt="n",xlab="", cex.lab = 1.5,ylab = "Fear", las = 1)# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(2,t.m, by = 3)], type = "b", col = "red", pch = 16)
abline(h =lavPredict(msst@lavaanres)[i,(m+1)], xpd = FALSE)
#abline(h = mean(as.numeric(PlosOne_W[i, seq(2,t.m, by = 3)]), na.rm = T), col="green")
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")
legend(31.5, 5, legend = c("Observed", "Predicted"), col = c( "red", "blue"),
       lty = 1, pch = 16, lwd = 2, cex = 1)


# plot sad observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m] * parameterEstimates(msst@lavaanres)[seq(3, t.m, by = 3), 5], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xlab="Measurement occasion", cex.lab = 1.5, ylab = "Sad", las = 1)# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(3,t.m, by = 3)], type = "b", pch= 16, col = "red")
abline(h =lavPredict(msst@lavaanres)[i,(m+1)], xpd = FALSE)
#abline(h = mean(as.numeric(PlosOne_W[i, seq(3,t.m, by = 3)]), na.rm = T), col="green")
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")


# plot variance components by variable.
var.comp <- as.data.frame(msst@lsttheory)

par(mfrow= c(3,1),mar= c(0,5,0,2), oma= c(4,0,1,6), xpd = NA)
matplot(var.comp[seq(1, t.m, by= 3),], type = "b", ylim = c(0, 1), col = c("black", "blue", "red"),
        ylab = "Worry", xaxt="n",xlab="", cex.lab = 1.5, lty = 1, pch =16, las = 1) # worry variance components

matplot(var.comp[seq(2, t.m, by= 3),], type = "b", ylim = c(0, 1), col = c("black", "blue", "red"),
        ylab = "Fear",xaxt="n",xlab="", cex.lab = 1.5, lty = 1, pch =16, las = 1) # fear variance components
legend(32, 0.75, legend = c("Rel", "Spe", "Con"), col = c("black", "blue", "red"),
       lty = 1, lwd = 2, pch = 16, cex = 1)

matplot(var.comp[seq(3, t.m, by= 3),], type = "b",ylim = c(0, 1), col = c("black", "blue", "red"),
        ylab = "Sad", xlab = "Measurement occasion", cex.lab = 1.5, 
        lty = 1, pch =16, las = 1) # sad variance components


# 3.0 Fit multistate-singletrait with Mplus ----
# To work on this section, the complete version of Mplus should be installed.

#  Prepare data and syntax

m <- 30 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

# h1: congenericity in states and congenericity in traits

msst.file.name <- paste0("msst_h1_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                  equiv.assumption = list(tau = "cong", theta = "cong"),
                                  scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h1 <- proc.time() - time0
rm(time0)

# h2: congenericity in states and essential equivalence in traits

msst.file.name <- paste0("msst_h2_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "cong", theta = "ess"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h2 <- proc.time() - time0
rm(time0)

# h3: congenericity in states and equivalence in traits

msst.file.name <- paste0("msst_h3_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "cong", theta = "equi"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h3 <- proc.time() - time0
rm(time0)

# h4: essential equivalence in states and congenericity in traits

msst.file.name <- paste0("msst_h4_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "ess", theta = "cong"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h4 <- proc.time() - time0
rm(time0)

# h5: essential equivalence in states and essential equivalence in traits

msst.file.name <- paste0("msst_h5_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "ess", theta = "ess"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h5 <- proc.time() - time0
rm(time0)

# h6: essential equivalence in states and equivalence in traits

msst.file.name <- paste0("msst_h6_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "ess", theta = "equi"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h6 <- proc.time() - time0
rm(time0)

# h7: equivalence in states and congenericity in traits

msst.file.name <- paste0("msst_h7_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "equi", theta = "cong"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h7 <- proc.time() - time0
rm(time0)

# h8: equivalence in states and essential equivalence in traits

msst.file.name <- paste0("msst_h8_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "equi", theta = "ess"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h8 <- proc.time() - time0
rm(time0)

# h9: equivalence in states and equivalence in traits

msst.file.name <- paste0("msst_h9_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",msst.file.name,".dat"), inpfile = T)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

msst_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                   equiv.assumption = list(tau = "equi", theta = "equi"),
                                   scale.invariance = list(lait0 = F, lait1 = F, lat0 = F, lat1 = F))

write(iter, paste0("Mplus_files/",msst.file.name,".inp"), append = T) # Write Analysis specifications
write(msst_syntax, paste0("Mplus_files/",msst.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",msst.file.name,".inp"))
time.msst.mplus.h9 <- proc.time() - time0
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

# h1: congenericity in the occasions and congenericity in the trait

tso.file.name <- paste0("tso3b_h1_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",tso.file.name,".dat"), inpfile = T)

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "cong", theta = "cong"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)

iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",tso.file.name,".inp"), append = T) 
write(tso_syntax, paste0("Mplus_files/",tso.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"))
time.tso.mplus.h1 <- proc.time() - time0
rm(time0)

# h2: congenericity in the occasions and equivalence in the trait

tso.file.name <- paste0("tso3b_h2_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",tso.file.name,".dat"), inpfile = T)

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "cong", theta = "equi"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",tso.file.name,".inp"), append = T) 
write(tso_syntax, paste0("Mplus_files/",tso.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"))
time.tso.mplus.h2 <- proc.time() - time0
rm(time0)

# h3: equivalence in the occasions and congenericity in the trait

tso.file.name <- paste0("tso3b_h3_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",tso.file.name,".dat"), inpfile = T)

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "equi", theta = "cong"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",tso.file.name,".inp"), append = T) 
write(tso_syntax, paste0("Mplus_files/",tso.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"))
time.tso.mplus.h3 <- proc.time() - time0
rm(time0)

# h4: equivalence in the occasions and equivalence in the trait

tso.file.name <- paste0("tso3b_h4_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",tso.file.name,".dat"), inpfile = T)

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "equi", theta = "equi"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",tso.file.name,".inp"), append = T) 
write(tso_syntax, paste0("Mplus_files/",tso.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"))
time.tso.mplus.h4 <- proc.time() - time0
rm(time0)

# h5: congenericity in the occasions and equivalence in the trait + homogeneity 
#     of autoregressive effects 

tso.file.name <- paste0("tso3b_h5_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",tso.file.name,".dat"), inpfile = T)

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "cong", theta = "equi"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = TRUE)
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",tso.file.name,".inp"), append = T) 
write(tso_syntax, paste0("Mplus_files/",tso.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"))
time.tso.mplus.h5 <- proc.time() - time0
rm(time0)

# h6: congenericity in the occasions and equivalence in the trait + scale 
#     invariance 

tso.file.name <- paste0("tso3b_h6_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",tso.file.name,".dat"), inpfile = T)

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "cong", theta = "equi"),
                                 scale.invariance = list(int = TRUE, lambda = TRUE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",tso.file.name,".inp"), append = T) 
write(tso_syntax, paste0("Mplus_files/",tso.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",tso.file.name,".inp"))
time.tso.mplus.h6 <- proc.time() - time0
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

# h1: Least retricted CUTS #! Does not really work, some loadings are extremely large.

cuts.file.name <- paste0("cuts_m-1_h1_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = FALSE, lambda = FALSE),
                                   state.trait.invariance = FALSE,
                                   fixed.method.loadings = FALSE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)  
write(cuts_syntax, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",cuts.file.name,".inp"))
time.cuts.mplus.h1 <- proc.time() - time0
rm(time0)

# h2: Cuts with all method loadings fixed to 1

cuts.file.name <- paste0("cuts_m-1_h2_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = FALSE, lambda = FALSE),
                                   state.trait.invariance = FALSE,
                                   fixed.method.loadings = TRUE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)  
write(cuts_syntax, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",cuts.file.name,".inp"))
time.cuts.mplus.h2 <- proc.time() - time0
rm(time0)

# h3: Cuts + weak factorial invariance across time for the CT and the CS 

cuts.file.name <- paste0("cuts_m-1_h3_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = TRUE, lambda = TRUE),
                                   state.trait.invariance = FALSE,
                                   fixed.method.loadings = TRUE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)  
write(cuts_syntax, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",cuts.file.name,".inp"))
time.cuts.mplus.h3 <- proc.time() - time0
rm(time0)

# h4: Cuts + weak factorial invariance across traits and states

cuts.file.name <- paste0("cuts_m-1_h4_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = FALSE, lambda = FALSE),
                                   state.trait.invariance = TRUE,
                                   fixed.method.loadings = TRUE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )
iter <- "
ANALYSIS:
H1iterations=30000;" # increase H1 iterations

write(iter, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)  
write(cuts_syntax, paste0("Mplus_files/",cuts.file.name,".inp"), append = T)

# run model in Mplus
time0 <- proc.time()
runModels(paste0(getwd(),"/Mplus_files/",cuts.file.name,".inp"))
time.cuts.mplus.h4 <- proc.time() - time0
rm(time0)


# Read Mplus output into R
fit.cuts <- readModels(paste0(getwd(),"/Mplus_files/",cuts.file.name,".out"))
fit.cuts$parameters$unstandardized # parameter estimates


# 6.0 mplus to lavaan ----
# In this section, we use the function mplus2lavaan (lavaan package) to verify 
# that both programs lead to the same solutions for the TSO and the CUTS models.
# This section is buggy, mplus2lavaan function is failing in reading the data 
# from the mplus input file, and also it seems it is not converting mplus syntax
# successfully into lavaan syntax. Instead, we used the mplus2lavaan.modelSyntax
# that only translate Mplus syntax in lavaan syntax.

# 6.1 fitting the TSO with lavaan ----

m <- 30 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

# h1: congenericity in the occasions and congenericity in the trait

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3a",
                                 equiv.assumption = list(occ = "cong", theta = "cong"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
tso.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", tso_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.tso.lavaan.h1 <- sem(model = tso.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.tso.lavaan.h1 <- proc.time() - time0
rm(time0)

save(fit.tso.lavaan.h1, file = paste0("lavaan_files/tso3a_h1_m",m, ".R"))


# h2: congenericity in the occasions and equivalence in the trait

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3a",
                                 equiv.assumption = list(occ = "cong", theta = "equi"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
tso.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", tso_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.tso.lavaan.h2 <- sem(model = tso.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.tso.lavaan.h2 <- proc.time() - time0
rm(time0)

save(fit.tso.lavaan.h2, file = paste0("lavaan_files/tso3a_h2_m",m, ".R"))

# h3: equivalence in the occasions and congenericity in the trait

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3a",
                                 equiv.assumption = list(occ = "equi", theta = "cong"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
tso.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", tso_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.tso.lavaan.h3 <- sem(model = tso.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.tso.lavaan.h3 <- proc.time() - time0
rm(time0)

save(fit.tso.lavaan.h3, file = paste0("lavaan_files/tso3a_h3_m",m, ".R"))

# h4: equivalence in the occasions and equivalence in the trait

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3a",
                                 equiv.assumption = list(occ = "equi", theta = "equi"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
tso.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", tso_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.tso.lavaan.h4 <- sem(model = tso.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.tso.lavaan.h4 <- proc.time() - time0
rm(time0)

save(fit.tso.lavaan.h4, file = paste0("lavaan_files/tso3a_h4_m",m, ".R"))

# load analyses
#figure 3a
for(i in 1:4){
  load(file = paste0("lavaan_files/tso3a_h", i, "_m",m, ".R"))
}

#figure 3b
for(i in 1:4){
  load(file = paste0("lavaan_files/tso3b_h", i, "_m",m, ".R"))
}

summary(fit.tso.lavaan.h4, estimates = FALSE, fit.measures = TRUE)

anova(fit.tso.lavaan.h1,fit.tso.lavaan.h4)

# h5: congenericity in the occasions and equivalence in the trait + homogeneity 
#     of autoregressive effects 

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "cong", theta = "equi"),
                                 scale.invariance = list(int = FALSE, lambda = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = TRUE)
tso.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", tso_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.tso.lavaan.h5 <- sem(model = tso.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.tso.lavaan.h5 <- proc.time() - time0
rm(time0)

save(fit.tso.lavaan.h5, file = paste0("lavaan_files/tso3b_h5_m",m, ".R"))

# h6: congenericity in the occasions and equivalence in the trait + scale 
#     invariance 

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
                                 equiv.assumption = list(occ = "cong", theta = "equi"),
                                 scale.invariance = list(int = TRUE, lambda = TRUE),
                                 homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                 autoregressive.homogeneity = FALSE)
tso.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", tso_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.tso.lavaan.h6 <- sem(model = tso.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.tso.lavaan.h6 <- proc.time() - time0
rm(time0)

save(fit.tso.lavaan.h6, file = paste0("lavaan_files/tso3b_h6_m",m, ".R"))


summary(fit.tso.lavaan.h6, estimates= FALSE, fit.measures = TRUE)

anova(fit.tso.lavaan.h2, fit.tso.lavaan.h6)

# 6.2 fitting the CUTS with lavaan -----

# The next warning messages were obtained when 60 measurements were used.
#Warning messages:
#   1: In lav_data_full(data = data, group = group, cluster = cluster,  :
#     lavaan WARNING: small number of observations (nobs < nvar)
#     nobs = 129 nvar = 180
#   2: In lav_object_post_check(object) :
#     lavaan WARNING: some estimated ov variances are negative

m <- 30 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

# h1: Least retricted CUTS

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = FALSE, lambda = FALSE),
                                   state.trait.invariance = FALSE,
                                   fixed.method.loadings = FALSE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )

cuts.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", cuts_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.cuts.lavaan.h1 <- sem(model = cuts.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.cuts.lavaan.h1 <- proc.time() - time0
rm(time0)

save(fit.cuts.lavaan.h1, file = paste0("lavaan_files/cuts_m-1_h1_m",m, ".R"))

# h2: Cuts with all method loadings fixed to 1

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = FALSE, lambda = FALSE),
                                   state.trait.invariance = FALSE,
                                   fixed.method.loadings = TRUE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )

cuts.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", cuts_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.cuts.lavaan.h2 <- sem(model = cuts.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.cuts.lavaan.h2 <- proc.time() - time0
rm(time0)

save(fit.cuts.lavaan.h2, file = paste0("lavaan_files/cuts_m-1_h2_m",m, ".R"))

# h3: Cuts + weak factorial invariance across time for the CT and the CS

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = TRUE, lambda = TRUE),
                                   state.trait.invariance = FALSE,
                                   fixed.method.loadings = TRUE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )

cuts.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", cuts_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.cuts.lavaan.h3 <- sem(model = cuts.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.cuts.lavaan.h3 <- proc.time() - time0
rm(time0)

save(fit.cuts.lavaan.h3, file = paste0("lavaan_files/cuts_m-1_h3_m",m, ".R"))

# h4: Cuts + weak factorial invariance across traits and states

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "m-1",
                                   scale.invariance = list(int = FALSE, lambda = FALSE),
                                   state.trait.invariance = TRUE,
                                   fixed.method.loadings = TRUE,
                                   homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                   fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) )

cuts.lavaan.syntax <- mplus2lavaan.modelSyntax(gsub("MODEL:", "", cuts_syntax))

# Fitting Cuts with lavaan
time0 <- proc.time()
fit.cuts.lavaan.h4 <- sem(model = cuts.lavaan.syntax, data = PlosOne_Wr, missing = "ml")
time.cuts.lavaan.h4 <- proc.time() - time0
rm(time0)

save(fit.cuts.lavaan.h4, file = paste0("lavaan_files/cuts_m-1_h4_m",m, ".R"))


# load analyses
# Orthogonal methods
for(i in 1:4){
  load(file = paste0("lavaan_files/cuts_h", i, "_m",m, ".R"))
}
rm(i)
# m-1
for(i in 1:4){
  load(file = paste0("lavaan_files/cuts_m-1_h", i, "_m",m, ".R"))
}
rm(i)

summary(fit.cuts.lavaan.h4, estimates = FALSE, fit.measures = TRUE)

anova(fit.cuts.lavaan.h2, fit.cuts.lavaan.h4)
 
# 7.0 Computing and plotting variance components of TSO and CUTS ----
# 7.1 Selected TSO model ----
m <- 30
t.m <- 3*m
tso  <- fit.tso.lavaan.h2

# 7.1.1 compute factor scores
trait.scores <- lavPredict(tso)[,(m+1):(m+3)] #scores trait
occ.scores <- lavPredict(tso)[,1:m] #scores occasion variables

# 7.1.2 Extract parameters from the lavaan object
ltrait <- parameterEstimates(tso)[ (t.m+1):(2*t.m),"est"] #loadings traits
locc <- parameterEstimates(tso)[ 1:t.m,"est"] #loadings occasion

tso.ins <- parameterEstimates(tso)[(2*t.m+1):(3*t.m), "est"] #intercepts

beta <- parameterEstimates(tso)[(3*t.m+((m+3)*(m+2)/2)+m-2):(3*t.m+((m+3)*(m+2)/2)+2*m-4),"est"]


err_var <- parameterEstimates(tso)[(3*t.m+((m+3)*(m+2)/2)+2*m-3):(4*t.m+((m+3)*(m+2)/2)+2*m-4), "est"]
occ_var <- parameterEstimates(tso)[(4*t.m+((m+3)*(m+2)/2)+2*m-3):(4*t.m+((m+3)*(m+2)/2)+3*m-4), "est"]
trait_var <- parameterEstimates(tso)[(4*t.m+((m+3)*(m+2)/2)+3*m-3):(4*t.m+((m+3)*(m+2)/2)+3*m-1), "est"]

# 7.1.3 Compute variance coeffcients

beta2 <- beta*beta
betaXocc <- rep(NA, (m-1))
# This loop compute the variance of each occasion variable as a linear 
#   combination of the previous state residuals
for(i in 1:(m-1)){
  beta2_prod <- rep(NA, i)
  for(j in 1:i){
    beta2_prod[j] <- prod(beta2[j:i])
  }
  betaXocc[i] <- sum(beta2_prod*occ_var[1:i])
}
rm(i, j, beta2_prod)

betaXocc <- c(0, betaXocc)
#total variance
tso_tot_var <- (ltrait^2)*trait_var + (locc^2)*rep(betaXocc, each = 3) + 
  (locc^2)*rep(occ_var, each = 3) + err_var
#reliability
tso_rel <- ((ltrait^2)*trait_var + (locc^2)*rep(betaXocc, each = 3) + 
              (locc^2)*rep(occ_var, each = 3))/tso_tot_var
#Consistency
tso_con <- ((ltrait^2)*trait_var + (locc^2)*rep(betaXocc, each = 3))/tso_tot_var
#Trait predictability
tso_pred <- ((ltrait^2)*trait_var)/tso_tot_var
#Trait unpredictability
tso_upred <- ((locc^2)*rep(betaXocc, each = 3))/tso_tot_var
#Occasion specificity
tso_spe <- ((locc^2)*rep(occ_var, each = 3))/tso_tot_var

# 7.1.4 Compute predicted values based on the factor scores
tso_pred_worry <- matrix(ltrait[seq(1, 90, by = 3)],129,m, byrow = TRUE)*trait.scores[,1] +
  matrix(locc[seq(1, 90, by = 3)],129,m, byrow = TRUE)*occ.scores +
  tso.ins[seq(1, 90, by = 3)]

tso_pred_fear <- matrix(ltrait[seq(2, 90, by = 3)],129,m, byrow = TRUE)*trait.scores[,2] +
  matrix(locc[seq(2, 90, by = 3)],129,m, byrow = TRUE)*occ.scores +
  tso.ins[seq(2, 90, by = 3)]

tso_pred_sad <- matrix(ltrait[seq(3, 90, by = 3)],129,m, byrow = TRUE)*trait.scores[,3] +
  matrix(locc[seq(3, 90, by = 3)],129,m, byrow = TRUE)*occ.scores +
  tso.ins[seq(3, 90, by = 3)]

# 7.2.5 Plots
i <- 53 # select person to plot
# plot worry observed vs latent state variables person i
par(mfrow= c(3,1),mar= c(0,5,0,2), oma= c(4,0,1,8), xpd = NA)

plot(tso_pred_worry[i,], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xaxt="n",xlab="", cex.lab = 1.5, ylab = "Worry", las = 1)# traceplot of the latent state scores of one person
points(1:m, PlosOne_W[i, seq(1,t.m, by = 3)], type = "b", col = "red", pch = 16)
abline(h =trait.scores[i,1], xpd = FALSE) #plot latent trait factor score
#abline(h = mean(as.numeric(PlosOne_W[i, seq(1,t.m, by = 3)]), na.rm = T), 
#      col="green", xpd = TRUE) # plot mean of worry of person i
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange", 
#      xpd = TRUE) # plot combined mean of all observed variables of person i

# plot fear observed vs latent state variables person i
plot(tso_pred_fear[i,], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xaxt="n",xlab="", cex.lab = 1.5,ylab = "Fear", las = 1)# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(2,t.m, by = 3)], type = "b", col = "red", pch = 16)
abline(h =trait.scores[i,2], xpd = FALSE)
#abline(h = mean(as.numeric(PlosOne_W[i, seq(2,t.m, by = 3)]), na.rm = T), col="green")
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")
legend(31.5, 5, legend = c("Observed", "Predicted"), col = c( "red", "blue"),
       lty = 1, pch = 16, lwd = 2, cex = 1)


# plot sad observed vs latent state variables person i
plot(tso_pred_sad[i,], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xlab="Measurement occasion", cex.lab = 1.5, ylab = "Sad", las = 1)# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(3,t.m, by = 3)], type = "b", pch= 16, col = "red")
abline(h =trait.scores[i,3], xpd = FALSE)
#abline(h = mean(as.numeric(PlosOne_W[i, seq(3,t.m, by = 3)]), na.rm = T), col="green")
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")


# plot variance components by variable.
var.comp <- cbind(tso_rel, tso_pred, tso_upred, tso_spe)

par(mfrow= c(3,1),mar= c(0,5,0,2), oma= c(4,0,1,6), xpd = NA)
matplot(var.comp[seq(1, t.m, by= 3),], type = "b", ylim = c(0, 1), col = c("black", "red", "green", "blue"),
        ylab = "Worry", xaxt="n",xlab="", cex.lab = 1.5, lty = 1, pch =16, las = 1) # worry variance components

matplot(var.comp[seq(2, t.m, by= 3),], type = "b", ylim = c(0, 1), col = c("black", "red", "green", "blue"),
        ylab = "Fear",xaxt="n",xlab="", cex.lab = 1.5, lty = 1, pch =16, las = 1) # fear variance components
legend(32, 0.75, legend = c("Rel", "Pred", "UPred", "Spe"), col = c("black", "red", "green", "blue"),
       lty = 1, lwd = 2, pch = 16, cex = 1)

matplot(var.comp[seq(3, t.m, by= 3),], type = "b",ylim = c(0, 1), col = c("black", "red", "green", "blue"),
        ylab = "Sad", xlab = "Measurement occasion", cex.lab = 1.5, 
        lty = 1, pch =16, las = 1) # sad variance components


# 7.2 Selected CUTS model with om ----

m <- 30
t.m <- 3*m
cuts <- fit.cuts.lavaan.h2

# 7.2.1 Compute factor scores
CT.scores <- lavPredict(cuts)[,m+1] #scores common trait
CS.scores <- lavPredict(cuts)[,1:m] #scores common states
UT.scores <- lavPredict(cuts)[,(m+2):(m+4)] #scores unique traits

# 7.2.2 Extract parameters from the lavaan object
lct <- parameterEstimates(cuts)[ (t.m+1):(2*t.m),"est"] #loadings common traits
lut <- parameterEstimates(cuts)[(2*t.m+1):(3*t.m), "est"] #loadings unique traits
lcs <- parameterEstimates(cuts)[ 1:t.m,"est"] #loadings common states

cuts.ins <- parameterEstimates(cuts)[(3*t.m+1):(4*t.m), "est"] #intercepts 
 
us_var <- parameterEstimates(cuts)[ (4*t.m+((m+4)*(m+3)/2)+1):(5*t.m+((m+4)*(m+3)/2)),"est"] # variance unique states
cs_var <- parameterEstimates(cuts)[ (5*t.m+((m+4)*(m+3)/2)+1):(5*t.m+((m+4)*(m+3)/2)+m),"est"] # variance common states
ut_var <- parameterEstimates(cuts)[ (5*t.m+((m+4)*(m+3)/2)+m+1):(5*t.m+((m+4)*(m+3)/2)+m+3),"est"] # variance unique traits
ct_var <- parameterEstimates(cuts)[ (5*t.m+((m+4)*(m+3)/2)+m+4),"est"]# variance common trait

#7.2.3 Compute variance coefficients
#Total variance
cuts_tot_var <- (lct^2)*ct_var + (lut^2)*ut_var + (lcs^2)*rep(cs_var, each = 3) + us_var
#Reliability
cuts_rel <-  ((lct^2)*ct_var + (lut^2)*ut_var + (lcs^2)*rep(cs_var, each = 3)) / cuts_tot_var
#Total Consistency
cuts_tcon <- ((lct^2)*ct_var + (lut^2)*ut_var) / cuts_tot_var
#Occasion specificity
cuts_spe <- ((lcs^2)*rep(cs_var, each = 3)) / cuts_tot_var
#Common consistency
cuts_ccon <- ((lct^2)*ct_var) / cuts_tot_var
#Unique consistency
cuts_ucon <- ((lut^2)*ut_var) / cuts_tot_var

# 7.2.4 Compute predictes scores based on the factor scores
cuts_pred_worry <- matrix(lct[seq(1, 90, by = 3)],129,m, byrow = TRUE)*CT.scores +
  matrix(lut[seq(1, 90, by = 3)],129,m, byrow = TRUE)*UT.scores[,1]+
  matrix(lcs[seq(1, 90, by = 3)],129,m, byrow = TRUE)*CS.scores +
  cuts.ins[seq(1, 90, by = 3)]

cuts_pred_fear <- matrix(lct[seq(2, 90, by = 3)],129,m, byrow = TRUE)*CT.scores +
  matrix(lut[seq(2, 90, by = 3)],129,m, byrow = TRUE)*UT.scores[,2]+
  matrix(lcs[seq(2, 90, by = 3)],129,m, byrow = TRUE)*CS.scores +
  cuts.ins[seq(2, 90, by = 3)]

cuts_pred_sad <- matrix(lct[seq(3, 90, by = 3)],129,m, byrow = TRUE)*CT.scores +
  matrix(lut[seq(3, 90, by = 3)],129,m, byrow = TRUE)*UT.scores[,3]+
  matrix(lcs[seq(3, 90, by = 3)],129,m, byrow = TRUE)*CS.scores +
  cuts.ins[seq(3, 90, by = 3)]

# 7.2.5 Plots
i <- 53 # select person to plot
# plot worry observed vs latent state variables person i
par(mfrow= c(3,1),mar= c(0,5,0,2), oma= c(4,0,1,8), xpd = NA)

plot(cuts_pred_worry[i,], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xaxt="n",xlab="", cex.lab = 1.5, ylab = "Worry", las = 1)# traceplot of the latent state scores of one person
points(1:m, PlosOne_W[i, seq(1,t.m, by = 3)], type = "b", col = "red", pch = 16)
abline(h =CT.scores[i], xpd = FALSE) #plot latent trait factor score
#abline(h = mean(as.numeric(PlosOne_W[i, seq(1,t.m, by = 3)]), na.rm = T), 
#      col="green", xpd = TRUE) # plot mean of worry of person i
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange", 
#      xpd = TRUE) # plot combined mean of all observed variables of person i

# plot fear observed vs latent state variables person i
plot(cuts_pred_fear[i,], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xaxt="n",xlab="", cex.lab = 1.5,ylab = "Fear", las = 1)# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(2,t.m, by = 3)], type = "b", col = "red", pch = 16)
abline(h =CT.scores[i], xpd = FALSE)
#abline(h = mean(as.numeric(PlosOne_W[i, seq(2,t.m, by = 3)]), na.rm = T), col="green")
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")
legend(31.5, 5, legend = c("Observed", "Predicted"), col = c( "red", "blue"),
       lty = 1, pch = 16, lwd = 2, cex = 1)


# plot sad observed vs latent state variables person i
plot(cuts_pred_sad[i,], 
     type = "b", pch = 16, ylim = c(0,7.5), col = "blue",
     xlab="Measurement occasion", cex.lab = 1.5, ylab = "Sad", las = 1)# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(3,t.m, by = 3)], type = "b", pch= 16, col = "red")
abline(h =CT.scores[i], xpd = FALSE)
#abline(h = mean(as.numeric(PlosOne_W[i, seq(3,t.m, by = 3)]), na.rm = T), col="green")
#abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")


# plot variance components by variable.
var.comp <- cbind(cuts_rel, cuts_ccon, cuts_ucon, cuts_spe)

par(mfrow= c(3,1),mar= c(0,5,0,2), oma= c(4,0,1,6), xpd = NA)
matplot(var.comp[seq(1, t.m, by= 3),], type = "b", ylim = c(0, 1), col = c("black", "red", "green", "blue"),
        ylab = "Worry", xaxt="n",xlab="", cex.lab = 1.5, lty = 1, pch =16, las = 1) # worry variance components

matplot(var.comp[seq(2, t.m, by= 3),], type = "b", ylim = c(0, 1), col = c("black", "red", "green", "blue"),
        ylab = "Fear",xaxt="n",xlab="", cex.lab = 1.5, lty = 1, pch =16, las = 1) # fear variance components
legend(32, 0.75, legend = c("Rel", "CCon", "UCon", "Spe"), col = c("black", "red", "green", "blue"),
       lty = 1, lwd = 2, pch = 16, cex = 1)

matplot(var.comp[seq(3, t.m, by= 3),], type = "b",ylim = c(0, 1), col = c("black", "red", "green", "blue"),
        ylab = "Sad", xlab = "Measurement occasion", cex.lab = 1.5, 
        lty = 1, pch =16, las = 1) # sad variance components

  


