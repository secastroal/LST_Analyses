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
#load( file = paste0("lavaan_files/msst_h9_m",m, ".R"))



# print lsttheory object returns the variance components: reliability, consistency, and 
# occasion specifity   
summary(msst.h3@lavaanres, fit.measures = TRUE)

anova(msst.h1@lavaanres,msst.h2@lavaanres)
anova(msst.h1@lavaanres,msst.h3@lavaanres)
anova(msst.h2@lavaanres,msst.h3@lavaanres)


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
(sum(is.na(PlosOne_Wr)))/(dim(PlosOne_Wr)[1]*dim(PlosOne_Wr)[2])

i <- 123 # select person to plot
# plot worry observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m], 
     type = "l", ylim = c(0,7.5), col = "blue")# traceplot of the latent state scores of one person
points(1:m, PlosOne_W[i, seq(1,t.m, by = 3)], type = "b", col = "red")
abline(h =lavPredict(msst@lavaanres)[i,(m+1)]) #plot latent trait factor score
abline(h = mean(as.numeric(PlosOne_W[i, seq(1,t.m, by = 3)]), na.rm = T), col="green") # plot mean of worry of person i
abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange") # plot combined mean of all observed variables of person i

# plot fear observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m] * parameterEstimates(msst@lavaanres)[seq(2, t.m, by = 3), 5], 
     type = "l", ylim = c(0,7.5), col = "blue")# traceplot of the latent state scores (multiplied by loadings) of one person
points(1:m, PlosOne_W[i, seq(2,t.m, by = 3)], type = "b", col = "red")
abline(h =lavPredict(msst@lavaanres)[i,(m+1)])
abline(h = mean(as.numeric(PlosOne_W[i, seq(2,t.m, by = 3)]), na.rm = T), col="green")
abline(h = mean(as.numeric(PlosOne_W[i, 1:t.m]), na.rm = T), col="orange")


# plot sad observed vs latent state variables person i
plot(lavPredict(msst@lavaanres)[i,1:m] * parameterEstimates(msst@lavaanres)[seq(3, t.m, by = 3), 5], 
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

cuts.file.name <- paste0("cuts_h1_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
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

cuts.file.name <- paste0("cuts_h2_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
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

cuts.file.name <- paste0("cuts_h3_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
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

cuts.file.name <- paste0("cuts_h4_m", m)

prepareMplusData(PlosOne_Wr,paste0("Mplus_files/",cuts.file.name,".dat"), inpfile = T)

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
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

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
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

save(fit.tso.lavaan.h1, file = paste0("lavaan_files/tso3b_h1_m",m, ".R"))


# h2: congenericity in the occasions and equivalence in the trait

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
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

save(fit.tso.lavaan.h2, file = paste0("lavaan_files/tso3b_h2_m",m, ".R"))

# h3: equivalence in the occasions and congenericity in the trait

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
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

save(fit.tso.lavaan.h3, file = paste0("lavaan_files/tso3b_h3_m",m, ".R"))

# h4: equivalence in the occasions and equivalence in the trait

tso_syntax <- write.tso.to.Mplus(PlosOne_Wr, nocc = m, figure = "3b",
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

save(fit.tso.lavaan.h4, file = paste0("lavaan_files/tso3b_h4_m",m, ".R"))


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
                                   method.trait = "om",
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

save(fit.cuts.lavaan.h1, file = paste0("lavaan_files/cuts_h1_m",m, ".R"))

# h2: Cuts with all method loadings fixed to 1

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
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

save(fit.cuts.lavaan.h2, file = paste0("lavaan_files/cuts_h2_m",m, ".R"))

# h3: Cuts + weak factorial invariance across time for the CT and the CS

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
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

save(fit.cuts.lavaan.h3, file = paste0("lavaan_files/cuts_h3_m",m, ".R"))

# h4: Cuts + weak factorial invariance across traits and states

cuts_syntax <- write.cuts.to.Mplus(PlosOne_Wr, nstate = m, 
                                   method.trait = "om",
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

save(fit.cuts.lavaan.h4, file = paste0("lavaan_files/cuts_h4_m",m, ".R"))


load( file = paste0("lavaan_files/cuts_h2_m",m, ".R"))



anova(fit.cuts.lavaan.h1, fit.cuts.lavaan.h2)




fit.cuts.lavaan.h1@test



summary(fit.cuts.lavaan.h4, fit.measures = T, estimates = F)

anova(fit.cuts.lavaan.h1, fit.cuts.lavaan.h2)

parameterEstimates(fit.cuts.lavaan.h4)[1:100,]

