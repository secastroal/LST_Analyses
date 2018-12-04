# Contents
# 0.0 Prepare environment
# 1.0 Read data
# 1.1 Selection of measurements and variables, turn data into wide format
# 2.0 Select number of measurements and delete full NAs cases
# 3.0 CUTS vs MLCUTS
# 4.0 MSST vs MLMSST
# 5.0 TSO vs MLTSO


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
source("mlmsst_to_Mplus.R")
source("mltso_to_Mplus.R")
source("mlcuts_to_Mplus.R")


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

# Missing values in neur are eliminated

for(i in 1:130){
  if(length(unique(PlosOne_C[((60*(i-1))+1):(60*i), "neur"]))==1){
    PlosOne_C[((60*(i-1))+1):(60*i), "neur"] <- unique(PlosOne_C[((60*(i-1))+1):(60*i), "neur"])[1]
  }else{
    PlosOne_C[((60*(i-1))+1):(60*i), "neur"] <- unique(PlosOne_C[((60*(i-1))+1):(60*i), "neur"])[-which(is.na(unique(PlosOne_C[((60*(i-1))+1):(60*i), "neur"])))]
  }
}
rm(i)

# Remove cases due to complete NAs

pos.cases.NA <- which(tapply(rowSums(is.na(PlosOne_C[,6:11])), PlosOne_C$subjn, sum, na.rm = TRUE)==360)
cases.NA <- unique(PlosOne_C$subjn)[pos.cases.NA]
PlosOne_C <- PlosOne_C[-which(PlosOne_C$subjn %in% cases.NA),]

rm(pos.cases.NA, cases.NA)

# Only the variables related to negative affect are selected to try LST models.
PlosOne_W <- reshape(PlosOne_C[, c(1, 8:10, 13)], v.names=c("worry", "fear", "sad"),
                     timevar = "time", idvar="subjn", direction="wide")
# The variable subject number is discarded.
PlosOne_W <- PlosOne_W[, -1]

# Change variable names to delete "." given that this is not allowed in Mplus syntax.
names(PlosOne_W) <- gsub("\\.", "", names(PlosOne_W))

# 2.0 Select number of measurements and delete full NAs cases ----

m <- 10 # number of measurement occasions
t.m <- m*3 # total number of variables
PlosOne_Cr <- PlosOne_C[PlosOne_C$time <= m, c("subjn","time", "worry", "fear", "sad") ]
PlosOne_Wr <- PlosOne_W[, 1:(t.m)]

if(m != 60){
  pos.cases.NA <- which(tapply(rowSums(is.na(PlosOne_Cr[,3:5])), PlosOne_Cr$subjn, sum, na.rm = TRUE)==3*m)
  cases.NA <- unique(PlosOne_Cr$subjn)[pos.cases.NA]
  if(length(pos.cases.NA) != 0){PlosOne_Cr <- PlosOne_Cr[-which(PlosOne_Cr$subjn %in% cases.NA),]}
  if(length(pos.cases.NA) != 0){PlosOne_Wr <- PlosOne_Wr[-pos.cases.NA,]}
  rm(pos.cases.NA, cases.NA)
}


# 3.0 CUTS vs MLCUTS ----

file.name <- paste0("mlcuts_m", m)

prepareMplusData(PlosOne_Cr,paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)


prepareMplusData(PlosOne_C[, c("subjn", "time", "worry", "fear", "sad", "cheer", "relax", "neur")],
                 paste0("ML_Mplus_files/","PlosOneAllVar",".dat"), inpfile = T)


analysis_syntax <- "USEVAR = worry fear sad;
CLUSTER = subjn;

ANALYSIS:
TYPE = TWOLEVEL;
ESTIMATOR = ML;
ITERATIONS = 500000;" # increase H1 iterations

ml_syntax <- write.mlcuts.to.Mplus(PlosOne_Cr[, -1])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))

rm(file.name, analysis_syntax, ml_syntax)

file.name <- paste0("cuts_m", m)

prepareMplusData(PlosOne_Wr,paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

analysis_syntax <- "ANALYSIS:
H1ITERATIONS = 5000;" # increase H1 iterations


model_syntax <- write.cuts.to.Mplus(PlosOne_Wr,  nstate = m,
                            method.trait = "om",
                            scale.invariance = list(int = TRUE, lambda = TRUE),
                            state.trait.invariance = FALSE,
                            fixed.method.loadings = TRUE,
                            homocedasticity.assumption = list(error = TRUE, cs.red = TRUE, ut.red = FALSE))


write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(model_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))

rm(file.name, analysis_syntax, model_syntax)

# 4.0 MSST vs MLMSST ----

file.name <- paste0("mlmsst_m", m,"_time")

prepareMplusData(PlosOne_Cr,paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

analysis_syntax <- "USEVAR = worry fear sad;
CLUSTER = subjn;

ANALYSIS:
TYPE = TWOLEVEL;
ESTIMATOR = ML;
ITERATIONS = 500000;" # increase H1 iterations

ml_syntax <- write.mlmsst.to.Mplus(PlosOne_Cr[, -1])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))

rm(file.name, analysis_syntax, ml_syntax)

file.name <- paste0("msst_m", m, "_time")

prepareMplusData(PlosOne_Wr,paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

analysis_syntax <- "ANALYSIS:
H1ITERATIONS = 5000;" # increase H1 iterations


model_syntax <- write.msst.to.Mplus(PlosOne_Wr, neta = m, ntheta = 1, 
                                    equiv.assumption = list(tau = "cong", theta = "equi"),
                                    scale.invariance = list(lait0 = TRUE, lait1 = TRUE, lat0 = TRUE, lat1 = TRUE),
                                    homocedasticity.assumption = list(error = TRUE, state.red = TRUE))
  

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(model_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))

rm(file.name, analysis_syntax, model_syntax)

# 5.0 TSO vs MLTSO ----

file.name <- paste0("mltso_m", m)

prepareMplusData(PlosOne_Cr,paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

analysis_syntax <- "USEVAR = worry fear sad;
CLUSTER = subjn;

ANALYSIS:
TYPE = TWOLEVEL;
ESTIMATOR = ML;
ITERATIONS = 500000;" # increase H1 iterations

ml_syntax <- write.mltso.to.Mplus(PlosOne_Cr[, -1])

write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(ml_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))

rm(file.name, analysis_syntax, ml_syntax)

file.name <- paste0("tso_m", m)

prepareMplusData(PlosOne_Wr,paste0("ML_Mplus_files/",file.name,".dat"), inpfile = T)

analysis_syntax <- "ANALYSIS:
H1ITERATIONS = 5000;" # increase H1 iterations


model_syntax <- write.tso.to.Mplus(PlosOne_Wr, m, figure = "3b",
                                   equiv.assumption = list(occ = "cong", theta = "equi"),
                                   scale.invariance = list(int = TRUE, lambda = TRUE),
                                   homocedasticity.assumption = list(error = TRUE, occ.red = TRUE),
                                   autoregressive.homogeneity = TRUE)


write(analysis_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T) # Write Analysis specifications
write(model_syntax, paste0("ML_Mplus_files/",file.name,".inp"), append = T)

runModels(paste0(getwd(),"/ML_Mplus_files/",file.name,".inp"))

rm(file.name, analysis_syntax, model_syntax)

# End
