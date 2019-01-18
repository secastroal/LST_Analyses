# Contents
# 0.0 Prepare environment


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
library(invgamma)
library(MASS)
source("R/write.Mplus.options.R")
source("R/write.cuts.to.Mplus.R")
source("R/write.mlcuts.to.Mplus.R")
source("R/write.msst.to.Mplus.R")
source("R/write.mlmsst.to.Mplus.R")
source("R/write.tso.to.Mplus.R")
source("R/write.mltso.to.Mplus.R")
