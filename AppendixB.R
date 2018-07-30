## Load packages and files ----
# First, the packages MplusAutomation and RCurl have to be loaded. MplusAutomation
# allows to use Mplus within R and RCurl allows to source online files. The latest
# is used to source the functions that write the Mplus syntax for the different 
# models.

# Either load the packages or install and load the packages. 
if(!require(MplusAutomation)) {install.packages("MplusAutomation"); require(installr)} 
if(!require(RCurl)) {install.packages("RCurl"); require(installr)} 

# Source the three functions to write the Mplus syntax of the msst model, the 
# tso model or the cuts model.

eval(parse(text = getURL("https://raw.githubusercontent.com/secastroal/LST_Analyses/master/msst_to_Mplus.R", 
                         ssl.verifypeer = FALSE)))

eval(parse(text = getURL("https://raw.githubusercontent.com/secastroal/LST_Analyses/master/tso_to_Mplus.R", 
                         ssl.verifypeer = FALSE)))

eval(parse(text = getURL("https://raw.githubusercontent.com/secastroal/LST_Analyses/master/cuts_to_Mplus.R", 
                         ssl.verifypeer = FALSE)))

## Data ----
# To perform the analyses, the data has to be in wide format (each repeated 
# measure of each variable has to be in one column). If the data is loaded into
# R in long format use the function "reshape" to turn it into wide format. 
# Change the names of the variables to be of 5 characters or less.
# Delete all the variables (e.g., ID), that are not needed for 
# the analysis.

# Example of how to read the data.
data <- read.table("data.txt", header = TRUE)

# How many measurement occasions are there?
n.measures <- 10

## Multistate-singletrait analysis in Mplus ----

# Write data and input file for Mplus.
prepareMplusData(data, filename = "msst.dat", inpfile = T)

# Create msst Mplus syntax.
# neta defines the number of measurement occasions.
# equiv.assumption defines the equivalence assumptions made on the states (tau) 
# and the trait (theta). It can be "equi" for equivalence, "ess" for essential
# equivalence, or "cong" for  congenericity.
# More arguments and its options can be found in 
# "https://raw.githubusercontent.com/secastroal/LST_Analyses/master/msst_to_Mplus.R"
msst_syntax <- write.msst.to.Mplus(data, neta = n.measures, 
                                   equiv.assumption = list(tau = "cong", theta = "cong"))

# Paste msst Mplus syntax in the Mplus input file.
write(msst_syntax, file = "msst.inp", append = T)

# run msst in Mplus
runModels("msst.inp")

# read msst Mplus output
msst <- readModels("msst.out")

# get parameter estimates
msst$parameters$unstandardized

## TSO analysis in Mplus ----

# Write data and input file for Mplus.
prepareMplusData(data, filename = "tso.dat", inpfile = T)

# Create tso Mplus syntax.
# nocc defines the number of measurement occasions.
# figure can be either "3a" or "3b". If figure is "3a", the latent state
# variables are modeled. If figure is "3b", the latent state variables 
# are not modeled. 
# More arguments and its options can be found in 
# "https://raw.githubusercontent.com/secastroal/LST_Analyses/master/tso_to_Mplus.R"
tso_syntax <- write.tso.to.Mplus(data, nocc = n.measures, figure = "3b")

# Paste tso Mplus syntax in the Mplus input file.
write(tso_syntax, file = "tso.inp", append = T)

# run tso in Mplus
runModels("tso.inp")

# read msst Mplus output
tso <- readModels("tso.out")

# get parameter estimates
tso$parameters$unstandardized

## CUTS analysis in Mplus ----

# Write data and input file for Mplus.
prepareMplusData(data, filename = "cuts.dat", inpfile = T)

# Create cuts Mplus syntax.
# nstates defines the number of measurement occasions.
# More arguments and its options can be found in 
# "https://raw.githubusercontent.com/secastroal/LST_Analyses/master/cuts_to_Mplus.R"
cuts_syntax <- write.cuts.to.Mplus(data, nstate = n.measures)

# Paste cuts Mplus syntax in the Mplus input file.
write(cuts_syntax, "cuts.inp", append = T)

# run cuts in Mplus
runModels("cuts.inp")

# read cuts Mplus output
cuts <- readModels("cuts.out")

# get parameter estimates
cuts$parameters$unstandardized
