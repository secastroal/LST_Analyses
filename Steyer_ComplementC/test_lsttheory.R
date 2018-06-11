library(lsttheory)


########### dataset: CTT model software appendix ###########

head(taucong)

m1 <- lsttheory(1,0,data=taucong,
                equiv.assumption=list(tau="cong", theta="cong"))
m1
summary(m1@lavaanres)
print(m1)
is(m1@lavaanres)
parameterEstimates(m1@lavaanres) # estimates + se + CI bounds
parameterTable(m1@lavaanres)
standardizedSolution(m1@lavaanres)
m1@lavaansyntax
cor(lavPredict(m1@lavaanres),eta)

modelCTT_supC <- "
eta =~ y1 + la21 *y2 + la31 *y3
eta ~ NA*1
y1 ~ 0*1

eta ~~ vareta * eta
y1 ~~ vareps1 *y1
y2 ~~ vareps2 *y2
y3 ~~ vareps3 *y3

rely1 := vareta / ( vareta+vareps1 )
rely2 := la21^2* vareta / ( la21^2* vareta+vareps2 )
rely3 := la31^2* vareta / ( la31^2* vareta+vareps3 )
"
m1_supC <- sem(model = modelCTT_supC, data = taucong)
parameterEstimates(m1_supC)
parameterTable(m1_supC)
semPaths(m1_supC, style="lisrel", intercepts=T, layout="tree2", rotation=4, nCharNodes=4, nCharEdges=4, optimizeLatRes=F, residScale=10)
m1@lavaansyntax

############## dataset: multistate model software appendix ###############

head(multistate)

m1 <- lsttheory(2,1,data=multistate,
                equiv.assumption=list(tau="cong", theta="ess"))
 m1
 cat(m1@lavaansyntax)
 summary(m1@lavaanres)
 parTable(m1@lavaanres)



######## dataset: multistate model (older version with 3 time points) ########
head(multistate02)

m1 <- lsttheory(3,0,data=multistate02,
                equiv.assumption=list(tau="ess", theta="cong"))
summary(m1@lavaanres)


################ dataset: multitrait-multistate model ###################

head(multitraitmultistate)
m2 <- lsttheory(4,2,multitraitmultistate)

m2
cat(m2@lavaansyntax)
summary(m2@lavaanres)
parTable(m2@lavaanres)
m1@lstmodel


m3 <- lsttheory(4,0,multitraitmultistate)
m3
cat(m3@lavaansyntax)
summary(m3@lavaanres)
parTable(m3@lavaanres)


m4 <- lsttheory(1,0,multitraitmultistate, 
                equiv.assumption=list(tau="equi", theta="cong"), 
                scale.invariance=list(lait0=FALSE, lait1=TRUE, lat0=TRUE, lat1=TRUE))

m4
cat(m4@lavaansyntax)
summary(m4@lavaanres)
parTable(m4@lavaanres)


m5 <- lsttheory(4,2,multitraitmultistate, 
                equiv.assumption=list(tau="cong", theta="cong"), 
                scale.invariance=list(lait0=FALSE, lait1=FALSE, 
                                      lat0=FALSE, lat1=FALSE))
m5
cat(m5@lavaansyntax)
summary(m5@lavaanres)
parTable(m5@lavaanres)



m1 <- lsttheory(4, 2, multitraitmultistate, 
                equiv.assumption=list(tau="cong", theta="cong"), 
                scale.invariance=list(lait0=TRUE,lait1=TRUE,lat0=TRUE,lat1=TRUE))
cat(m1@lavaansyntax)



################# semPlot ##########################

# 
library(semPlot)
semPaths(m1@lavaanres, style="lisrel", intercepts=F, layout="tree2", rotation=4, nCharNodes=4, nCharEdges=4, optimizeLatRes=F, residScale=10)
 
semPaths(m5@lavaanres, style="lisrel", intercepts=F, layout="tree2", rotation=4, nCharNodes=4, nCharEdges=4, optimizeLatRes=F, residScale=10)
 