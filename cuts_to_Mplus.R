
# Function to write Mplus syntax for the common unique trait state model model given a data.frame ----

# Several possible constraints are mentioned in page 51 (Hamaker, Schuurman, & Zijlmans, 2017)
# weak factorial invariance across trait and state, 
# weak factorial invariance over time for common state and/or common traits, 
# homocedasticity of variances over time
# Furthermore unique traits/states and common traits/states are assumed 
#     to have mean 0. Although, this constraint is not explicit in the Mplus syntax
#     available in the suplementary material.

# dt is the data.frame. It is supposed to be a longitudinal data in wide format.
#     Only the variables that are measured repeatedly are included in the data set, 
#     ID or sex variables should not be included.
# nstate is the number of common state variables. In other words, it is the number
#     of measurement occasions.
# method.trait defines whether there is a unique trait to each variable (orthogonal methods "om" Default)
#     or whether there is the number of variables minus one unique traits ("m-1").
#     This restriction might be neccesary when "om" don't converge (see suplementary
#     material).
# scale.invariance define whether the measurement model is invariant over time or not.
#     It can be define for intercepts (int) and loadings. 
#     If scale.invariance$lambda is TRUE, scale invariance of the 
#     loadings is assumed for both common states and common traits.
# state.trait.invariance define whether the trait structure is equal to the state 
#     structure or not. This means that the loading of the common trait on the variable 
#     Yit is equal to the loading of the common state on the variable Yit.
# homocedasticity.assumption define whether the variance components are invariant
#     over time. The variances that can be constraint are the variances of the unique
#     states (error), the variances of the common states (cs.red), and the variances
#     of the unique traits (ut.red).

write.cuts.to.Mplus <- function(dt, nstate,
                                method.trait = "om",
                                scale.invariance = list(int = FALSE, lambda = FALSE),
                                state.trait.invariance = FALSE,
                                fixed.method.loadings = FALSE,
                                homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) ){

nobs <- dim(dt)[2]

#labels
ind <- paste0(rep(1:(nobs/nstate) ,nstate),
             rep(1:nstate, each=(nobs/nstate)))


#lambdas for states and traits
lcs <- paste0("lcs", ind) #lambda common states
lct <- paste0("lct", ind) #lambda common trait
lut <- paste0("lut", ind) #lambda unique traits
# intercepts for observed variables
ins <- paste0("ins", ind)


# Common states names
CS <- paste0("CS", 1:nstate) # Common state

# Unique and Common traits names

CT <- "CT"

UT <- paste0("UT", 1:(nobs/nstate))

# observed variables names
obs <- names(dt)

# variances and residuals' labels

o_red <- paste0("err", ind)
cs_var <- paste0("cs_var", 1:nstate)
ut_var <- paste0("ut_var", 1:(nobs/nstate))
ct_var <- "ct_var"


# computed parameters labels
# In the cuts model there are not variances components defined such as
# reliability or consistency. However, based on the definition of these variance
# components in the LST, similar variances components could be defined.

#rel <- paste0("rel", ind)
#pre <- paste0("pre", ind)
#upre <- paste0("upre", ind)
#spe <- paste0("spe", ind)


# if-else statements to define constraints (parameters equalities)
fixed_s <- seq(1, nobs, by = nobs/nstate)

lcs[fixed_s] <- 1
lct[fixed_s] <- 1
lut[1:(nobs/nstate)]<-1

m.obs <- obs

if(method.trait == "m-1"){
  UT <- UT[-1]
  ut_var <- ut_var[-1]
  m.obs <- obs[-(seq(1, nobs, by = nobs/nstate))]
  lut <- lut[-(seq(1, nobs, by = nobs/nstate))]
}

# define invariance of intercepts

if(scale.invariance$int){
  ins <- rep(ins[1:(nobs/nstate)], nstate)
}

# define invariance over time of lambdas for states and traits

if(scale.invariance$lambda){
  lcs <- rep(lcs[1:(nobs/nstate)], nstate)
  lct <- rep(lct[1:(nobs/nstate)], nstate)
}

# define state trait invariance

if(state.trait.invariance){
  lcs <- lct
}

# fixed mehtod loadings

if(fixed.method.loadings){
  lut <- rep(1, length(m.obs))
}

# homogeneity of CS, US and UT variances

if(homocedasticity.assumption$error){
  o_red <- rep(o_red[1:(nobs/nstate)], nstate)
}

if(homocedasticity.assumption$cs.red){
  cs_var <- rep(cs_var[1], nstate)
}

if(homocedasticity.assumption$ut.red){
  ut_var <- rep(ut_var[1], nobs/nstate)
}


# write equations for the common states
cs_syntax <- rep(NA, nstate)
for(i in 1:nstate){
  ms <- ((i-1)*(nobs/nstate) + 1):(i*(nobs/nstate))
  paths <- rep(NA, length(ms) )
  ix <- which(lcs[1:length(ms)]=="1")
  
    paths[ix] <- paste0(obs[ms[ix]],"@",lcs[ms[ix]])
    paths[-ix] <- paste0(obs[ms[-ix]]," (",lcs[ms[-ix]],")")
  
  paths <- paste(paths, collapse = "\n")
  cs_syntax[i] <- paste(CS[i], "BY", paths, ";", sep = " ")
}
rm(i, ms, paths, ix)


# write equations for the common and unique trait

ct_syntax <- NA
paths <- rep(NA, nobs )
ix <- which(lct=="1")

paths[ix] <- paste0(obs[ix],"@",lct[ix])
paths[-ix] <- paste0(obs[-ix]," (",lct[-ix],")")

paths <- paste(paths, collapse = "\n")
ct_syntax <- paste(CT, "BY", paths, ";", sep = " ")
rm(paths, ix)

ut_syntax <- rep(NA, length(UT))
  
  for(i in 1:length(UT)){
    paths <- rep(NA, nstate )
    ms <- seq(i, length(m.obs), by = length(m.obs)/nstate)
    
    ix <- which(lut[ms]=="1")
    
    paths[ix] <- paste0(m.obs[ms[ix]],"@",lut[ms[ix]])
    paths[-ix] <- paste0(m.obs[ms[-ix]]," (",lut[ms[-ix]],")")
    
    paths <- paste(paths, collapse = "\n")
    ut_syntax[i] <- paste(UT[i], "BY", paths, ";", sep = " ")
  }
  rm(i, paths)
  




# write intercetps of observed variables
intercepts_syntax <- paste0("[", obs, "*]", " (",ins,")", " ;", collapse = "\n")



# fix means of latent variables to 0. This is optional, it is not explicitly
# in the example syntaxes of the suplementary material.

cs_means_syntax <- ""
ut_means_syntax <- ""
ct_means_syntax <- ""

if(fixed.means$cs){cs_means_syntax <- paste0("[", CS, "@0] ;", collapse = "\n")}

if(fixed.means$ut){ut_means_syntax <- paste0("[", UT, "@0] ;", collapse = "\n")}

if(fixed.means$ct){ct_means_syntax <- paste0("[", CT, "@0] ;", collapse = "\n")}


# write 0 correlations among occasions  and traits
latent.v <- c(CS, UT, CT)

cor_syntax <- rep(NA, length(latent.v)-1)
for(i in 1:(length(latent.v)-1)){
  
  cor_syntax[i] <- paste0(latent.v[i], " WITH ",
                          paste0(latent.v[-(1:i)], "@0", collapse = "\n"), " ;")
}
rm(i)


# write variance components

var_syntax <- paste(paste0(obs, " (", o_red, ") ;", collapse = "\n"),
                     paste0(CS, " (", cs_var, ") ;", collapse = "\n"),
                     paste0(UT, " (", ut_var, ") ;", collapse = "\n"),
                     paste0(CT, " (", ct_var, ") ;", collapse =  "\n"),
                     sep = "\n")

# define computed parameters
#def_syntax <- paste(paste0("NEW(", rel, ") ;", collapse = "\n" ),
#                    paste0("NEW(", con, ") ;", collapse = "\n" ),
#                    paste0("NEW(", spe, ") ;", collapse = "\n" ),
#                    sep = "\n")

#lat_r <- rep(lat, each = nobs/neta)

#total_var <- paste0(lat_r, "^2*",las, "^2*", 
#                     rep(t_var, each = neta/ntheta), "+",
#                     las, "^2*", rep(s_red, each = nobs/neta), 
#                    "+", o_red)

#trait_var <- paste0(lat_r, "^2*",las, "^2*", 
#                    rep(t_var, each = neta/ntheta))

#state_var <- paste0(las, "^2*", rep(s_red, each = nobs/neta))


#rel_syntax <- paste0(rel, "= (", trait_var, "+", state_var,
#                     ")/", "\n(", total_var, ") ;", collapse = "\n")

#spe_syntax <- paste0(spe, "= (", state_var,
#                     ")/", "\n(", total_var, ") ;", collapse = "\n")

#con_syntax <- paste0(con, "= (", trait_var,
#                     ")/", "\n(", total_var, ") ;", collapse = "\n")



complete_syntax <- paste("\nMODEL:",
                        "\n! common states",
                        paste0(cs_syntax, collapse = "\n"),
                        "\n! common trait",
                        paste0(ct_syntax, collapse = "\n"),
                        "\n! unique trait",
                        paste0(ut_syntax, collapse = "\n"),
                        "\n! Intercepts and means",
                        paste0(intercepts_syntax, collapse = "\n"),
                        paste0(cs_means_syntax, collapse = "\n"),
                        paste0(ut_means_syntax, collapse = "\n"),
                        paste0(ct_means_syntax, collapse = "\n"),
                        "\n! Set 0 covariances",
                        paste0(cor_syntax, collapse = "\n"),
                        "\n! Variance components",
                        var_syntax,
                        #"\nMODEL CONSTRAINT:",
                        #def_syntax,
                        #"\n! Reliabily",
                        #rel_syntax,
                        #"\n! Occasion specificity",
                        #spe_syntax,
                        #"\n! Consistency",
                        #con_syntax,
                        sep = "\n")

return(complete_syntax)

}


