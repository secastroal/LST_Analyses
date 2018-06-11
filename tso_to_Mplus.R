
# Function to write Mplus syntax for the extended trait-state-occasion model given a data.frame ----

# This function writes a syntax of the extended TSO model based on the figure 3B
# of Eid et al. (2017)

# Restrictions homogeneity of autoregressive effects, homogeneity of trait loadings (equivalence assumption),
# homogeneity of state residual variances, homogenity of error variances.

# dt is the data.frame. It is supposed to be a longitudinal data in wide format.
#     Only the variables that are measured repeatedly are included in the data set, 
#     ID or sex variables should not be included.
# nocc is the number of latent occasion variables or the number of measurement 
#     occasions.
# trait.indicators is a logical value to define whether the model will have
#     a trait indicator for each variable (TRUE) or there will be a trait variable
#     associated to all the observed variables (FALSE). The lastest option is not 
#     presented by Eid et al. (2017), hence, it is uncertain whether it would be right.
# equiv.assumption define if the model assumes tau-equivalence, equivalence or 
#     tau-congenericity in relation to the latent occasions (occ) and the latent traits (theta).
#     Use "equi" or "cong". The essential tau-equivalence assumption (as in LST)
#     was ommited because given the way we are modelling the TSO model, it is unclear 
#     whether the intercepts are associated to the trait indicators or the latent 
#     occasions. Also, in this case, these assumptions are just modifying whether
#     the loadings are fixed to 1 or not, they do not affect the intercepts as they
#      do in LST theory.
# scale.invariance define whether the measurement model is invariant over time or not.
#     It can be define for the intercepts and for the loadings of the latent occasions.
#     it does not apply for the loadings of the traits.
# homocedasticity.assumption define whether the variances of the latent occasions (occ.red)
#     and the measurement errors are invariant over time.
# autoregressive.homogeneity is a logical value that defines whether the auto-
#     regressive effects are constant over time or not.

write.tso.to.Mplus <- function(dt, nocc, trait.indicators = TRUE,
                                equiv.assumption = list(occ = "cong", theta = "equi"),
                                scale.invariance = list(int = FALSE, lambda = FALSE),
                                homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                                autoregressive.homogeneity = FALSE){


nobs <- dim(dt)[2]

#labels
ind <- paste0(rep(1:(nobs/nocc) ,nocc),
             rep(1:nocc, each=(nobs/nocc)))
#lambdas for occasions and traits
las <- paste0("las", ind)
lat <- paste0("lat", ind)

# intercepts for observed variables
ins <- paste0("ins", ind)


#latent occasion variables' names
occ <- paste0("occ", 1:nocc)

# trait and trait indicators names
if(trait.indicators){
  theta <- paste0("theta", 1:(nobs/nocc))}else{
    theta <- "theta"
}


# observed variables names
obs <- names(dt)

# variances and residuals' labels

o_red <- paste0("err", ind)
occ_var <- paste0("oc_var", 1:nocc)

if(trait.indicators){
  t_var <- paste0("t_var", 1:(nobs/nocc))}else{
    t_var <- "t_var"
  }

# autoregressive efects labels

beta <- paste0("beta", 1:(nocc-1))

# computed parameters labels
# In the TSO the consistency is divided in predictability by trait and 
#     unpredictability by trait. These components can be computed with Mplus
#     or afterwards. If, it is needed I can make the syntax to compute them 
#     with Mplus, for now it does not seem necessary

#rel <- paste0("rel", ind)
#pre <- paste0("pre", ind)
#upre <- paste0("upre", ind)
#spe <- paste0("spe", ind)


# if-else statements to define constraints (parameters equalities)
fixed_s <- seq(1, nobs, by = nobs/nocc)

if(equiv.assumption$occ == "equi"){
  las <- rep(1, nobs)
  
}

if(equiv.assumption$occ == "cong"){
  las[fixed_s] <- 1
  if(scale.invariance$lambda){
    las <- rep(las[1:(nobs/nocc)], nocc)
  }
}


if(equiv.assumption$theta == "equi"){
  lat <- rep(1, nobs)
}


if(equiv.assumption$theta == "cong"){
  lat[1] <- "1"
  if(trait.indicators){
    lat[1:(nobs/nocc)] <- "1"
  }
}

# define invariance of intercepts

if(scale.invariance$int){
  if(trait.indicators){
    ins <- rep(ins[1:(nobs/nocc)], nocc)
  }else{
    ins <- rep(ins[1], nobs)
  }
  
}

# homogeneity of autoregresive effects

if(autoregressive.homogeneity){
  beta <- rep(beta[1], nocc-1)
}

# homogeneity of errors and occasion variances

if(homocedasticity.assumption$error){
  o_red <- rep(o_red[1:(nobs/nocc)], nocc)
}

if(homocedasticity.assumption$occ.red){
  occ_var <- rep(occ_var[1], nocc)
}


# write equations for the latent occasions
occ_syntax <- rep(NA, nocc)
for(i in 1:nocc){
  ms <- ((i-1)*(nobs/nocc) + 1):(i*(nobs/nocc))
  paths <- rep(NA, length(ms) )
  ix <- which(las[1:length(ms)]=="1")
  
    paths[ix] <- paste0(obs[ms[ix]],"@",las[ms[ix]])
    paths[-ix] <- paste0(obs[ms[-ix]]," (",las[ms[-ix]],")")
  
  paths <- paste(paths, collapse = "\n")
  occ_syntax[i] <- paste(occ[i], "BY", paths, ";", sep = " ")
}
rm(i, ms, paths, ix)


# write equations for the latent trait or trait indicators

trait_syntax <- NA
paths <- rep(NA, nobs )
ix <- which(lat=="1")

paths[ix] <- paste0(obs[ix],"@",lat[ix])
paths[-ix] <- paste0(obs[-ix]," (",lat[-ix],")")

paths <- paste(paths, collapse = "\n")
trait_syntax <- paste(theta, "BY", paths, ";", sep = " ")
rm(paths, ix)
if(trait.indicators){
  trait_syntax <- rep(NA, nobs/nocc)
  
  for(i in 1:(nobs/nocc)){
    paths <- rep(NA, nocc )
    ix <- which(lat[seq(i, nobs, by = nobs/nocc)]=="1")
    
    paths[ix] <- paste0(obs[seq(i, nobs, by = nobs/nocc)[ix]],"@",lat[seq(i, nobs, by = nobs/nocc)[ix]])
    paths[-ix] <- paste0(obs[seq(i, nobs, by = nobs/nocc)[-ix]]," (",lat[seq(i, nobs, by = nobs/nocc)[-ix]],")")
    
    paths <- paste(paths, collapse = "\n")
    trait_syntax[i] <- paste(theta[i], "BY", paths, ";", sep = " ")
  }
  rm(i, paths, ix)
  
}



# write intercetps of observed variables
intercepts_syntax <- paste0("[", obs, "*]", " (",ins,")", " ;", collapse = "\n")



# write means of trait indicators and occasions

occ_mean_syntax <- paste0("[", occ, "@0] ;", collapse = "\n" )

trait_mean_syntax <- paste0("[", theta, "] ;", collapse = "\n")

# write 0 correlations among occasions  and traits
latent.v <- c(occ, theta)

cor_syntax <- rep(NA, nocc)
for(i in 1:nocc){
  
  cor_syntax[i] <- paste0(latent.v[i], " WITH ",
                          paste0(latent.v[-(1:i)], "@0", collapse = "\n"), " ;")
}
rm(i)


# write autoregressive effects

beta_syntax <- paste0(occ[2:nocc], " ON ", occ[1:(nocc-1)], " (", beta, ") ;", 
                      collapse = "\n")


# write variance components

var_syntax <- paste(paste0(obs, " (", o_red, ") ;", collapse = "\n"),
                     paste0(occ, " (", occ_var, ") ;", collapse = "\n"),
                     paste0(theta, " (", t_var, ") ;", collapse =  "\n"),
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
                        "\n! Latent occasion variables",
                        paste0(occ_syntax, collapse = "\n"),
                        "\n! Latent trait variable or trait indicators",
                        paste0(trait_syntax, collapse = "\n"),
                        "\n! Intercepts and means",
                        paste0(intercepts_syntax, collapse = "\n"),
                        paste0(occ_mean_syntax, collapse = "\n"),
                        paste0(trait_mean_syntax, collapse = "\n"),
                        "\n! Set 0 correlations",
                        paste0(cor_syntax, collapse = "\n"),
                        "\n! autoregressive effects",
                        paste0(beta_syntax, collapse = "\n"),
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