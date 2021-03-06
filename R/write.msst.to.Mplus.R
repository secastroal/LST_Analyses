
# Function to write Mplus syntax for the multistate-singletrait model given a data.frame----


# dt is the data.frame. It is supposed to be a longitudinal data in wide format.
#     Only the variables that are measured repeatedly are included in the data set, 
#     ID or sex variables should not be included.
# neta is the number of latent state variables
# ntheta is the number of latent trait variables
# equiv.assumption are the assumptions tau-equivalence, essential-tau-equivalence, 
#     and tau-congenericity in relation to the latent states (tau) and the latent traits (theta).
#     use "equi", "ess", and "cong" 
# scale.invariance define whether the measurement model is invariant over time or not.
#     It can be define for both intercepts and loadings of the latent states and of the
#     latent traits.

write.msst.to.Mplus <- function (dt, neta, ntheta = 1, 
                                 equiv.assumption = list(tau = "cong", theta = "equi"), 
                                 scale.invariance = list(lait0 = FALSE, lait1 = FALSE, lat0 = FALSE, lat1 = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, state.red = FALSE),
                                 second.order.trait = TRUE)
{
nobs <- dim(dt)[2]  
#labels
ind_s <- paste0(rep(1:(nobs/neta) ,neta),
             rep(1:neta, each=(nobs/neta)))
ind_t <- paste0(1:neta, rep(1,neta))
if(!second.order.trait){
  ind_t <- ind_s
}

#lambdas for states and traits
las <- paste0("las", ind_s)
lat <- paste0("lat", ind_t)

# intercepts for observed and state variables
ins <- paste0("ins", ind_s)
int <- paste0("int", ind_t)
if(!second.order.trait){
  ins <- NULL
}


#latent variables' names
eta <- paste0(rep("eta", neta), 1:neta)

theta <- paste0(rep("theta", ntheta), 1:ntheta)

# observed variables names
obs <- names(dt)

# variances and residuals' labels

o_red <- paste0("o_err", ind_s)
s_red <- paste0("s_var", ind_t)
if(!second.order.trait){
  s_red <- paste0("s_var", 1:neta)
}
t_var <- paste0("t_var", 1:ntheta)

# computed parameters labels

rel <- paste0("rel", ind_s)
con <- paste0("con", ind_s)
spe <- paste0("spe", ind_s)

# if-else statements to define constraints (parameters equalities)
fixed_s <- seq(1, nobs, by = nobs/neta)
fixed_t <- seq(1, neta, by = neta/ntheta)
if(!second.order.trait){
  fixed_t <- seq(1, nobs, by = nobs/ntheta)
}

if(equiv.assumption$tau == "equi"){
  ins <- rep(0, nobs)
  las <- rep(1, nobs)
  if(!second.order.trait){
    ins <- NULL
  }
  
}

if(equiv.assumption$tau == "ess"){
  ins[fixed_s] <- 0
  las <- rep(1, nobs)
  if(scale.invariance$lait0){
    ins <- rep(ins[1:(nobs/neta)], neta)
  }
  if(!second.order.trait){
    stop('When fitting non second order LST models, the equivalence assumption "ess" is not possible')
  }
}

if(equiv.assumption$tau == "cong"){
  ins[fixed_s] <- 0
  las[fixed_s] <- 1
  if(scale.invariance$lait0){
    ins <- rep(ins[1:(nobs/neta)], neta)
  }
  if(scale.invariance$lait1){
    las <- rep(las[1:(nobs/neta)], neta)
  }
  if(!second.order.trait){
    ins <- NULL
  }
}


if(equiv.assumption$theta == "equi"){
  int <- rep(0, neta)
  lat <- rep(1, neta)
  if(!second.order.trait){
    int <- rep(0, nobs)
    lat <- rep(1, nobs)
  }
}

if(equiv.assumption$theta == "ess"){
  int[fixed_t] <- 0
  lat <- rep(1, neta)
  if(scale.invariance$lat0){
    int <- rep(int[1:(neta/ntheta)], ntheta)
    if(!second.order.trait){
      int <- rep(int[1:(nobs/neta)], neta)
    }
  }
}

if(equiv.assumption$theta == "cong"){
  int[fixed_t] <- 0
  lat[fixed_t] <- 1
  if(scale.invariance$lat0){
    int <- rep(int[1:(neta/ntheta)], ntheta)
    if(!second.order.trait){
      int <- rep(int[1:(nobs/neta)], neta)
    }
  }
  if(scale.invariance$lat1){
    lat <- rep(lat[1:(neta/ntheta)], ntheta)
      if(!second.order.trait){
        lat <- rep(lat[1:(nobs/neta)], neta)
      }
  }
}

if(homocedasticity.assumption$error){
  o_red <- rep(o_red[1:(nobs/neta)], neta)
}


if(homocedasticity.assumption$state.red){
  s_red <- rep(s_red[1:(neta/ntheta)], ntheta)
  if(!second.order.trait){
    s_red <- rep(s_red[1], neta)
  }
}

# write equations for the latent states
eta_syntax <- rep(NA, neta)
for(i in 1:neta){
  ms <- ((i-1)*(nobs/neta) + 1):(i*(nobs/neta))
  paths <- rep(NA, length(ms) )
  ix <- which(las[1:length(ms)]=="1")
  
    paths[ix] <- paste0(obs[ms[ix]],"@",las[ms[ix]])
    paths[-ix] <- paste0(obs[ms[-ix]]," (",las[ms[-ix]],")")
  
  paths <- paste(paths, collapse = "\n")
  eta_syntax[i] <- paste(eta[i], "BY", paths, ";", sep = " ")
}
rm(i, ms, paths, ix)


# write equations for the latent trait

theta_syntax <- rep(NA, ntheta)
for(i in 1:ntheta){
  ms <- ((i-1)*(neta/ntheta) + 1):(i*(neta/ntheta))
  paths <- rep(NA, length(ms) )
  ix <- which(lat[1:length(ms)]=="1")
  
  paths[ix] <- paste0(eta[ms[ix]],"@",lat[ms[ix]])
  paths[-ix] <- paste0(eta[ms[-ix]]," (",lat[ms[-ix]],")")
  
  paths <- paste(paths, collapse = "\n")
  theta_syntax[i] <- paste(theta[i], "BY", paths, ";", sep = " ")
}
rm(i, ms, paths, ix)

if(!second.order.trait){
  theta_syntax <- rep(NA, ntheta)
  for(i in 1:ntheta){
    ms <- ((i-1)*(nobs/ntheta) + 1):(i*(nobs/ntheta))
    paths <- rep(NA, length(ms) )
    ix <- which(lat[1:length(ms)]=="1")
    
    paths[ix] <- paste0(obs[ms[ix]],"@",lat[ms[ix]])
    paths[-ix] <- paste0(obs[ms[-ix]]," (",lat[ms[-ix]],")")
    
    paths <- paste(paths, collapse = "\n")
    theta_syntax[i] <- paste(theta[i], "BY", paths, ";", sep = " ")
  }
  rm(i, ms, paths, ix)
}

# write intercetps for the latent states
eta_syntax_int <- rep(NA, neta)
for(i in 1:neta){
  ms <- ((i-1)*(nobs/neta) + 1):(i*(nobs/neta))
  paths <- rep(NA, length(ms) )
  ix <- which(ins[1:length(ms)]=="0")
  
  paths[ix] <- paste0( "[", obs[ms[ix]],"@",ins[ms[ix]], "]", " ;")
  paths[-ix] <- paste0("[", obs[ms[-ix]], "*]", " (",ins[ms[-ix]],")", " ;")
  
  eta_syntax_int[i] <- paste(paths, collapse = "\n") 
}
rm(i, ms, paths, ix)

if(!second.order.trait){
  eta_syntax_int <- NULL
}


# write intercepts for the latent trait

theta_syntax_int <- rep(NA, ntheta)
for(i in 1:ntheta){
  ms <- ((i-1)*(neta/ntheta) + 1):(i*(neta/ntheta))
  paths <- rep(NA, length(ms) )
  ix <- which(int[1:length(ms)]=="0")
  
  paths[ix] <- paste0("[", eta[ms[ix]],"@",int[ms[ix]], "]", " ;")
  paths[-ix] <- paste0("[", eta[ms[-ix]], "*]", " (",int[ms[-ix]],")", " ;")
  
  theta_syntax_int[i] <- paste(paths, collapse = "\n")
}
rm(i, ms, paths, ix)

if(!second.order.trait){
  theta_syntax_int <- rep(NA, ntheta)
  for(i in 1:ntheta){
    ms <- ((i-1)*(nobs/ntheta) + 1):(i*(nobs/ntheta))
    paths <- rep(NA, length(ms) )
    ix <- which(int[1:length(ms)]=="0")
    
    paths[ix] <- paste0("[", obs[ms[ix]],"@",int[ms[ix]], "]", " ;")
    paths[-ix] <- paste0("[", obs[ms[-ix]], "*]", " (",int[ms[-ix]],")", " ;")
    
    theta_syntax_int[i] <- paste(paths, collapse = "\n")
  }
  rm(i, ms, paths, ix)
  
}


theta_syntax_means <- paste0("[", theta, "*] ;", collapse = "\n")

# write variance components

var_syntax <- paste(paste0(obs, "(", o_red, ") ;", collapse = "\n"),
                     paste0(eta, "(", s_red, ") ;", collapse = "\n"),
                     paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                     sep = "\n")

# write 0 correlations among occasions  and traits
if(!second.order.trait){latent.v <- c(eta, theta)

cor_syntax <- rep(NA, length(latent.v)-1)
for(i in 1:(length(latent.v)-1)){
  
  cor_syntax[i] <- paste0(latent.v[i], " WITH ",
                          paste0(latent.v[-(1:i)], "@0", collapse = "\n"), " ;")
}
rm(i)
}

# define computed parameters
def_syntax <- paste(paste0("NEW(", rel, ") ;", collapse = "\n" ),
                    paste0("NEW(", con, ") ;", collapse = "\n" ),
                    paste0("NEW(", spe, ") ;", collapse = "\n" ),
                    sep = "\n")

lat_r <- rep(lat, each = nobs/neta)
if(!second.order.trait){
  lat_r <- lat
}

total_var <- paste0(lat_r, "^2*",las, "^2*", 
                     rep(t_var, each = neta/ntheta), "+",
                     las, "^2*", rep(s_red, each = nobs/neta), 
                    "+", o_red)

trait_var <- paste0(lat_r, "^2*",las, "^2*", 
                    rep(t_var, each = neta/ntheta))

state_var <- paste0(las, "^2*", rep(s_red, each = nobs/neta))


rel_syntax <- paste0(rel, "= (", trait_var, "+", state_var,
                     ")/", "\n(", total_var, ") ;", collapse = "\n")

spe_syntax <- paste0(spe, "= (", state_var,
                     ")/", "\n(", total_var, ") ;", collapse = "\n")

con_syntax <- paste0(con, "= (", trait_var,
                     ")/", "\n(", total_var, ") ;", collapse = "\n")



complete_syntax <- paste("\nMODEL:",
                        "\n! Latent state variables",
                        paste0(eta_syntax, collapse = "\n"),
                        "\n! Latent trait variables",
                        paste0(theta_syntax, collapse = "\n"),
                        "\n! Intercepts and means",
                        paste0(eta_syntax_int, collapse = "\n"),
                        paste0(theta_syntax_int, collapse = "\n"),
                        paste0(theta_syntax_means, collapse = "\n"),
                        "\n! Variance components",
                        var_syntax,
                        "\nMODEL CONSTRAINT:",
                        def_syntax,
                        "\n! Reliabily",
                        rel_syntax,
                        "\n! Occasion specificity",
                        spe_syntax,
                        "\n! Consistency",
                        con_syntax,
                        sep = "\n")

if(!second.order.trait){
  complete_syntax <- paste("\nMODEL:",
                           "\n! Latent state variables",
                           paste0(eta_syntax, collapse = "\n"),
                           "\n! Latent trait variables",
                           paste0(theta_syntax, collapse = "\n"),
                           "\n! Intercepts and means",
                           paste0(eta_syntax_int, collapse = "\n"),
                           paste0(theta_syntax_int, collapse = "\n"),
                           paste0(theta_syntax_means, collapse = "\n"),
                           "\n! Variance components",
                           var_syntax,
                           "\n! Set zero covariances",
                           paste0(cor_syntax, collapse = "\n"),
                           sep = "\n")
}
return(complete_syntax)

}
