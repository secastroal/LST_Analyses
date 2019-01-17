write.mlmsst.to.Mplus <- function(dt, scale.equivalence = list(tau = FALSE, theta = FALSE) ){
  nobs <- dim(dt)[2]
  #observed variables
  obs <- names(dt)
  #latent state residual
  eta <- "eta"
  #latent trait variable
  theta <- "theta"
  #latent trait indicator variables
  theta_ind <- paste0("theta", 1:nobs)
  #loadings
  las <- paste0("las", 1:nobs) #lambda latent states residuals
  lat <- paste0("lat", 1:nobs) #lambda latent trait
  int <- paste0("int", 1:nobs) #intercepts trait
  
  #variance components
  
  o_red <- paste0("o_err", 1:nobs)
  s_var <- "s_var"
  t_var <- "t_var"
  
  # if-else statements to define constraints (parameters equalities)
  las[1] <- 1
  
  lat[1] <- 1
  
  int[1] <- 0
  
  # scale equivalence
  
  if(scale.equivalence$tau){
    las <- 1
    }
  
  if(scale.equivalence$theta){
    lat <- 1
    int <- 0
  }
  
  #syntanxes
  
  theta_ind_syntax <- paste(theta_ind, " BY ", paste0(obs, "@1 ;"), collapse = "\n")
  
  paths <- rep(NA, nobs )
  ix <- which(lat=="1")
  
  paths[ix] <- paste0(theta_ind[ix],"@",lat[ix])
  paths[-ix] <- paste0(theta_ind[-ix]," (",lat[-ix],")")
  
  paths <- paste(paths, collapse = "\n")
  theta_syntax <- paste(theta, "BY", paths, ";", sep = " ")
  rm(paths, ix)
  
  
  paths <- rep(NA, nobs )
  ix <- which(las=="1")
  
  paths[ix] <- paste0(obs[ix],"@",las[ix])
  paths[-ix] <- paste0(obs[-ix]," (",las[-ix],")")
  
  paths <- paste(paths, collapse = "\n")
  eta_syntax <- paste(eta, "BY", paths, ";", sep = " ")
  rm(paths, ix)
  
  paths <- rep(NA, nobs )
  ix <- which(int=="0")
  
  paths[ix] <- paste0("[", theta_ind[ix],"@",int[ix], "] ;")
  paths[-ix] <- paste0("[", theta_ind[-ix],"*] (",int[-ix],") ;")
  
  theta_int_syntax <- paste(paths, collapse = "\n")
  rm(paths, ix)
  
  theta_ind_int_syntax <- paste0("[", obs, "@0] ;", collapse = "\n" )
  
  theta_mean_syntax <- paste0("[", theta, "*] ;", collapse = "\n")
  
  
  var_within_syntax <- paste(paste0(obs, "(", o_red, ") ;", collapse = "\n"),
                      paste0(eta, "(", s_var, ") ;", collapse = "\n"),
                      sep = "\n")
  
  
  var_between_syntax <- paste(paste0(c(obs, theta_ind), "@0;", collapse = "\n"),
                              paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                              sep = "\n")
  
  
  
  complete_syntax <- paste("\nMODEL:",
                           "\n%WITHIN%",
                           "\n! Latent state residual",
                           paste0(eta_syntax, collapse = "\n"),
                           "\n! Latent errors and latent state variances",
                           var_within_syntax,
                           "\n%BETWEEN%",
                           "\n! Indicator trait variables",
                           theta_ind_syntax,
                           "\n! Latent trait variable",
                           paste0(theta_syntax, collapse = "\n"),
                           "\n! Means and Intercepts",
                           theta_ind_int_syntax,
                           paste0(theta_int_syntax, collapse = "\n"),
                           theta_mean_syntax,
                           "\n! Fixing level-2 residual variances at 0 and latent trait variance",
                           var_between_syntax,
                           sep = "\n")
  
  return(complete_syntax)
  
}
