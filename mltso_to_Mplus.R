write.mltso.to.Mplus <- function(dt, scale.equivalence = list(tau = FALSE) ){
  nobs <- dim(dt)[2]
  #observed variables
  obs <- names(dt)
  #Latent state residual
  eta <- "eta"
  #Latent trait indicator variables
  theta <- paste0("theta", 1:nobs)
  #loadings
  las <- rep(" ", nobs) #lambda latent states residuals
  
  #variance components
  
  o_red <- paste0("o_err", 1:nobs)
  s_var <- "s_var"
  t_var <- paste0("t_var", 1:nobs)
  
  # if-else statements to define constraints (parameters equalities)
  las[1] <- 1
  
  
  # scale equivalence
  
  if(scale.equivalence$tau){
    las <- 1
  }
  
  #syntanxes
  
  theta_syntax <- paste0(theta, " BY ", obs, "@1 ;", collapse = "\n")
  
  
  paths <- rep(NA, nobs )
  ix <- which(las=="1")
  
  paths[ix] <- paste0(obs[ix],"@",las[ix])
  paths[-ix] <- paste0(obs[-ix],las[-ix])
  
  paths <- paste(paths, collapse = "\n")
  eta_syntax <- paste(eta, "BY", paths, " (&1) ;", sep = " ")
  rm(paths, ix)
  
  
  
  theta_mean_syntax <- paste0("[", theta, "@0] ;", collapse = "\n")
  
  theta_int_syntax <- paste0("[", obs, "*] ;", collapse = "\n")
  
  
  var_within_syntax <- paste(paste0(obs, "(", o_red, ") ;", collapse = "\n"),
                             paste0(eta, "(", s_var, ") ;", collapse = "\n"),
                             sep = "\n")
  
  
  var_between_syntax <- paste(paste0(obs, "@0.001;", collapse = "\n"),
                              paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                              sep = "\n")
  
  beta_syntax <- paste(eta, "ON", paste0(eta, "&1"), ";", sep = " ")
  
  
  
  complete_syntax <- paste("\nMODEL:",
                           "\n%WITHIN%",
                           "\n! Latent state residual",
                           paste0(eta_syntax, collapse = "\n"),
                           "\n! Autoregressive effect on the latent state residual",
                           beta_syntax,
                           "\n! Latent errors and latent state variances",
                           var_within_syntax,
                           "\n%BETWEEN%",
                           "\n! Latent trait-indicator variables",
                           theta_syntax,
                           "\n! Intercepts",
                           theta_int_syntax,
                           "\n! Fixing means of latent trait-indicator variables at 0",
                           theta_mean_syntax,
                           "\n! Fixing level-2 variances at 0 and latent trait-indicator variances",
                           var_between_syntax,
                           sep = "\n")
  
  return(complete_syntax)
  
}
