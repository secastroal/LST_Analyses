write.mltso.to.Mplus <- function(dt, scale.equivalence = list(tau = FALSE) ){
  nobs <- dim(dt)[2]
  #observed variables
  obs <- names(dt)
  #common state
  eta <- "eta"
  #common trait
  theta <- paste0("theta", 1:nobs)
  #loadings
  las <- paste0("las", 1:nobs) #lambda common states
  
  #variance components
  
  o_red <- paste0("o_err", 1:nobs)
  s_var <- "s_var"
  t_var <- "t_var"
  
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
  paths[-ix] <- paste0(obs[-ix]," (",las[-ix],")")
  
  paths <- paste(paths, collapse = "\n")
  eta_syntax <- paste(eta, "BY", paths, " (&1) ;", sep = " ")
  rm(paths, ix)
  
  
  
  theta_int_syntax <- paste0("[", obs, c("@0] ;", rep("] ;", nobs-1)), collapse = "\n")
  
  theta_mean_syntax <- paste0("[", theta, "*] ;", collapse = "\n")
  
  
  var_within_syntax <- paste(paste0(obs, "(", o_red, ") ;", collapse = "\n"),
                             paste0(eta, "(", s_var, ") ;", collapse = "\n"),
                             sep = "\n")
  
  
  var_between_syntax <- paste(paste0(obs, "@0;", collapse = "\n"),
                              paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                              sep = "\n")
  
  beta_syntax <- paste(eta, "ON", paste0(eta, "&1"), ";", sep = " ")
  
  
  
  complete_syntax <- paste("\nMODEL:",
                           "\n%WITHIN%",
                           "\n! common states",
                           paste0(eta_syntax, collapse = "\n"),
                           beta_syntax,
                           var_within_syntax,
                           "\n%BETWEEN%",
                           "\n! common trait",
                           theta_syntax,
                           theta_int_syntax,
                           theta_mean_syntax,
                           var_between_syntax,
                           sep = "\n")
  
  return(complete_syntax)
  
}
