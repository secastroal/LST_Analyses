write.mlmsst.to.Mplus <- function(dt, scale.equivalence = list(tau = FALSE, theta = FALSE) ){
  nobs <- dim(dt)[2]
  #observed variables
  obs <- names(dt)
  #common state
  eta <- "eta"
  #common trait
  theta <- "theta"
  #loadings
  las <- paste0("las", 1:nobs) #lambda common states
  lat <- paste0("lat", 1:nobs) #lambda common trait
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
  
  paths <- rep(NA, nobs )
  ix <- which(lat=="1")
  
  paths[ix] <- paste0(obs[ix],"@",lat[ix])
  paths[-ix] <- paste0(obs[-ix]," (",lat[-ix],")")
  
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
  
  paths[ix] <- paste0("[", obs[ix],"@",int[ix], "] ;")
  paths[-ix] <- paste0("[", obs[-ix],"*] (",int[-ix],") ;")
  
  theta_int_syntax <- paste(paths, collapse = "\n")
  rm(paths, ix)
  
  theta_mean_syntax <- paste0("[", theta, "*] ;", collapse = "\n")
  
  
  var_within_syntax <- paste(paste0(obs, "(", o_red, ") ;", collapse = "\n"),
                      paste0(eta, "(", s_var, ") ;", collapse = "\n"),
                      sep = "\n")
  
  
  var_between_syntax <- paste(paste0(obs, "@0;", collapse = "\n"),
                              paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                              sep = "\n")
  
  
  
  complete_syntax <- paste("\nMODEL:",
                           "\n%WITHIN%",
                           "\n! common states",
                           paste0(eta_syntax, collapse = "\n"),
                           var_within_syntax,
                           "\n%BETWEEN%",
                           "\n! common trait",
                           paste0(theta_syntax, collapse = "\n"),
                           paste0(theta_int_syntax, collapse = "\n"),
                           theta_mean_syntax,
                           var_between_syntax,
                           sep = "\n")
  
  return(complete_syntax)
  
}
