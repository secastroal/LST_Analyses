write.mlcuts.to.Mplus <- function(dt, state.trait.invariance = FALSE ){
  nobs <- dim(dt)[2]
  #observed variables
  obs <- names(dt)
  #common state
  cs <- "CS"
  #common trait
  ct <- "CT"
  #loadings
  lcs <- paste0("lcs", 1:nobs) #lambda common states
  lct <- paste0("lct", 1:nobs) #lambda common trait
  
  # if-else statements to define constraints (parameters equalities)
  lcs[1] <- 1
  
  lct[1] <- 1
  
  #state-trait invariance
  
  if(state.trait.invariance){
    lcs <- lct
  }
  
  #syntanxes
  
  paths <- rep(NA, nobs )
  ix <- which(lct=="1")
  
  paths[ix] <- paste0(obs[ix],"@",lct[ix])
  paths[-ix] <- paste0(obs[-ix]," (",lct[-ix],")")
  
  paths <- paste(paths, collapse = "\n")
  ct_syntax <- paste(ct, "BY", paths, ";", sep = " ")
  rm(paths, ix)
  
  
  paths <- rep(NA, nobs )
  ix <- which(lcs=="1")
  
  paths[ix] <- paste0(obs[ix],"@",lcs[ix])
  paths[-ix] <- paste0(obs[-ix]," (",lcs[-ix],")")
  
  paths <- paste(paths, collapse = "\n")
  cs_syntax <- paste(cs, "BY", paths, ";", sep = " ")
  rm(paths, ix)
  
  complete_syntax <- paste("\nMODEL:",
                           "\n%WITHIN%",
                           "\n! common states",
                           paste0(cs_syntax, collapse = "\n"),
                           "\n%BETWEEN%",
                           "\n! common trait",
                           paste0(ct_syntax, collapse = "\n"),
                           sep = "\n")
  
  return(complete_syntax)
  
}
