# Functions to compute the variance coefficients of the different models with time varying parameters.
cuts.var.coeff.tv <- function(I, nT, within.parameters, between.parameters){
  
  CT.var <- matrix(between.parameters$CT.var, ncol = nT, nrow = I)
  UT.var <- matrix(between.parameters$UT.var, ncol = nT, nrow = I)
  CS.var <- matrix(within.parameters$CS.var, ncol = nT, nrow = I, byrow = TRUE)
  US.var <- t(within.parameters$US.var)
  
  #multiply variances by loadings
  
  CT.var <- t(between.parameters$loadings ^ 2) * CT.var
  UT.var <- t(between.parameters$loadings.UT ^ 2) * UT.var
  CS.var <- t(within.parameters$loadings ^ 2) * CS.var
  
  # Compute variance components
  
  total.variance <- CT.var + UT.var + CS.var + UT.var
  rel <- (CT.var + UT.var + CS.var) / total.variance
  spe <- (CS.var) / total.variance
  tcon <- (CT.var + UT.var) / total.variance
  ccon <- (CT.var) / total.variance
  ucon <- (UT.var) / total.variance
  
  out <- data.frame(rbind(ccon, ucon, tcon, spe, rel))
  row.names(out) <- paste0(rep(c("ccon", "ucon", "tcon", "spe", "rel"), each = I), 
                           1:I)
  names(out) <- paste0("time", 1:nT)
  
  return(out)
}

msst.var.coeff.tv <- function(I, nT, within.parameters, between.parameters){
  
  trait.var <- matrix(between.parameters$trait.var, ncol = nT, nrow = I)
  state.var <- matrix(within.parameters$state.var, ncol = nT, nrow = I, byrow = TRUE)
  error.var <- t(within.parameters$error.var)
  
  #multiply variances by loadings
  
  trait.var <- t(between.parameters$loadings ^ 2) * trait.var
  state.var <- t(within.parameters$loadings ^ 2) * state.var
  
  # Compute variance components
  
  total.variance <- trait.var + state.var + error.var
  rel <- (trait.var + state.var) / total.variance
  con <- (trait.var) / total.variance
  spe <- (state.var) / total.variance
  
  out <- data.frame(rbind(con, spe, rel))
  row.names(out) <- paste0(rep(c("con","spe", "rel"), each = I), 
                           1:I)
  names(out) <- paste0("time", 1:nT)
  
  return(out)
}

tso.var.coeff.tv <- function(I, nT, within.parameters, between.parameters){
  
  trait.var <- matrix(between.parameters$trait.ind.var, ncol = nT, nrow = I)
  state.var <- matrix(within.parameters$state.var, ncol = nT, nrow = I, byrow = TRUE)
  error.var <- t(within.parameters$error.var)
  
  occ.var <- rep(NA, nT - 1) #create components of the explained variance due to the autoregressive effect
  ar.effect_2 <- within.parameters$ar.effect ^ 2
  
  # This loop computes the variance of each occasion variable as a linear 
  #   combination of the previous state residuals
  for(i in 1:(nT-1)){
    ar_2_prod <- rep(NA, i)
    for(j in 1:i){
      ar_2_prod[j] <- prod(ar.effect_2[j:i])
    }
    occ.var[i] <- sum(ar_2_prod*within.parameters$state.var[1:i])
  }
  rm(i, j, ar_2_prod)
  
  occ.var <- c(0, occ.var)
  
  
  #multiply variances by loadings
  
  trait.var <- t(between.parameters$loadings ^ 2) * trait.var
  state.var <- t(within.parameters$loadings ^ 2) * state.var
  occ.var <- t((within.parameters$loadings ^ 2) * occ.var)
  
  #compute variance components
  total.variance <- trait.var + occ.var + state.var + error.var
  
  rel <- (trait.var + occ.var + state.var) / total.variance
  con <- (trait.var + occ.var) / total.variance
  pred <- (trait.var) / total.variance
  upred <- (occ.var) / total.variance
  spe <- (state.var) / total.variance
  
  out <- data.frame(rbind( pred, upred, con, spe, rel))
  row.names(out) <- paste0(rep(c( "pred", "upred","con","spe", "rel"), each = I), 1:I)
  names(out) <- paste0("time", 1:nT)
  
  return(out)
}