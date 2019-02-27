# Functions to compute the variance coefficients of the different models with time varying parameters.

# cuts coefficients are based on eq. 2, 4, 5, 9, 10, and 11 of Geiser and Lockhart (2012)
cuts.var.coeff.tv <- function(I, nT, within.parameters, between.parameters){
  
  CT.var <- matrix(between.parameters$CT.var, ncol = nT, nrow = I)
  UT.var <- matrix(between.parameters$UT.var, ncol = nT, nrow = I)
  CS.var <- matrix(within.parameters$CS.var, ncol = nT, nrow = I, byrow = TRUE)
  US.var <- t(within.parameters$US.var)
  
  #multiply variances by loadings / compute weighted variances
  
  CT.var.w <- t(between.parameters$loadings ^ 2) * CT.var
  UT.var.w <- t(between.parameters$loadings.UT ^ 2) * UT.var
  CS.var.w <- t(within.parameters$loadings ^ 2) * CS.var
  
  # Compute variance components
  
  total.variance <- CT.var.w + UT.var.w + CS.var.w + UT.var # total variance
  rel <- (CT.var.w + UT.var.w + CS.var.w) / total.variance # reliability
  spe <- (CS.var.w) / total.variance # occasion specificity
  tcon <- (CT.var.w + UT.var.w) / total.variance # total consistency
  ccon <- (CT.var.w) / total.variance # common consistency
  ucon <- (UT.var.w) / total.variance # unique consistency
  
  out <- data.frame(rbind(ccon, ucon, tcon, spe, rel))
  row.names(out) <- paste0(rep(c("ccon", "ucon", "tcon", "spe", "rel"), each = I), 
                           1:I)
  names(out) <- paste0("time", 1:nT)
  
  return(out)
}

# msst coefficients are based on eq. 2, 3, 4, and 5 of Geiser and Lockhart (2012)
msst.var.coeff.tv <- function(I, nT, within.parameters, between.parameters){
  
  trait.var <- matrix(between.parameters$trait.var, ncol = nT, nrow = I)
  state.var <- matrix(within.parameters$state.var, ncol = nT, nrow = I, byrow = TRUE)
  error.var <- t(within.parameters$error.var)
  
  #multiply variances by loadings / compute weighted variances
  
  trait.var.w <- t(between.parameters$loadings ^ 2) * trait.var
  state.var.w <- t(within.parameters$loadings ^ 2) * state.var
  
  # Compute variance components
  
  total.variance <- trait.var.w + state.var.w + error.var # total variance
  rel <- (trait.var.w + state.var.w) / total.variance # reliability
  con <- (trait.var.w) / total.variance # consistency
  spe <- (state.var.w) / total.variance # occasion specificity
  
  out <- data.frame(rbind(con, spe, rel))
  row.names(out) <- paste0(rep(c("con","spe", "rel"), each = I), 
                           1:I)
  names(out) <- paste0("time", 1:nT)
  
  return(out)
}

# tso coefficients are based on eq. 20 - 24, of Eid et al.(2017). However, I use the total variance instead
# of the "reliable" variance as the denominator.
tso.var.coeff.tv <- function(I, nT, within.parameters, between.parameters){
  
  trait.var <- matrix(between.parameters$trait.ind.var, ncol = nT, nrow = I)
  state.var <- matrix(within.parameters$state.var, ncol = nT, nrow = I, byrow = TRUE)
  error.var <- t(within.parameters$error.var)
  
  occ.var <- rep(NA, nT - 1) #create components of the explained variance due to the autoregressive effect
  ar.effect_2 <- within.parameters$ar.effect ^ 2
  
  # This loop computes the variance of each occasion variable as a linear 
  #   combination of the previous state residuals variances weigthed by the product 
  #   of the autoregressive effects
  for(i in 1:(nT-1)){
    ar_2_prod <- rep(NA, i)
    for(j in 1:i){
      ar_2_prod[j] <- prod(ar.effect_2[j:i])
    }
    occ.var[i] <- sum(ar_2_prod*within.parameters$state.var[1:i])
  }
  rm(i, j, ar_2_prod)
  
  occ.var <- c(0, occ.var)
  
  
  #multiply variances by loadings /  compute weighted variances
  
  trait.var.w <- t(between.parameters$loadings ^ 2) * trait.var
  state.var.w <- t(within.parameters$loadings ^ 2) * state.var
  occ.var.w <- t((within.parameters$loadings ^ 2) * occ.var)
  
  #compute variance components
  total.variance <- trait.var.w + occ.var.w + state.var.w + error.var # total variance
  rel <- (trait.var.w + occ.var.w + state.var.w) / total.variance # reliability
  con <- (trait.var.w + occ.var.w) / total.variance # consistency
  pred <- (trait.var.w) / total.variance # predictability by trait
  upred <- (occ.var.w) / total.variance # unpredictability by trait
  spe <- (state.var.w) / total.variance # occasion specificity
  
  out <- data.frame(rbind( pred, upred, con, spe, rel))
  row.names(out) <- paste0(rep(c( "pred", "upred","con","spe", "rel"), each = I), 1:I)
  names(out) <- paste0("time", 1:nT)
  
  return(out)
}