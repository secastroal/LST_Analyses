# Functions to compute the variance coefficients of the different models.

# cuts coefficients are based on eq. 2, 4, 5, 9, 10, and 11 of Geiser and Lockhart (2012)
cuts.var.coeff <- function(within.parameters, between.parameters){
  # Total variances. There are not loadings weigthing the variances of the unique traits because they are 
  # equal to 1 when asusming time-invariant parameters.
  total.variance <- (between.parameters$loadings ^ 2) * between.parameters$CT.var + between.parameters$UT.var +
    (within.parameters$loadings ^ 2) * within.parameters$CS.var + within.parameters$US.var
  # Reliability
  rel <- ((between.parameters$loadings ^ 2) * between.parameters$CT.var + between.parameters$UT.var +
            (within.parameters$loadings ^ 2) * within.parameters$CS.var) / total.variance
  # Occasion specificity
  spe <- ((within.parameters$loadings ^ 2) * within.parameters$CS.var) / total.variance
  # Total consistency
  tcon <- ((between.parameters$loadings ^ 2) * between.parameters$CT.var + between.parameters$UT.var) / total.variance
  # Common consistency
  ccon <- ((between.parameters$loadings ^ 2) * between.parameters$CT.var) / total.variance
  # Unique consistency
  ucon <- (between.parameters$UT.var) / total.variance
  
  out <- data.frame(c(ccon, ucon, tcon, spe, rel))
  row.names(out) <- paste0(rep(c("ccon", "ucon", "tcon", "spe", "rel"), each = length(within.parameters$loadings)), 
                           1:length(within.parameters$loadings))
  names(out) <- "Var.Coeff"
  
  return(out)
}

# msst coefficients are based on eq. 2, 3, 4, and 5 of Geiser and Lockhart (2012)
msst.var.coeff <- function(within.parameters, between.parameters){
  # Total variance
  total.variance <- (between.parameters$loadings ^ 2) * between.parameters$trait.var + 
    (within.parameters$loadings ^ 2) * within.parameters$state.var + within.parameters$error.var
  # Reliability
  rel <- ((between.parameters$loadings ^ 2) * between.parameters$trait.var + 
    (within.parameters$loadings ^ 2) * within.parameters$state.var) / total.variance
  # Consistency
  con <- ((between.parameters$loadings ^ 2) * between.parameters$trait.var) / total.variance
  # Occasion specificity
  spe <- ((within.parameters$loadings ^ 2) * within.parameters$state.var) / total.variance
  
  out <- data.frame(c(con, spe, rel))
  row.names(out) <- paste0(rep(c("con","spe", "rel"), each = length(within.parameters$loadings)), 
                           1:length(within.parameters$loadings))
  names(out) <- "Var.Coeff"
  
  return(out)
}

# tso coefficients are based on eq. 20 - 24, of Eid et al.(2017). However, I use the total variance instead
# of the "reliable" variance as the denominator.
tso.var.coeff <- function(I, nT, within.parameters, between.parameters){
  
  trait.var <- matrix(between.parameters$trait.ind.var, ncol = nT, nrow = I)
  state.var <- matrix(within.parameters$state.var, ncol = nT, nrow = I)
  error.var <- matrix(within.parameters$error.var, ncol = nT, nrow = I)
  
  occ.var <- rep(NA, nT - 1) #create components of the explained variance due to the autoregressive effect
  ar.effect_2 <- within.parameters$ar.effect ^ 2
  # As parameters are time invariant, the terms for each occasion can be simplified as the product of the
  # state residual variance and the geometric serie of the autoregressive effect at each specific occasion.
  for(i in 1:(nT-1)){
    occ.var[i] <- (ar.effect_2 * (1 - ar.effect_2 ^ i))/(1 - ar.effect_2) #geometric serie of the autoregressive effect at each time point
  }
  
  occ.var <- c(0, occ.var) * within.parameters$state.var #product of the geometric series and the state residual variance
  
  #multiply variances by loadings / compute weighted variances
  
  state.var.w <- (within.parameters$loadings ^ 2) * state.var
  occ.var.w <- (within.parameters$loadings ^ 2) %o% occ.var
  
  #compute variance components
  total.variance <- trait.var + occ.var.w + state.var.w + error.var # total variance
  
  rel <- (trait.var + occ.var.w + state.var.w) / total.variance # reliability
  con <- (trait.var + occ.var.w) / total.variance # consistency
  pred <- (trait.var) / total.variance # predictability by trait
  upred <- (occ.var.w) / total.variance # unpredictability by trait
  spe <- (state.var.w) / total.variance # occasion specificity
  
  out <- data.frame(rbind( pred, upred, con, spe, rel))
  row.names(out) <- paste0(rep(c( "pred", "upred","con","spe", "rel"), each = I), 1:I)
  names(out) <- paste0("time", 1:nT)
  
  return(out)
}