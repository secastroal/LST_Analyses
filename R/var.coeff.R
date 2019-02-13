# Functions to compute the variance coefficients of the different models.
cuts.var.coeff <- function(within.parameters, between.parameters){
  total.variance <- (between.parameters$loadings ^ 2) * between.parameters$CT.var + between.parameters$UT.var +
    (within.parameters$loadings ^ 2) * within.parameters$CS.var + within.parameters$US.var
  rel <- ((between.parameters$loadings ^ 2) * between.parameters$CT.var + between.parameters$UT.var +
            (within.parameters$loadings ^ 2) * within.parameters$CS.var) / total.variance
  spe <- ((within.parameters$loadings ^ 2) * within.parameters$CS.var) / total.variance
  tcon <- ((between.parameters$loadings ^ 2) * between.parameters$CT.var + between.parameters$UT.var) / total.variance
  ccon <- ((between.parameters$loadings ^ 2) * between.parameters$CT.var) / total.variance
  ucon <- (between.parameters$UT.var) / total.variance
  
  out <- data.frame(c(ccon, ucon, tcon, spe, rel))
  row.names(out) <- paste0(rep(c("ccon", "ucon", "tcon", "spe", "rel"), each = length(within.parameters$loadings)), 1:length(within.parameters$loadings))
  names(out) <- "Var.Coeff"
  
  return(out)
}

msst.var.coeff <- function(within.parameters, between.parameters){
  total.variance <- (between.parameters$loadings ^ 2) * between.parameters$trait.var + 
    (within.parameters$loadings ^ 2) * within.parameters$state.var + within.parameters$error.var
  rel <- ((between.parameters$loadings ^ 2) * between.parameters$trait.var + 
    (within.parameters$loadings ^ 2) * within.parameters$state.var) / total.variance
  con <- ((between.parameters$loadings ^ 2) * between.parameters$trait.var) / total.variance
  spe <- ((within.parameters$loadings ^ 2) * within.parameters$state.var) / total.variance
  
  out <- data.frame(c(con, spe, rel))
  row.names(out) <- paste0(rep(c("con","spe", "rel"), each = length(within.parameters$loadings)), 1:length(within.parameters$loadings))
  names(out) <- "Var.Coeff"
  
  return(out)
}

tso.var.coeff <- function(I, nT, within.parameters, between.parameters){
  
  trait.var <- matrix(between.parameters$trait.ind.var, ncol = nT, nrow = I)
  state.var <- matrix(within.parameters$state.var, ncol = nT, nrow = I)
  error.var <- matrix(within.parameters$error.var, ncol = nT, nrow = I)
  
  occ.var <- rep(NA, nT - 1) #create components of the explained variance due to the autoregressive effect
  ar.effect_2 <- within.parameters$ar.effect ^ 2
  for(i in 1:(nT-1)){
    occ.var[i] <- (ar.effect_2 * (1 - ar.effect_2 ^ i))/(1 - ar.effect_2)
  }
  
  occ.var <- c(0, occ.var) * within.parameters$state.var
  
  #multiply variances by loadings
  
  state.var <- (within.parameters$loadings ^ 2) * state.var
  occ.var <- (within.parameters$loadings ^ 2) %o% occ.var
  
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