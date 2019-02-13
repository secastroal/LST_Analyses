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
  
  trait.var <- rep(between.parameters$trait.ind.var, times = nT)
  state.var <- rep(within.parameters$state.var, times = nT * I)
  error.var <- rep(within.parameters$error.var, times = nT)
  
  
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