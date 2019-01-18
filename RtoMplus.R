## This document contains functions that facilitate writing Mplus sintaxes of three LST models and their
##    respective multilevel version.

# Contents
# 1.0 Mplus_analysis_syntax
# 2.0 write.msst.to.Mplus
# 3.0 write.cuts.to.Mplus
# 4.0 write.tso.to.Mplus
# 5.0 write.mlmsst.to.Mplus
# 6.0 write.mlcuts.to.Mplus
# 7.0 write.mltso.to.Mplus

# 1.0 Mplus_analysis_syntax ----

Mplus_analysis_syntax <- function(usevariables, cluster = NULL, analysis_type = c("GENERAL", "TWOLEVEL"),
                                  estimator = c("ML", "MLR", "BAYES"), iterations, h1iterations = NULL, 
                                  processors = NULL){
  usevariables <- paste("USEVAR", "=", paste(usevariables, collapse = "\n"), ";", sep = " ")
  if(!is.null(cluster)){
    cluster <- paste("CLUSTER", "=", cluster, ";", sep = " ")
  }
  analysis_type <- paste("TYPE", "=", analysis_type, ";", sep = " ")
  estimator <- paste("ESTIMATOR", "=", estimator, ";", sep = " ")
  if(estimator != "BAYES"){
    iterations <- paste("ITERATIONS", "=", iterations, ";", sep = " ")
    h1iterations <- paste("H1ITERATIONS", "=", h1iterations, ";", sep = " ")
    processors <- NULL
  }else{
    iterations <- paste("BITERATIONS", "=", paste0("(", iterations, ")"), ";", sep = " " )
    h1iterations <- NULL
    processors <- paste("PROCESSORS", "=", processors, sep = " ")
  }
  
  analysis_syntax <- paste(usevariables,
                           cluster,
                           "\nANALYSIS:",
                           analysis_type,
                           estimator,
                           iterations, 
                           h1iterations,
                           processors,
                           sep = "\n")
  cat(analysis_syntax)
  return(analysis_syntax)
}

# 2.0 write.msst.to.Mplus ----

# Function to write Mplus syntax for the multistate-singletrait model given a data.frame

# dt is the data.frame. It is supposed to be a longitudinal data in wide format.
#     Only the variables that are measured repeatedly are included in the data set, 
#     ID or sex variables should not be included.
# neta is the number of latent state variables
# ntheta is the number of latent trait variables
# equiv.assumption are the assumptions tau-equivalence, essential-tau-equivalence, 
#     and tau-congenericity in relation to the latent states (tau) and the latent traits (theta).
#     use "equi", "ess", and "cong" 
# scale.invariance define whether the measurement model is invariant over time or not.
#     It can be define for both intercepts and loadings of the latent states and of the
#     latent traits.

write.msst.to.Mplus <- function (dt, neta, ntheta = 1, 
                                 equiv.assumption = list(tau = "cong", theta = "equi"), 
                                 scale.invariance = list(lait0 = FALSE, lait1 = FALSE, lat0 = FALSE, lat1 = FALSE),
                                 homocedasticity.assumption = list(error = FALSE, state.red = FALSE),
                                 second.order.trait = TRUE)
{
  nobs <- dim(dt)[2]  
  #labels
  ind_s <- paste0(rep(1:(nobs/neta) ,neta),
                  rep(1:neta, each=(nobs/neta)))
  ind_t <- paste0(1:neta, rep(1,neta))
  if(!second.order.trait){
    ind_t <- ind_s
  }
  
  #lambdas for states and traits
  las <- paste0("las", ind_s)
  lat <- paste0("lat", ind_t)
  
  # intercepts for observed and state variables
  ins <- paste0("ins", ind_s)
  int <- paste0("int", ind_t)
  if(!second.order.trait){
    int <- NULL
  }
  
  
  #latent variables' names
  eta <- paste0(rep("eta", neta), 1:neta)
  
  theta <- paste0(rep("theta", ntheta), 1:ntheta)
  
  # observed variables names
  obs <- names(dt)
  
  # variances and residuals' labels
  
  o_red <- paste0("o_err", ind_s)
  s_red <- paste0("s_var", ind_t)
  t_var <- paste0("t_var", 1:ntheta)
  
  # computed parameters labels
  
  rel <- paste0("rel", ind_s)
  con <- paste0("con", ind_s)
  spe <- paste0("spe", ind_s)
  
  # if-else statements to define constraints (parameters equalities)
  fixed_s <- seq(1, nobs, by = nobs/neta)
  fixed_t <- seq(1, neta, by = neta/ntheta)
  if(!second.order.trait){
    fixed_t <- fixed_s
  }
  
  if(equiv.assumption$tau == "equi"){
    ins <- rep(0, nobs)
    las <- rep(1, nobs)
    
  }
  
  if(equiv.assumption$tau == "ess"){
    ins[fixed_s] <- 0
    las <- rep(1, nobs)
    if(scale.invariance$lait0){
      ins <- rep(ins[1:(nobs/neta)], neta)
    }
  }
  
  if(equiv.assumption$tau == "cong"){
    ins[fixed_s] <- 0
    las[fixed_s] <- 1
    if(scale.invariance$lait0){
      ins <- rep(ins[1:(nobs/neta)], neta)
    }
    if(scale.invariance$lait1){
      las <- rep(las[1:(nobs/neta)], neta)
    }
  }
  
  
  if(equiv.assumption$theta == "equi"){
    int <- rep(0, neta)
    lat <- rep(1, neta)
    if(!second.order.trait){
      int <- NULL
      lat <- rep(1, nobs)
    }
  }
  
  if(equiv.assumption$theta == "ess"){
    int[fixed_t] <- 0
    lat <- rep(1, neta)
    if(scale.invariance$lat0){
      int <- rep(int[1:(neta/ntheta)], ntheta)
    }
    if(!second.order.trait){
      stop('When fitting non second order LST models, the equivalence assumption "ess" is not possible')
    }
  }
  
  if(equiv.assumption$theta == "cong"){
    int[fixed_t] <- 0
    lat[fixed_t] <- 1
    if(scale.invariance$lat0){
      int <- rep(int[1:(neta/ntheta)], ntheta)
    }
    if(scale.invariance$lat1){
      lat <- rep(lat[1:(neta/ntheta)], ntheta)
      if(!second.order.trait){
        lat <- rep(lat[1:(nobs/neta)], neta)
      }
    }
    if(!second.order.trait){
      int <- NULL
    }
  }
  
  if(homocedasticity.assumption$error){
    o_red <- rep(o_red[1:(nobs/neta)], neta)
  }
  
  
  if(homocedasticity.assumption$state.red){
    s_red <- rep(s_red[1:(neta/ntheta)], ntheta)
  }
  
  # write equations for the latent states
  eta_syntax <- rep(NA, neta)
  for(i in 1:neta){
    ms <- ((i-1)*(nobs/neta) + 1):(i*(nobs/neta))
    paths <- rep(NA, length(ms) )
    ix <- which(las[1:length(ms)]=="1")
    
    paths[ix] <- paste0(obs[ms[ix]],"@",las[ms[ix]])
    paths[-ix] <- paste0(obs[ms[-ix]]," (",las[ms[-ix]],")")
    
    paths <- paste(paths, collapse = "\n")
    eta_syntax[i] <- paste(eta[i], "BY", paths, ";", sep = " ")
  }
  rm(i, ms, paths, ix)
  
  
  # write equations for the latent trait
  
  theta_syntax <- rep(NA, ntheta)
  for(i in 1:ntheta){
    ms <- ((i-1)*(neta/ntheta) + 1):(i*(neta/ntheta))
    paths <- rep(NA, length(ms) )
    ix <- which(lat[1:length(ms)]=="1")
    
    paths[ix] <- paste0(eta[ms[ix]],"@",lat[ms[ix]])
    paths[-ix] <- paste0(eta[ms[-ix]]," (",lat[ms[-ix]],")")
    
    paths <- paste(paths, collapse = "\n")
    theta_syntax[i] <- paste(theta[i], "BY", paths, ";", sep = " ")
  }
  rm(i, ms, paths, ix)
  
  if(!second.order.trait){
    theta_syntax <- rep(NA, ntheta)
    for(i in 1:ntheta){
      ms <- ((i-1)*(nobs/ntheta) + 1):(i*(nobs/ntheta))
      paths <- rep(NA, length(ms) )
      ix <- which(lat[1:length(ms)]=="1")
      
      paths[ix] <- paste0(obs[ms[ix]],"@",lat[ms[ix]])
      paths[-ix] <- paste0(obs[ms[-ix]]," (",lat[ms[-ix]],")")
      
      paths <- paste(paths, collapse = "\n")
      theta_syntax[i] <- paste(theta[i], "BY", paths, ";", sep = " ")
    }
    rm(i, ms, paths, ix)
  }
  
  # write intercetps for the latent states
  eta_syntax_int <- rep(NA, neta)
  for(i in 1:neta){
    ms <- ((i-1)*(nobs/neta) + 1):(i*(nobs/neta))
    paths <- rep(NA, length(ms) )
    ix <- which(ins[1:length(ms)]=="0")
    
    paths[ix] <- paste0( "[", obs[ms[ix]],"@",ins[ms[ix]], "]", " ;")
    paths[-ix] <- paste0("[", obs[ms[-ix]], "*]", " (",ins[ms[-ix]],")", " ;")
    
    eta_syntax_int[i] <- paste(paths, collapse = "\n") 
  }
  rm(i, ms, paths, ix)
  
  
  # write intercepts for the latent trait
  
  theta_syntax_int <- rep(NA, ntheta)
  for(i in 1:ntheta){
    ms <- ((i-1)*(neta/ntheta) + 1):(i*(neta/ntheta))
    paths <- rep(NA, length(ms) )
    ix <- which(int[1:length(ms)]=="0")
    
    paths[ix] <- paste0("[", eta[ms[ix]],"@",int[ms[ix]], "]", " ;")
    paths[-ix] <- paste0("[", eta[ms[-ix]], "*]", " (",int[ms[-ix]],")", " ;")
    
    theta_syntax_int[i] <- paste(paths, collapse = "\n")
  }
  rm(i, ms, paths, ix)
  
  if(!second.order.trait){
    theta_syntax_int <- NULL
  }
  
  theta_syntax_means <- paste0("[", theta, "*] ;", collapse = "\n")
  
  # write variance components
  
  var_syntax <- paste(paste0(obs, "(", o_red, ") ;", collapse = "\n"),
                      paste0(eta, "(", s_red, ") ;", collapse = "\n"),
                      paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                      sep = "\n")
  
  # define computed parameters
  def_syntax <- paste(paste0("NEW(", rel, ") ;", collapse = "\n" ),
                      paste0("NEW(", con, ") ;", collapse = "\n" ),
                      paste0("NEW(", spe, ") ;", collapse = "\n" ),
                      sep = "\n")
  
  lat_r <- rep(lat, each = nobs/neta)
  if(!second.order.trait){
    lat_r <- lat
  }
  
  total_var <- paste0(lat_r, "^2*",las, "^2*", 
                      rep(t_var, each = neta/ntheta), "+",
                      las, "^2*", rep(s_red, each = nobs/neta), 
                      "+", o_red)
  
  trait_var <- paste0(lat_r, "^2*",las, "^2*", 
                      rep(t_var, each = neta/ntheta))
  
  state_var <- paste0(las, "^2*", rep(s_red, each = nobs/neta))
  
  
  rel_syntax <- paste0(rel, "= (", trait_var, "+", state_var,
                       ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  spe_syntax <- paste0(spe, "= (", state_var,
                       ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  con_syntax <- paste0(con, "= (", trait_var,
                       ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  
  
  complete_syntax <- paste("\nMODEL:",
                           "\n! Latent state variables",
                           paste0(eta_syntax, collapse = "\n"),
                           "\n! Latent trait variables",
                           paste0(theta_syntax, collapse = "\n"),
                           "\n! Intercepts and means",
                           paste0(eta_syntax_int, collapse = "\n"),
                           paste0(theta_syntax_int, collapse = "\n"),
                           paste0(theta_syntax_means, collapse = "\n"),
                           "\n! Variance components",
                           var_syntax,
                           "\nMODEL CONSTRAINT:",
                           def_syntax,
                           "\n! Reliabily",
                           rel_syntax,
                           "\n! Occasion specificity",
                           spe_syntax,
                           "\n! Consistency",
                           con_syntax,
                           sep = "\n")
  return(complete_syntax)
  
}


# 3.0 write.cuts.to.Mplus ----

# Function to write Mplus syntax for the common unique trait state model model given a data.frame

# Several possible constraints are mentioned in page 51 (Hamaker, Schuurman, & Zijlmans, 2017)
# weak factorial invariance across trait and state, 
# weak factorial invariance over time for common state and/or common traits, 
# homocedasticity of variances over time
# Furthermore unique traits/states and common traits/states are assumed 
#     to have mean 0. Although, this constraint is not explicit in the Mplus syntax
#     available in the suplementary material.

# dt is the data.frame. It is supposed to be a longitudinal data in wide format.
#     Only the variables that are measured repeatedly are included in the data set, 
#     ID or sex variables should not be included.
# nstate is the number of common state variables. In other words, it is the number
#     of measurement occasions.
# method.trait defines whether there is a unique trait to each variable (orthogonal methods "om" Default)
#     or whether there is the number of variables minus one unique traits ("m-1").
#     This restriction might be neccesary when "om" don't converge (see suplementary
#     material).
# scale.invariance define whether the measurement model is invariant over time or not.
#     It can be define for intercepts (int) and loadings. 
#     If scale.invariance$lambda is TRUE, scale invariance of the 
#     loadings is assumed for both common states and common traits.
# state.trait.invariance define whether the trait structure is equal to the state 
#     structure or not. This means that the loading of the common trait on the variable 
#     Yit is equal to the loading of the common state on the variable Yit.
# homocedasticity.assumption define whether the variance components are invariant
#     over time. The variances that can be constraint are the variances of the unique
#     states (error), the variances of the common states (cs.red), and the variances
#     of the unique traits (ut.red).

write.cuts.to.Mplus <- function(dt, nstate,
                                method.trait = "om",
                                scale.invariance = list(int = FALSE, lambda = FALSE),
                                state.trait.invariance = FALSE,
                                fixed.method.loadings = TRUE,
                                homocedasticity.assumption = list(error = FALSE, cs.red = FALSE, ut.red = FALSE),
                                fixed.means = list(cs = FALSE, ut = FALSE, ct = FALSE) ){
  
  nobs <- dim(dt)[2]
  
  #labels
  ind <- paste0(rep(1:(nobs/nstate) ,nstate),
                rep(1:nstate, each=(nobs/nstate)))
  
  
  #lambdas for states and traits
  lcs <- paste0("lcs", ind) #lambda common states
  lct <- paste0("lct", ind) #lambda common trait
  lut <- paste0("lut", ind) #lambda unique traits
  # intercepts for observed variables
  ins <- paste0("ins", ind)
  
  
  # Common states names
  CS <- paste0("CS", 1:nstate) # Common state
  
  # Unique and Common traits names
  
  CT <- "CT"
  
  UT <- paste0("UT", 1:(nobs/nstate))
  
  # observed variables names
  obs <- names(dt)
  
  # variances and residuals' labels
  
  o_red <- paste0("err", ind)
  cs_var <- paste0("cs_var", 1:nstate)
  ut_var <- paste0("ut_var", 1:(nobs/nstate))
  ct_var <- "ct_var"
  
  
  # computed parameters labels
  # In the cuts model there are not variances components defined such as
  # reliability or consistency. However, based on the definition of these variance
  # components in the LST, similar variances components could be defined.
  
  #rel <- paste0("rel", ind)
  #pre <- paste0("pre", ind)
  #upre <- paste0("upre", ind)
  #spe <- paste0("spe", ind)
  
  
  # if-else statements to define constraints (parameters equalities)
  fixed_s <- seq(1, nobs, by = nobs/nstate)
  
  lcs[fixed_s] <- 1
  lct[fixed_s] <- 1
  lut[1:(nobs/nstate)]<-1
  
  m.obs <- obs
  
  if(method.trait == "m-1"){
    UT <- UT[-1]
    ut_var <- ut_var[-1]
    m.obs <- obs[-(seq(1, nobs, by = nobs/nstate))]
    lut <- lut[-(seq(1, nobs, by = nobs/nstate))]
  }
  
  # define invariance of intercepts
  
  if(scale.invariance$int){
    ins <- rep(ins[1:(nobs/nstate)], nstate)
  }
  
  # define invariance over time of lambdas for states and traits
  
  if(scale.invariance$lambda){
    lcs <- rep(lcs[1:(nobs/nstate)], nstate)
    lct <- rep(lct[1:(nobs/nstate)], nstate)
  }
  
  # define state trait invariance
  
  if(state.trait.invariance){
    lcs <- lct
  }
  
  # fixed mehtod loadings
  
  if(fixed.method.loadings){
    lut <- rep(1, length(m.obs))
  }
  
  # homogeneity of CS, US and UT variances
  
  if(homocedasticity.assumption$error){
    o_red <- rep(o_red[1:(nobs/nstate)], nstate)
  }
  
  if(homocedasticity.assumption$cs.red){
    cs_var <- rep(cs_var[1], nstate)
  }
  
  if(homocedasticity.assumption$ut.red){
    ut_var <- rep(ut_var[1], nobs/nstate)
  }
  
  
  # write equations for the common states
  cs_syntax <- rep(NA, nstate)
  for(i in 1:nstate){
    ms <- ((i-1)*(nobs/nstate) + 1):(i*(nobs/nstate))
    paths <- rep(NA, length(ms) )
    ix <- which(lcs[1:length(ms)]=="1")
    
    paths[ix] <- paste0(obs[ms[ix]],"@",lcs[ms[ix]])
    paths[-ix] <- paste0(obs[ms[-ix]]," (",lcs[ms[-ix]],")")
    
    paths <- paste(paths, collapse = "\n")
    cs_syntax[i] <- paste(CS[i], "BY", paths, ";", sep = " ")
  }
  rm(i, ms, paths, ix)
  
  
  # write equations for the common and unique trait
  
  ct_syntax <- NA
  paths <- rep(NA, nobs )
  ix <- which(lct=="1")
  
  paths[ix] <- paste0(obs[ix],"@",lct[ix])
  paths[-ix] <- paste0(obs[-ix]," (",lct[-ix],")")
  
  paths <- paste(paths, collapse = "\n")
  ct_syntax <- paste(CT, "BY", paths, ";", sep = " ")
  rm(paths, ix)
  
  ut_syntax <- rep(NA, length(UT))
  
  for(i in 1:length(UT)){
    paths <- rep(NA, nstate )
    ms <- seq(i, length(m.obs), by = length(m.obs)/nstate)
    
    ix <- which(lut[ms]=="1")
    
    paths[ix] <- paste0(m.obs[ms[ix]],"@",lut[ms[ix]])
    paths[-ix] <- paste0(m.obs[ms[-ix]]," (",lut[ms[-ix]],")")
    
    paths <- paste(paths, collapse = "\n")
    ut_syntax[i] <- paste(UT[i], "BY", paths, ";", sep = " ")
  }
  rm(i, paths)
  
  
  
  
  
  # write intercetps of observed variables
  intercepts_syntax <- paste0("[", obs, "*]", " (",ins,")", " ;", collapse = "\n")
  
  
  
  # fix means of latent variables to 0. This is optional, it is not explicitly
  # in the example syntaxes of the suplementary material.
  
  cs_means_syntax <- ""
  ut_means_syntax <- ""
  ct_means_syntax <- ""
  
  if(fixed.means$cs){cs_means_syntax <- paste0("[", CS, "@0] ;", collapse = "\n")}
  
  if(fixed.means$ut){ut_means_syntax <- paste0("[", UT, "@0] ;", collapse = "\n")}
  
  if(fixed.means$ct){ct_means_syntax <- paste0("[", CT, "@0] ;", collapse = "\n")}
  
  
  # write 0 correlations among occasions  and traits
  latent.v <- c(CS, UT, CT)
  
  cor_syntax <- rep(NA, length(latent.v)-1)
  for(i in 1:(length(latent.v)-1)){
    
    cor_syntax[i] <- paste0(latent.v[i], " WITH ",
                            paste0(latent.v[-(1:i)], "@0", collapse = "\n"), " ;")
  }
  rm(i)
  
  
  # write variance components
  
  var_syntax <- paste(paste0(obs, " (", o_red, ") ;", collapse = "\n"),
                      paste0(CS, " (", cs_var, ") ;", collapse = "\n"),
                      paste0(UT, " (", ut_var, ") ;", collapse = "\n"),
                      paste0(CT, " (", ct_var, ") ;", collapse =  "\n"),
                      sep = "\n")
  
  # define computed parameters
  #def_syntax <- paste(paste0("NEW(", rel, ") ;", collapse = "\n" ),
  #                    paste0("NEW(", con, ") ;", collapse = "\n" ),
  #                    paste0("NEW(", spe, ") ;", collapse = "\n" ),
  #                    sep = "\n")
  
  #lat_r <- rep(lat, each = nobs/neta)
  
  #total_var <- paste0(lat_r, "^2*",las, "^2*", 
  #                     rep(t_var, each = neta/ntheta), "+",
  #                     las, "^2*", rep(s_red, each = nobs/neta), 
  #                    "+", o_red)
  
  #trait_var <- paste0(lat_r, "^2*",las, "^2*", 
  #                    rep(t_var, each = neta/ntheta))
  
  #state_var <- paste0(las, "^2*", rep(s_red, each = nobs/neta))
  
  
  #rel_syntax <- paste0(rel, "= (", trait_var, "+", state_var,
  #                     ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  #spe_syntax <- paste0(spe, "= (", state_var,
  #                     ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  #con_syntax <- paste0(con, "= (", trait_var,
  #                     ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  
  
  complete_syntax <- paste("\nMODEL:",
                           "\n! common states",
                           paste0(cs_syntax, collapse = "\n"),
                           "\n! common trait",
                           paste0(ct_syntax, collapse = "\n"),
                           "\n! unique trait",
                           paste0(ut_syntax, collapse = "\n"),
                           "\n! Intercepts and means",
                           paste0(intercepts_syntax, collapse = "\n"),
                           paste0(cs_means_syntax, collapse = "\n"),
                           paste0(ut_means_syntax, collapse = "\n"),
                           paste0(ct_means_syntax, collapse = "\n"),
                           "\n! Set 0 covariances",
                           paste0(cor_syntax, collapse = "\n"),
                           "\n! Variance components",
                           var_syntax,
                           #"\nMODEL CONSTRAINT:",
                           #def_syntax,
                           #"\n! Reliabily",
                           #rel_syntax,
                           #"\n! Occasion specificity",
                           #spe_syntax,
                           #"\n! Consistency",
                           #con_syntax,
                           sep = "\n")
  
  return(complete_syntax)
  
}

# 4.0 write.tso.to.Mplus ----

# Function to write Mplus syntax for the extended trait-state-occasion model given a data.frame

# This function writes a syntax of the extended TSO model based on the figure 3B
# of Eid et al. (2017)

# Restrictions homogeneity of autoregressive effects, homogeneity of trait loadings (equivalence assumption),
# homogeneity of state residual variances, homogenity of error variances.

# dt is the data.frame. It is supposed to be a longitudinal data in wide format.
#     Only the variables that are measured repeatedly are included in the data set, 
#     ID or sex variables should not be included.
# nocc is the number of latent occasion variables or the number of measurement 
#     occasions.
# figure define if the TSO model to fit is the one of figure 3a or the one of
#     figure 3b in Eid et al. (2017) 
# equiv.assumption define if the model assumes tau-equivalence, equivalence or 
#     tau-congenericity in relation to the latent occasions (occ) and the latent traits (theta).
#     Use "equi" or "cong". The essential tau-equivalence assumption (as in LST)
#     was ommited because given the way we are modelling the TSO model, it is unclear 
#     whether the intercepts are associated to the trait indicators or the latent 
#     occasions. Also, in this case, these assumptions are just modifying whether
#     the loadings are fixed to 1 or not, they do not affect the intercepts as they
#      do in LST theory.
# scale.invariance define whether the measurement model is invariant over time or not.
#     It can be define for the intercepts and for the loadings of the latent occasions.
#     it does not apply for the loadings of the traits.
# homocedasticity.assumption define whether the variances of the latent occasions (occ.red)
#     and the measurement errors are invariant over time.
# autoregressive.homogeneity is a logical value that defines whether the auto-
#     regressive effects are constant over time or not.

write.tso.to.Mplus <- function(dt, nocc, figure = c("3a", "3b"),
                               equiv.assumption = list(occ = "cong", theta = "equi"),
                               scale.invariance = list(int = FALSE, lambda = FALSE),
                               homocedasticity.assumption = list(error = FALSE, occ.red = FALSE),
                               autoregressive.homogeneity = FALSE){
  
  
  nobs <- dim(dt)[2]
  
  #labels
  ind <- paste0(rep(1:(nobs/nocc) ,nocc),
                rep(1:nocc, each=(nobs/nocc)))
  #lambdas for occasions and traits
  las <- paste0("las", ind)
  lat <- paste0("lat", ind)
  if(figure == "3a"){
    lat <- paste0("lat", 1:nocc)
  }
  
  # intercepts for observed variables
  ins <- paste0("ins", ind)
  
  
  #latent occasion variables' names
  occ <- paste0("occ", 1:nocc)
  if(figure == "3a"){
    eta <- paste0("eta", 1:nocc)
  }
  
  # trait and trait indicators names
  if(figure == "3b"){
    theta <- paste0("theta", 1:(nobs/nocc))}else{
      theta <- "theta"
    }
  
  
  # observed variables names
  obs <- names(dt)
  
  # variances and residuals' labels
  
  o_red <- paste0("err", ind)
  occ_var <- paste0("oc_var", 1:nocc)
  
  if(figure == "3b"){
    t_var <- paste0("t_var", 1:(nobs/nocc))}else{
      t_var <- "t_var"
    }
  
  # autoregressive efects labels
  
  beta <- paste0("beta", 1:(nocc-1))
  
  # computed parameters labels
  # In the TSO the consistency is divided in predictability by trait and 
  #     unpredictability by trait. These components can be computed with Mplus
  #     or afterwards. If, it is needed I can make the syntax to compute them 
  #     with Mplus, for now it does not seem necessary
  
  #rel <- paste0("rel", ind)
  #pre <- paste0("pre", ind)
  #upre <- paste0("upre", ind)
  #spe <- paste0("spe", ind)
  
  
  # if-else statements to define constraints (parameters equalities)
  fixed_s <- seq(1, nobs, by = nobs/nocc)
  
  if(equiv.assumption$occ == "equi"){
    las <- rep(1, nobs)
    
  }
  
  if(equiv.assumption$occ == "cong"){
    las[fixed_s] <- 1
    if(scale.invariance$lambda){
      las <- rep(las[1:(nobs/nocc)], nocc)
    }
  }
  
  
  if(equiv.assumption$theta == "equi"){
    lat <- rep(1, length(lat))
  }
  
  
  if(equiv.assumption$theta == "cong"){
    if(figure == "3b"){
      lat[1:(nobs/nocc)] <- "1"
    }else{lat[1] <- "1"}
  }
  
  # define invariance of intercepts
  
  if(scale.invariance$int){
    if(figure == "3b"){
      ins <- rep(ins[1:(nobs/nocc)], nocc)
    }else{
      ins <- rep(ins[1], nobs)
    }
    
  }
  
  # homogeneity of autoregresive effects
  
  if(autoregressive.homogeneity){
    beta <- rep(beta[1], nocc-1)
  }
  
  # homogeneity of errors and occasion variances
  
  if(homocedasticity.assumption$error){
    o_red <- rep(o_red[1:(nobs/nocc)], nocc)
  }
  
  if(homocedasticity.assumption$occ.red){
    occ_var <- rep(occ_var[1], nocc)
  }
  
  
  # write equations for the latent occasions and latent states
  if(figure == "3b"){
    occ_syntax <- rep(NA, nocc)
    for(i in 1:nocc){
      ms <- ((i-1)*(nobs/nocc) + 1):(i*(nobs/nocc))
      paths <- rep(NA, length(ms) )
      ix <- which(las[1:length(ms)]=="1")
      
      paths[ix] <- paste0(obs[ms[ix]],"@",las[ms[ix]])
      paths[-ix] <- paste0(obs[ms[-ix]]," (",las[ms[-ix]],")")
      
      paths <- paste(paths, collapse = "\n")
      occ_syntax[i] <- paste(occ[i], "BY", paths, ";", sep = " ")
    }
    rm(i, ms, paths, ix)
  }
  
  if(figure == "3a"){
    occ_syntax_a <- rep(NA, nocc)
    for(i in 1:nocc){
      ms <- ((i-1)*(nobs/nocc) + 1):(i*(nobs/nocc))
      paths <- rep(NA, length(ms) )
      ix <- which(las[1:length(ms)]=="1")
      
      paths[ix] <- paste0(obs[ms[ix]],"@",las[ms[ix]])
      paths[-ix] <- paste0(obs[ms[-ix]]," (",las[ms[-ix]],")")
      
      paths <- paste(paths, collapse = "\n")
      occ_syntax_a[i] <- paste(eta[i], "BY", paths, ";", sep = " ")
    }
    rm(i, ms, paths, ix)
    
    occ_syntax_b <- paste(occ,"BY", paste0(eta, "@1"), ";", sep = " ")
    
    occ_syntax <- c(occ_syntax_a, occ_syntax_b)
  }
  
  
  
  # write equations for the latent trait or trait indicators
  if(figure == "3a"){
    trait_syntax <- NA
    paths <- rep(NA, nocc )
    ix <- which(lat=="1")
    
    paths[ix] <- paste0(eta[ix],"@",lat[ix])
    paths[-ix] <- paste0(eta[-ix]," (",lat[-ix],")")
    
    paths <- paste(paths, collapse = "\n")
    trait_syntax <- paste(theta, "BY", paths, ";", sep = " ")
    rm(paths, ix)
  }
  if(figure == "3b"){
    trait_syntax <- rep(NA, nobs/nocc)
    
    for(i in 1:(nobs/nocc)){
      paths <- rep(NA, nocc )
      ix <- which(lat[seq(i, nobs, by = nobs/nocc)]=="1")
      
      paths[ix] <- paste0(obs[seq(i, nobs, by = nobs/nocc)[ix]],"@",lat[seq(i, nobs, by = nobs/nocc)[ix]])
      paths[-ix] <- paste0(obs[seq(i, nobs, by = nobs/nocc)[-ix]]," (",lat[seq(i, nobs, by = nobs/nocc)[-ix]],")")
      
      paths <- paste(paths, collapse = "\n")
      trait_syntax[i] <- paste(theta[i], "BY", paths, ";", sep = " ")
    }
    rm(i, paths, ix)
    
  }
  
  
  
  # write intercetps of observed variables
  intercepts_syntax <- paste0("[", obs, "*]", " (",ins,")", " ;", collapse = "\n")
  
  
  
  # write means of trait indicators and occasions
  
  occ_mean_syntax <- paste0("[", occ, "@0] ;", collapse = "\n" )
  
  # When the means of the traits are included,the model is not identified.
  #trait_mean_syntax <- paste0("[", theta, "* ] ;", collapse = "\n")
  
  # write 0 correlations among occasions  and traits
  latent.v <- c(occ, theta)
  
  cor_syntax <- rep(NA, nocc)
  for(i in 1:nocc){
    
    cor_syntax[i] <- paste0(latent.v[i], " WITH ",
                            paste0(latent.v[-(1:i)], "@0", collapse = "\n"), " ;")
  }
  rm(i)
  
  
  # write autoregressive effects
  
  beta_syntax <- paste0(occ[2:nocc], " ON ", occ[1:(nocc-1)], " (", beta, ") ;", 
                        collapse = "\n")
  
  
  # write variance components
  if(figure == "3a"){
    var_syntax <- paste(paste0(obs, " (", o_red, ") ;", collapse = "\n"),
                        paste0(occ, " (", occ_var, ") ;", collapse = "\n"),
                        paste0(eta, "@0 ;", collapse = "\n"),
                        paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                        sep = "\n")
  }
  
  
  if(figure == "3b"){
    var_syntax <- paste(paste0(obs, " (", o_red, ") ;", collapse = "\n"),
                        paste0(occ, " (", occ_var, ") ;", collapse = "\n"),
                        paste0(theta, "(", t_var, ") ;", collapse =  "\n"),
                        sep = "\n")
  }
  
  # define computed parameters
  #def_syntax <- paste(paste0("NEW(", rel, ") ;", collapse = "\n" ),
  #                    paste0("NEW(", con, ") ;", collapse = "\n" ),
  #                    paste0("NEW(", spe, ") ;", collapse = "\n" ),
  #                    sep = "\n")
  
  #lat_r <- rep(lat, each = nobs/neta)
  
  #total_var <- paste0(lat_r, "^2*",las, "^2*", 
  #                     rep(t_var, each = neta/ntheta), "+",
  #                     las, "^2*", rep(s_red, each = nobs/neta), 
  #                    "+", o_red)
  
  #trait_var <- paste0(lat_r, "^2*",las, "^2*", 
  #                    rep(t_var, each = neta/ntheta))
  
  #state_var <- paste0(las, "^2*", rep(s_red, each = nobs/neta))
  
  
  #rel_syntax <- paste0(rel, "= (", trait_var, "+", state_var,
  #                     ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  #spe_syntax <- paste0(spe, "= (", state_var,
  #                     ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  #con_syntax <- paste0(con, "= (", trait_var,
  #                     ")/", "\n(", total_var, ") ;", collapse = "\n")
  
  
  complete_syntax <- paste("\nMODEL:",
                           "\n! Latent occasion variables",
                           paste0(occ_syntax, collapse = "\n"),
                           "\n! Latent trait variable or trait indicators",
                           paste0(trait_syntax, collapse = "\n"),
                           "\n! Intercepts and means",
                           paste0(intercepts_syntax, collapse = "\n"),
                           paste0(occ_mean_syntax, collapse = "\n"),
                           #paste0(trait_mean_syntax, collapse = "\n"),
                           "\n! Set 0 correlations",
                           paste0(cor_syntax, collapse = "\n"),
                           "\n! autoregressive effects",
                           paste0(beta_syntax, collapse = "\n"),
                           "\n! Variance components",
                           var_syntax,
                           #"\nMODEL CONSTRAINT:",
                           #def_syntax,
                           #"\n! Reliabily",
                           #rel_syntax,
                           #"\n! Occasion specificity",
                           #spe_syntax,
                           #"\n! Consistency",
                           #con_syntax,
                           sep = "\n")
  
  return(complete_syntax)
  
}

# 5.0 write.mlmsst.to.Mplus ----

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

# 6.0 write.mlcuts.to.Mplus ----

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

# 7.0 write.mltso.to.Mplus ----

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


