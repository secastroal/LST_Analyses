# Function to personalize some analysis options in the Mplus syntax ----

write.Mplus.options <- function(usevariables, cluster = NULL, analysis_type = c("GENERAL", "TWOLEVEL"),
                                  estimator = c("ML", "MLR", "BAYES"), iterations, h1iterations = 50000, 
                                  processors = 4, chains = 4, thin = 10){
  usevariables_syntax <- paste("USEVAR", "=", paste(usevariables, collapse = "\n"), ";", sep = " ")
  if(!is.null(cluster)){
    cluster_syntax <- paste("CLUSTER", "=", cluster, ";", sep = " ")
  }else{
    cluster_syntax <- NULL
  }
  analysis_type_syntax <- paste("TYPE", "=", analysis_type, ";", sep = " ")
  estimator_syntax <- paste("ESTIMATOR", "=", estimator, ";", sep = " ")
  if(estimator == "BAYES"){
    iterations_syntax <- paste("BITERATIONS", "=", paste0("(", iterations, ")"), ";", sep = " " )
    processors_syntax <- paste("PROCESSORS", "=", processors, ";", sep = " ")
    chains_syntax <- paste("CHAINS", "=", chains, ";", sep = " ")
    thin_syntax <- paste("THIN", "=", thin, ";", sep = " ")
    
    analysis_syntax <- paste(usevariables_syntax,
                             cluster_syntax,
                             "\nANALYSIS:",
                             analysis_type_syntax,
                             estimator_syntax,
                             iterations_syntax, 
                             processors_syntax,
                             chains_syntax,
                             thin_syntax,
                             sep = "\n")
    
  }
  if(estimator != "BAYES"){
    iterations_syntax <- paste("ITERATIONS", "=", iterations, ";", sep = " ")
    h1iterations_syntax <- paste("H1ITERATIONS", "=", h1iterations, ";", sep = " ")
    
    analysis_syntax <- paste(usevariables_syntax,
                             cluster_syntax,
                             "\nANALYSIS:",
                             analysis_type_syntax,
                             estimator_syntax,
                             iterations_syntax, 
                             h1iterations_syntax,
                             sep = "\n")
  }
  
  
  cat(analysis_syntax)
  return(analysis_syntax)
}

