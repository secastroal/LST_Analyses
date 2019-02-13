# Function to personalize some analysis options in the Mplus syntax ----

write.Mplus.options <- function(usevariables, cluster = NULL, analysis_type = c("GENERAL", "TWOLEVEL"),
                                  estimator = c("ML", "MLR", "BAYES"), iterations, h1iterations = NULL, 
                                  processors = NULL){
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
    h1iterations_syntax <- NULL
    processors_syntax <- paste("PROCESSORS", "=", processors, ";", sep = " ")
  }
  if(estimator != "BAYES"){
    iterations_syntax <- paste("ITERATIONS", "=", iterations, ";", sep = " ")
    h1iterations_syntax <- paste("H1ITERATIONS", "=", h1iterations, ";", sep = " ")
    processors_syntax <- NULL
  }
  
  analysis_syntax <- paste(usevariables_syntax,
                           cluster_syntax,
                           "\nANALYSIS:",
                           analysis_type_syntax,
                           estimator_syntax,
                           iterations_syntax, 
                           h1iterations_syntax,
                           processors_syntax,
                           sep = "\n")
  cat(analysis_syntax)
  return(analysis_syntax)
}

