# Function to verify whether the model estimation in Mplus was terminated normally without any warning or error message.
# The input is an mplus.model object returned by the function MplusAutomation::readModels. 

check.mplus <- function(mplus.model){
  modelstring <- "THE MODEL ESTIMATION TERMINATED NORMALLY"
  file.path <- grep(as.character(mplus.model$summaries["Filename"]), list.files(pattern = ".out", recursive = TRUE), 
                    value = TRUE)
  if(any(grep(modelstring, readLines(file.path)))){
    if(length(mplus.model$warnings) == 0 && length(mplus.model$errors) == 0){
      return("Ok")
    }else{return("Warnings and Errors")}
  }else{
    return("Non-Convergence")}
}

