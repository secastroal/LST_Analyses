#' Run Mplus Models
#'
#' This function runs a group of Mplus models (.inp files) located within a
#' single directory or nested within subdirectories.
#'
#' @param target the directory containing Mplus input files (.inp) to run
#'   OR the single .inp file to be run. May be a full path, relative path,
#'   or a filename within the working directory.
#'   Defaults to the current working directory. Example: \dQuote{C:/Users/Michael/Mplus Runs}
#' @param recursive optional. If \code{TRUE}, run all models nested in subdirectories
#'   within \code{directory}. Defaults to \code{FALSE}. Not relevant if \code{target} is a single file.
#' @param filefilter a Perl regular expression (PCRE-compatible) specifying particular input
#'   files to be run within \code{directory}. See \code{regex} or \url{http://www.pcre.org/pcre.txt}
#'   for details about regular expression syntax. Not relevant if \code{target} is a single file.
#' @param showOutput optional. If \code{TRUE}, show estimation output (TECH8)
#'   in the R console. Note that if run within Rgui, output will display within R,
#'   but if run via Rterm, a separate window will appear during estimation.
#' @param replaceOutfile optional. Currently supports three settings: \dQuote{always}, which
#'   runs all models, regardless of whether an output file for the model exists; \dQuote{never},
#'   which does not run any model that has an existing output file; and \dQuote{modifiedDate}, which
#'   only runs a model if the modified date for the input file is more recent than
#'   the output file modified date (implying there have been updates to the model).
#' @param logFile optional. If non-null, specifies a file (and optionally, directory)
#'   that records the settings passed into the function and the models run (or skipped)
#'   during the run.
#' @param Mplus_command optional. N.B.: No need to pass this parameter for most users (has intelligent
#'   defaults). Allows the user to specify the name/path of the Mplus executable to be used for
#'   running models. This covers situations where Mplus is not in the system's path,
#'   or where one wants to test different versions of the Mplus program.
#' @param killOnFail optional. Windows only for now. If \code{TRUE}, kill all processes named mplus.exe when
#'   \code{runModels} does not terminate normally. Defaults to \code{TRUE}.
#' @param local_tmpdir optional. Linux/Mac for now. If \code{TRUE}, set the TMPDIR environment variable to the
#'   location of the .inp file prior to execution. This is useful in Monte Carlo studies where many instances of Mplus
#'   may run in parallel and we wish to avoid collisions in temporary files among processes.
#'
#' @return None. Function is used for its side effects (running models).
#' @author Michael Hallquist
#' @seealso \code{\link{runModels_Interactive}}
#' @keywords interface
#' @export
#' @examples
#' \dontrun{
#'   runModels("C:/Users/Michael/Mplus Runs", recursive=TRUE, showOutput=TRUE,
#'     replaceOutfile="modifiedDate", logFile="MH_RunLog.txt",
#'     Mplus_command="C:\\Users\\Michael\\Mplus Install\\Mplus51.exe")
#' }
runModels_2 <- function(target=getwd(), recursive=FALSE, filefilter = NULL, showOutput=FALSE,
                      replaceOutfile="always", logFile="Mplus Run Models.log", Mplus_command="Mplus", killOnFail=TRUE, local_tmpdir=FALSE,
                      timeout = 0) {
  
  stopifnot(replaceOutfile %in% c("always", "never", "modifiedDate"))
  
  #TODO: would be good to come back and make this more versatile, supporting a vector target
  if (length(target) > 1L) { stop("target for runModels must be a single file or single directory")}
  
  #retain working directory and reset at end of run
  #need to set here to ensure that logTarget initialization below is within target directory, not getwd()
  curdir <- getwd()
  
  #check whether target is a directory or a single file
  if (grepl(".*\\.inp?$", target, perl=TRUE)) { #tolerate .in and .inp files
    directory <- dirname(target)
    filelist <- basename(target)
    
    if (!is.null(filefilter)) { warning("Using runModels() with a single .inp target ignores filefilter") }
    if (!file.exists(target)) { stop("runModels cannot locate target file: ", target) }
    setwd(directory)
    
    #look for .out file of the same name to handle skipping existing files downstream.
    if (file.exists(outtest <- sub("\\.inp?$", ".out", filelist, perl=TRUE))) { filelist <- c(filelist, outtest)}
    
  } else {
    ## remove trailing slash, which generates file.exists error on windows: https://bugs.r-project.org/bugzilla/show_bug.cgi?id=14721
    directory <- sub("(\\\\|/)?$", "", target, perl=TRUE)
    ## however, trailing slash is needed if at root of a drive on Windows
    if (.Platform$OS.type == "windows" && grepl("^[a-zA-Z]:$", directory)) {
      directory <- paste0(directory, "/")
    }
    
    if (!file.exists(directory)) { stop("runModels cannot change to directory: ", directory) }
    setwd(directory)
    
    #list files in the current directory
    filelist <- list.files(recursive=recursive, pattern=filefilter)
  }
  
  normalComplete <- FALSE
  
  #if log file requested, then open file connection for writing
  if (!is.null(logFile)) {
    logTarget <- file(description = logFile, open = "wt", blocking = TRUE)
    writeLines(c(paste("------Begin Mplus Model Run: ", format(Sys.time(), "%d%b%Y %H:%M:%S"), "------", sep=""),
                 paste("Target directory: ", directory, sep=""),
                 "Run options:",
                 paste("\tRecursive (run models in subdirectories):", as.character(recursive)),
                 paste("\tShow output on console:", as.character(showOutput)),
                 paste("\tReplace existing outfile:", replaceOutfile),
                 "------"
    ), con=logTarget)
    #need to flush after each writeLines to keep the text file current.
    flush(logTarget)
  }
  
  isLogOpen <- function() {
    #if null is passed as the log file, then it is by definition not open (non-existent)
    if (is.null(logFile)) return(FALSE)
    
    connections <- data.frame(showConnections(all = FALSE))
    if (length(grep(MplusAutomation:::splitFilePath(logFile)$filename, connections$description, ignore.case=TRUE)) > 0) return(TRUE)
    else return(FALSE)
  }
  
  #if the function gets interrupted (e.g., the user presses escape), kill the Mplus process (doesn't happen automatically).
  exitRun <- function() {
    deleteOnKill <- TRUE #whether to delete unfinished output
    
    if (normalComplete == FALSE && isLogOpen()) {
      writeLines("Run terminated abnormally", logTarget)
      flush(logTarget)
    }
    
    #create a data frame consisting of the process names and pids
    #uses str split on the output of wmic to extract name and pid columns
    #depends on windows tools
    if (.Platform$OS.type == "windows" && normalComplete == FALSE && killOnFail == TRUE) {
      processList <- ldply(strsplit(shell("wmic process get caption, processid", intern=TRUE), split="\\s+", perl=TRUE),
                           function(element) {
                             return(data.frame(procname=element[1], pid=element[2], stringsAsFactors = FALSE))
                           })
      
      if(length(grep("mplus.exe", processList$procname, ignore.case=TRUE)) > 0) {
        if(isLogOpen()) {
          writeLines("Killing wayward Mplus processes", logTarget)
          flush(logTarget)
        }
        shell("taskkill /f /im mplus.exe")
        
        #if the process is currently running and we kill it, then the output and gph files will be incomplete.
        #in general, it would be good to delete these.
        if(deleteOnKill == TRUE) {
          noExtension <- substr(absFilename, length(absFilename) - 4, length(absFilename))
          outDelete <- paste(noExtension, ".out", sep="")
          gphDelete <- paste(noExtension, ".gph", sep="")
          if (file.exists(outDelete)) {
            unlink(outDelete)
            if(isLogOpen()) {
              writeLines(paste("Deleting unfinished output file:", outDelete), logTarget)
              flush(logTarget)
            }
          }
          if (file.exists(gphDelete)) {
            unlink(gphDelete)
            if(isLogOpen()) {
              writeLines(paste("Deleting unfinished graph file:", gphDelete), logTarget)
              flush(logTarget)
            }
          }
        }
      }
    }
    
    #close logfile
    if (isLogOpen()) { close(logTarget) }
    
    #reset working directory
    setwd(curdir)
  } #end exitRun definition
  
  on.exit(exitRun())
  
  #select only .inp files using grep
  inpfiles <- filelist[grep(".*\\.inp?$", filelist, ignore.case=TRUE)] #tolerate .in and .inp files
  outfiles <- filelist[grep(".*\\.out$", filelist, ignore.case=TRUE)]
  
  if(length(inpfiles) < 1) stop("No Mplus input files detected in the target directory: ", directory)
  
  dropOutExtensions <- sapply(outfiles, function(x) {
    if (nchar(x) >= 4) return(tolower(substr(x, 1, (nchar(x)-4))))
  })
  
  for (i in 1:length(inpfiles)) {
    
    if (!replaceOutfile == "always") {
      #if there is a match in the outfiles for this input file, then decide whether to skip
      if (tolower(sub("\\.inp?$", "", inpfiles[i], perl=TRUE)) %in% dropOutExtensions) {
        
        if (replaceOutfile == "modifiedDate") {
          #if check date is true, then the output file must exist and it must be
          #older than the input file to re-run
          inpmtime <- file.info(inpfiles[i])$mtime
          
          #need to locate the exact outfile match
          matchPos <- grep(tolower(substr(inpfiles[i], 1, (nchar(inpfiles[i]) - 4))), dropOutExtensions)
          if (length(matchPos) < 1) warning("Could not locate matching outfile")
          outmtime <- file.info(outfiles[matchPos[1]])$mtime
          
          if (inpmtime <= outmtime) {
            if (isLogOpen()) {
              writeLines(paste("Skipping model because output file is newer than input file:", inpfiles[i]), logTarget)
              flush(logTarget)
            }
            next
          }
        }
        else if (replaceOutfile == "never"){
          if (isLogOpen()) {
            writeLines(paste("Skipping model because output file already exists for:", inpfiles[i]), logTarget)
            flush(logTarget)
          }
          next
        }
      }
    }
    
    #split input file into one element per directory (e.g., "dir1/dir2/mytest.inp" becomes a 3-element vector
    #code adapted from R.utils filePath command
    inputSplit <- MplusAutomation:::splitFilePath(inpfiles[i])
    if (is.na(inputSplit$directory)) dirtocd <- directory
    else dirtocd <- file.path(directory, inputSplit$directory)
    
    #the absolute path to the file combines the passed-in directory with the subdirectories
    #identified in the case of recursive=T
    
    absFilename <- file.path(directory, inpfiles[i])
    
    #UPDATE 21Oct2011: Since mplus has been released for linux, don't default to wine.
    if (.Platform$OS.type == "unix" && Mplus_command == "Mplus") {
      if (Sys.info()["sysname"] == "Darwin") Mplus_command <- "/Applications/Mplus/mplus"
      else Mplus_command <- "mplus" #linux is case sensitive
    }
    
    #navigate to working directory in DOS using cd command so that Mplus finds the appropriate files (support rel paths)
    #switched over to use relative filename because of problems in Mplus via Wine misinterpreting absolute paths due to forward slashes.
    #25Jul2012: Quote Mplus_command in case it's in a path with spaces.
    
    #create batch files to be able to set a time limit
    #cm file
    #cm.syntax <- paste0("@ECHO OFF\n", paste("\"", Mplus_command, "\" \"", inputSplit$filename, "\"", sep=""), 
     #                    "\ntaskkill /im Mplus.exe /f\ntaskkill /im cmd.exe /f\nEXIT")
    #writeLines(cm.syntax, con = "cm.bat")
    
    #run file
    #run.syntax <- paste0("@ECHO OFF\ntimeout /t ", timeout, 
     #                   "\ntaskkill /im Mplus.exe /f\ntaskkill /im cmd.exe /f\nEXIT")
    #writeLines(run.syntax, con = "run.bat")
    
    #process file
    #process.syntax <- "@ECHO OFF \nstart /b cm.bat \nstart cmd.exe /c run.bat \nEXIT"
    #writeLines(process.syntax, con = "process.bat")
    
    #rm(cm.syntax, run.syntax, process.syntax)
    
    #command <- paste("cd \"", dirtocd, "\" && \"process.bat\"", sep="")
    command <- paste("cd \"", dirtocd, "\" && \"", Mplus_command, 
                     "\" \"", inputSplit$filename, "\"", sep = "")
    
    #allow for divergence if the package is being run in Linux (Mplus via wine)
    if (.Platform$OS.type == "windows") {
      #Given that Mplus is a Windows program, should generally automatically remap / as \ for Windows
      #remap forward slashes to backslashes
      command <- chartr("/", "\\", command)
      
      #code swiped from shell because shell didn't support suppressing output
      shellcommand <- Sys.getenv("COMSPEC")
      flag <- "/c"
      
      #assemble full command from the shell call, flag, and command
      command <- paste(shellcommand, flag, command)
    }
    else if (.Platform$OS.type == "unix") {
      #allow for unix-specific customization here
    }
    
    if (isLogOpen()) {
      writeLines(paste("Currently running model:", inputSplit$filename), logTarget)
      flush(logTarget)
    }
    
    #run the model
    cat("\nRunning model:", inputSplit$filename, "\n")
    cat("System command:", command, "\n")
    
    #unix system command does not have show.output.on.console or invisible parameters
    if (.Platform$OS.type == "windows")	{
      system(command, show.output.on.console = showOutput, invisible=(!showOutput), wait=TRUE, timeout = timeout)
    } else {
      if(showOutput) stdout.value = ""
      else stdout.value = NULL
      #need to switch to each directory, then run Mplus within using just the filename
      oldwd <- getwd()
      setwd(dirtocd)
      if (local_tmpdir) { Sys.setenv(TMPDIR=dirtocd) } #define TMPDIR local to the .inp file to execute
      exitCode <- system2(Mplus_command, args=c(shQuote(inputSplit$filename)), stdout=stdout.value, wait=TRUE, timeout = timeout)
      if (exitCode > 0L) {
        return("timeout")
        warning("Mplus returned error code: ", exitCode, ", for model: ", inputSplit$filename, "\n")
      }
      setwd(oldwd)
    }
  }
  
  if (isLogOpen()) {
    writeLines(c("", paste("------End Mplus Model Run: ", format(Sys.time(), "%d%b%Y %H:%M:%S"), "------", sep="")), logTarget)
    flush(logTarget)
  }
  normalComplete <- TRUE #if we made it here, then all completed normally
  
  #exitRun will fire here if all successful.
}

