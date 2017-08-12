Renamer <-R6Class(
  classname = "Renamer",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,editable=FALSE){
      super$initialize("Renamer",editable,list(arg1=atacProc))
      if(!is.null(atacProc)){
        private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
        private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
        regexProcName<-sprintf("(fastq|fq|%s)",atacProc$getProcName())
      }else{
        regexProcName<-"(fastq|fq)"
      }

      if(!is.null(fastqInput1)){
        private$paramlist[["fastqInput1"]] <- fastqInput1;
      }
      if(!is.null(fastqInput2)){
        private$paramlist[["fastqInput2"]] <- fastqInput2;
      }

      if(is.null(private$paramlist[["fastqInput2"]])){
        private$singleEnd<-TRUE
      }else{
        private$singleEnd<-FALSE
      }

      if(is.null(fastqOutput1)){
        if(!is.null(private$paramlist[["fastqInput1"]])){
          prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],regexProcName)
          private$paramlist[["fastqOutput1"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".fq"))
        }
      }else{
        private$paramlist[["fastqOutput1"]] <- fastqOutput1;
      }
      if(is.null(fastqOutput2)){
        if(!is.null(private$paramlist[["fastqInput2"]])){
          prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput2"]],regexProcName)
          private$paramlist[["fastqOutput2"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".fq"));
        }
      }else{
        private$paramlist[["fastqOutput2"]] <- fastqOutput2;
      }

      private$paramValidation()

    }
  ),
  private = list(
    processing = function(){
      private$writeLog(paste0("processing file:"))
      private$writeLog(sprintf("source:%s",private$paramlist[["fastqInput1"]]))
      private$writeLog(sprintf("destination:%s",private$paramlist[["fastqOutput1"]]))
      .renamer_call(private$paramlist[["fastqInput1"]],private$paramlist[["fastqOutput1"]])
      if(!is.null(private$paramlist[["fastqInput2"]])){
        private$writeLog(paste0("processing file:"))
        private$writeLog(sprintf("source:%s",private$paramlist[["fastqInput2"]]))
        private$writeLog(sprintf("destination:%s",private$paramlist[["fastqOutput2"]]))
        .renamer_call(private$paramlist[["fastqInput2"]],private$paramlist[["fastqOutput2"]])
       }
    },
    checkRequireParam = function(){
      if(is.null(private$paramlist[["fastqInput1"]])){
        stop("fastqInput1 is required.")
      }

    },
    checkAllPath = function(){
        private$checkFileExist(private$paramlist[["fastqInput1"]]);
        private$checkFileExist(private$paramlist[["fastqInput2"]]);
        private$checkFileCreatable(private$paramlist[["fastqOutput1"]]);
        private$checkFileCreatable(private$paramlist[["fastqOutput2"]]);
    }
  )

)


atacRenamer <- function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL){
    atacproc <- Renamer$new(atacProc, fastqOutput1, fastqOutput2,fastqInput1,fastqInput2)
    atacproc$process()
    return(atacproc)
}

