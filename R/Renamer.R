Renamer <-R6Class(
  classname = "Renamer",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL, interleave = FALSE, editable=FALSE){
      super$initialize("Renamer",editable,list(arg1=atacProc))
      if(!is.null(atacProc)){
        private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
        private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
        regexProcName<-sprintf("(fastq|fq|%s)",atacProc$getProcName())
        private$paramlist[["interleave"]] <- atacProc$getParam("interleave")
      }else{
        regexProcName<-"(fastq|fq)"
        private$paramlist[["interleave"]] <- interleave
        if(is.null(fastqInput2)){
            private$singleEnd<-TRUE
        }else{
            private$singleEnd<-FALSE
        }
      }

      if(!is.null(fastqInput1)){
        private$paramlist[["fastqInput1"]] <- fastqInput1;
      }
      if(!is.null(fastqInput2)){
        private$paramlist[["fastqInput2"]] <- fastqInput2;
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
        if(private$singleEnd||private$paramlist[["interleave"]]){
            private$singleCall(1)
        }else if(.obtainConfigure("threads")>=2){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s",private$paramlist[["fastqInput2"]]))
            private$writeLog(sprintf("destination:%s",private$paramlist[["fastqOutput2"]]))
            cl<-makeCluster(2)
            parLapply(cl = cl,X = 1:2,fun = private$singleCall)
            stopCluster(cl)
        }else{
            private$singleCall(1)
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s",private$paramlist[["fastqInput2"]]))
            private$writeLog(sprintf("destination:%s",private$paramlist[["fastqOutput2"]]))
            private$singleCall(2)
        }
    },
    singleCall = function(number){
        if(number==1){
            .renamer_call(inputFile = private$paramlist[["fastqInput1"]],
                          outputFile = private$paramlist[["fastqOutput1"]],
                          interleave = private$paramlist[["interleave"]])
        }else if(number==2){
            .renamer_call(inputFile = private$paramlist[["fastqInput2"]],
                          outputFile = private$paramlist[["fastqOutput2"]], 
                          interleave = private$paramlist[["interleave"]])
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


atacRenamer <- function(atacProc, 
                        fastqOutput1=NULL,
                        fastqOutput2=NULL,
                        fastqInput1=NULL, 
                        fastqInput2=NULL, 
                        interleave = FALSE){
    atacproc <- Renamer$new(atacProc = atacProc,
                            fastqOutput1 = fastqOutput1,
                            fastqOutput2 = fastqOutput2,
                            fastqInput1 = fastqInput1,
                            fastqInput2 = fastqInput2)
    atacproc$process()
    return(atacproc)
}

