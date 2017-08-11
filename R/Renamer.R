Renamer <-R6Class(
  classname = "Renamer",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,editable=FALSE){
      super$initialize("Renamer",editable,list(arg1=atacProc))
      if(!is.null(atacProc)){
        private$paramlist[["fastqInput1"]] <- atacProc$getParam("fastqOutput1");
        private$paramlist[["fastqInput2"]] <- atacProc$getParam("fastqOutput2");
        regexProcName<-paste0("|",atacProc$getProcName())
      }else{
        regexProcName<-""
      }

      if(!is.null(fastqInput1)){
        private$paramlist[["fastqInput1"]] <- fastqInput1;
      }
      if(!is.null(fastqInput2)){
        private$paramlist[["fastqInput2"]] <- fastqInput2;
      }

      if(is.null(fastqOutput1)){
        prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],sprintf("(fastq|fq%s)",regexProcName))
        private$paramlist[["fastqOutput1"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".fq"))
      }else{
        private$paramlist[["fastqOutput1"]] <- fastqOutput1;
      }
      if(is.null(fastqOutput2)){
        if(!is.null(private$paramlist[["fastqInput2"]])){
          prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput2"]],sprintf("(fastq|fq%s)",regexProcName))
          private$paramlist[["fastqOutput2"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".fq"));
        }
      }else{
        private$paramlist[["fastqOutput2"]] <- fastqOutput2;
        private$singleEnd<-FALSE
      }

      private$checkFileExist(private$paramlist[["fastqInput1"]]);
      private$checkFileExist(private$paramlist[["fastqInput2"]]);

      private$checkFileCreatable(private$paramlist[["fastqOutput1"]]);
      private$checkFileCreatable(private$paramlist[["fastqOutput2"]]);

      private$checkRequireParam()
    },
    processing = function(){
      if(!super$processing()){
          return()
      }
      .renamer_call(private$paramlist[["fastqInput1"]],private$paramlist[["fastqOutput1"]])
      if(!is.null(private$paramlist[["fastqInput2"]])){
        .renamer_call(private$paramlist[["fastqInput2"]],private$paramlist[["fastqOutput2"]])
      }
      private$setFinish()
    },
    setResultParam = function(fastqOutput1, fastqOutput2=NULL){
      super$setResultParam();
      private$paramlist[["fastqOutput1"]] <- fastqOutput1
      private$paramlist[["fastqOutput2"]] <- fastqOutput2
    }
  ),
  private = list(
    checkRequireParam = function(){
      if(private$editable){
        return();
      }
      if(is.null(private$paramlist[["fastqInput1"]])){
        stop("fastqInput1 is required.")
      }

    }
  )

)


atacRenamer <- function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL){
    atacproc <- Renamer$new(atacProc, fastqOutput1, fastqOutput2,fastqInput1,fastqInput2)
    atacproc$processing()
    return(atacproc)
}

