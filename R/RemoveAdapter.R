RemoveAdapter <-R6Class(
  classname = "RemoveAdapter",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc,adapter1=NULL,adapter2=NULL,fastqOutput1=NULL,reportPrefix=NULL,
                          fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,editable=FALSE){
      super$initialize("RemoveAdapter",editable,list(arg1=atacProc))
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
      if(is.null(reportPrefix)){
        if(!is.null(private$paramlist[["fastqInput1"]])){
          prefix<-private$getBasenamePrefix(private$paramlist[["fastqInput1"]],regexProcName)
          private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
        }
      }else{
        private$paramlist[["reportPrefix"]] <- reportPrefix;
      }

      private$paramlist[["adapter1"]] <- adapter1
      private$paramlist[["adapter2"]] <- adapter2

      private$paramValidation()


    }
  ),
  private = list(
      processing = function(){
          if(private$singleEnd){
              private$writeLog("begin to remove adapter")
              private$writeLog("source:")
              private$writeLog(private$paramlist[["fastqInput1"]])
              private$writeLog(paste0("Adapter1:",private$paramlist[["adapter1"]]))
              private$writeLog("Destination:")
              private$writeLog(private$paramlist[["fastqOutput1"]])
              private$writeLog(private$paramlist[["reportPrefix"]])
              private$writeLog(paste0("Threads:",.obtainConfigure("threads")))
              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=private$paramlist[["adapter1"]],
                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
                                    threads=.obtainConfigure("threads"))
          }else{
              if(is.null(private$paramlist[["adapter1"]])){
                  private$writeLog("begin to find adapter")
                  adapters<-.identify_adapters_call(private$paramlist[["fastqInput1"]],
                                                    private$paramlist[["fastqInput2"]],threads=getConfigure("threads"))
                  private$paramlist[["adapter1"]] <- adapters$adapter1
                  private$paramlist[["adapter2"]] <- adapters$adapter2
              }
              private$writeLog("begin to remove adapter")
              private$writeLog("source:")
              private$writeLog(private$paramlist[["fastqInput1"]])
              private$writeLog(private$paramlist[["fastqInput2"]])
              private$writeLog(paste0("Adapter1:",private$paramlist[["adapter1"]]))
              private$writeLog(paste0("Adapter2:",private$paramlist[["adapter2"]]))
              private$writeLog("Destination:")
              private$writeLog(private$paramlist[["fastqOutput1"]])
              private$writeLog(private$paramlist[["fastqOutput2"]])
              private$writeLog(private$paramlist[["reportPrefix"]])
              private$writeLog(paste0("Threads:",.obtainConfigure("threads")))
              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=private$paramlist[["adapter1"]],
                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
                                    inputFile2=private$paramlist[["fastqInput2"]],adapter2=private$paramlist[["adapter2"]],
                                    outputFile2 = private$paramlist[["fastqOutput2"]],threads=.obtainConfigure("threads"))
          }
      },
    checkRequireParam = function(){
      if(is.null(private$paramlist[["fastqInput1"]])){
        stop("'fastqInput1' is required.")
      }
      if(private$singleEnd&&is.null(private$paramlist[["adapter1"]])){
        stop("Parameter 'adapter1' is requied for single end sequencing data.")
      }
    },
    checkAllPath = function(){
       private$checkFileExist(private$paramlist[["fastqInput1"]]);
       private$checkFileExist(private$paramlist[["fastqInput2"]]);
       private$checkFileCreatable(private$paramlist[["fastqOutput1"]]);
       private$checkFileCreatable(private$paramlist[["fastqOutput2"]]);
       private$checkPathExist(private$paramlist[["reportPrefix"]]);

    }
  )


)



atacRemoveAdapter <- function(atacProc,adapter1=NULL,adapter2=NULL,fastqOutput1=NULL,reportPrefix=NULL,
                              fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL){
    removeAdapter <- RemoveAdapter$new(atacProc,adapter1,adapter2,fastqOutput1,reportPrefix,
                                       fastqOutput2,fastqInput1, fastqInput2)
    removeAdapter$process()
    return(removeAdapter)
}

