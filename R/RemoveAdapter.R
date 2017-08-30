RemoveAdapter <-R6Class(
  classname = "RemoveAdapter",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc,adapter1=NULL,adapter2=NULL,fastqOutput1=NULL,reportPrefix=NULL,
                          fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,paramList="default",findParamList="default",editable=FALSE){
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

      if(.obtainConfigure("threads")>1){
          private$paramlist[["paramList"]]<-c("--threads",as.character(.obtainConfigure("threads")))
      }
      if(!is.null(paramList)&&sum(paramList!="default")>0){

          rejectp<-"--file1|--adapter1|--output1|--file2|--adapter2|--output2|--threads|--basename"
          private$checkParam(paramList,rejectp)
          private$paramlist[["paramList"]]<-c(paramList,private$paramlist[["paramList"]])
      }


      if(.obtainConfigure("threads")>1){
          private$paramlist[["findParamList"]]<-c("--threads",as.character(.obtainConfigure("threads")))
      }
      if(!is.null(findParamList)&&sum(findParamList!="default")>0){
          rejectp<-"--file1|--file2|--threads|--identify-adapters|--basename"
          private$checkParam(findParamList,rejectp)
          private$paramlist[["findParamList"]]<-c(findParamList,private$paramlist[["findParamList"]])
      }

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
              private$writeLog(paste0("other parameters:",.obtainConfigure("threads")))
              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=private$paramlist[["adapter1"]],
                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
                                    paramlist=private$paramlist[["paramList"]])
          }else{
              adapter1<-private$paramlist[["adapter1"]]
              adapter2<-private$paramlist[["adapter2"]]
              if(is.null(private$paramlist[["adapter1"]])){
                  private$writeLog("begin to find adapter")
                  if(length(private$paramlist[["findParamList"]])>0){
                      adapters<-identify_adapters(file1 = private$paramlist[["fastqInput1"]],
                                                  file2 = private$paramlist[["fastqInput2"]],
                                                  paste(private$paramlist[["findParamList"]],collapse = " "),
                                                  basename = private$paramlist[["reportPrefix"]], overwrite=TRUE)
                  }else{
                      adapters<-identify_adapters(file1 = private$paramlist[["fastqInput1"]],
                                                  file2 = private$paramlist[["fastqInput2"]],
                                                  basename = private$paramlist[["reportPrefix"]],overwrite=TRUE)
                  }
                  
                  adapter1 <- adapters[1]
                  adapter2 <- adapters[2]
              }
              private$writeLog("begin to remove adapter")
              private$writeLog("source:")
              private$writeLog(private$paramlist[["fastqInput1"]])
              private$writeLog(private$paramlist[["fastqInput2"]])
              private$writeLog(paste0("Adapter1:",adapter1))
              private$writeLog(paste0("Adapter2:",adapter2))
              private$writeLog("Destination:")
              private$writeLog(private$paramlist[["fastqOutput1"]])
              private$writeLog(private$paramlist[["fastqOutput2"]])
              private$writeLog(private$paramlist[["reportPrefix"]])
              private$writeLog(paste0("Threads:",.obtainConfigure("threads")))
              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=adapter1,
                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
                                    inputFile2=private$paramlist[["fastqInput2"]],adapter2=adapter2,
                                    outputFile2 = private$paramlist[["fastqOutput2"]],paramlist=private$paramlist[["paramList"]])
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
                              fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,paramList="default",findParamList="default"){
    removeAdapter <- RemoveAdapter$new(atacProc,adapter1,adapter2,fastqOutput1,reportPrefix,
                                       fastqOutput2,fastqInput1, fastqInput2,paramList,findParamList)
    removeAdapter$process()
    return(removeAdapter)
}

