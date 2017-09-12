RemoveAdapter <-R6Class(
  classname = "RemoveAdapter",
  inherit = BaseProc,
  public = list(
    initialize = function(atacProc,adapter1=NULL,adapter2=NULL,fastqOutput1=NULL,reportPrefix=NULL,
                          fastqOutput2=NULL,fastqInput1=NULL, fastqInput2=NULL,interleave=FALSE,paramList="default",findParamList="default",editable=FALSE){
      super$initialize("RemoveAdapter",editable,list(arg1=atacProc))
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

      
      if(!is.null(paramList)&&sum(paramList!="default")>0){

          rejectp<-"--file1|--adapter1|--output1|--file2|--adapter2|--output2|--threads|--basename"
          private$checkParam(paramList,rejectp)
          private$paramlist[["paramList"]]<-c(paramList,private$paramlist[["paramList"]])
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
          if(.obtainConfigure("threads")>1){
              threadparam<-c("--threads",as.character(.obtainConfigure("threads")))
          }else{
              threadparam<-NULL
          }
          findParamList <- paste(c(private$paramlist[["findParamList"]],threadparam),collapse = " ")
          paramList <- paste(c(private$paramlist[["paramList"]],threadparam), collapse = " ")
          if(private$singleEnd){
              private$writeLog("begin to remove adapter")
              private$writeLog("source:")
              private$writeLog(private$paramlist[["fastqInput1"]])
              private$writeLog(paste0("Adapter1:",private$paramlist[["adapter1"]]))
              private$writeLog("Destination:")
              private$writeLog(private$paramlist[["fastqOutput1"]])
              private$writeLog(private$paramlist[["reportPrefix"]])
              private$writeLog(paste0("other parameters:",.obtainConfigure("threads")))
#              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=private$paramlist[["adapter1"]],
#                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
#                                    paramlist=private$paramlist[["paramList"]])
              if(length(paramList)>0){
                  remove_adapters(file1 = private$paramlist[["fastqInput1"]], 
                                  paramList,
                                  adapter1 = private$paramlist[["adapter1"]], 
                                  output1 = private$paramlist[["fastqOutput1"]],
                                  basename = private$paramlist[["reportPrefix"]],
                                  interleaved = FALSE,
                                  overwrite = TRUE)
              }else{
                  remove_adapters(file1 = private$paramlist[["fastqInput1"]], 
                                  adapter1 = private$paramlist[["adapter1"]], 
                                  output1 = private$paramlist[["fastqOutput1"]],
                                  basename = private$paramlist[["reportPrefix"]],
                                  interleaved = FALSE,
                                  overwrite = TRUE)
              }
              
          }else if(private$paramlist[["interleave"]]){
              adapter1<-private$paramlist[["adapter1"]]
              adapter2<-private$paramlist[["adapter2"]]
              if(is.null(private$paramlist[["adapter1"]])){
                  private$writeLog("begin to find adapter")
                  if(length(findParamList)>0){
                      adapters<-identify_adapters(file1 = private$paramlist[["fastqInput1"]],
                                                  file2 = NULL,
                                                  findParamList,
                                                  basename = private$paramlist[["reportPrefix"]], overwrite=TRUE)
                  }else{
                      adapters<-identify_adapters(file1 = private$paramlist[["fastqInput1"]],
                                                  file2 = NULL,
                                                  basename = private$paramlist[["reportPrefix"]],overwrite=TRUE)
                  }
                  
                  adapter1 <- adapters[1]
                  adapter2 <- adapters[2]
              }
              if(length(paramList)>0){
                  remove_adapters(file1 = private$paramlist[["fastqInput1"]], 
                                  paramList,
                                  adapter1 = adapter1, 
                                  output1 = private$paramlist[["fastqOutput1"]],
                                  adapter2 = adapter2, 
                                  basename = private$paramlist[["reportPrefix"]],
                                  interleaved = TRUE,
                                  overwrite = TRUE)
              }else{
                  remove_adapters(file1 = private$paramlist[["fastqInput1"]], 
                                  adapter1 = adapter1, 
                                  output1 = private$paramlist[["fastqOutput1"]],
                                  adapter2 = adapter2, 
                                  basename = private$paramlist[["reportPrefix"]],
                                  interleaved = TRUE,
                                  overwrite = TRUE)
              }
 
              
          }else{
              adapter1<-private$paramlist[["adapter1"]]
              adapter2<-private$paramlist[["adapter2"]]
              if(is.null(private$paramlist[["adapter1"]])){
                  private$writeLog("begin to find adapter")
                  if(length(findParamList)>0){
                      adapters<-identify_adapters(file1 = private$paramlist[["fastqInput1"]],
                                                  file2 = private$paramlist[["fastqInput2"]],
                                                  findParamList,
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
#              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=adapter1,
#                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
#                                    inputFile2=private$paramlist[["fastqInput2"]],adapter2=adapter2,
#                                    outputFile2 = private$paramlist[["fastqOutput2"]],paramlist=private$paramlist[["paramList"]])
              if(length(paramList)>0){
                  remove_adapters(file1 = private$paramlist[["fastqInput1"]], 
                                  paramList,
                                  adapter1 = adapter1, 
                                  output1 = private$paramlist[["fastqOutput1"]],
                                  file2 = private$paramlist[["fastqInput2"]], 
                                  adapter2 = adapter2, 
                                  output2 = private$paramlist[["fastqOutput2"]],
                                  basename = private$paramlist[["reportPrefix"]],
                                  interleaved = FALSE,
                                  overwrite = TRUE)
              }else{
                  remove_adapters(file1 = private$paramlist[["fastqInput1"]], 
                                  adapter1 = adapter1, 
                                  output1 = private$paramlist[["fastqOutput1"]],
                                  file2 = private$paramlist[["fastqInput2"]], 
                                  adapter2 = adapter2, 
                                  output2 = private$paramlist[["fastqOutput2"]],
                                  basename = private$paramlist[["reportPrefix"]],
                                  interleaved = FALSE,
                                  overwrite = TRUE)
                  
              }
              
          }
      },
    checkRequireParam = function(){
      if(is.null(private$paramlist[["fastqInput1"]])){
        stop("'fastqInput1' is required.")
      }
      if(private$singleEnd&&is.null(private$paramlist[["adapter1"]])){
        stop("Parameter 'adapter1' is requied for single end sequencing data.")
      }
        if(private$paramlist[["interleave"]]&&private$singleEnd){
            stop("Single end data should not be interleave")
        }
    },
    checkAllPath = function(){
       private$checkFileExist(private$paramlist[["fastqInput1"]]);
       private$checkFileExist(private$paramlist[["fastqInput2"]]);
       private$checkFileCreatable(private$paramlist[["fastqOutput1"]]);
       private$checkFileCreatable(private$paramlist[["fastqOutput2"]]);
       private$checkPathExist(private$paramlist[["reportPrefix"]]);

    },
    getReportValImp = function(item){
        if(sum(item == c("adapter1","adapter2"))>0){
            adapter<-readLines(paste0(private$paramlist[["reportPrefix"]],".",item))
            return(adapter[1])
        }
        if(item == "settings"){
            tblist <- private$getTopic("\\[Adapter trimming\\]")
            splitlist <- strsplit(tblist,": ")
            lst <- list()
            for(i in 1:length(tblist)){
                lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
            }
            return(lst)
            
        }
        if(item == "statistics"){
            tblist <- private$getTopic("\\[Trimming statistics\\]")
            splitlist <- strsplit(tblist,": ")
            lst <- list()
            for(i in 1:length(tblist)){
                lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
            }
            return(lst)
        }
        if(item == "distribution"){
            tblist <- private$getTopic("\\[Length distribution\\]")
            splitlist <- strsplit(tblist,"\t")
            colkey <- splitlist[[1]]
            tbdt <- NULL
            for(i in 2:length(tblist)){
                tbdt <- c(tbdt,splitlist[[i]])
            }
            tbdt<-as.integer(tbdt)
            df<-as.data.frame(matrix(tbdt,length(tblist)-1,4,TRUE))
            colnames(df) <- colkey
            return(df)
        }
        stop(paste0(item," is not an item of report value."))
    },
    getReportItemsImp = function(){
        return(c("adapter1","adapter2","settings","statistics","distribution"))
    },
    getTopic = function(topic){
        setLine<-readLines(paste0(private$paramlist[["reportPrefix"]],".settings"))
        itemstarts<-grep("\\[",setLine)
        
        itemstart<-grep(topic,setLine)
        itemsendidx<-which(itemstarts == itemstart) + 1
        if(itemsendidx>length(itemstarts)){
            itemend <- length(setLine)
        }else{
            itemend <- itemstarts[itemsendidx]
            itemend <- itemend - 3
        }   
        
        return(setLine[(itemstart+1):itemend])
    }
  )


)



atacRemoveAdapter <- function(atacProc,adapter1=NULL,adapter2=NULL,
                              fastqOutput1=NULL,reportPrefix=NULL,
                              fastqOutput2=NULL,fastqInput1=NULL, 
                              fastqInput2=NULL,interleave=FALSE,
                              paramList="default",findParamList="default"){
    removeAdapter <- RemoveAdapter$new(atacProc = atacProc,
                                       adapter1 = adapter1,adapter2 = adapter2,
                                       fastqOutput1 = fastqOutput1, reportPrefix = reportPrefix,
                                       fastqOutput2 = fastqOutput2, fastqInput1 = fastqInput1, 
                                       fastqInput2 = fastqInput2, interleave = interleave,
                                       paramList = paramList,findParamList = findParamList)
    removeAdapter$process()
    return(removeAdapter)
}

