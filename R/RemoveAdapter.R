RemoveAdapter <-R6Class(
    classname = "RemoveAdapter",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, adapter1=NULL,adapter2=NULL,fastqOutput1=NULL, reportPrefix=NULL,
                              fastqOutput2=NULL, fastqInput1=NULL, fastqInput2=NULL,interleave=FALSE,
                              threads=NULL, paramList= NULL,findParamList=NULL,editable=FALSE){
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
            
            private$paramlist[["paramList"]]<-""
            if(!is.null(paramList)){
                paramList<-trimws(as.character(paramList))
                paramList<-paste(paramList,collapse = " ")
                paramList <- strsplit(paramList,"\\s+")[[1]]
                if(length(paramList)>0){
                    rejectp<-"--file1|--adapter1|--output1|--file2|--adapter2|--output2|--threads|--basename"
                    private$checkParam(paramList,rejectp)
                    private$paramlist[["paramList"]]<-paramList
                }
            }
            
            
            private$paramlist[["findParamList"]] <- ""         
            if(!is.null(findParamList)){
                findParamList<-trimws(as.character(findParamList))
                findParamList<-paste(findParamList,collapse = " ")
                findParamList <- strsplit(findParamList,"\\s+")[[1]]  
                if(length(paramList)>0){
                    rejectp<-"--file1|--file2|--threads|--identify-adapters|--basename"
                    private$checkParam(findParamList,rejectp)
                    private$paramlist[["findParamList"]]<-findParamList
                }
            }
            
            if(!is.null(threads)){
                private$paramlist[["threads"]] <- as.integer(threads)
            }
            
            private$paramValidation()
            
            
        }
    ),
    private = list(
        processing = function(){
            if(!is.null(private$paramlist[["threads"]])){
                if(private$paramlist[["threads"]]>1){
                    threadparam<-c("--threads",as.character(private$paramlist[["threads"]]))
                }
            }else if(.obtainConfigure("threads")>1){
                threadparam<-c("--threads",as.character(.obtainConfigure("threads")))
            }else{
                threadparam<-NULL
            }
            findParamList <- paste(c(threadparam, private$paramlist[["findParamList"]]),collapse = " ")
            paramList <- paste(c(threadparam, private$paramlist[["paramList"]]), collapse = " ")
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
            if(item == "adapterslist"){
                tblist <- private$getTopic("\\[Adapter sequences\\]")
                splitlist <- strsplit(tblist,": ")
                if(private$singleEnd){
                    return(data.frame(adapter=c("adapter for single end data"),sequence=c(splitlist[[1]][2])))
                }else{
                    return(data.frame(adapter=c("adapter for paired end data mate 1","adapter for paired end data mate 2"),
                                      sequence=c(splitlist[[1]][2],splitlist[[2]][2])))
                    #return(list(adapter1=splitlist[[1]][2],
                    #            adapter2=splitlist[[2]][2]))
                }
            }
            if(item == "adapters"){
                tblist <- private$getTopic("\\[Adapter sequences\\]")
                splitlist <- strsplit(tblist,": ")
                if(private$singleEnd){
                    return(private$listToFrame(list(adpater1=splitlist[[1]][2])))
                }else{
                    return(private$listToFrame(list(adapter1=splitlist[[1]][2],
                                                    adapter2=splitlist[[2]][2])))
                }
            }
            if(sum(item == c("adapter1","adapter2"))>0){
                tblist <- private$getTopic("\\[Adapter sequences\\]")
                splitlist <- strsplit(tblist,": ")
                if(item == "adapter1"){
                    return(splitlist[[1]][2])
                }
                if(item == "adapter2"){
                    return(splitlist[[2]][2])
                }
            }
            if(item == "settingslist"){
                tblist <- private$getTopic("\\[Adapter trimming\\]")
                splitlist <- strsplit(tblist,": ")
                lst <- list()
                for(i in 1:length(tblist)){
                    lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
                }
                return(lst)
            }
            if(item == "settings"){
                tblist <- private$getTopic("\\[Adapter trimming\\]")
                splitlist <- strsplit(tblist,": ")
                lst <- list()
                for(i in 1:length(tblist)){
                    lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
                }
                return(private$listToFrame(lst))
                
            }
            if(item == "statisticslist"){
                tblist <- private$getTopic("\\[Trimming statistics\\]")
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
                return(private$listToFrame(lst))
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
                if(self$isSingleEnd()){
                    colsize<-4
                }else{
                    colsize<-6
                }
                df<-as.data.frame(matrix(tbdt,length(tblist)-1,colsize,TRUE))
                colnames(df) <- colkey
                return(df)
            }
            stop(paste0(item," is not an item of report value."))
        },
        getReportItemsImp = function(){
            if(self$isSingleEnd()){
                return(c("adapters","adapterslist","adapter1","settings","statistics","settingslist","statisticslist","distribution"))
            }else{
                return(c("adapters","adapterslist","adapter1","adapter2","settings","statistics","settingslist","statisticslist","distribution"))
            }
            
        },
        getTopic = function(topic){
            setLine<-readLines(paste0(private$paramlist[["reportPrefix"]],".settings"))
            itemstarts<-grep("\\]$",setLine)
            
            itemstart<-grep(topic,setLine)
            itemsendidx<-which(itemstarts == itemstart) + 1
            if(itemsendidx>length(itemstarts)){
                itemend <- length(setLine)
            }else{
                itemend <- itemstarts[itemsendidx]
                itemend <- itemend - 3
            }   
            
            return(setLine[(itemstart+1):itemend])
        },
        listToFrame = function(a){
            return(data.frame(Item=names(a),Value=as.character(a)))
        }
    )
    
    
)


#' @name atacRemoveAdapter
#' @aliases atacRemoveAdapter
#' @aliases removeAdapter
#' @title Use AdapterRemoval to remove adapters
#' @description 
#' Use AdapterRemoval to remove adapters
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacRenamer}} 
#' \code{\link{renamer}}
#' \code{\link{atacUnzipAndMerge}} 
#' \code{\link{unzipAndMerge}}
#' @param fastqInput1 \code{Character} vector. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file paths with #1 mates paired
#' with file paths in fastqInput2
#' And it can also be interleaved file paths when argument
#' interleaved=\code{TRUE}
#' @param fastqInput2 \code{Character} vector. It contains file paths with #2
#' mates paired with file paths in fastqInput1
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param adapter1 \code{Character}. It is an adapter sequence for file1. 
#' For single end data, it is requied.
#' @param adapter2 \code{Character}. It is an adapter sequence for file2.
#' @param fastqOutput1 \code{Character}. The trimmed mate1 reads output file
#' path for fastqInput2. Defualt:
#' basename.pair1.truncated (paired-end),
#' basename.truncated (single-end), or
#' basename.paired.truncated (interleaved)
#' @param fastqOutput2 \code{Character}. The trimmed mate2 reads output file
#' path for fastqInput2. Default:
#' BASENAME.pair2.truncated (only used in PE mode, but not if
#' --interleaved-output is enabled)
#' @param interleave \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param paramList Additional arguments to be passed on to the binaries 
#' for removing adapter. See below for details.
#' @param findParamList Additional arguments to be passed on to the binaries
#' for identifying adapter. See below for details.
#' @param reportPrefix \code{Character}. The prefix of report files path. 
#' Default: generate from known parameters
#' @details The parameter related to input and output file path
#' will be automatically 
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{fregLenDistr} instead.
#' You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8") with white space
#' splited just like command line,
#' or put them in \code{Character} vector(e.g. c("--threads","8")).
#' Note that some arguments(
#' "--file1","--file2","--adapter1","--adapter2","--output1","--output2",
#' "--basename","--interleaved","thread") to the
#' paramList and findParamList are invalid if they are already handled as explicit
#' function arguments. See the output of
#' \code{adapterremoval_usage()} for details about available parameters.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacRenamer}} 
#' \code{\link{renamer}}
#' \code{\link{atacUnzipAndMerge}} 
#' \code{\link{unzipAndMerge}}
#' \code{\link{atacBowtie2Mapping}} 
#' @examples 
#' library(magrittr)
#' td <- tempdir()
#' setConfigure("tmpdir",td)
#' 
#' # Identify adapters
#' prefix<-system.file(package="ATACFlow", "extdata", "uzmg")
#' (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
#' (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
#' 
#' reads_merged_1 <- file.path(td,"reads1.fastq")
#' reads_merged_2 <- file.path(td,"reads2.fastq")
#' atacproc <- 
#' atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2) %>%
#' atacRenamer %>% atacRemoveAdapter
#' 
#' dir(td) 
#' 


#' @rdname atacRemoveAdapter
#' @export 
atacRemoveAdapter <- function(atacProc,adapter1=NULL,adapter2=NULL,
                              fastqOutput1=NULL,reportPrefix=NULL,
                              fastqOutput2=NULL,fastqInput1=NULL, 
                              fastqInput2=NULL,interleave=FALSE,
                              paramList= NULL,findParamList=NULL){
    atacproc <- RemoveAdapter$new(atacProc = atacProc,
                                       adapter1 = adapter1,adapter2 = adapter2,
                                       fastqOutput1 = fastqOutput1, reportPrefix = reportPrefix,
                                       fastqOutput2 = fastqOutput2, fastqInput1 = fastqInput1, 
                                       fastqInput2 = fastqInput2, interleave = interleave,
                                       paramList = paramList,findParamList = findParamList)
    atacproc$process()
    invisible(atacproc)
}
#' @rdname atacRemoveAdapter
#' @export 
removeAdapter <- function(fastqInput1, fastqInput2=NULL,
                              adapter1=NULL,adapter2=NULL,
                              fastqOutput1=NULL,reportPrefix=NULL,
                              fastqOutput2=NULL, 
                              interleave=FALSE,
                              paramList = NULL,findParamList = NULL){
    atacproc <- RemoveAdapter$new(atacProc = NULL,
                                       adapter1 = adapter1,adapter2 = adapter2,
                                       fastqOutput1 = fastqOutput1, reportPrefix = reportPrefix,
                                       fastqOutput2 = fastqOutput2, fastqInput1 = fastqInput1, 
                                       fastqInput2 = fastqInput2, interleave = interleave,
                                       paramList = paramList,findParamList = findParamList)
    atacproc$process()
    invisible(atacproc)
}

