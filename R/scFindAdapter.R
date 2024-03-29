setClass(Class = "SCFindAdapter",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "SCFindAdapter",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        fastqInput1 <- allparam[["fastqInput1"]]
        fastqInput2 <- allparam[["fastqInput2"]]
        paramList <- allparam[["paramList"]]
        findParamList <- allparam[["findParamList"]]
        reportPrefix <- allparam[["reportPrefix"]]
        threads <- allparam[["threads"]]
        
        if(length(prevSteps) >0){
            prevSteps <- prevSteps[[1]]
            input(.Object)[["fastqInput1"]] <- output(prevSteps)[["fastqOutput1"]]
            input(.Object)[["fastqInput2"]] <- output(prevSteps)[["fastqOutput2"]]
        }
        
        if(!is.null(fastqInput1)){
            input(.Object)[["fastqInput1"]] <- fastqInput1;
        }
        
        if(!is.null(fastqInput2)){
            input(.Object)[["fastqInput2"]] <- fastqInput2;
        }
        
        if(is.null(reportPrefix)){
            param(.Object)$reportPrefix <- getStepWorkDir(.Object,"find")
        }else{
            param(.Object)$reportPrefix <- reportPrefix
        }
        
        output(.Object)$reportPrefix_adapter1_Output <- paste0(param(.Object)$reportPrefix,".adapter1")
        output(.Object)$reportPrefix_adapter2_Output <- paste0(param(.Object)$reportPrefix,".adapter2")
        
        
        param(.Object)$findParamList<-""
        if(!is.null(findParamList)){
            findParamList<-trimws(as.character(findParamList))
            findParamList<-paste(findParamList,collapse = " ")
            findParamList <- strsplit(findParamList,"\\s+")[[1]]
            if(length(paramList)>0){
                rejectp<-"--file1|--file2|--threads|--identify-adapters|--basename"
                checkParam(.Object,findParamList,rejectp)
                param(.Object)$findParamList<-findParamList
            }
        }
        
        stopifnot(is.numeric(threads))
        param(.Object)$threads <- threads
        
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "SCFindAdapter",
    definition = function(.Object,...){
        threadparam <- NULL
        if(!is.null(param(.Object)$threads)){
            if(param(.Object)$threads>1){
                threadparam<-c("--threads",as.character(param(.Object)$threads))
            }
        }else if(getThreads()>1){
            threadparam<-c("--threads",as.character(getThreads()))
        }else{
            threadparam<-NULL
        }
        findParamList <- paste(c(threadparam, param(.Object)$findParamList),collapse = " ")
            writeLog(.Object,"begin to find adapter")
            if(length(findParamList)>0){
                print(input(.Object)[["fastqInput1"]])
                print(input(.Object)[["fastqInput2"]])
                print(findParamList)
                print(param(.Object)[["reportPrefix"]])
                adapters<-identify_adapters(file1 = input(.Object)[["fastqInput1"]],
                                            file2 = input(.Object)[["fastqInput2"]],
                                            findParamList,
                                            basename = param(.Object)[["reportPrefix"]], overwrite=TRUE)
            }else{
                adapters<-identify_adapters(file1 = output(.Object)[["fastqInput1"]],
                                            file2 = output(.Object)[["fastqInput2"]],
                                            basename = param(.Object)[["reportPrefix"]], overwrite=TRUE)
            }
            property(.Object)[["adapter1"]] <- adapters[1]
            property(.Object)[["adapter2"]] <- adapters[2]
        
       
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "SCFindAdapter",
    definition = function(.Object, ...){
        
        report(.Object)[["adapter1"]] <- readLines(output(.Object)$reportPrefix_adapter1_Output)
        report(.Object)[["adapter2"]] <- readLines(output(.Object)$reportPrefix_adapter2_Output)
        .Object
    }
)



#' @name SCFindAdapter
#' @title Use AdapterRemoval to identify adapters
#' @description
#' Use AdapterRemoval to identify adapters for paired end data
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
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
#' @param interleave \code{Logical}. Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param findParamList Additional arguments to be passed on to the binaries
#' for identifying adapter. See below for details.
#' @param reportPrefix \code{Character}. The prefix of report files path.
#' Default: generate from known parameters
#' @param threads The number of threads used in this step.
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{scFindAdapter} instead.
#' You can put all aditional
#' arguments in one \code{Character}(e.g. "--threads 8") with white space
#' splited just like command line,
#' or put them in \code{Character} vector(e.g. c("--threads","8")).
#' Note that some arguments(
#' "--file1","--file2","--adapter1","--adapter2","--output1","--output2",
#' "--basename","--interleaved","thread") to the
#' findParamList are invalid if they are already handled as explicit
#' function arguments. See the output of
#' \code{adapterremoval_usage()} for details about available parameters.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
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
#' setTmpDir(td)
#'
#' # Identify adapters
#' prefix<-system.file(package="esATAC", "extdata", "uzmg")
#' (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
#' (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
#'
#' reads_merged_1 <- file.path(td,"reads1.fastq")
#' reads_merged_2 <- file.path(td,"reads2.fastq")
#' atacproc <-
#' atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2) %>%
#' atacRenamer %>% atacSCFindAdapter
#'
#' dir(td)
#' @importFrom Rbowtie2 identify_adapters



setGeneric("atacSCFindAdapter",function(atacProc,fastqInput1=NULL,
                                        fastqInput2=NULL, reportPrefix = NULL,
                                        findParamList=NULL, threads = getThreads(), ...) standardGeneric("atacSCFindAdapter"))
#' @rdname SCFindAdapter
#' @aliases atacSCFindAdapter
#' @export
setMethod(
    f = "atacSCFindAdapter",
    signature = "ATACProc",
    definition = function(atacProc,fastqInput1=NULL,
                          fastqInput2=NULL, reportPrefix = NULL, 
                          findParamList=NULL, threads = getThreads(), ...){
        allpara <- c(list(Class = "SCFindAdapter", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname SCFindAdapter
#' @aliases scFindAdapter
#' @export
scFindAdapter <- function(fastqInput1, fastqInput2 = NULL, reportPrefix = NULL,
                           findParamList = NULL, threads = getThreads(), ...){
    allpara <- c(list(Class = "SCFindAdapter", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}

