setClass(Class = "RemoveAdapter",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "RemoveAdapter",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        fastqInput1 <- allparam[["fastqInput1"]]
        fastqInput2 <- allparam[["fastqInput2"]]
        adapter1 <- allparam[["adapter1"]]
        adapter2 <- allparam[["adapter2"]]
        fastqOutput1 <- allparam[["fastqOutput1"]]
        fastqOutput2 <- allparam[["fastqOutput2"]]
        reportPrefix <- allparam[["reportPrefix"]]
        interleave <- allparam[["interleave"]]
        paramList <- allparam[["paramList"]]
        threads <- allparam[["threads"]]


        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                fastqSteps <- prevSteps[[1]]
                fastqSteps<-c(unlist(fastqSteps),list())
                fastqStep <- fastqSteps[[length(fastqSteps)]]
                param(.Object)[["interleave"]] <- property(fastqStep)[["interleave"]]
                param(.Object)[["singleEnd"]] <- property(fastqStep)[["singleEnd"]]
                input(.Object)[["fastqInput1"]] <- output(fastqStep)[["fastqOutput1"]]
                if(!param(.Object)[["interleave"]] && !is.null(input(.Object)[["fastqInput1"]])){
                    input(.Object)[["fastqInput2"]] <- output(fastqStep)[["fastqOutput2"]]
                }
                print("++++++++++")
                print(input(.Object))
                print("++++++++++")
            }
            if(!is.null(prevSteps[[2]])){
                findAdapterStep <- prevSteps[[2]]
                param(.Object)$adapter1 <- property(findAdapterStep)[["adapter1"]]
                param(.Object)$adapter2 <- property(findAdapterStep)[["adapter2"]]
            }
        }else{
            param(.Object)[["interleave"]] <- interleave
            property(.Object)[["interleave"]] <- interleave
            if(is.null(fastqInput2)){
                property(.Object)[["singleEnd"]] <- TRUE
                property(.Object)[["singleEnd"]] <- TRUE
            }else{
                property(.Object)[["singleEnd"]] <- FALSE
                property(.Object)[["singleEnd"]] <- FALSE
            }
        }


        print(fastqInput1)
        if(!is.null(fastqInput1)){
            input(.Object)[["fastqInput1"]] <- fastqInput1;
        }

        if(!is.null(fastqInput2)){
            input(.Object)[["fastqInput2"]] <- fastqInput2;
        }


        print("++++++++++")
        print(input(.Object))
        print("++++++++++")
        if(is.null(fastqOutput1)){
            output(.Object)[["fastqOutput1"]] <- getAutoPath(.Object, input(.Object)[["fastqInput1"]], "fq|fastq", "fq")
        }else{
            output(.Object)[["fastqOutput1"]] <- fastqOutput1;
        }
        if(is.null(fastqOutput2)){
            if(!is.null(input(.Object)[["fastqInput2"]])){
                output(.Object)[["fastqOutput2"]] <- getAutoPath(.Object, input(.Object)[["fastqInput2"]], "fq|fastq", "fq")
            }
        }else{
            output(.Object)[["fastqOutput2"]] <- fastqOutput2;
        }
        if(is.null(reportPrefix)){
            param(.Object)$reportPrefix <- getStepWorkDir(.Object,"adrm.report")
        }else{
            param(.Object)$reportPrefix <- reportPrefix
        }

        if(!is.null(adapter1)){
            param(.Object)$adapter1 <- adapter1
        }
        if(!is.null(adapter2)){
            param(.Object)$adapter1 <- adapter1
        }
        
       
        param(.Object)$paramList<-""
        if(!is.null(paramList)){
            paramList<-trimws(as.character(paramList))
            paramList<-paste(paramList,collapse = " ")
            paramList <- strsplit(paramList,"\\s+")[[1]]
            if(length(paramList)>0){
                rejectp<-"--file1|--adapter1|--output1|--file2|--adapter2|--output2|--threads|--basename"
                checkParam(paramList,rejectp)
                param(.Object)$paramList<-paramList
            }
        }

       

        stopifnot(is.numeric(threads))
        param(.Object)$threads <- threads

        .Object
    }
)

setMethod(
    f = "processing",
    signature = "RemoveAdapter",
    definition = function(.Object,...){
        if(!is.null(param(.Object)$threads)){
            if(param(.Object)$threads>1){
                threadparam<-c("--threads",as.character(param(.Object)$threads))
            }
        }else if(.obtainConfigure("threads")>1){
            threadparam<-c("--threads",as.character(.obtainConfigure("threads")))
        }else{
            threadparam<-NULL
        }
        paramList <- paste(c(threadparam, param(.Object)$paramList), collapse = " ")
        if(param(.Object)$singleEnd){
            writeLog(.Object,"begin to remove adapter")
            writeLog(.Object,"source:",showMsg=FALSE)
            writeLog(.Object,input(.Object)[["fastqInput1"]],showMsg=FALSE)
            writeLog(.Object,paste0("Adapter1:",param(.Object)$adapter1))
            writeLog(.Object,"Destination:",showMsg=FALSE)
            writeLog(.Object,output(.Object)[["fastqOutput1"]],showMsg=FALSE)
            writeLog(.Object,param(.Object)$reportPrefix,showMsg=FALSE)
            writeLog(.Object,paste0("other parameters:",.obtainConfigure("threads")),showMsg=FALSE)
            #              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=private$paramlist[["adapter1"]],
            #                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
            #                                    paramlist=private$paramlist[["paramList"]])
            if(length(paramList)>0){
                remove_adapters(file1 = input(.Object)[["fastqInput1"]],
                                paramList,
                                adapter1 = output(.Object)$adapter1,
                                output1 = output(.Object)[["fastqOutput1"]],
                                basename = param(.Object)$reportPrefix,
                                interleaved = FALSE,
                                overwrite = TRUE)
            }else{
                remove_adapters(file1 = input(.Object)[["fastqInput1"]],
                                adapter1 = output(.Object)[["adapter1"]],
                                output1 = output(.Object)[["fastqOutput1"]],
                                basename = param(.Object)[["reportPrefix"]],
                                interleaved = FALSE,
                                overwrite = TRUE)
            }

        }else if(param(.Object)[["interleave"]]){
            adapter1<-param(.Object)[["adapter1"]]
            adapter2<-param(.Object)[["adapter2"]]
            if(length(paramList)>0){
                remove_adapters(file1 = input(.Object)[["fastqInput1"]],
                                paramList,
                                adapter1 = adapter1,
                                output1 = output(.Object)[["fastqOutput1"]],
                                adapter2 = adapter2,
                                basename = param(.Object)[["reportPrefix"]],
                                interleaved = TRUE,
                                overwrite = TRUE)
            }else{
                remove_adapters(file1 = input(.Object)[["fastqInput1"]],
                                adapter1 = adapter1,
                                output1 = output(.Object)[["fastqOutput1"]],
                                adapter2 = adapter2,
                                basename = param(.Object)[["reportPrefix"]],
                                interleaved = TRUE,
                                overwrite = TRUE)
            }


        }else{
            adapter1<-param(.Object)[["adapter1"]]
            adapter2<-param(.Object)[["adapter2"]]
            writeLog(.Object,"begin to remove adapter")
            writeLog(.Object,"source:",showMsg=FALSE)
            writeLog(.Object,input(.Object)[["fastqInput1"]],showMsg=FALSE)
            writeLog(.Object,input(.Object)[["fastqInput2"]],showMsg=FALSE)
            writeLog(.Object,paste0("Adapter1:",adapter1))
            writeLog(.Object,paste0("Adapter2:",adapter2))
            writeLog(.Object,"Destination:",showMsg=FALSE)
            writeLog(.Object,output(.Object)[["fastqOutput1"]],showMsg=FALSE)
            writeLog(.Object,output(.Object)[["fastqOutput2"]],showMsg=FALSE)
            writeLog(.Object,param(.Object)[["reportPrefix"]],showMsg=FALSE)
            writeLog(.Object,paste0("Threads:",.obtainConfigure("threads")),showMsg=FALSE)
            #              .remove_adapters_call(inputFile1=private$paramlist[["fastqInput1"]],adapter1=adapter1,
            #                                    outputFile1 = private$paramlist[["fastqOutput1"]],basename = private$paramlist[["reportPrefix"]],
            #                                    inputFile2=private$paramlist[["fastqInput2"]],adapter2=adapter2,
            #                                    outputFile2 = private$paramlist[["fastqOutput2"]],paramlist=private$paramlist[["paramList"]])
            if(length(paramList)>0){
                remove_adapters(file1 = input(.Object)[["fastqInput1"]],
                                paramList,
                                adapter1 = adapter1,
                                output1 = output(.Object)[["fastqOutput1"]],
                                file2 = input(.Object)[["fastqInput2"]],
                                adapter2 = adapter2,
                                output2 = output(.Object)[["fastqOutput2"]],
                                basename = param(.Object)[["reportPrefix"]],
                                interleaved = FALSE,
                                overwrite = TRUE)
            }else{
                remove_adapters(file1 = input(.Object)[["fastqInput1"]],
                                adapter1 = adapter1,
                                output1 = output(.Object)[["fastqOutput1"]],
                                file2 = input(.Object)[["fastqInput2"]],
                                adapter2 = adapter2,
                                output2 = output(.Object)[["fastqOutput2"]],
                                basename = param(.Object)[["reportPrefix"]],
                                interleaved = FALSE,
                                overwrite = TRUE)
            }

        }
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "RemoveAdapter",
    definition = function(.Object,...){
        if(is.null(input(.Object)[["fastqInput1"]])){
            stop("'fastqInput1' is required.")
        }
        if(is.null(param(.Object)[["adapter1"]])){
            stop("Parameter 'adapter1' is requied")
        }
        if(is.null(param(.Object)[["adapter2"]])){
            stop("Parameter 'adapter2' is requied")
        }
        if(param(.Object)[["interleave"]]&&param(.Object)$singleEnd){
            stop("Single end data should not be interleave")
        }
    }
)


setMethod(
    f = "genReport",
    signature = "RemoveAdapter",
    definition = function(.Object,...){
        #adapterslist
        tblist <- getTopic(.Object,"\\[Adapter sequences\\]")
        splitlist <- strsplit(tblist,": ")
        if(.Object@singleEnd){
            report(.Object)[["adapterslist"]] <- (data.frame(adapter=c("adapter for single end data"),sequence=c(splitlist[[1]][2])))
        }else{
            report(.Object)[["adapterslist"]] <- (data.frame(adapter=c("adapter for paired end data mate 1","adapter for paired end data mate 2"),
                                                             sequence=c(splitlist[[1]][2],splitlist[[2]][2])))
            #return(list(adapter1=splitlist[[1]][2],
            #            adapter2=splitlist[[2]][2]))
        }   

        # adapters
        tblist <- getTopic(.Object,"\\[Adapter sequences\\]")
        splitlist <- strsplit(tblist,": ")
        if(.Object@singleEnd){
            report(.Object)[["adapters"]] <- (listToFrame(list(adpater1=splitlist[[1]][2])))
        }else{
            report(.Object)[["adapters"]] <- (listToFrame(list(adapter1=splitlist[[1]][2],
                                                               adapter2=splitlist[[2]][2])))
        }
            
        # adapter1 adapter2
        tblist <- getTopic(.Object,"\\[Adapter sequences\\]")
        splitlist <- strsplit(tblist,": ")
        report(.Object)[["adapter1"]] <- (splitlist[[1]][2])
        report(.Object)[["adapter2"]] <- (splitlist[[2]][2])

        #settingslist
        tblist <- getTopic(.Object,"\\[Adapter trimming\\]")
        splitlist <- strsplit(tblist,": ")
        lst <- list()
        for(i in 1:length(tblist)){
            lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
        }
        report(.Object)[["settingslist"]] <- (lst)
        
        # settings
        tblist <- getTopic(.Object,"\\[Adapter trimming\\]")
        splitlist <- strsplit(tblist,": ")
        lst <- list()
        for(i in 1:length(tblist)){
            lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
        }
        report(.Object)[["settings"]] <- (listToFrame(lst))
       
        #statisticslist
        tblist <- getTopic(.Object,"\\[Trimming statistics\\]")
        splitlist <- strsplit(tblist,": ")
        lst <- list()
        for(i in 1:length(tblist)){
            lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
        }
        report(.Object)[["statisticslist"]] <- (lst)
        
        #statistics
        tblist <- getTopic(.Object,"\\[Trimming statistics\\]")
        splitlist <- strsplit(tblist,": ")
        lst <- list()
        for(i in 1:length(tblist)){
            lst[[splitlist[[i]][1]]]<-splitlist[[i]][2]
        }
        report(.Object)[["statistics"]] <- (listToFrame(lst))
        
        #distribution
        tblist <- getTopic(.Object,"\\[Length distribution\\]")
        splitlist <- strsplit(tblist,"\t")
        colkey <- splitlist[[1]]
        tbdt <- NULL
        for(i in 2:length(tblist)){
            tbdt <- c(tbdt,splitlist[[i]])
        }
        tbdt<-as.integer(tbdt)
        if(isSingleEnd(.Object)){
            colsize<-4
        }else{
            colsize<-6
        }
        df<-as.data.frame(matrix(tbdt,length(tblist)-1,colsize,TRUE))
        colnames(df) <- colkey
        report(.Object)[["distribution"]] <- (df)
        .Object
    }
)

setGeneric(
    name = "getTopic",
    def = function(.Object,topic,...){
        standardGeneric("getTopic")
    }
)

setMethod(
    f = "getTopic",
    signature = "RemoveAdapter",
    definition = function(.Object, topic,...){
        setLine<-readLines(paste0(param(.Object)[["reportPrefix"]],".settings"))
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
    }
)



listToFrame <- function(a){
    return(data.frame(Item=names(a),Value=as.character(a)))
}




#' @name RemoveAdapter
#' @title Use AdapterRemoval to remove adapters
#' @description
#' Use AdapterRemoval to remove adapters
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
#' @param threads \code{Numeric}. The max threads allowed to be used by this step.
#' Default: getThreads()
#' @param paramList Additional arguments to be passed on to the binaries
#' for removing adapter. See below for details.
#' @param findParamList Additional arguments to be passed on to the binaries
#' for identifying adapter. See below for details.
#' @param reportPrefix \code{Character}. The prefix of report files path.
#' Default: generate from known parameters
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{removeAdapter} instead.
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
#' atacRenamer %>% atacFindAdapter %>% atacRemoveAdapter
#'
#' dir(td)
#' @importFrom Rbowtie2 identify_adapters
#' @importFrom Rbowtie2 remove_adapters



setGeneric("atacRemoveAdapter",function(atacProc,adapter1=NULL,adapter2=NULL,
                                        fastqOutput1=NULL,reportPrefix=NULL,
                                        fastqOutput2=NULL,fastqInput1=NULL,
                                        fastqInput2=NULL,interleave=FALSE, threads = getThreads(),
                                        paramList= NULL,findParamList=NULL, ...) standardGeneric("atacRemoveAdapter"))
#' @rdname RemoveAdapter
#' @aliases atacRemoveAdapter
#' @export
setMethod(
    f = "atacRemoveAdapter",
    signature = "ATACProc",
    definition = function(atacProc,adapter1=NULL,adapter2=NULL,
                          fastqOutput1=NULL,reportPrefix=NULL,
                          fastqOutput2=NULL,fastqInput1=NULL,
                          fastqInput2=NULL,interleave=FALSE, threads = getThreads(),
                          paramList= NULL,findParamList=NULL, ...){
        allpara <- c(list(Class = "RemoveAdapter", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname RemoveAdapter
#' @aliases removeAdapter
#' @export
removeAdapter <- function(fastqInput1, fastqInput2,
                              adapter1, adapter2,
                              fastqOutput1=NULL,reportPrefix=NULL,
                              fastqOutput2=NULL,
                              interleave=FALSE, threads = getThreads(),
                              paramList = NULL,findParamList = NULL, ...){
    allpara <- c(list(Class = "RemoveAdapter", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}

