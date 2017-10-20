setClass(Class = "Renamer",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "Renamer",
    definition = function(.Object,atacProc,..., fastqOutput1=NULL, fastqOutput2=NULL,
                          fastqInput1=NULL, fastqInput2=NULL, 
                          interleave = FALSE, threads = NULL, editable=FALSE){
        .Object <- init(.Object,"Renamer",editable,list(arg1=atacProc))
        if(!is.null(atacProc)){
            .Object@paramlist[["fastqInput1"]] <- getParam(atacProc,"fastqOutput1");
            .Object@paramlist[["fastqInput2"]] <- getParam(atacProc,"fastqOutput2");
            regexProcName<-sprintf("(fastq|fq|%s)",getProcName(atacProc))
            .Object@paramlist[["interleave"]] <- getParam(atacProc,"interleave")
        }else{
            regexProcName<-"(fastq|fq)"
            .Object@paramlist[["interleave"]] <- interleave
            if(is.null(fastqInput2)){
                .Object@singleEnd<-TRUE
            }else{
                .Object@singleEnd<-FALSE
            }
        }
        
        if(!is.null(fastqInput1)){
            .Object@paramlist[["fastqInput1"]] <- fastqInput1;
        }
        if(!is.null(fastqInput2)){
            .Object@paramlist[["fastqInput2"]] <- fastqInput2;
        }
        
        
        if(is.null(fastqOutput1)){
            if(!is.null(.Object@paramlist[["fastqInput1"]])){
                prefix<-getBasenamePrefix(.Object,.Object@paramlist[["fastqInput1"]],regexProcName)
                .Object@paramlist[["fastqOutput1"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".fq"))
            }
        }else{
            .Object@paramlist[["fastqOutput1"]] <- fastqOutput1;
        }
        if(is.null(fastqOutput2)){
            if(!is.null(.Object@paramlist[["fastqInput2"]])){
                prefix<-getBasenamePrefix(.Object,.Object@paramlist[["fastqInput2"]],regexProcName)
                .Object@paramlist[["fastqOutput2"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".fq"));
            }
        }else{
            .Object@paramlist[["fastqOutput2"]] <- fastqOutput2;
        }
        if(!is.null(threads)){
            .Object@paramlist[["threads"]] <- as.integer(threads)
        }
        
        paramValidation(.Object)
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "Renamer",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["fastqInput1"]]))
        .Object <- writeLog(.Object,sprintf("destination:%s",.Object@paramlist[["fastqOutput1"]]))
        threads <- .obtainConfigure("threads")
        if(!is.null(.Object@paramlist[["threads"]])){
            threads <- .Object@paramlist[["threads"]]
        }
        if(.Object@singleEnd||.Object@paramlist[["interleave"]]){
            singleCall(number=1,.Object=.Object)
        }else if(threads>=2){
            .Object <- writeLog(.Object,paste0("processing file:"))
            .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["fastqInput2"]]))
            .Object <- writeLog(.Object,sprintf("destination:%s",.Object@paramlist[["fastqOutput2"]]))
            cl<-makeCluster(2)
            parLapply(cl = cl,X = 1:2,fun = singleCall, .Object=.Object)
            stopCluster(cl)
        }else{
            singleCall(1,.Object=.Object)
            .Object <- writeLog(.Object,paste0("processing file:"))
            .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["fastqInput2"]]))
            .Object <- writeLog(.Object,sprintf("destination:%s",.Object@paramlist[["fastqOutput2"]]))
            singleCall(2)
        }
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "Renamer",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["fastqInput1"]])){
            stop("fastqInput1 is required.")
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "Renamer",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["fastqInput1"]]);
        checkFileExist(.Object,.Object@paramlist[["fastqInput2"]]);
        checkFileCreatable(.Object,.Object@paramlist[["fastqOutput1"]]);
        checkFileCreatable(.Object,.Object@paramlist[["fastqOutput2"]]);
    }
)

singleCall<-function(number,.Object){
    if(number==1){
        .renamer_call(inputFile = .Object@paramlist[["fastqInput1"]],
                      outputFile = .Object@paramlist[["fastqOutput1"]],
                      interleave = .Object@paramlist[["interleave"]])
    }else if(number==2){
        .renamer_call(inputFile = .Object@paramlist[["fastqInput2"]],
                      outputFile = .Object@paramlist[["fastqOutput2"]], 
                      interleave = .Object@paramlist[["interleave"]])
    }
}



Renamer <-R6Class(
  classname = "Renamer",
  inherit = ATACProc,
  public = list(
    initialize = function(atacProc, fastqOutput1=NULL, fastqOutput2=NULL,
                          fastqInput1=NULL, fastqInput2=NULL, 
                          interleave = FALSE, threads = NULL, editable=FALSE){
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
        if(!is.null(threads)){
            private$paramlist[["threads"]] <- as.integer(threads)
        }

      private$paramValidation()

    }
  ),
  private = list(
    processing = function(){
        private$writeLog(paste0("processing file:"))
        private$writeLog(sprintf("source:%s",private$paramlist[["fastqInput1"]]))
        private$writeLog(sprintf("destination:%s",private$paramlist[["fastqOutput1"]]))
        threads <- .obtainConfigure("threads")
        if(!is.null(private$paramlist[["threads"]])){
             threads <- private$paramlist[["threads"]]
        }
        if(private$singleEnd||private$paramlist[["interleave"]]){
            private$singleCall(1)
        }else if(threads>=2){
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

#' @name atacRenamer
#' @aliases atacRenamer
#' @aliases renamer
#' @title Rename reads name in fastq 
#' @description 
#' Rename reads name in fastq with increasing integer  
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacUnzipAndMerge}} 
#' \code{\link{unzipAndMerge}}
#' @param fastqInput1 \code{Character} scalar. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file path with #1 mates paired
#' with file path in file2
#' And it can also be interleaved file paths when argument
#' interleave=\code{TRUE}
#' @param fastqInput2 \code{Character} scalar. It contains file path with #2
#' mates paired with file paths in fastqInput1
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param fastqOutput1 \code{Character} scalar. 
#' The output file path of renamed fastqInput1.
#' @param fastqOutput2 \code{Character} scalar. 
#' The output file path of renamed fastqInput2.
#' @param interleave \code{Character} scalar. 
#' Set \code{TRUE} when files are
#' interleaved paired-end sequencing data.
#' @param threads \code{Integer} scalar. 
#' The threads will be created in this process. default: 1
#' @details The parameter related to input and output file path
#' will be automatically 
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{renamer} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacUnzipAndMerge}} 
#' \code{\link{unzipAndMerge}}
#' \code{\link{atacQCReport}} 
#' \code{\link{atacRemoveAdapter}} 
#' @examples 
#' library(magrittr)
#' td <- tempdir()
#' setConfigure("tmpdir",td)
#' 
#' # Identify adapters
#' prefix<-system.file(package="ATACpipe", "extdata", "uzmg")
#' (reads_1 <-file.path(prefix,"m1",dir(file.path(prefix,"m1"))))
#' (reads_2 <-file.path(prefix,"m2",dir(file.path(prefix,"m2"))))
#' 
#' reads_merged_1 <- file.path(td,"reads1.fastq")
#' reads_merged_2 <- file.path(td,"reads2.fastq")
#' atacproc <- 
#' atacUnzipAndMerge(fastqInput1 = reads_1,fastqInput2 = reads_2) %>%
#' atacRenamer
#' 
#' dir(td) 
#' 


#' @rdname atacRenamer
#' @exportMethod atacRenamer
setGeneric("atacRenamer",function(atacProc,fastqOutput1=NULL,
                                  fastqOutput2=NULL,
                                  fastqInput1=NULL, 
                                  fastqInput2=NULL, 
                                  interleave = FALSE) standardGeneric("atacRenamer")) 
setMethod(
    f = "atacRenamer",
    signature = "ATACProc",
    definition = function(atacProc, 
             fastqOutput1=NULL,
             fastqOutput2=NULL,
             fastqInput1=NULL, 
             fastqInput2=NULL, 
             interleave = FALSE){
        atacproc <- new(
            "Renamer",
            atacProc = atacProc,
            fastqOutput1 = fastqOutput1,
            fastqOutput2 = fastqOutput2,
            fastqInput1 = fastqInput1,
            fastqInput2 = fastqInput2)
        atacproc <- process(atacproc, atacProc)
        invisible(atacproc)
    }
)
#' @rdname atacRenamer
#' @export 
renamer <- function(fastqInput1=NULL, 
                    fastqInput2=NULL, 
                    fastqOutput1=NULL,
                    fastqOutput2=NULL,
                    interleave = FALSE,
                    threads = NULL){
    atacproc <- Renamer$new(atacProc = NULL,
                            fastqOutput1 = fastqOutput1,
                            fastqOutput2 = fastqOutput2,
                            fastqInput1 = fastqInput1,
                            fastqInput2 = fastqInput2,
                            threads = NULL)
    atacproc$process()
    invisible(atacproc)
}

