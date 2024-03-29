setClass(Class = "SCRenamer",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "SCRenamer",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        fastqOutput1 <- allparam[["fastqOutput1"]]
        fastqOutput2 <- allparam[["fastqOutput2"]]
        fastqBarcodeInput <- allparam[["fastqBarcodeInput"]]
        fastqInput1 <- allparam[["fastqInput1"]]
        fastqInput2 <- allparam[["fastqInput2"]]
        threads <- allparam[["threads"]]
        if(length(prevSteps)>0){
            atacProc <- prevSteps[[1]]
            input(.Object)[["fastqInput1"]] <- output(atacProc)[["fastqOutput1"]]
            input(.Object)[["fastqBarcodeInput"]] <- output(atacProc)[["fastqBarcodeOutput"]]
            input(.Object)[["fastqInput2"]] <- output(atacProc)[["fastqOutput2"]]
        }else{
            param(.Object)[["interleave"]] <- interleave
            if(is.null(fastqInput2)){
                property(.Object)[["interleave"]] <-TRUE
            }else{
                property(.Object)[["interleave"]]<-FALSE
            }
        }

        if(!is.null(fastqInput1)){
            input(.Object)[["fastqInput1"]] <- fastqInput1;
        }
        if(!is.null(fastqInput2)){
            input(.Object)[["fastqInput2"]] <- fastqInput2;
        }
        if(!is.null(fastqBarcodeInput)){
            input(.Object)[["fastqBarcodeInput"]] <- fastqBarcodeInput;
        }


        if(is.null(fastqOutput1)){
            if(!is.null(input(.Object)[["fastqInput1"]])){
                output(.Object)[["fastqOutput1"]] <- getAutoPath(.Object,originPath = input(.Object)[["fastqInput1"]],"fastq|fq","fq")
            }
        }else{
            output(.Object)[["fastqOutput1"]] <- fastqOutput1
        }
        if(is.null(fastqOutput2)){
            if(!is.null(input(.Object)[["fastqInput2"]])){
                output(.Object)[["fastqOutput2"]] <- getAutoPath(.Object,originPath = input(.Object)[["fastqInput2"]],"fastq|fq","fq")
            }
        }else{
            output(.Object)[["fastqOutput2"]] <- fastqOutput2
        }
        stopifnot(is.numeric(threads))
        param(.Object)[["threads"]] <- threads
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "SCRenamer",
    definition = function(.Object,...){
        writeLog(.Object,paste0("processing file:"))
        writeLog(.Object,sprintf("source:%s",input(.Object)[["fastqInput1"]]))
        writeLog(.Object,sprintf("destination:%s",output(.Object)[["fastqOutput1"]]))
        threads <- param(.Object)[["threads"]]
      
        if(threads>=2){
            writeLog(.Object,paste0("processing file:"))
            writeLog(.Object,sprintf("source:%s",input(.Object)[["fastqInput2"]]))
            writeLog(.Object,sprintf("destination:%s",output(.Object)[["fastqOutput2"]]))
            cl <- makeCluster(2)
            parLapply(cl = cl, X = 1:2, fun = scSingleCall,.renamer_call=.sc_renamer_call, .Object=.Object)
            stopCluster(cl)
#            mclapply(X = 1:2,FUN = singleCall,.renamer_call=.renamer_call, .Object=.Object, mc.cores=2)
        }else{
            scSingleCall(1,.renamer_call=.sc_renamer_call,.Object=.Object)
            writeLog(.Object,paste0("processing file:"))
            writeLog(.Object,sprintf("source:%s",input(.Object)[["fastqInput2"]]))
            writeLog(.Object,sprintf("destination:%s",output(.Object)[["fastqOutput2"]]))
            scSingleCall(2,.renamer_call=.sc_renamer_call,.Object=.Object)
        }
        .Object
    }
)


rename_fq <- function(inputfq,barcodefq,outputfq){
    f <- file(inputfq,'r')
    b <- file(barcodefq, 'r')
    if(file.exists(outputfq)){
        file.remove(outputfq)
    }
    line_size <- 10000000
    while(TRUE){
        lines <- readLines(f, n=line_size)
        blines <- readLines(b, n=line_size)
        if(length(lines)==0){
            break
        }
        sz <- length(lines)
        lines[seq(1,min(line_size,sz),4)] <- paste0('@',blines[seq(2,min(line_size,sz),4)],':',substring(lines[seq(1,min(line_size,sz),4)],2))
        write(lines,file = outputfq, append = TRUE, sep = "\n")
    }
    close(f)
    close(b)
}

scSingleCall<-function(number,.renamer_call,.Object){
    if(number==1){
        rename_fq(inputfq=input(.Object)[["fastqInput1"]],
              barcodefq=input(.Object)[["fastqBarcodeInput"]],
              outputfq=output(.Object)[["fastqOutput1"]])
    }else if(number==2){
        rename_fq(inputfq=input(.Object)[["fastqInput2"]],
              barcodefq=input(.Object)[["fastqBarcodeInput"]],
              outputfq=output(.Object)[["fastqOutput2"]])
    }
}

setMethod(
  f = "genReport",
  signature = "SCRenamer",
  definition = function(.Object, ...){
    .Object
  }
)






#' @title Rename reads name in fastq
#' @description
#' Rename reads name in fastq with increasing integer
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacUnzipAndMerge}}
#' \code{\link{unzipAndMerge}}
#' @param fastqInput1 \code{Character} scalar. For single-end sequencing,
#' it contains sequence file paths.
#' For paired-end sequencing, it can be file path with #1 mates paired
#' with file path in file2
#' And it can also be interleaved file paths when argument
#' interleave=\code{TRUE}
#' @param fastqBarcodeInput \code{Character} scalar. It contains file path with barcode.
#' @param fastqInput2 \code{Character} scalar. It contains file path with #2
#' mates paired with file paths in fastqInput1
#' For single-end sequencing files and interleaved paired-end sequencing
#' files(argument interleaved=\code{TRUE}),
#' it must be \code{NULL}.
#' @param fastqOutput1 \code{Character} scalar.
#' The output file path of renamed fastqInput1.
#' @param fastqOutput2 \code{Character} scalar.
#' The output file path of renamed fastqInput2.
#' @param threads \code{Integer} scalar.
#' The threads will be created in this process. default: 1
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{renamer} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacUnzipAndMerge}}
#' \code{\link{unzipAndMerge}}
#' \code{\link{atacQCReport}}
#' \code{\link{atacRemoveAdapter}}
#' @examples
#' 
#' ignoreCheck() # warnning: run this for fast test only
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
#' atacRenamer
#'
#' dir(td)
#'
#' @importFrom parallel makeCluster parLapply stopCluster mclapply

#' @name SCRenamer


setGeneric("atacSCRenamer",function(atacProc,fastqOutput1=NULL,
                                  fastqOutput2=NULL,
                                  fastqInput1=NULL,
                                  fastqBarcodeInput = NULL,
                                  fastqInput2=NULL,
                                  interleave = FALSE, 
                                  threads = getThreads(), 
                                  ...) standardGeneric("atacSCRenamer"))

#' @rdname SCRenamer
#' @aliases atacSCRenamer
#' @export
setMethod(
    f = "atacSCRenamer",
    signature = "ATACProc",
    definition = function(atacProc,
             fastqOutput1=NULL,
             fastqOutput2=NULL,
             fastqInput1=NULL,
             fastqBarcodeInput = NULL,
             fastqInput2=NULL,
             interleave = FALSE, 
             threads = getThreads(), ...){
        allpara <- c(list(Class = "SCRenamer", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname SCRenamer
#' @aliases scRenamer
#' @export
scRenamer <- function(fastqInput1=NULL,
                    fastqInput2=NULL,
                    fastqOutput1=NULL,
                    fastqBarcodeInput = NULL,
                    fastqOutput2=NULL,
                    interleave = FALSE,
                    threads = getThreads(), ...){
    allpara <- c(list(Class = "SCRenamer", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}

