setClass(Class = "Rsortbam",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "Rsortbam",
    definition = function(.Object,prevSteps = list(), ...){
        atacProc <- NULL
        if(length(prevSteps)>0){
            atacProc <- prevSteps[[1]]
        }
        allparam <- list(...)
        bamInput <- allparam[["bamInput"]]
        bamOutput <- allparam[["bamOutput"]]
        # necessary parameters
        if(!is.null(atacProc)){ # atacproc from SamToBam
            input(.Object)[["bamInput"]] <- getParam(atacProc, "bamOutput")
        }else if(is.null(atacProc)){ # input
            input(.Object)[["bamInput"]] <- bamInput
        }
        # unnecessary parameters
        if(is.null(bamOutput)){
            param(.Object)[["bamOutput_tmp"]] <- getAutoPath(.Object, input(.Object)[["bamInput"]],"bam", "")
            output(.Object)[["bamOutput"]] <- getAutoPath(.Object, input(.Object)[["bamInput"]],"bam", ".bam")
        }else{
            bamOutput <- addFileSuffix(bamOutput,".bam")
            param(.Object)[["bamOutput_tmp"]] <- substring(bamOutput,first = 1,last = nchar(bamOutput) - 4)
            output(.Object)[["bamOutput"]] <- bamOutput
            
        }
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "Rsortbam",
    definition = function(.Object,...){
        Rsamtools::sortBam(file = input(.Object)[["bamInput"]], destination = param(.Object)[["bamOutput_tmp"]])
        Rsamtools::indexBam(files = output(.Object)[["bamOutput"]])
        .Object
    }
)



setMethod(
    f = "genReport",
    signature = "Rsortbam",
    definition = function(.Object, ...){
        .Object
    }
)

#' @name Rsortbam
#' @title Sort bam file and rebuild bai index.
#' @description
#' Sort bamfile and build index.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSam2Bam}}.
#' @param bamInput \code{Character} scalar.
#' Input bam file path.
#' @param bamOutput \code{Character} scalar.
#' Output bam file path.
#' @param ... Additional arguments, currently unused.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(Rsamtools)
#' ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' bamsort(bamInput = ex1_file)
#'
#' @seealso
#' \code{\link{atacSam2Bam}}
#' \code{\link{atacBam2Bed}}


setGeneric("atacBamSort",function(atacProc,
                                  bamInput = NULL, bamOutput = NULL, ...) standardGeneric("atacBamSort"))

#' @rdname Rsortbam
#' @aliases atacBamSort
#' @export
setMethod(
    f = "atacBamSort",
    signature = "ATACProc",
    definition = function(atacProc,
                          bamInput = NULL, bamOutput = NULL, ...){
        allpara <- c(list(Class = "Rsortbam", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname Rsortbam
#' @aliases bamsort
#' @export
bamsort <- function(bamInput = NULL, bamOutput = NULL, ...){
    allpara <- c(list(Class = "Rsortbam", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
