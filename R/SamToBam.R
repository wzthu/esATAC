setClass(Class = "SamToBam",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "SamToBam",
    definition = function(.Object,prevSteps = list(),...){
        atacProc <- NULL
        if(length(prevSteps)>0){
            atacProc <- prevSteps[[1]]
        }
        allparam <- list(...)
        samInput <- allparam[["samInput"]]
        bamOutput <- allparam[["bamOutput"]]
        isSort <- allparam[["isSort"]]
        
        # necessary parameters
        if((!is.null(atacProc)) ){
            input(.Object)[["samInput"]] <- getParam(atacProc, "samOutput")
        }else if(is.null(atacProc)){ # input
            input(.Object)[["samInput"]] <- samInput
        }
        # unnecessary parameters
        if(is.null(bamOutput)){
            param(.Object)[["name_tmp"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","" )
        
            output(.Object)[["bamOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","bam")
            if(isSort){
                output(.Object)[["baiOutput"]] <- getAutoPath(.Object, input(.Object)[["samInput"]],"sam|SAM","bam.bai")
            }
        }else{
            bamOutput <- addFileSuffix(bamOutput,".bam")
            param(.Object)[["name_tmp"]] <- substring(bamOutput,first = 1,last = nchar(bamOutput) - 4)
            output(.Object)[["bamOutput"]] <- bamOutput
            if(isSort){
                output(.Object)[["baiOutput"]] <- addFileSuffix(bamOutput,".bai")
            }
            
        }
        param(.Object)[['isSort']] <- isSort
        .Object

    }
)




setMethod(
    f = "processing",
    signature = "SamToBam",
    definition = function(.Object,...){
        
        Rsamtools::asBam(file = input(.Object)[["samInput"]], 
                         destination = param(.Object)[["name_tmp"]],
                         overwrite = TRUE, indexDestination = param(.Object)[['isSort']])
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "SamToBam",
    definition = function(.Object, ...){
        .Object
    }
)


#' @name SamToBam
#' @title Convert sam format to bam format.
#' @description
#' This function convert a sam file into a bam file.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacBowtie2Mapping}}.
#' @param samInput \code{Character} scalar.
#' Sam file input path.
#' @param bamOutput \code{Character} scalar.
#' Bam file output path. If ignored, bed file will be put in the same path as
#' the sam file.
#' @param isSort \code{Logical} scalar.
#' Sort bam.
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{bamToBed} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' sam_bz <- system.file("extdata", "Example.sam.bz2", package="esATAC")
#' sam_path <- as.vector(bunzip2(filename = sam_bz,
#' destname = file.path(getwd(), "Example.sam"),
#' ext="bz2", FUN=bzfile, remove = FALSE))
#' sam2bam(samInput = sam_path)
#'
#' @seealso
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{atacBam2Bed}}
#' \code{\link{atacBamSort}}


setGeneric("atacSam2Bam",function(atacProc,
                                  samInput = NULL, bamOutput = NULL, isSort=TRUE, ...) standardGeneric("atacSam2Bam"))
#' @rdname SamToBam
#' @aliases atacSam2Bam
#' @export
setMethod(
    f = "atacSam2Bam",
    signature = "ATACProc",
    definition = function(atacProc,
                          samInput = NULL, bamOutput = NULL, isSort=TRUE, ...){
        allpara <- c(list(Class = "SamToBam", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname SamToBam
#' @aliases sam2bam
#' @export
sam2bam <- function(samInput, bamOutput = NULL, isSort=TRUE, ...){
    allpara <- c(list(Class = "SamToBam", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
