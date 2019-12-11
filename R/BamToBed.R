setClass(Class = "BamToBed",
         contains = "ATACProc"
)

setMethod(
    f = "init",
    signature = "BamToBed",
    definition = function(.Object, prevSteps = list(),...){
        allparam <- list(...)
        bamInput <- allparam[["bamInput"]]
        bedOutput <- allparam[["bedOutput"]]
       
        if(length(prevSteps) > 0){
            prevSteps <- prevSteps[[1]]
            input(.Object)[["bamInput"]] <- output(prevSteps)[["bamOutput"]]
        }else{
            input(.Object)[["bamInput"]] <- bamInput
        }

        if(is.null(bedOutput)){
            output(.Object)[["bedOutput"]] <- getAutoPath(.Object, input(.Object)[["bamInput"]],"bam","bed")
        }else{
            output(.Object)[["bedOutput"]] <- addFileSuffix(bedOutput, ".bed")
        }

        .Object
    } # definition end
) # setMethod initialize end



setMethod(
    f = "processing",
    signature = "BamToBed",
    definition = function(.Object, ...){
        rtracklayer::export(con = output(.Object)[["bedOutput"]],
                            object = rtracklayer::import(con = input(.Object)[["bamInput"]], format = "bam", paired = TRUE, use.names = TRUE),
                            format = "bed")
        .Object
    }
)





#' @name BamToBed
#' @title Convert bam format to bed format.
#' @description
#' This function convert a bam file into a bed file.
#' Note:bed file is 0-based.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacBamSort}},
#' \code{\link{atacSam2Bam}}.
#' @param bamInput \code{Character} scalar.
#' Bam file input path.
#' @param bedOutput \code{Character} scalar.
#' Bed file output path. If ignored, bed file will be put in the same path as
#' the bam file.
#' @param ... Additional arguments, currently unused.
#' @details The bam file wiil be automatically obtained from
#' object(\code{atacProc}) or input by hand. Output can be ignored.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(Rsamtools)
#' ex1_file <- system.file("extdata", "ex1.bam", package="Rsamtools")
#' bam2bed(bamInput = ex1_file)
#'
#' @seealso
#' \code{\link{atacBamSort}}
#' \code{\link{atacSam2Bam}}
#' @importFrom rtracklayer export



setGeneric("atacBam2Bed", function(atacProc, bamInput = NULL, bedOutput = NULL, ...) standardGeneric("atacBam2Bed"))

#' @rdname BamToBed
#' @aliases atacBam2Bed
#' @export
setMethod(
    f = "atacBam2Bed",
    signature = "ATACProc",

    definition = function(atacProc, bamInput = NULL, bedOutput = NULL, ...){
        allpara <- c(list(Class = "BamToBed", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)

#' @rdname BamToBed
#' @aliases bam2bed
#' @export

bam2bed <- function(bamInput, bedOutput = NULL, ...){
    allpara <- c(list(Class = "BamToBed", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
