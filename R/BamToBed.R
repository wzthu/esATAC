setClass(Class = "BamToBed",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "BamToBed",
    definition = function(.Object, atacProc, ..., bamInput = NULL, bedOutput = NULL, editable = FALSE){
        .Object <- init(.Object,"BamToBed",editable,list(arg1 = atacProc))

        if(!is.null(atacProc)){
            .Object@paramlist[["bamInput"]] <- getParam(atacProc, "bamOutput")
            regexProcName <- sprintf("(bam|%s)", getProcName(atacProc))
        }else{
            .Object@paramlist[["bamInput"]] <- bamInput
            regexProcName <- "(bam)"
        }

        if(is.null(bedOutput)){
            prefix <- getBasenamePrefix(.Object, .Object@paramlist[["bamInput"]], regexProcName)
            .Object@paramlist[["bedOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                          paste(prefix, ".", getProcName(.Object), ".bed", sep = ""))
        }else{
            name_split <- unlist(base::strsplit(x = bedOutput, split = ".", fixed = TRUE))
            suffix <- tail(name_split, 1)
            if(suffix == "bed"){
                .Object@paramlist[["bedOutput"]] <- bedOutput
            }else{
                .Object@paramlist[["bedOutput"]] <- paste(bedOutput, ".bed", sep = "")
            }
        }
        # parameter check and return
        paramValidation(.Object)
        .Object
    } # definition end
) # setMethod initialize end



setMethod(
    f = "processing",
    signature = "BamToBed",
    definition = function(.Object, ...){
        .Object <- writeLog(.Object, paste0("processing file:"))
        .Object <- writeLog(.Object, sprintf("Bam Input Source:%s", .Object@paramlist[["bamInput"]]))
        .Object <- writeLog(.Object, sprintf("Bed Output Destination:%s", .Object@paramlist[["bedOutput"]]))
        rtracklayer::export(con = .Object@paramlist[["bedOutput"]],
                            object = rtracklayer::import(con = .Object@paramlist[["bamInput"]], format = "bam", paired = TRUE, use.names = TRUE),
                            format = "bed")
        .Object
    }
)




setMethod(
    f = "checkRequireParam",
    signature = "BamToBed",
    definition = function(.Object, ...){
        if(is.null(.Object@paramlist[["bamInput"]])){
            stop("bamInput is required.")
        }
    }
)

setMethod(
    f = "checkAllPath",
    signature = "BamToBed",
    definition = function(.Object, ...){
        checkFileExist(.Object, .Object@paramlist[["bamInput"]])
        checkFileCreatable(.Object, .Object@paramlist[["bedOutput"]])
    }
)


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
#' \link[rtracklayer]{import}
#' \link[rtracklayer]{export}
#' @importFrom rtracklayer export
#' @name atacBam2Bed
#' @export
#' @docType methods
#' @rdname atacBam2Bed-methods
setGeneric("atacBam2Bed", function(atacProc, bamInput = NULL, bedOutput = NULL) standardGeneric("atacBam2Bed"))

#' @rdname atacBam2Bed-methods
#' @aliases atacBam2Bed
setMethod(
    f = "atacBam2Bed",
    signature = "ATACProc",
    definition = function(atacProc, bamInput = NULL, bedOutput = NULL){
        atacproc <- new("BamToBed", atacProc = atacProc, bamInput = bamInput, bedOutput = bedOutput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname atacBam2Bed-methods
#' @export
bam2bed <- function(bamInput, bedOutput = NULL){
    atacproc <- new("BamToBed", atacProc = NULL, bamInput = bamInput, bedOutput = bedOutput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
