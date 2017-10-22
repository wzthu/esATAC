setClass(Class = "Rsortbam",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "Rsortbam",
    definition = function(.Object,atacProc, ..., bamInput = NULL, bamOutput = NULL, editable=FALSE){
        .Object <- init(.Object,"Rsortbam",editable,list(arg1=atacProc))
        # necessary parameters
        if(!is.null(atacProc)){ # atacproc from SamToBam
            .Object@paramlist[["bamInput"]] <- getParam(atacProc, "bamOutput")
            regexProcName <- sprintf("(bam|%s)", getProcName(atacProc))
        }else if(is.null(atacProc)){ # input
            .Object@paramlist[["bamInput"]] <- bamInput
            regexProcName <- "(bam)"
        }
        # unnecessary parameters
        if(is.null(bamOutput)){
            prefix <- getBasenamePrefix(.Object, .Object@paramlist[["bamInput"]], regexProcName)
            .Object@paramlist[["bamOutput_tmp"]] <- file.path(.obtainConfigure("tmpdir"),
                                                              paste0(prefix, ".", getProcName(.Object)))
            .Object@paramlist[["bamOutput"]] <- paste0(.Object@paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
        }else{
            name_split <- unlist(base::strsplit(x = bamOutput, split = ".", fixed = TRUE))
            suffix <- tail(name_split, 1)
            if(suffix == "bam"){
                .Object@paramlist[["bamOutput_tmp"]] <- paste0(name_split[-length(name_split)], collapse = ".")
                .Object@paramlist[["bamOutput"]] <- bamOutput
            }else{
                .Object@paramlist[["bamOutput_tmp"]] <- bamOutput
                .Object@paramlist[["bamOutput"]] <- paste0(.Object@paramlist[["bamOutput_tmp"]], ".bam", collapse = "")
            }
        }
        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "Rsortbam",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["bamInput"]]))
        .Object <- writeLog(.Object,sprintf("destination:%s",.Object@paramlist[["bamOutput"]]))
        Rsamtools::sortBam(file = .Object@paramlist[["bamInput"]], destination = .Object@paramlist[["bamOutput_tmp"]])
        Rsamtools::indexBam(files = .Object@paramlist[["bamOutput"]])
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "Rsortbam",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bamInput"]])){
            stop("Parameter bamInput or atacProc is required!")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "Rsortbam",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bamInput"]])
        checkPathExist(.Object,.Object@paramlist[["bamOutput"]])
    }
)


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
#' \link[Rsamtools]{sortBam}
#' \link[Rsamtools]{indexBam}
#' \code{\link{atacSam2Bam}}
#' \code{\link{atacBam2Bed}}

#' @name atacBamSort
#' @export
#' @docType methods
#' @rdname atacBamSort-methods
setGeneric("atacBamSort",function(atacProc = NULL,
                                  bamInput = NULL, bamOutput = NULL) standardGeneric("atacBamSort"))

#' @rdname atacBamSort-methods
#' @aliases atacBamSort
setMethod(
    f = "atacBamSort",
    signature = "ATACProc",
    definition = function(atacProc = NULL,
                          bamInput = NULL, bamOutput = NULL){
        atacproc <- new(
            "Rsortbam",
            atacProc = atacProc,
            bamInput = bamInput,
            bamOutput = bamOutput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)

#' @rdname atacBamSort-methods
#' @export
bamsort <- function(bamInput = NULL, bamOutput = NULL){
    atacproc <- new(
        "Rsortbam",
        atacProc = NULL,
        bamInput = bamInput,
        bamOutput = bamOutput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
