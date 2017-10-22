setClass(Class = "SamToBam",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "SamToBam",
    definition = function(.Object,atacProc, ..., samInput = NULL, bamOutput = NULL, editable = FALSE){
        .Object <- init(.Object,"SamToBam",editable,list(arg1=atacProc))
        # necessary parameters
        if((!is.null(atacProc)) ){
            .Object@paramlist[["samInput"]] <- getParam(atacProc, "samOutput")
            regexProcName <- sprintf("(sam|%s)", getProcName(atacProc))
        }else if(is.null(atacProc)){ # input
            .Object@paramlist[["samInput"]] <- samInput
            regexProcName <- "(sam)"
        }
        # unnecessary parameters
        if(is.null(bamOutput)){
            prefix <- getBasenamePrefix(.Object, .Object@paramlist[["samInput"]], regexProcName)
            .Object@paramlist[["name_tmp"]] <- file.path(.obtainConfigure("tmpdir"),
                                                         paste0(prefix, ".", getProcName(.Object)))
            .Object@paramlist[["bamOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                          paste0(prefix, ".", getProcName(.Object), ".bam"))
            .Object@paramlist[["baiOutput"]] <- file.path(.obtainConfigure("tmpdir"),
                                                          paste0(prefix, ".", getProcName(.Object), ".bam.bai"))
        }else{
            name_split <- unlist(base::strsplit(x = bamOutput, split = ".", fixed = TRUE))
            suffix <- tail(name_split, 1)
            if(suffix == "bam"){
                .Object@paramlist[["name_tmp"]] <- paste0(name_split[-length(name_split)], collapse = ".")
                .Object@paramlist[["bamOutput"]] <- bamInput
                .Object@paramlist[["baiOutput"]] <- paste0(.Object@paramlist[["bamOutput"]], ".bai", collapse = "")
            }else{
                .Object@paramlist[["name_tmp"]] <- bamInput
                .Object@paramlist[["bamOutput"]] <- paste0(.Object@paramlist[["name_tmp"]], ".bam", collapse = "")
                .Object@paramlist[["baiOutput"]] <- paste0(.Object@paramlist[["bamOutput"]], ".bai", collapse = "")
            }

        }
        paramValidation(.Object)
        .Object

    }
)




setMethod(
    f = "processing",
    signature = "SamToBam",
    definition = function(.Object,...){
        .Object <- writeLog(.Object,paste0("processing file:"))
        .Object <- writeLog(.Object,sprintf("source:%s",.Object@paramlist[["samInput"]]))
        .Object <- writeLog(.Object,sprintf("Bam destination:%s",.Object@paramlist[["bamOutput"]]))
        .Object <- writeLog(.Object,sprintf("Bai destination:%s",.Object@paramlist[["baiOutput"]]))
        Rsamtools::asBam(file = .Object@paramlist[["samInput"]], destination = .Object@paramlist[["name_tmp"]])
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "SamToBam",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["samInput"]])){
            stop("Parameter samInput is required!")
        }
    }
)

setMethod(
    f = "checkAllPath",
    signature = "SamToBam",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["samInput"]])
        checkPathExist(.Object,.Object@paramlist[["bamOutput"]])
    }
)



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
#' @details The sam file wiil be automatically obtained from
#' object(\code{atacProc}) or input by hand. bamOutput can be ignored.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for
#' downstream analysis.
#' @author Wei Zhang
#' @examples
#'
#' library(R.utils)
#' sam_bz <- system.file("extdata", "Example.sam.bz2", package="ATACpipe")
#' sam_path <- as.vector(bunzip2(filename = sam_bz,
#' destname = file.path(getwd(), "Example.sam"),
#' ext="bz2", FUN=bzfile, remove = FALSE))
#' sam2bam(samInput = sam_path)
#'
#' @seealso
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{atacBam2Bed}}
#' \code{\link{atacBamSort}}

#' @name atacSam2Bam
#' @export
#' @docType methods
#' @rdname atacSam2Bam-methods
setGeneric("atacSam2Bam",function(atacProc = NULL,
                                  samInput = NULL, bamOutput = NULL) standardGeneric("atacSam2Bam"))
#' @rdname atacSam2Bam-methods
#' @aliases atacSam2Bam
setMethod(
    f = "atacSam2Bam",
    signature = "ATACProc",
    definition = function(atacProc = NULL,
                          samInput = NULL, bamOutput = NULL){
        atacproc <- new(
            "SamToBam",
            atacProc = atacProc,
            samInput = samInput,
            bamOutput = bamOutput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacSam2Bam-methods
#' @export
sam2bam <- function(samInput, bamOutput = NULL){
    atacproc <- new(
        "SamToBam",
        atacProc = NULL,
        samInput = samInput,
        bamOutput = bamOutput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
