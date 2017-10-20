setClass(Class = "CutSitePre",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "CutSitePre",
    definition = function(.Object, atacProc, ..., bedInput = NULL, csOutput.dir = NULL,
                          prefix = NULL, editable = FALSE){
        .Object <- init(.Object, "CutSitePre", editable, list(arg1 = atacProc))

        if(!is.null(atacProc)){
            .Object@paramlist[["bedInput"]] <- getParam(atacProc, "bedOutput")
            regexProcName <- sprintf("(bed|%s)", getProcName(atacProc))
        }else{
            .Object@paramlist[["bedInput"]] <- bedInput
            regexProcName <- "(bed)"
        }

        if(is.null(csOutput.dir)){
            dir_name <- getBasenamePrefix(.Object, .Object@paramlist[["bedInput"]], regexProcName)
            csOutput_dir <- file.path(.obtainConfigure("tmpdir"),
                                      dir_name)
            dir.create(csOutput_dir)
        }else{
            csOutput_dir <- csOutput.dir
        }

        if(is.null(prefix)){
            .Object@paramlist[["csfile.dir"]] <- paste(csOutput_dir, "/", "Cutsite", sep = "")
        }else{
            .Object@paramlist[["csfile.dir"]] <- paste(csOutput_dir, "/", prefix, sep = "")
        }

        # may use inthe future
        # .Object@paramlist[["csfile.rds"]] <- paste(.Object@paramlist[["csfile.dir"]], ".rds", sep = "")

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "CutSitePre",
    definition = function(.Object,...){
        .Object <- writeLog(.Object, paste0("processing file:"))
        .Object <- writeLog(.Object, sprintf("source:%s", .Object@paramlist[["bedInput"]]))
        .Object <- writeLog(.Object, sprintf("destination:%s", .Object@paramlist[["csfile.dir"]]))
        file_path_data <- .CutSite_call(InputFile = .Object@paramlist[["bedInput"]], OutputFile = .Object@paramlist[["csfile.dir"]])
        # may use in the future
        # saveRDS(object = file_path_data, file = .Object@paramlist[["csfile.rds"]])
        .Object
    }

)


setMethod(
    f = "checkRequireParam",
    signature = "CutSitePre",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bedInput"]])){
            stop("Parameter bedInput is required!")
        }
    }
)


setMethod(
    f = "checkAllPath",
    signature = "CutSitePre",
    definition = function(.Object,...){
        checkFileExist(.Object, .Object@paramlist[["bedInput"]])
        checkPathExist(.Object, .Object@paramlist[["csfile.dir"]])
    }
)


CutSitePre <- R6::R6Class(
    classname = "CutSitePre",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, bedInput = NULL, csOutput.dir = NULL,
                              prefix = NULL, editable = FALSE){
            super$initialize("CutSitePre", editable, list(arg1 = atacProc))

            # necessary parameters
            if(!is.null(atacProc)){ # get parameter from class BamToBed or SamToBed
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput")
                regexProcName <- sprintf("(bed|%s)", atacProc$getProcName())
            }else{
                private$paramlist[["bedInput"]] <- bedInput
                regexProcName <- "(bed)"
            }
            # unnecessary parameters
            if(is.null(csOutput.dir)){
                dir_name <- private$getBasenamePrefix(private$paramlist[["bedInput"]], regexProcName)
                csOutput_dir <- file.path(.obtainConfigure("tmpdir"),
                                          dir_name)
                dir.create(csOutput_dir)
            }else{
                csOutput_dir <- csOutput.dir
            }
            if(is.null(prefix)){
                private$paramlist[["csfile.dir"]] <- paste0(csOutput_dir, "/", "Cutsite", collapse = "")
            }else{
                private$paramlist[["csfile.dir"]] <- paste0(csOutput_dir, "/", prefix, collapse = "")
            }
            # may use inthe future
            # private$paramlist[["csfile.rds"]] <- paste0(private$paramlist[["csfile.dir"]], ".rds", collapse = "")
            # parameter check
            private$paramValidation()
        } # initialization end

    ), # public end


    private = list(
        processing = function(){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s", private$paramlist[["bedInput"]]))
            private$writeLog(sprintf("destination:%s", private$paramlist[["csfile.dir"]]))
            file_path_data <- .CutSite_call(InputFile = private$paramlist[["bedInput"]], OutputFile = private$paramlist[["csfile.dir"]])
            # may use in the future
            # saveRDS(object = file_path_data, file = private$paramlist[["csfile.rds"]])
        }, # processing end

        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop("Parameter bedInput is required!")
            }
        }, # checkRequireParam end

        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]])
            private$checkPathExist(private$paramlist[["csfile.dir"]])
        } # checkAllPath end

    ) # private end

) # class end


#' @name atacExtractCutSite
#' @aliases atacExtractCutSite
#' @aliases extractcutsite
#' @title Extract ATAC-seq cutting site from bed file.
#' @description
#' Extract cutting site from ATAC-seq fangment bed file
#' (from \code{\link{atacSamToBed}}).
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}.
#' @param bedInput \code{Character} scalar.
#' Input bed file path, must be merged bed file(a line is a fragment). The
#' input file should be UCSC bed format(0-based).
#' @param csOutput.dir \code{Character} scalar.
#' The output path, an empty folder would be great. Default: a folder in the
#' same path as input bed file.
#' @param prefix \code{Character} scalar.
#' Output file name prefix, e.g. prefix_chr*.bed, default "Cutsite".
#' @details In ATAC-seq data, every line in merged bed file
#' (from \code{\link{atacSamToBed}}, the first 3 column is chr, start, end)
#' means a DNA fragment, the cutting site is start+1 and end, this function
#' extract and sort this information for the next step
#' (\code{\link{atacCutSiteCount}}).
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream
#' analysis.
#' @author Wei Zhang
#'
#' @examples
#'
#' library(R.utils)
#' fra_path <- system.file("extdata", "chr20.50000.bed.bz2", package="ATACpipe")
#' frag <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20.50000.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' extractcutsite(bedInput = frag, prefix = "ATAC")
#'
#' @seealso
#' \code{\link{atacCutSiteCount}}


#' @rdname atacExtractCutSite
#' @exportMethod atacExtractCutSite
setGeneric("atacExtractCutSite",
           function(atacProc, bedInput = NULL, csOutput.dir = NULL, prefix = NULL) standardGeneric("atacExtractCutSite"))

setMethod(
    f = "atacExtractCutSite",
    signature = "ATACProc",
    definition = function(atacProc, bedInput = NULL, csOutput.dir = NULL, prefix = NULL){
        atacproc <- new("CutSitePre", atacProc = atacProc, bedInput = bedInput, csOutput.dir = csOutput.dir, prefix = prefix)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)


#' @rdname atacExtractCutSite
#' @export
extractcutsite <- function(bedInput, csOutput.dir = NULL, prefix = NULL){
    atacproc <- new("CutSitePre", atacProc = NULL, bedInput = bedInput, csOutput.dir = csOutput.dir, prefix = prefix)
    atacproc <- process(atacproc)
    invisible(atacproc)
}

