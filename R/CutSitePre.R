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
            # parameter check
            private$paramValidation()
        } # initialization end

    ), # public end


    private = list(
        processing = function(){
            private$writeLog(paste0("processing file:"))
            private$writeLog(sprintf("source:%s", private$paramlist[["bedInput"]]))
            private$writeLog(sprintf("destination:%s", private$paramlist[["csfile.dir"]]))
            Cut_Site_Number <- .CutSite_call(InputFile = private$paramlist[["bedInput"]], OutputFile = private$paramlist[["csfile.dir"]])
            print(sprintf("Total Cut Site Number is: %d", Cut_Site_Number))
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
#' (from \code{\link{atacSam2Bed}}).
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSam2Bed}}.
#' @param bedInput \code{Character} scalar.
#' Input bed file path, must be merged bed file(a line is a fragment). The
#' input file should be UCSC bed format(0-based).
#' @param csOutput.dir \code{Character} scalar.
#' The output path, an empty folder would be great. Default: a folder in the
#' same path as input bed file.
#' @param prefix \code{Character} scalar.
#' Output file name prefix, e.g. prefix_chr*.bed, default "Cutsite".
#' @details In ATAC-seq data, every line in merged bed file
#' (from \code{\link{atacSam2Bed}}, the first 3 column is chr, start, end)
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
#' fra_path <- system.file("extdata", "chr20.50000.bed.bz2", package="ATACFlow")
#' frag <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20.50000.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' extractcutsite(bedInput = frag, prefix = "ATAC")
#'
#' @seealso
#' \code{\link{atacCutSiteCount}}

#' @rdname atacExtractCutSite
#' @export
atacExtractCutSite <- function(atacProc, bedInput = NULL, csOutput.dir = NULL, prefix = NULL){
    tmp <- CutSitePre$new(atacProc, bedInput, csOutput.dir, prefix)
    tmp$process()
    invisible(tmp)
}

#' @rdname atacExtractCutSite
#' @export
extractcutsite <- function(bedInput, csOutput.dir = NULL, prefix = NULL){
    tmp <- CutSitePre$new(atacProc = NULL, bedInput, csOutput.dir, prefix)
    tmp$process()
    invisible(tmp)
}

