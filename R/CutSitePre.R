setClass(Class = "CutSitePre",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "CutSitePre",
    definition = function(.Object, prevSteps = list(), ...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        csOutput.dir <- allparam[["csOutput.dir"]]
        prefix <- allparam[["prefix"]]
        snp.info <- allparam[["snp.info"]]

        atacProc <- NULL
        if(length(prevSteps) > 0){
            atacProc <- prevSteps[[1]]
        }
        
        if(!is.null(atacProc)){
            input(.Object)[["bedInput"]] <- getParam(atacProc, "bedOutput")
        }else{
            input(.Object)[["bedInput"]] <- bedInput
            regexProcName <- "(bed)"
        }

        if(is.null(csOutput.dir)){
 
            output(.Object)[["csOutput.dir"]] <- getStepWorkDir(.Object, "csOutput")
            
        }else{
            output(.Object)[["csOutput.dir"]] <- csOutput.dir
        }

        if(is.null(prefix)){
            param(.Object)[["csfile.dir"]] <- file.path(output(.Object)[["csOutput.dir"]], "Cutsite")
        }else{
            param(.Object)[["csfile.dir"]] <-  file.path(output(.Object)[["csOutput.dir"]], prefix)
        }


        
        
        # may use inthe future
        # .Object@paramlist[["csfile.rds"]] <- paste(.Object@paramlist[["csfile.dir"]], ".rds", sep = "")

        .Object
    }
)


setMethod(
    f = "processing",
    signature = "CutSitePre",
    definition = function(.Object,...){
       
        dir.create(output(.Object)[["csOutput.dir"]])
       # dir.create(output(.Object)[["csfile.dir"]])

        
        file_path_data <- .CutSite_call(InputFile = input(.Object)[["bedInput"]], OutputFile = param(.Object)[["csfile.dir"]])
        # may use in the future
        # saveRDS(object = file_path_data, file = .Object@paramlist[["csfile.rds"]])
        .Object
    }

)

setMethod(
  f = "genReport",
  signature = "CutSitePre",
  definition = function(.Object, ...){
    .Object
  }
)


#' @name CutSitePre
#' @title Extract ATAC-seq cutting site from bed file.
#' @description
#' Extract cutting site from ATAC-seq fangment bed file
#' (from \code{\link{atacSamToBed}}).
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
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
#' @param ... Additional arguments, currently unused.
#' @details In ATAC-seq data, every line in merged bed file
#' (from \code{\link{atacSamToBed}}, the first 3 column is chr, start, end)
#' means a DNA fragment, the cutting site is start+1 and end, this function
#' extract and sort this information for the next step
#' (\code{\link{atacCutSiteCount}}).
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream
#' analysis.
#' @author Wei Zhang
#'
#' @examples
#'
#' library(R.utils)
#' fra_path <- system.file("extdata", "chr20.50000.bed.bz2", package="esATAC")
#' frag <- as.vector(bunzip2(filename = fra_path,
#' destname = file.path(getwd(), "chr20.50000.bed"),
#' ext="bz2", FUN=bzfile, overwrite=TRUE, remove = FALSE))
#' extractcutsite(bedInput = frag, prefix = "ATAC")
#'
#' @seealso
#' \code{\link{atacCutSiteCount}}


setGeneric("atacExtractCutSite",
           function(atacProc, bedInput = NULL, csOutput.dir = NULL, prefix = NULL, ...) standardGeneric("atacExtractCutSite"))

#' @rdname CutSitePre
#' @aliases atacExtractCutSite
#' @export
setMethod(
    f = "atacExtractCutSite",
    signature = "ATACProc",
    definition = function(atacProc, bedInput = NULL, csOutput.dir = NULL, prefix = NULL, ...){
        allpara <- c(list(Class = "CutSitePre", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)


#' @rdname CutSitePre
#' @aliases extractcutsite
#' @export

extractcutsite <- function(bedInput, csOutput.dir = NULL, prefix = NULL, ...){
    allpara <- c(list(Class = "CutSitePre", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}

