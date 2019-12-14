setClass(Class = "BedToBigWig",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "BedToBigWig",
    definition = function(.Object,prevSteps = list(), ...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        bsgenome <- allparam[["bsgenome"]]
        bwOutput <- allparam[["bwOutput"]]
        toWig <- allparam[["toWig"]]
        

        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc<-c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["bedInput"]] <- output(atacProc)[["bedOutput"]]
            }
        }
       

        if(!is.null(bedInput)){
            input(.Object)[["bedInput"]] <- bedInput;
        }
        
        if(toWig){
            sfx<-"wig"
        }else{
            sfx<-"bw"
        }

        if(is.null(bwOutput)){
            if(!is.null(input(.Object)[["bedInput"]])){
                output(.Object)[["bwOutput"]] <- getAutoPath(.Object, input(.Object)[["bedInput"]], "BED|Bed|bed", sfx)
            }
        }else{
            output(.Object)[["bwOutput"]] <- bwOutput;
        }

        param(.Object)[["bsgenome"]] <- bsgenome

        param(.Object)[["toWig"]] <- toWig

        .Object
    }
)

setMethod(
    f = "processing",
    signature = "BedToBigWig",
    definition = function(.Object,...){
        if(is.null(param(.Object)[["bsgenome"]])){
            genome <- seqinfo(getRefRc("bsgenome"))
        }else{
            genome <- seqinfo(param(.Object)[["bsgenome"]])
        }
#        bedranges <- import(input(.Object)[["bedInput"]], genome = genome)
       bedranges <- import(input(.Object)[["bedInput"]])
        cov <- coverage(bedranges)
        ans <- GRanges(cov)
        ans <- subset(ans, score > 0)
        if(param(.Object)[["toWig"]]){
            export.wig(as(ans, "UCSCData"),output(.Object)[["bwOutput"]])
        }else{
            export.bw(ans,output(.Object)[["bwOutput"]])
        }
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "BedToBigWig",
    definition = function(.Object, ...){
        .Object
    }
)


#' @name BedToBigWig
#' @title generate BigWig file from BED file
#' @description
#' This function is used to generate BigWig file
#' from BED reads file.
#' The BigWig file can be shown reads coverage on genome browser.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}},
#' \code{\link{atacBedUtils}}.
#' @param bedInput \code{Character} scalar.
#' Bed file input path.
#' @param bsgenome \code{BSGenome} object scalar.
#' BSGenome object for specific species.
#' @param bwOutput \code{Character} scalar.
#' BigWig file output path.
#' @param toWig \code{Logical} scalar.
#' @param ... Additional arguments, currently unused.
#' Save as wig file instead of binary BigWig file
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{bedToBigWig} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacBedUtils}}
#' \code{\link{bedUtils}}
#'
#' @examples
#' library(R.utils)
#' td <- tempdir()
#' setTmpDir(td)
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' \dontrun{
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' bedToBigWig(bedfile, BSgenome.Hsapiens.UCSC.hg19)
#'
#' dir(td)
#' }
#' @importFrom rtracklayer export.wig export.bw


setGeneric("atacBedToBigWig",function(atacProc, bedInput = NULL,
                                      bsgenome = NULL, bwOutput = NULL,
                                      toWig = FALSE, ...) standardGeneric("atacBedToBigWig"))


#' @rdname BedToBigWig
#' @aliases atacBedToBigWig
#' @export
setMethod(
    f = "atacBedToBigWig",
    signature = "ATACProc",
    definition = function(atacProc, bedInput = NULL,
                          bsgenome = NULL, bwOutput = NULL,
                          toWig = FALSE, ...){
        allpara <- c(list(Class = "BedToBigWig", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname BedToBigWig
#' @aliases bedToBigWig
#' @export
bedToBigWig <- function(bedInput, bsgenome = NULL, bwOutput = NULL, toWig = FALSE, ...){
    allpara <- c(list(Class = "BedToBigWig", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
