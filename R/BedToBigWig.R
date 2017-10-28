setClass(Class = "BedToBigWig",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "BedToBigWig",
    definition = function(.Object,atacProc, ..., bedInput = NULL, bsgenome = NULL, bwOutput = NULL, toWig = FALSE, editable=FALSE){
        .Object <- init(.Object,"BedToBigWig",editable,list(arg1=atacProc))

        if(!is.null(atacProc)){
            .Object@paramlist[["bedInput"]] <- getParam(atacProc, "bedOutput");
            regexProcName<-sprintf("(BED|bed|Bed|%s)",getProcName(atacProc))
        }else{
            regexProcName<-"(BED|bed|Bed)"
        }

        if(!is.null(bedInput)){
            .Object@paramlist[["bedInput"]] <- bedInput;
        }
        if(toWig){
            sfx<-".wig"
        }else{
            sfx<-".bw"
        }

        if(is.null(bwOutput)){
            if(!is.null(.Object@paramlist[["bedInput"]])){
                prefix <- getBasenamePrefix(.Object, .Object@paramlist[["bedInput"]],regexProcName)
                .Object@paramlist[["bwOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),sfx))
            }
        }else{
            .Object@paramlist[["bwOutput"]] <- bwOutput;
        }

        .Object@paramlist[["bsgenome"]] <- bsgenome

        .Object@paramlist[["toWig"]] <- toWig

        paramValidation(.Object)
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "BedToBigWig",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bsgenome"]])){
            genome <- seqinfo(.obtainConfigure("bsgenome"))
        }else{
            genome <- seqinfo(.Object@paramlist[["bsgenome"]])
        }
        bedranges <- import(.Object@paramlist[["bedInput"]], genome = genome)
        cov <- coverage(bedranges)
        ans <- GRanges(cov)
        ans <- subset(ans, score > 0)
        if(.Object@paramlist[["toWig"]]){
            export.wig(ans,.Object@paramlist[["bwOutput"]])
        }else{
            export.bw(ans,.Object@paramlist[["bwOutput"]])
        }
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "BedToBigWig",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bedInput"]])){
            stop(paste("bedInput is requied"));
        }
    }
)

setMethod(
    f = "checkAllPath",
    signature = "BedToBigWig",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bedInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["bwOutput"]]);
    }
)


#' @name atacBedToBigWig
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
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{bedToBigWig} instead.
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
#' options(atacConf=setConfigure("tmpdir",td))
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#'
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' bedToBigWig(bedfile, BSgenome.Hsapiens.UCSC.hg19)
#'
#' dir(td)


#' @name atacBedToBigWig
#' @export
#' @docType methods
#' @rdname atacBedToBigWig-methods
setGeneric("atacBedToBigWig",function(atacProc, bedInput = NULL,
                                      bsgenome = NULL, bwOutput = NULL,
                                      toWig = FALSE, ...) standardGeneric("atacBedToBigWig"))


#' @rdname atacBedToBigWig-methods
#' @aliases atacBedToBigWig
setMethod(
    f = "atacBedToBigWig",
    signature = "ATACProc",
    definition = function(atacProc, bedInput = NULL,
                          bsgenome = NULL, bwOutput = NULL,
                          toWig = FALSE, ...){

        atacproc <- new(
            "BedToBigWig",
            atacProc = atacProc,
            bedInput = bedInput,
            bsgenome = bsgenome,
            bwOutput = bwOutput,
            toWig = toWig)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacBedToBigWig-methods
#' @export
bedToBigWig <- function(bedInput, bsgenome = NULL, bwOutput = NULL, toWig = FALSE, ...){
    atacproc <- new(
        "BedToBigWig",
        atacProc = NULL,
        bedInput = bedInput,
        bsgenome = bsgenome,
        bwOutput = bwOutput,
        toWig = toWig)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
