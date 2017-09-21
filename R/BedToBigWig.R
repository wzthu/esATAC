BedToBigWig <- R6Class(
    classname = "BedToBigWig",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, bedInput = NULL, bwOutput = NULL, toWig = FALSE, editable=FALSE){
            super$initialize("BedToBigWig",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            
            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }
            if(toWig){
                sfx<-".wig"
            }else{
                sfx<-".bw"
            }
            
            if(is.null(bwOutput)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["bwOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),sfx))
                }
            }else{
                private$paramlist[["bwOutput"]] <- bwOutput;
            }
            
            private$paramlist[["toWig"]] <- toWig
            private$paramValidation()
        }
    ),
    
    private = list(
        processing = function(){
            genome <- seqinfo(.obtainConfigure("bsgenome"))
            bedranges <- import(private$paramlist[["bedInput"]], genome = genome)
            cov <- coverage(bedranges)
            ans <- GRanges(cov)
            ans <- subset(ans, score > 0)
            if(private$paramlist[["toWig"]]){
                export.wig(ans,private$paramlist[["bwOutput"]])
            }else{
                export.bw(ans,private$paramlist[["bwOutput"]])
            }
            
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop(paste("bedInput is requied"));
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileCreatable(private$paramlist[["bwOutput"]]);
        }
    )
    
    
)


#' @name atacBedToBigWig
#' @aliases atacBedToBigWig
#' @aliases bedToBigWig
#' @title generate BigWig file from BED file
#' @description 
#' This function is used to generate BigWig file 
#' from BED reads file.
#' The BigWig file can be shown reads coverage on genome browser.
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}, 
#' \code{\link{atacBedUtils}}.
#' @param bedInput \code{Character} scalar. 
#' Bed file input path. 
#' @param bwOutput \code{Character} scalar. 
#' BigWig file output path.
#' @param toWig \code{Logical} scalar.
#' Save as wig file instead of binary BigWig file
#' @details The parameter related to input and output file path
#' will be automatically 
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{bedToBigWig} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacSamToBed}} 
#' \code{\link{atacBedUtils}}

#' @rdname atacBedToBigWig
#' @export 
atacBedToBigWig <- function(atacProc, bedInput = NULL, bwOutput = NULL, toWig = FALSE){
    atacproc <- BedToBigWig$new(atacProc = atacProc, bedInput = bedInput, bwOutput = bwOutput, toWig = toWig)
    atacproc$process()
    invisible(atacproc)
}
#' @rdname atacBedToBigWig
#' @export
bedToBigWig <- function(bedInput, bwOutput = NULL, toWig = FALSE){
    atacproc <- BedToBigWig$new(atacProc = NULL, bedInput = bedInput, bwOutput = bwOutput, toWig = toWig)
    atacproc$process()
    invisible(atacproc)
}
