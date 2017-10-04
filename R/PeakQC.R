PeakQC <-R6Class(
    classname = "PeakQC",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, reportOutput=NULL,bsgenome = NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),bedInput = NULL,editable=FALSE){
            super$initialize("PeakQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            qcbedInput <- qcbedInput[1]
            if(qcbedInput == "DHS"){
                private$paramlist[["qcbedInput"]]<-.obtainConfigure("DHS");
                private$fixtag = "DHS"
            }else if(qcbedInput == "blacklist"){
                private$paramlist[["qcbedInput"]]<-.obtainConfigure("blacklist");
                private$fixtag = "blacklist"
            }else{
                private$paramlist[["qcbedInput"]]<-qcbedInput;
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }

            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report.txt"))
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }

            private$paramlist[["bsgenome"]] <- bsgenome
            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            if(is.null(private$paramlist[["bsgenome"]])){
                genome <- seqinfo(.obtainConfigure("bsgenome"))
            }else{
                genome <- seqinfo(private$paramlist[["bsgenome"]])
            }

            inputbed <- import(con = private$paramlist[["bedInput"]], genome = genome,format = "bed")


            qcbedInput<-import(con = private$paramlist[["qcbedInput"]], genome = genome,format = "bed")



            qcval=list();

            qcval[["totalInput"]]<-length(inputbed)
            qcval[["qcbedInput"]]<-length(subsetByOverlaps(inputbed, qcbedInput,ignore.strand = TRUE))
            qcval[["qcbedRate"]]<-qcval[["qcbedInput"]]/qcval[["totalInput"]]

            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop("bedInput is required.")
            }
            if(is.null(private$paramlist[["qcbedInput"]])){
                stop("qcbedInput is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileExist(private$paramlist[["qcbedInput"]]);
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            if(item == "report"){
                cqcval<-as.character(qcval)
                cqcval[3]<-sprintf("%.2f",as.numeric(cqcval[[3]]))
                if(private$fixtag=="DHS"){
                    return(data.frame(Item=c("Total peaks","DHS regions", "Ratio"),
                                      Value=cqcval))
                }else if(private$fixtag=="blacklist"){
                    return(data.frame(Item=c("Total peaks","Blacklist regions", "Ratio"),Value=cqcval))
                }else{
                    return(data.frame(Item=names(qcval),Value=as.character(qcval)))
                }

            }else{
                return(qcval[[item]])
            }
        },
        getReportItemsImp = function(){
            return(c("report","totalInput","qcbedInput","qcbedRate"))
        },
        fixtag="User file"
    )


)

#' @name atacPeakQC
#' @aliases atacPeakQC
#' @aliases peakQC
#' @title Quality control for peak overlap
#' @description
#' These functions are used to generate fregment distribution plot.
#' The fourier transform of fregment distribution will be calculated.
#' Strength distribution around period at 10.4bp and 180bp
#' will be shown in another two plots.
#' @param atacProc \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}},
#' \code{\link{atacBedUtils}}.
#' @param reportOutput \code{Character} scalar.
#' The report file path.
#' @param bsgenome \code{BSGenome} object scalar.
#' BSGenome object for specific species.
#' @param qcbedInput \code{Character} scalar.
#' It can be "DHS","blacklist" or
#' Other quality control BED file input path.
#' @param bedInput \code{Character} scalar.
#' BED file input path for quality control.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{peakQC} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}
#'
#' @examples
#' library(R.utils)
#' library(magrittr)
#' td <- tempdir()
#' setConfigure("tmpdir",td)
#'
#' bedbzfile <- system.file(package="ATACFlow", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' blacklistfile <- system.file(package="ATACFlow", "extdata", "hg19.blacklist.bed")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' bedUtils(bedInput = bedfile,maxFregLen = 100, chrFilterList = NULL) %>%
#' atacPeakCalling %>% atacPeakQC(qcbedInput = blacklistfile, bsgenome = BSgenome.Hsapiens.UCSC.hg19)
#' dir(td)
#' @rdname atacPeakQC
#' @export
atacPeakQC<-function(atacProc, bsgenome = NULL, reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"), bedInput = NULL){
    atacproc<-PeakQC$new(atacProc, bsgenome = bsgenome, reportOutput=reportOutput,qcbedInput = qcbedInput,bedInput = bedInput,editable=FALSE)
    atacproc$process()
    invisible(atacproc)
}
#' @rdname atacPeakQC
#' @export
peakQC<-function(bedInput, bsgenome = NULL, reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed")){
    atacproc<-PeakQC$new(atacProc=NULL, bsgenome = bsgenome, reportOutput=reportOutput,qcbedInput = qcbedInput,bedInput = bedInput,editable=FALSE)
    atacproc$process()
    invisible(atacproc)
}
