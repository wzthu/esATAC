setClass(Class = "PeakQC",
         contains = "ATACProc",
         slots = list(fixtag = "character"),
         prototype = list(fixtag = "User file")
)


setMethod(
    f = "initialize",
    signature = "PeakQC",
    definition = function(.Object,atacProc, reportOutput=NULL,
                          bsgenome = NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),
                          bedInput = NULL,editable=FALSE){
        .Object <- init(.Object,"PeakQC",editable,list(arg1=atacProc))
        if(!is.null(atacProc)){
            .Object@paramlist[["bedInput"]] <- getParam(atacProc, "bedOutput");
            regexProcName<-sprintf("(BED|bed|Bed|%s)",getProcName(atacProc))
        }else{
            regexProcName<-"(BED|bed|Bed)"
        }
        qcbedInput <- qcbedInput[1]
        if(qcbedInput == "DHS"){
            .Object@paramlist[["qcbedInput"]]<-.obtainConfigure("DHS");
            .Object@fixtag = "DHS"
        }else if(qcbedInput == "blacklist"){
            .Object@paramlist[["qcbedInput"]]<-.obtainConfigure("blacklist");
            .Object@fixtag = "blacklist"
        }else{
            .Object@paramlist[["qcbedInput"]]<-qcbedInput;
        }

        if(!is.null(bedInput)){
            .Object@paramlist[["bedInput"]] <- bedInput;
        }

        if(is.null(reportOutput)){
            if(!is.null(.Object@paramlist[["bedInput"]])){
                prefix <- getBasenamePrefix(.Object, .Object@paramlist[["bedInput"]], regexProcName)
                .Object@paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".report.txt"))
            }
        }else{
            .Object@paramlist[["reportOutput"]] <- reportOutput;
        }

        .Object@paramlist[["bsgenome"]] <- bsgenome

        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "PeakQC",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bsgenome"]])){
            genome <- seqinfo(.obtainConfigure("bsgenome"))
        }else{
            genome <- seqinfo(.Object@paramlist[["bsgenome"]])
        }

        inputbed <- import(con = .Object@paramlist[["bedInput"]], genome = genome,format = "bed")


        qcbedInput<-import(con = .Object@paramlist[["qcbedInput"]], genome = genome,format = "bed")



        qcval=list();

        qcval[["totalInput"]]<-length(inputbed)
        qcval[["qcbedInput"]]<-length(subsetByOverlaps(inputbed, qcbedInput,ignore.strand = TRUE))
        qcval[["qcbedRate"]]<-qcval[["qcbedInput"]]/qcval[["totalInput"]]

        write.table(as.data.frame(qcval),file = .Object@paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)

        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "PeakQC",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bedInput"]])){
            stop("bedInput is required.")
        }
        if(is.null(.Object@paramlist[["qcbedInput"]])){
            stop("qcbedInput is required.")
        }
    }
)

setMethod(
    f = "checkAllPath",
    signature = "PeakQC",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bedInput"]]);
        checkFileExist(.Object,.Object@paramlist[["qcbedInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["reportOutput"]]);
    }
)

setMethod(
    f = "getReportValImp",
    signature = "PeakQC",
    definition = function(.Object, item){
        qcval <- as.list(read.table(file= .Object@paramlist[["reportOutput"]],header=TRUE))
        if(item == "report"){
            cqcval<-as.character(qcval)
            cqcval[3]<-sprintf("%.2f",as.numeric(cqcval[[3]]))
            if(.Object@fixtag=="DHS"){
                return(data.frame(Item=c("Total peaks","DHS regions", "Ratio"),
                                  Value=cqcval))
            }else if(.Object@fixtag=="blacklist"){
                return(data.frame(Item=c("Total peaks","Blacklist regions", "Ratio"),Value=cqcval))
            }else{
                return(data.frame(Item=names(qcval),Value=as.character(qcval)))
            }

        }else{
            return(qcval[[item]])
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "PeakQC",
    definition = function(.Object){
        return(c("report","totalInput","qcbedInput","qcbedRate"))
    }
)



#' @name atacPeakQC
#' @title Quality control for peak overlap
#' @description
#' These functions are used to generate fregment distribution plot.
#' The fourier transform of fregment distribution will be calculated.
#' Strength distribution around period at 10.4bp and 180bp
#' will be shown in another two plots.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
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
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{peakQC} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}
#'
#' @examples
#' library(R.utils)
#' library(magrittr)
#' td <- tempdir()
#' options(atacConf=setConfigure("tmpdir",td))
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#' blacklistfile <- system.file(package="esATAC", "extdata", "hg19.blacklist.bed")
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' bedUtils(bedInput = bedfile,maxFregLen = 100, chrFilterList = NULL) %>%
#' atacPeakCalling %>% atacPeakQC(qcbedInput = blacklistfile, bsgenome = BSgenome.Hsapiens.UCSC.hg19)
#' dir(td)
#' @name  atacPeakQC
#' @export
#' @docType methods
#' @rdname atacPeakQC-methods
setGeneric("atacPeakQC",function(atacProc, bsgenome = NULL,
                                 reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),
                                 bedInput = NULL, ...) standardGeneric("atacPeakQC"))
#' @rdname atacPeakQC-methods
#' @aliases atacPeakQC
setMethod(
    f = "atacPeakQC",
    signature = "ATACProc",
    definition = function(atacProc, bsgenome = NULL,
                          reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"),
                          bedInput = NULL, ...){
        atacproc <- new(
            "PeakQC",
            atacProc = atacProc,
            bsgenome = bsgenome,
            reportOutput = reportOutput,
            qcbedInput = qcbedInput,
            bedInput = bedInput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacPeakQC-methods
#' @export
peakQC<-function(bedInput, bsgenome = NULL, reportOutput=NULL,qcbedInput = c("DHS","blacklist","path/to/bed"), ...){
    atacproc <- new(
        "PeakQC",
        atacProc = NULL,
        bsgenome = bsgenome,
        reportOutput = reportOutput,
        qcbedInput = qcbedInput,
        bedInput = bedInput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
