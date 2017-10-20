setClass(Class = "FRiPQC",
         contains = "ATACProc"
)


setMethod(
    f = "initialize",
    signature = "FRiPQC",
    definition = function(.Object,atacProcReads,atacProcPeak,bsgenome = NULL,
                          reportOutput=NULL,readsBedInput=NULL,
                          peakBedInput=NULL,editable=FALSE){
        .Object <- init(.Object,"FripQC",editable,list(arg1 = atacProcReads, arg2 = atacProcPeak))
        if(!is.null(atacProcReads)){
            .Object@paramlist[["readsBedInput"]] <- getParam(atacProcReads, "bedOutput");
            #.Object@paramlist[["readsCount"]] <- atacProcReads$getParam("readsCount")
            regexProcName<-sprintf("(BED|bed|Bed|%s)",getProcName(atacProcReads))
        }else{
            regexProcName<-"(BED|bed|Bed)"
        }
        if(!is.null(atacProcPeak)){
            .Object@paramlist[["peakBedInput"]] <- getParam(atacProcPeak, "bedOutput");
        }

        if(!is.null(readsBedInput)){
            .Object@paramlist[["readsBedInput"]] <- readsBedInput
        }
        if(!is.null(peakBedInput)){
            .Object@paramlist[["peakBedInput"]] <- peakBedInput
        }

        if(is.null(reportOutput)){
            if(!is.null(.Object@paramlist[["peakBedInput"]])){
                prefix<-getBasenamePrefix(.Object, .Object@paramlist[["peakBedInput"]],regexProcName)
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
    signature = "FRiPQC",
    definition = function(.Object,...){
        qcval=list();
        if(is.null(.Object@paramlist[["bsgenome"]])){
            genome <- seqinfo(.obtainConfigure("bsgenome"))
        }else{
            genome <- seqinfo(.Object@paramlist[["bsgenome"]])
        }
        message("load reads")
        gr_a <- rtracklayer::import(con = .Object@paramlist[["readsBedInput"]], genome = genome,format = "bed")
        message("load peak")
        gr_b <- rtracklayer::import(con = .Object@paramlist[["peakBedInput"]], genome = genome,format = "bed")
        message("get overlap number")
        qcval[["peakReads"]]<-length(subsetByOverlaps(gr_a, gr_b, ignore.strand = TRUE))
        message("finish overlap")
        qcval[["totalReads"]]<-length(gr_a)
        qcval[["totalPeaks"]]<-length(gr_b)
        qcval[["FRiP"]]<-qcval[["peakReads"]]/qcval[["totalReads"]]
        #unlink(paste0(.Object@paramlist[["reportPrefix"]],".tmp"))
        ####.Object@paramlist[["qcval"]]<-qcval
        write.table(as.data.frame(qcval),file = .Object@paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "FRiPQC",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["readsBedInput"]])){
            stop("readsBedInput is required.")
        }
        if(is.null(.Object@paramlist[["peakBedInput"]])){
            stop("peakBedInput is required.")
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "FRiPQC",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["readsBedInput"]]);
        checkFileExist(.Object,.Object@paramlist[["peakBedInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["reportOutput"]]);
    }
)



setMethod(
    f = "getReportValImp",
    signature = "FRiPQC",
    definition = function(.Object, item){
        qcval <- as.list(read.table(file= .Object@paramlist[["reportOutput"]],header=TRUE))
        if(item == "report"){
            cqcval<-as.character(qcval)
            cqcval[4] <- sprintf("%.2f",as.numeric(cqcval[4]))
            return(data.frame(Item=c("The number of reads in peak",
                                     "The number of total reads",
                                     "The number of total peaks",
                                     "FRiP"),Value=cqcval))
        }else{
            return(qcval[[item]])
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "FRiPQC",
    definition = function(.Object){
        return(c("report","peakReads","totalReads","totalPeaks","FRiP"))
    }
)














FRiPQC <-R6Class(
    classname = "FRiPQC",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProcReads,atacProcPeak,bsgenome = NULL,reportOutput=NULL,readsBedInput=NULL,peakBedInput=NULL,editable=FALSE){
            super$initialize("FRiPQC",editable,list(arg1=atacProcReads,arg2=atacProcPeak))
            if(!is.null(atacProcReads)){
                private$paramlist[["readsBedInput"]] <- atacProcReads$getParam("bedOutput");
                #private$paramlist[["readsCount"]] <- atacProcReads$getParam("readsCount")
                regexProcName<-sprintf("(BED|bed|Bed|%s)",atacProcReads$getProcName())
            }else{
                regexProcName<-"(BED|bed|Bed)"
            }
            if(!is.null(atacProcPeak)){
                private$paramlist[["peakBedInput"]] <- atacProcPeak$getParam("bedOutput");
            }

            if(!is.null(readsBedInput)){
                private$paramlist[["readsBedInput"]] <- readsBedInput
            }
            if(!is.null(peakBedInput)){
                private$paramlist[["peakBedInput"]] <- peakBedInput
            }

            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["peakBedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["peakBedInput"]],regexProcName)
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
            qcval=list();
            if(is.null(private$paramlist[["bsgenome"]])){
                genome <- seqinfo(.obtainConfigure("bsgenome"))
            }else{
                genome <- seqinfo(private$paramlist[["bsgenome"]])
            }
            message("load reads")
            gr_a <- rtracklayer::import(con = private$paramlist[["readsBedInput"]], genome = genome,format = "bed")
            message("load peak")
            gr_b <- rtracklayer::import(con = private$paramlist[["peakBedInput"]], genome = genome,format = "bed")
            message("get overlap number")
            qcval[["peakReads"]]<-length(subsetByOverlaps(gr_a, gr_b, ignore.strand = TRUE))
            message("finish overlap")
            qcval[["totalReads"]]<-length(gr_a)
            qcval[["totalPeaks"]]<-length(gr_b)
            qcval[["FRiP"]]<-qcval[["peakReads"]]/qcval[["totalReads"]]
            #unlink(paste0(private$paramlist[["reportPrefix"]],".tmp"))
            ####private$paramlist[["qcval"]]<-qcval
            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)

        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["readsBedInput"]])){
                stop("readsBedInput is required.")
            }
            if(is.null(private$paramlist[["peakBedInput"]])){
                stop("peakBedInput is required.")
            }

        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["readsBedInput"]]);
            private$checkFileExist(private$paramlist[["peakBedInput"]]);
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            if(item == "report"){
                cqcval<-as.character(qcval)
                cqcval[4] <- sprintf("%.2f",as.numeric(cqcval[4]))
                return(data.frame(Item=c("The number of reads in peak",
                                         "The number of total reads",
                                         "The number of total peaks",
                                         "FRiP"),Value=cqcval))
            }else{
                return(qcval[[item]])
            }
        },
        getReportItemsImp = function(){
            return(c("report","peakReads","totalReads","totalPeaks","FRiP"))
        }
    )


)

#' @name atacFripQC
#' @aliases atacFripQC
#' @aliases fripQC
#' @title Quality control for fraction of reads in peaks (FRiP)
#' @description
#' Calculate the fraction of reads falling within peak regions
#' @param atacProcReads \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacBedUtils}}
#' \code{\link{bedUtils}}
#' @param atacProcPeak \code{\link{ATACProc}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacPeakCalling}},
#' \code{\link{peakCalling}}.
#' @param bsgenome \code{BSGenome} object scalar.
#' BSGenome object for specific species.
#' @param reportOutput \code{Character} scalar.
#' The report file path
#' @param readsBedInput \code{Character} scalar.
#' Reads BED file for peak calling.
#' @param peakBedInput \code{Character} scalar.
#' Peaks BED file
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{fregLenDistr} instead.
#' @return An invisible \code{\link{fripQC}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}
#' @examples
#' library(R.utils)
#' library(magrittr)
#' td <- tempdir()
#' setConfigure("tmpdir",td)
#'
#' bedbzfile <- system.file(package="ATACpipe", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#'
#' readsProc<-bedUtils(bedInput = bedfile,maxFregLen = 100, chrFilterList = NULL)
#' peaksProc<- readsProc %>% atacPeakCalling
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' #atacFripQC(readsProc,peaksProc,bsgenome=BSgenome.Hsapiens.UCSC.hg19)
#'
#'
#' @rdname atacFripQC
#' @exportMethod atacFripQC
setGeneric("atacFripQC",function(atacProcReads,atacProcPeak,bsgenome = NULL,
                                  reportOutput=NULL,readsBedInput=NULL,
                                  peakBedInput=NULL) standardGeneric("atacFripQC"))
setMethod(
    f = "atacFripQC",
    signature = "ATACProc",
    definition = function(atacProcReads,atacProcPeak,bsgenome = NULL,
                          reportOutput=NULL,readsBedInput=NULL,
                          peakBedInput=NULL){
        atacproc <- new(
            "FRiPQC",
            atacProcReads = atacProcReads,
            atacProcPeak = atacProcPeak,
            bsgenome = bsgenome,
            reportOutput = reportOutput,
            readsBedInput = readsBedInput,
            peakBedInput = peakBedInput)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacFripQC
#' @export
fripQC<-function(readsBedInput, peakBedInput,bsgenome = NULL, reportOutput=NULL){
    atacproc <- new(
        "FRiPQC",
        atacProcReads = NULL,
        atacProcPeak = NULL,
        bsgenome = bsgenome,
        reportOutput = reportOutput,
        readsBedInput = readsBedInput,
        peakBedInput = peakBedInput)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
