setClass(Class = "LibComplexQC",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "LibComplexQC",
    definition = function(.Object,atacProc,...,reportOutput=NULL,samInput=NULL,
                          singleEnd = FALSE,subsampleSize=Inf,editable=FALSE){
        .Object <- init(.Object,"LibComplexQC",editable,list(arg1=atacProc))
        if(!is.null(atacProc)){
            .Object@paramlist[["samInput"]] <- getParam(atacProc,"samOutput");
            regexProcName<-sprintf("(SAM|Sam|sam|%s)",getProcName(atacProc))
        }else{
            regexProcName<-"(SAM|Sam|sam)"
            .Object@singleEnd<-singleEnd
        }
        if(!is.null(samInput)){
            .Object@paramlist[["samInput"]] <- samInput;
        }
        if(is.null(reportOutput)){
            if(!is.null(.Object@paramlist[["samInput"]])){
                prefix<-getBasenamePrefix(.Object,.Object@paramlist[["samInput"]],regexProcName)
                .Object@paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",getProcName(.Object),".report"))
            }
        }else{
            .Object@paramlist[["reportOutput"]] <- reportOutput;
        }

        if(is.infinite(subsampleSize)){
            .Object@paramlist[["subsample"]] <- FALSE;
            .Object@paramlist[["subsampleSize"]] <- 1e9;
        }else{
            .Object@paramlist[["subsample"]] <- TRUE;
            .Object@paramlist[["subsampleSize"]] <- subsampleSize;
        }
        paramValidation(.Object)
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "LibComplexQC",
    definition = function(.Object,...){
        if(!.Object@singleEnd){
            qcval0<-.sam2bed_merge_call(samfile = .Object@paramlist[["samInput"]], bedfile = paste0(.Object@paramlist[["reportOutput"]],".tmp"),
                                        posOffset = 0, negOffset = 0,sortBed = FALSE,
                                        uniqueBed = FALSE, filterList = NULL,minFregLen = 0,maxFregLen = 1000000,saveExtLen = FALSE ,downSample=.Object@paramlist[["subsampleSize"]])
        }else{
            qcval0<-.sam2bed_call(samfile = .Object@paramlist[["samInput"]], bedfile = paste0(.Object@paramlist[["reportOutput"]],".tmp"),
                                  posOffset = 0, negOffset = 0, sortBed = FALSE, uniqueBed = FALSE,  filterList = NULL,downSample=.Object@paramlist[["subsampleSize"]])
        }
        qcval<-.lib_complex_qc_call(bedfile=paste0(.Object@paramlist[["reportOutput"]],".tmp"), sortedBed=FALSE, max_reads=.Object@paramlist[["subsampleSize"]])
        qcval[["samTotal"]] <- qcval0[["total"]]
        qcval[["chrM"]] <- qcval0[["filted"]]
        qcval[["multimap"]] <- qcval0[["multimap"]]
        qcval[["nonMultimap"]] <-as.character(as.numeric(qcval0[["total"]])-as.numeric(qcval0[["multimap"]]))
        qcval[["NRF"]] <- as.numeric(qcval[["total"]])/
            (as.numeric(qcval[["nonMultimap"]]))

        unlink(paste0(.Object@paramlist[["reportOutput"]],".tmp"))
        print(qcval)
        print(.Object@paramlist[["reportOutput"]])
        write.table(as.data.frame(qcval),file = .Object@paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "LibComplexQC",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["samInput"]])){
            stop("samInput is required.")
        }
    }
)



setMethod(
    f = "checkAllPath",
    signature = "LibComplexQC",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["samInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["reportOutput"]]);
    }
)


setMethod(
    f = "getReportValImp",
    signature = "LibComplexQC",
    definition = function(.Object, item){
        qcval <- as.list(read.table(file= .Object@paramlist[["reportOutput"]],header=TRUE))
        if(item == "report"){
            showdf<-data.frame(
                Item = c(
                    "Total mapped reads (ratio of original reads)",
                    "Unique locations mapped uniquely by reads",
                    "Uniquely mappable reads",
                    "Non-Redundant Fraction (NRF)",
                    "Locations with only 1 reads mapping uniquely",
                    "Locations with only 2 reads mapping uniquely",
                    "PCR Bottlenecking Coefficients 1 (PBC1)",
                    "PCR Bottlenecking Coefficients 2 (PBC2)"),
                Value = c(
                    getVMShow(qcval[["samTotal"]],TRUE),
                    getVMShow(qcval[["total"]],TRUE),
                    getVMShow(qcval[["nonMultimap"]],TRUE),
                    sprintf("%.2f",qcval[["NRF"]]),
                    getVMShow(qcval[["one"]],TRUE),
                    getVMShow(qcval[["two"]],TRUE),
                    sprintf("%.2f",qcval[["PBC1"]]),
                    sprintf("%.2f",qcval[["PBC2"]])
                ),
                Reference = c("",
                              "",
                              "",
                              ">0.7",
                              "",
                              "",
                              ">0.7",
                              ">3"
                )
            )
            return(showdf)
            #return(data.frame(Item=names(qcval),Value=as.character(qcval)))
        }else{
            return(qcval[[item]])
        }
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "LibComplexQC",
    definition = function(.Object){
        return(c("report","NRF","PBC1","PBC2","one","two","total","reads","nonMultimap"))
    }
)

#' @name atacLibComplexQC
#' @title Quality control for library complexity
#' @description
#' The function calculate the nonredundant fraction of reads (NRF).
#' Its definition is number of distinct uniquely mapping reads (i.e. after removing duplicates) / Total number of reads.
#' The function also Calculate PCR Bottlenecking Coefficient 1 (PBC1) and
#' PCR Bottlenecking Coefficient 2 (PBC2).
#' PBC1=M1/M_DISTINCT and PBC2=M1/M2, where
#' M1: number of genomic locations where exactly one read maps uniquely,
#' M2: number of genomic locations where two reads map uniquely
#' M_DISTINCT: number of distinct genomic locations to which some read maps uniquely.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{bowtie2Mapping}}
#' @param reportOutput \code{Character} scalar.
#' The report file path
#' @param samInput \code{Character} scalar.
#' The SAM file input path.
#' @param singleEnd \code{Character} scalar.
#' Single end data if TRUE. Paired end data if FALSE.
#' @param subsampleSize \code{Integer} scalar.
#' Down sample reads if the number is less than total number
#' when \code{subsample} is TRUE
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{fregLenDistr} instead.
#' @return An invisible \code{\link{libComplexQC}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{bowtie2Mapping}}
#'
#' @examples
#' library(R.utils)
#' td <- tempdir()
#' options(atacConf=setConfigure("tmpdir",td))
#'
#' sambzfile <- system.file(package="esATAC", "extdata", "Example.sam.bz2")
#' samfile <- file.path(td,"Example.sam")
#' bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
#' atacproc<-libComplexQC(samInput = samfile)
#'
#' @name atacLibComplexQC
#' @export
#' @docType methods
#' @rdname atacLibComplexQC-methods
setGeneric("atacLibComplexQC",function(atacProc,reportOutput=NULL,samInput=NULL,
                                  singleEnd = FALSE,subsampleSize=Inf, ...) standardGeneric("atacLibComplexQC"))

#' @rdname atacLibComplexQC-methods
#' @aliases atacLibComplexQC
setMethod(
    f = "atacLibComplexQC",
    signature = "ATACProc",
    definition = function(atacProc,reportOutput=NULL,samInput=NULL,
                          singleEnd = FALSE,subsampleSize=Inf, ...){
        atacproc <- new(
            "LibComplexQC",
            atacProc = atacProc,
            reportOutput = reportOutput,
            samInput = samInput,
            singleEnd = singleEnd,
            subsampleSize = subsampleSize)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacLibComplexQC-methods
#' @export
libComplexQC<-function(samInput, reportOutput=NULL,singleEnd = FALSE,subsampleSize=Inf, ...){
    atacproc <- new(
        "LibComplexQC",
        atacProc = NULL,
        reportOutput = reportOutput,
        samInput = samInput,
        singleEnd = singleEnd,
        subsampleSize = subsampleSize)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
