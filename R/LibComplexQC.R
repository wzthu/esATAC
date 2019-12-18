setClass(Class = "LibComplexQC",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "LibComplexQC",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        samInput <- allparam[["samInput"]]
        reportOutput <- allparam[["reportOutput"]]
        singleEnd <- allparam[["singleEnd"]]
        subsampleSize <- allparam[["subsampleSize"]]
        
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc <- c(unlist(atacProc),list())
                atacProc <- atacProc[[length(atacProc)]]
                input(.Object)[["samInput"]] <- output(atacProc)[["samOutput"]]
                param(.Object)[["singleEnd"]] <- property(atacProc)[["singleEnd"]]
                singleEnd <- property(atacProc)[["singleEnd"]]
            }
        }else{
            param(.Object)[["singleEnd"]] <- singleEnd
            property(.Object)[["singleEnd"]] <- singleEnd
        }
        if(!is.null(samInput)){
            input(.Object)[["samInput"]] <- samInput;
        }
        if(is.null(reportOutput)){
            if(!is.null(input(.Object)[["samInput"]])){
                output(.Object)[["reportOutput"]] <- getAutoPath(.Object,input(.Object)[["samInput"]],"Sam|SAM|sam","report")
            }
        }else{
            output(.Object)[["reportOutput"]] <- reportOutput;
        }

        if(is.infinite(subsampleSize)){
            param(.Object)[["subsample"]] <- FALSE;
            param(.Object)[["subsampleSize"]] <- 1e9;
        }else{
            param(.Object)[["subsample"]] <- TRUE;
            param(.Object)[["subsampleSize"]] <- subsampleSize;
        }
        .Object
    }
)


setMethod(
    f = "processing",
    signature = "LibComplexQC",
    definition = function(.Object,...){
        if(!param(.Object)$singleEnd){
            qcval0<-.sam2bed_merge_call(samfile = input(.Object)[["samInput"]], bedfile = paste0(output(.Object)[["reportOutput"]],".tmp"),
                                        posOffset = 0, negOffset = 0,sortBed = FALSE,
                                        uniqueBed = FALSE, filterList = NULL,minFragLen = 0,maxFragLen = 1000000,saveExtLen = FALSE ,downSample=param(.Object)[["subsampleSize"]])
        }else{
            qcval0<-.sam2bed_call(samfile = input(.Object)[["samInput"]], bedfile = paste0(output(.Object)[["reportOutput"]],".tmp"),
                                  posOffset = 0, negOffset = 0, sortBed = FALSE, uniqueBed = FALSE,  filterList = NULL,downSample=param(.Object)[["subsampleSize"]])
        }
        qcval<-.lib_complex_qc_call(bedfile=paste0(output(.Object)[["reportOutput"]],".tmp"), sortedBed=FALSE, max_reads=param(.Object)[["subsampleSize"]])
        qcval[["samTotal"]] <- qcval0[["total"]]
        qcval[["chrM"]] <- qcval0[["filted"]]
        qcval[["multimap"]] <- qcval0[["multimap"]]
        qcval[["nonMultimap"]] <-as.character(as.numeric(qcval0[["total"]])-as.numeric(qcval0[["multimap"]]))
        qcval[["NRF"]] <- as.numeric(qcval[["total"]])/
            (as.numeric(qcval[["nonMultimap"]]))

        unlink(paste0(output(.Object)[["reportOutput"]],".tmp"))
        print(unlist(qcval))
        print(output(.Object)[["reportOutput"]])
        write.table(as.data.frame(qcval),file = output(.Object)[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)

        
#       
        .Object
    }
)

setMethod(
    f = "genReport",
    signature = "LibComplexQC",
    definition = function(.Object, ...){
        qcval <- as.list(read.table(file= output(.Object)[["reportOutput"]],header=TRUE))
        
        report(.Object)$table<-data.frame(
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
            ))
        for(n in names(qcval)){
            report(.Object)[[n]] <- qcval[[n]]
        }     
        
        .Object
    }
)

#' @name LibComplexQC
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
#' you can use \code{libComplexQC} instead.
#' @return An invisible \code{\link{libComplexQC}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacBowtie2Mapping}}
#' \code{\link{bowtie2Mapping}}
#'
#' @examples
#' library(R.utils)
#' td <- tempdir()
#' setTmpDir(td)
#'
#' sambzfile <- system.file(package="esATAC", "extdata", "Example.sam.bz2")
#' samfile <- file.path(td,"Example.sam")
#' bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
#' atacproc<-libComplexQC(samInput = samfile)
#'

setGeneric("atacLibComplexQC",function(atacProc,reportOutput=NULL,samInput=NULL,
                                  singleEnd = FALSE,subsampleSize=Inf, ...) standardGeneric("atacLibComplexQC"))

#' @rdname LibComplexQC
#' @aliases atacLibComplexQC
#' @export
setMethod(
    f = "atacLibComplexQC",
    signature = "ATACProc",
    definition = function(atacProc,reportOutput=NULL,samInput=NULL,
                          singleEnd = FALSE,subsampleSize=Inf, ...){
        allpara <- c(list(Class = "LibComplexQC", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname LibComplexQC
#' @aliases libComplexQC
#' @export
libComplexQC<-function(samInput, reportOutput=NULL,singleEnd = FALSE,subsampleSize=Inf, ...){
    allpara <- c(list(Class = "LibComplexQC", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
