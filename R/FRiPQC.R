setClass(Class = "FRiPQC",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "FRiPQC",
    definition = function(.Object, prevSteps, ...){
        allparam <- list(...)
        readsBedInput <- allparam[["readsBedInput"]]
        peakBedInput <- allparam[["peakBedInput"]]
        bsgenome <- allparam[["bsgenome"]]
        reportOutput <- allparam[["reportOutput"]]
        
        if(length(prevSteps) > 0){
            if(!is.null(prevSteps[[1]])){
                atacProc <- prevSteps[[1]]
                atacProc<-c(unlist(atacProc),list())
                atacProcReads <- atacProc[[length(atacProc)]]
                input(.Object)[["readsBedInput"]] <- output(atacProcReads)[["bedOutput"]]
            }
            if(!is.null(prevSteps[[2]])){
                atacProc <- prevSteps[[2]]
                atacProc<-c(unlist(atacProc),list())
                atacProcPeak <- atacProc[[length(atacProc)]]
                input(.Object)[["peakBedInput"]] <- output(atacProcPeak)[["bedOutput"]]
            }
        }
        

        if(!is.null(readsBedInput)){
            input(.Object)[["readsBedInput"]] <- readsBedInput
        }
        if(!is.null(peakBedInput)){
            input(.Object)[["peakBedInput"]] <- peakBedInput
        }

        if(is.null(reportOutput)){
            if(!is.null(input(.Object)[["peakBedInput"]])){
                output(.Object)[["reportOutput"]] <- getAutoPath(.Object, input(.Object)[["peakBedInput"]], "BED|Bed|bed" , "report.txt")
            }
        }else{
            output(.Object)[["reportOutput"]] <- reportOutput;
        }
        if(is.null(bsgenome)){
            param(.Object)[["bsgenome"]] <- getRefRc("bsgenome")
        }else{
            param(.Object)[["bsgenome"]] <- bsgenome
        }
        

        .Object
    }
)


setMethod(
    f = "processing",
    signature = "FRiPQC",
    definition = function(.Object,...){
        qcval=list();
        genome <- seqinfo(param(.Object)[["bsgenome"]])
        message("load reads")
#        gr_a <- rtracklayer::import(con = .Object@paramlist[["readsBedInput"]], genome = genome,format = "bed")
        gr_a <- rtracklayer::import(con = input(.Object)[["readsBedInput"]], format = "bed")
        message("load peak")
#        gr_b <- rtracklayer::import(con = .Object@paramlist[["peakBedInput"]], genome = genome,format = "bed")
        gr_b <- rtracklayer::import(con = input(.Object)[["peakBedInput"]],format = "bed")
        message("get overlap number")
        qcval[["peakReads"]]<-length(subsetByOverlaps(gr_a, gr_b, ignore.strand = TRUE))
        message("finish overlap")
        qcval[["totalReads"]]<-length(gr_a)
        qcval[["totalPeaks"]]<-length(gr_b)
        qcval[["FRiP"]]<-qcval[["peakReads"]]/qcval[["totalReads"]]
        #unlink(paste0(.Object@paramlist[["reportPrefix"]],".tmp"))
        ####.Object@paramlist[["qcval"]]<-qcval
        write.table(as.data.frame(qcval),file = output(.Object)[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        .Object
    }
)

setMethod(
    f = "checkRequireParam",
    signature = "FRiPQC",
    definition = function(.Object,...){
        if(is.null(input(.Object)[["readsBedInput"]])){
            stop("readsBedInput is required.")
        }
        if(is.null(input(.Object)[["peakBedInput"]])){
            stop("peakBedInput is required.")
        }
    }
)





setMethod(
    f = "genReport",
    signature = "FRiPQC",
    definition = function(.Object, ...){
        qcval <- as.list(read.table(file= output(.Object)[["reportOutput"]],header=TRUE))
        cqcval<-as.character(qcval)
        cqcval[4] <- sprintf("%.2f",as.numeric(cqcval[4]))
        report(.Object)$report <- (data.frame(Item=c("The number of reads in peak",
                                                     "The number of total reads",
                                                     "The number of total peaks",
                                                     "FRiP"),Value=cqcval))
        report(.Object) <- c(report(.Object), qcval)
        .Object
    }
)



#' @name FRiPQC
#' @title Quality control for fraction of reads in peaks (FRiP)
#' @description
#' Calculate the fraction of reads falling within peak regions
#' @param atacProcReads \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacBedUtils}}
#' \code{\link{bedUtils}}
#' @param atacProcPeak \code{\link{ATACProc-class}} object scalar.
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
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' or you can use \code{fripQC} instead.
#' @return An invisible \code{\link{fripQC}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacSamToBed}}
#' \code{\link{atacBedUtils}}
#' @examples
#' library(R.utils)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(magrittr)
#' td <- tempdir()
#' setTmpDir(td)
#'
#' bedbzfile <- system.file(package="esATAC", "extdata", "chr20.50000.bed.bz2")
#' bedfile <- file.path(td,"chr20.50000.bed")
#' bunzip2(bedbzfile,destname=bedfile,overwrite=TRUE,remove=FALSE)
#'
#' bedUtils(bedInput = bedfile,maxFragLen = 100, chrFilterList = NULL)  %>%
#' atacPeakCalling %>% atacFripQC(bsgenome=BSgenome.Hsapiens.UCSC.hg19)
#' 
#' 
#' dir(td)
#' 
#' 
#' 


setGeneric("atacFripQC",function(atacProc,bsgenome = NULL,
                                  reportOutput=NULL,readsBedInput=NULL,
                                  peakBedInput=NULL, ...) standardGeneric("atacFripQC"))


#' @rdname FRiPQC
#' @aliases atacFripQC
#' @export
setMethod(
    f = "atacFripQC",
    signature = "ATACProc",
    definition = function(atacProc,bsgenome = NULL,
                          reportOutput=NULL,readsBedInput=NULL,
                          peakBedInput=NULL, ...){
        allpara <- c(list(Class = "FRiPQC", prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)


#' @rdname FRiPQC
#' @aliases fripQC
#' @export
fripQC<-function(readsBedInput, peakBedInput,bsgenome = NULL, reportOutput=NULL, ...){
    allpara <- c(list(Class = "FRiPQC", prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}
