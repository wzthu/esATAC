setClass(Class = "BedUtils",
         contains = "ATACProc"
)

setMethod(
    f = "initialize",
    signature = "BedUtils",
    definition = function(.Object,atacProc, ..., bedInput = NULL, bedOutput = NULL, reportOutput = NULL, mergePair = FALSE, downSample = NULL,
                          posOffset = 0L, negOffset= 0L, chrFilterList= NULL,select = TRUE,
                          sortBed = TRUE, uniqueBed = TRUE, minFregLen = 0,maxFregLen = 2e9, editable=FALSE){
        .Object <- init(.Object,"BedUtils",editable,list(arg1=atacProc))
        if(!is.null(atacProc)){
            .Object@paramlist[["bedInput"]] <- getParam(atacProc, "bedOutput");
            regexProcName<-sprintf("(BED|Bed|bed|%s)",getProcName(atacProc))
        }else{
            regexProcName<-"(BED|Bed|bed)"
        }

        if(!is.null(bedInput)){
            .Object@paramlist[["bedInput"]] <- bedInput;
        }
        if(is.null(bedOutput)){
            if(!is.null(.Object@paramlist[["bedInput"]])){
                .Object@paramlist[["bedOutput"]] <- getAutoPath(.Object, .Object@paramlist[["bedInput"]],regexProcName,".bed")
            }
        }else{
            .Object@paramlist[["bedOutput"]] <- bedOutput;
        }

        if(is.null(reportOutput)){
            if(!is.null(.Object@paramlist[["bedInput"]])){

                .Object@paramlist[["reportOutput"]] <- getAutoPath(.Object, .Object@paramlist[["bedInput"]],regexProcName,".report")
            }
        }else{
            .Object@paramlist[["reportOutput"]] <- reportOutput;
        }

        .Object@paramlist[["mergePair"]] <- mergePair;
        if(is.null(downSample)){
            .Object@paramlist[["downSample"]]<-2e9
        }else{
            .Object@paramlist[["downSample"]]<-downSample
        }
        .Object@paramlist[["posOffset"]] <- posOffset;
        .Object@paramlist[["negOffset"]] <- negOffset;
        .Object@paramlist[["filterList"]] <- chrFilterList;
        .Object@paramlist[["select"]] <- select;
        .Object@paramlist[["sortBed"]] <- sortBed
        .Object@paramlist[["uniqueBed"]] <- uniqueBed
        .Object@paramlist[["minFregLen"]] <- minFregLen
        .Object@paramlist[["maxFregLen"]] <- maxFregLen

        paramValidation(.Object)
        .Object
    }
)

setMethod(
    f = "processing",
    signature = "BedUtils",
    definition = function(.Object,...){
        qcval<-.bedOprUtils_call(ibedfile = .Object@paramlist[["bedInput"]],
                                 obedfile = .Object@paramlist[["bedOutput"]],
                                 reportPrefix = "",
                                 mergePair = .Object@paramlist[["mergePair"]],
                                 downSample = .Object@paramlist[["downSample"]],
                                 posOffset = .Object@paramlist[["posOffset"]],
                                 negOffset = .Object@paramlist[["negOffset"]],
                                 sortBed = .Object@paramlist[["sortBed"]],
                                 uniqueBed = .Object@paramlist[["uniqueBed"]],
                                 minFregLen = .Object@paramlist[["minFregLen"]],
                                 maxFregLen = .Object@paramlist[["maxFregLen"]],
                                 filterList = .Object@paramlist[["filterList"]],
                                 select = .Object@paramlist[["select"]])
        write.table(as.data.frame(qcval),file = .Object@paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        .Object
    }
)


setMethod(
    f = "checkRequireParam",
    signature = "BedUtils",
    definition = function(.Object,...){
        if(is.null(.Object@paramlist[["bedInput"]])){
            stop(paste("bedInput is requied"));
        }
    }
)

setMethod(
    f = "checkAllPath",
    signature = "BedUtils",
    definition = function(.Object,...){
        checkFileExist(.Object,.Object@paramlist[["bedInput"]]);
        checkFileCreatable(.Object,.Object@paramlist[["bedOutput"]]);
    }
)


setMethod(
    f = "getReportValImp",
    signature = "BedUtils",
    definition = function(.Object, item){
        qcval <- as.list(read.table(file= .Object@paramlist[["reportOutput"]],header=TRUE))
        return(qcval[[item]])
    }
)


setMethod(
    f = "getReportItemsImp",
    signature = "BedUtils",
    definition = function(.Object){
        return(c("total","save","filted","extlen","unique"))
    }
)


#' @title process bed file with limit memory
#' @description
#' This function is used to
#' merge interleave paired end reads in bed,
#' downsample bed reads,
#' shift bed reads,
#' filter bed reads according to chromosome,
#' filter bed reads according to fregment size,
#' sort bed,
#' remove duplicates reads in bed.
#' @param atacProc \code{\link{ATACProc-class}} object scalar.
#' It has to be the return value of upstream process:
#' \code{\link{atacBam2Bed}}
#' \code{\link{bam2bed}}
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' @param bedInput \code{Character} scalar.
#' Bed file input path.
#' @param bedOutput \code{Character} scalar.
#' Bed file output path.
#' @param mergePair \code{Logical} scalar
#' Merge paired end interleave reads.
#' @param downSample \code{Integer} scalar
#' Down sample reads if the number is less than total number
#' @param reportOutput \code{Character} scalar.
#' Report output file path.
#' @param posOffset \code{Integer} scalar
#' The offset that positive strand reads will shift.
#' @param negOffset \code{Integer} scalar
#' The offset that negative strand reads will shift.
#' @param chrFilterList \code{Character} vector
#' The chromatin(or regex of chromatin) will be retain/discard
#' if \code{select} is TRUE/FALSE
#' @param select \code{Logical} scalar
#' The chromatin in \code{chrFilterList} will be retain if TRUE. default: FALSE
#' @param sortBed \code{Logical} scalar
#' Sort bed file in the order of chromatin, start, end
#' @param uniqueBed \code{Logical} scalar
#' Remove duplicates reads in bed if TRUE. default: FALSE
#' @param minFregLen \code{Integer} scalar
#' The minimum fregment size will be retained.
#' @param maxFregLen \code{Integer} scalar
#' The maximum fregment size will be retained.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' \code{atacProc} should be set \code{NULL}
#' or you can use \code{bedUtils} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacBam2Bed}}
#' \code{\link{bam2bed}}
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacFregLenDistr}}
#' \code{\link{atacExtractCutSite}}
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacTSSQC}}
#' \code{\link{atacBedToBigWig}}
#'
#' @examples
#' library(R.utils)
#' library(magrittr)
#' td <- tempdir()
#' options(atacConf=setConfigure("tmpdir",td))
#'
#' sambzfile <- system.file(package="ATACpipe", "extdata", "Example.sam.bz2")
#' samfile <- file.path(td,"Example.sam")
#' bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
#' atacproc<-samToBed(samInput = samfile) %>%
#' atacBedUtils(maxFregLen = 100, chrFilterList = NULL)
#'
#' @name atacBedUtils
#' @export
#' @docType methods
#' @rdname atacBedUtils-methods
setGeneric("atacBedUtils",function(atacProc, bedInput = NULL, bedOutput = NULL,  mergePair = FALSE, downSample = NULL,
                                  posOffset = 0L, negOffset= 0L, chrFilterList= c("chrM"),select = FALSE,
                                  sortBed = FALSE, uniqueBed = FALSE, minFregLen = 0,maxFregLen = 2e9) standardGeneric("atacBedUtils"))
#' @rdname atacBedUtils-methods
#' @aliases atacBedUtils
setMethod(
    f = "atacBedUtils",
    signature = "ATACProc",
    definition = function(atacProc, bedInput = NULL, bedOutput = NULL,  mergePair = FALSE, downSample = NULL,
                          posOffset = 0L, negOffset= 0L, chrFilterList= c("chrM"),select = FALSE,
                          sortBed = FALSE, uniqueBed = FALSE, minFregLen = 0,maxFregLen = 2e9){
        atacproc <- new(
            "BedUtils",
            atacProc = atacProc,
            bedInput = bedInput,
            bedOutput = bedOutput,
            mergePair = mergePair,
            downSample = downSample,
            posOffset = posOffset,
            negOffset = negOffset,
            chrFilterList = chrFilterList,
            select = select,
            sortBed = sortBed,
            uniqueBed = uniqueBed,
            minFregLen = minFregLen,
            maxFregLen = maxFregLen)
        atacproc <- process(atacproc)
        invisible(atacproc)
    }
)
#' @rdname atacBedUtils-methods
#' @export
bedUtils <- function(bedInput, bedOutput = NULL, mergePair = FALSE, downSample = NULL,reportOutput = NULL,
                         posOffset = 0L, negOffset= 0L, chrFilterList= c("chrM"),select = FALSE,
                         sortBed = FALSE, uniqueBed = FALSE, minFregLen = 0,maxFregLen = 2e9){
    atacproc <- new(
        "BedUtils",
        atacProc = NULL,
        bedInput = bedInput,
        bedOutput = bedOutput,
        mergePair = mergePair,
        downSample = downSample,
        posOffset = posOffset,
        negOffset = negOffset,
        chrFilterList = chrFilterList,
        select = select,
        sortBed = sortBed,
        uniqueBed = uniqueBed,
        minFregLen = minFregLen,
        maxFregLen = maxFregLen)
    atacproc <- process(atacproc)
    invisible(atacproc)
}
