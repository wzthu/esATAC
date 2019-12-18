setClass(Class = "BedUtils",
         contains = "ATACProc"
)


setMethod(
    f = "init",
    signature = "BedUtils",
    definition = function(.Object,prevSteps = list(),...){
        allparam <- list(...)
        bedInput <- allparam[["bedInput"]]
        bedOutput <- allparam[["bedOutput"]]
        mergePair <- allparam[["mergePair"]]
        downSample <- allparam[["downSample"]]
        reportOutput <- allparam[["reportOutput"]]
        posOffset <- allparam[["posOffset"]]
        negOffset <- allparam[["negOffset"]]
        chrFilterList <- allparam[["chrFilterList"]]
        select <- allparam[["select"]]
        sortBed <- allparam[["sortBed"]]
        uniqueBed <- allparam[["uniqueBed"]]
        minFragLen <- allparam[["minFragLen"]]
        maxFragLen <- allparam[["maxFragLen"]]
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
        if(is.null(bedOutput)){
            if(!is.null(input(.Object)[["bedInput"]])){
                output(.Object)[["bedOutput"]] <- 
                    getAutoPath(.Object, input(.Object)[["bedInput"]],"BED|bed|Bed","bed")
            }
        }else{
            output(.Object)[["bedOutput"]] <- bedOutput;
        }

        if(is.null(reportOutput)){
            if(!is.null(input(.Object)[["bedInput"]])){

                output(.Object)[["reportOutput"]] <- 
                    getAutoPath(.Object, input(.Object)[["bedInput"]],"BED|bed|Bed",".report")
            }
        }else{
            output(.Object)[["reportOutput"]] <- reportOutput;
        }

        param(.Object)[["mergePair"]] <- mergePair;
        if(is.null(downSample)){
            param(.Object)[["downSample"]]<-2e9
        }else{
            param(.Object)[["downSample"]]<-downSample
        }
        param(.Object)[["posOffset"]] <- posOffset;
        param(.Object)[["negOffset"]] <- negOffset;
        param(.Object)[["filterList"]] <- chrFilterList;
        param(.Object)[["select"]] <- select;
        param(.Object)[["sortBed"]] <- sortBed
        param(.Object)[["uniqueBed"]] <- uniqueBed
        param(.Object)[["minFragLen"]] <- minFragLen
        param(.Object)[["maxFragLen"]] <- maxFragLen

        .Object
    }
)

setMethod(
    f = "processing",
    signature = "BedUtils",
    definition = function(.Object,...){
        qcval<-.bedOprUtils_call(ibedfile = input(.Object)[["bedInput"]],
                                 obedfile = output(.Object)[["bedOutput"]],
                                 reportPrefix = "",
                                 mergePair = param(.Object)[["mergePair"]],
                                 downSample = param(.Object)[["downSample"]],
                                 posOffset = param(.Object)[["posOffset"]],
                                 negOffset = param(.Object)[["negOffset"]],
                                 sortBed = param(.Object)[["sortBed"]],
                                 uniqueBed = param(.Object)[["uniqueBed"]],
                                 minFragLen = param(.Object)[["minFragLen"]],
                                 maxFragLen = param(.Object)[["maxFragLen"]],
                                 filterList = param(.Object)[["filterList"]],
                                 select = param(.Object)[["select"]])
        write.table(as.data.frame(qcval),file = output(.Object)[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        
        .Object
    }
)


setMethod(
    f = "genReport",
    signature = "BedUtils",
    definition = function(.Object, ...){
        qcval <- as.list(read.table(file= output(.Object)[["reportOutput"]],header=TRUE))
        for(n in names(qcval)){
            report(.Object)[[n]] <- qcval[[n]]
        }  
        .Object
    }
)


#' @name  BedUtils
#' @title process bed file with limit memory
#' @description
#' This function is used to
#' merge interleave paired end reads in bed,
#' downsample bed reads,
#' shift bed reads,
#' filter bed reads according to chromosome,
#' filter bed reads according to fragment size,
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
#' @param minFragLen \code{Integer} scalar
#' The minimum fragment size will be retained.
#' @param maxFragLen \code{Integer} scalar
#' The maximum fragment size will be retained.
#' @param newStepType \code{Character} scalar.
#' New step type name for different default parameters.
#' @param ... Additional arguments, currently unused.
#' @details The parameter related to input and output file path
#' will be automatically
#' obtained from \code{\link{ATACProc-class}} object(\code{atacProc}) or
#' generated based on known parameters
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently,
#' you can use \code{bedUtils} instead.
#' @return An invisible \code{\link{ATACProc-class}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso
#' \code{\link{atacBam2Bed}}
#' \code{\link{bam2bed}}
#' \code{\link{atacSamToBed}}
#' \code{\link{samToBed}}
#' \code{\link{atacFragLenDistr}}
#' \code{\link{atacExtractCutSite}}
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacTSSQC}}
#' \code{\link{atacBedToBigWig}}
#'
#' @examples
#' library(R.utils)
#' library(magrittr)
#' td <- tempdir()
#' setTmpDir(td)
#'
#' sambzfile <- system.file(package="esATAC", "extdata", "Example.sam.bz2")
#' samfile <- file.path(td,"Example.sam")
#' bunzip2(sambzfile,destname=samfile,overwrite=TRUE,remove=FALSE)
#' atacproc<-samToBed(samInput = samfile) %>%
#' atacBedUtils(maxFragLen = 100, chrFilterList = NULL)
#'

setGeneric("atacBedUtils",function(atacProc, bedInput = NULL, 
                                   bedOutput = NULL,  mergePair = FALSE, 
                                   downSample = NULL,
                                  posOffset = 0L, negOffset= 0L, 
                                  chrFilterList= c("chrM"),select = FALSE,
                                  sortBed = FALSE, uniqueBed = FALSE, 
                                  minFragLen = 0,maxFragLen = 2e9,  
                                  newStepType = "BedUtils",...) 
    standardGeneric("atacBedUtils"))

#' @rdname BedUtils
#' @aliases atacBedUtils
#' @export
setMethod(
    f = "atacBedUtils",
    signature = "ATACProc",
    definition = function(atacProc, bedInput = NULL, 
                          bedOutput = NULL,  mergePair = FALSE, 
                          downSample = NULL,
                          posOffset = 0L, negOffset= 0L, 
                          chrFilterList = c("chrM"),select = FALSE,
                          sortBed = FALSE, uniqueBed = FALSE, 
                          minFragLen = 0,maxFragLen = 2e9, 
                          newStepType = "BedUtils", ...){
        allpara <- c(list(Class = regAttachedStep(newStepType,"BedUtils"), 
                          prevSteps = list(atacProc)),as.list(environment()),list(...))
        step <- do.call(new,allpara)
        invisible(step)
    }
)
#' @rdname BedUtils
#' @aliases bedUtils
#' @export
bedUtils <- function(bedInput, bedOutput = NULL, 
                     mergePair = FALSE, downSample = NULL,
                     reportOutput = NULL,
                     posOffset = 0L, negOffset= 0L, 
                     chrFilterList = c("chrM"),select = FALSE,
                     sortBed = FALSE, uniqueBed = FALSE, 
                     minFragLen = 0,maxFragLen = 2e9, 
                     newStepType = "BedUtils", ...){
    allpara <- c(list(Class = regAttachedStep(newStepType,"BedUtils"), 
                      prevSteps = list()),as.list(environment()),list(...))
    step <- do.call(new,allpara)
    invisible(step)
}




