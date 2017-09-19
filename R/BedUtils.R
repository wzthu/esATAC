BedUtils<-R6::R6Class(
    classname = "BedUtils",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc, bedInput = NULL, bedOutput = NULL, reportOutput = NULL, mergePair = FALSE, downSample = NULL,
                              posOffset = 0L, negOffset= 0L, chrFilterList= NULL,select = TRUE,
                               sortBed = TRUE, uniqueBed = TRUE, minFregLen = 0,maxFregLen = 2e9, editable=FALSE){
            super$initialize("BedUtils",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["bedInput"]] <- atacProc$getParam("bedOutput");
                regexProcName<-sprintf("(BED|Bed|bed|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(BED|Bed|bed)"
            }

            if(!is.null(bedInput)){
                private$paramlist[["bedInput"]] <- bedInput;
            }
            if(is.null(bedOutput)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    private$paramlist[["bedOutput"]] <- private$getAutoPath(private$paramlist[["bedInput"]],regexProcName,".bed")
                }
            }else{
                private$paramlist[["bedOutput"]] <- bedOutput;
            }
          
            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    
                    private$paramlist[["reportOutput"]] <- private$getAutoPath(private$paramlist[["bedInput"]],regexProcName,".report")
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }    

            



            private$paramlist[["mergePair"]] <- mergePair;
            if(is.null(downSample)){
                private$paramlist[["downSample"]]<-2e9
            }else{
                private$paramlist[["downSample"]]<-downSample
            }
            private$paramlist[["posOffset"]] <- posOffset;
            private$paramlist[["negOffset"]] <- negOffset;
            private$paramlist[["filterList"]] <- chrFilterList;
            private$paramlist[["select"]] <- select;
            private$paramlist[["sortBed"]] <- sortBed
            private$paramlist[["uniqueBed"]] <- uniqueBed
            private$paramlist[["minFregLen"]] <- minFregLen
            private$paramlist[["maxFregLen"]] <- maxFregLen



            private$paramValidation()
        }
    ),

    private = list(
        processing = function(){
            #reportOutput <- private$paramlist[["reportOutput"]]
            #if(is.null(reportOutput)){
            #    reportOutput<-""
            #}
            qcval<-.bedOprUtils_call(ibedfile = private$paramlist[["bedInput"]],
                                        obedfile = private$paramlist[["bedOutput"]],
                                        reportPrefix = "",
                                        mergePair = private$paramlist[["mergePair"]],
                                        downSample = private$paramlist[["downSample"]],
                                        posOffset = private$paramlist[["posOffset"]],
                                        negOffset = private$paramlist[["negOffset"]],
                                        sortBed = private$paramlist[["sortBed"]],
                                        uniqueBed = private$paramlist[["uniqueBed"]],
                                        minFregLen = private$paramlist[["minFregLen"]],
                                        maxFregLen = private$paramlist[["maxFregLen"]],
                                        filterList = private$paramlist[["filterList"]],
                                        select = private$paramlist[["select"]])
            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
            
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop(paste("bedInput is requied"));
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileCreatable(private$paramlist[["bedOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            return(qcval[[item]])
        },
        getReportItemsImp = function(){
            return(c("total","save","filted","extlen","unique"))
        }
    )


)



#' @name atacBedUtils
#' @aliases atacBedUtils
#' @aliases bedUtils
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
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacSamToBed}}, 
#' \code{\link{atacBamToBed}}.
#' @param bedInput \code{Character} scalar. 
#' Bed file input path. 
#' @param bedOutput \code{Character} scalar. 
#' Bed file output path.
#' @param mergePair \code{Logical} scalar
#' Merge paired end interleave reads.
#' @param downSample \code{Integer} scalar
#' Down sample reads if the number is less than total number
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
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{bedUtils} instead.
#' @return An invisible \code{\link{ATACProc}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacSamToBed}} 
#' \code{\link{atacBamToBed}}
#' \code{\link{atacPeakCalling}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}
#' \code{\link{atacBamToBed}}

#' @rdname atacBedUtils
#' @export 
atacBedUtils <- function(atacProc, bedInput = NULL, bedOutput = NULL,  mergePair = FALSE, downSample = NULL,
                         posOffset = 0L, negOffset= 0L, chrFilterList= c("chrM"),select = FALSE,
                         sortBed = FALSE, uniqueBed = FALSE, minFregLen = 0,maxFregLen = 2e9){
    atacproc <- BedUtils$new(atacProc, bedInput = bedInput, bedOutput = bedOutput, reportOutput = NULL, mergePair = mergePair, downSample = downSample,
                             posOffset = posOffset, negOffset= negOffset, chrFilterList= chrFilterList,select = select,
                             sortBed = sortBed, uniqueBed = uniqueBed, minFregLen = minFregLen,maxFregLen = maxFregLen)
    atacproc$process()
    invisible(atacproc)
}

#' @rdname atacBedUtils
#' @export 
bedUtils <- function(bedInput, bedOutput = NULL, mergePair = FALSE, downSample = NULL,reportOutput = NULL,
                         posOffset = 0L, negOffset= 0L, chrFilterList= c("chrM"),select = FALSE,
                         sortBed = FALSE, uniqueBed = FALSE, minFregLen = 0,maxFregLen = 2e9){
    atacproc <- BedUtils$new(atacProc = NULL, bedInput = bedInput, bedOutput = bedOutput, reportOutput = reportOutput, mergePair = mergePair, downSample = downSample,
                             posOffset = posOffset, negOffset= negOffset, chrFilterList= chrFilterList,select = select,
                             sortBed = sortBed, uniqueBed = uniqueBed, minFregLen = minFregLen,maxFregLen = maxFregLen)
    atacproc$process()
    invisible(atacproc)
}
