BedUtils<-R6::R6Class(
    classname = "BedUtils",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc, bedInput = NULL, bedOutput = NULL, report = TRUE, reportPrefix = NULL, mergePair = FALSE, downSample = NULL,
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
            if(report){
                if(is.null(reportPrefix)){
                    if(!is.null(private$paramlist[["bedInput"]])){

                        private$paramlist[["reportPrefix"]] <- private$getAutoPath(private$paramlist[["bedInput"]],regexProcName,".report")
                    }
                }else{
                    private$paramlist[["reportPrefix"]] <- reportPrefix;
                }

            }else{
                private$paramlist[["reportPrefix"]] <- ""
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
            .bedOprUtils_call(ibedfile = private$paramlist[["bedInput"]],
                                        obedfile = private$paramlist[["bedOutput"]],
                                        reportPrefix = private$paramlist[["reportPrefix"]],
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
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["bedInput"]])){
                stop(paste("bedInput is requied"));
            }
        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["bedInput"]]);
            private$checkFileCreatable(private$paramlist[["bedOutput"]]);
        }
    )


)



#' convert sam to bed
#' @param ATAC_obj obj returned from ATAC_mapping
#' @param bedInput sam file dir
#' @param bedfile bed file dir
#' @param readlen reads length
#' @export
atacBedUtils <- function(atacProc, bedInput = NULL, bedOutput = NULL, report = TRUE, reportPrefix = NULL, mergePair = FALSE, downSample = NULL,
                         posOffset = 0L, negOffset= 0L, chrFilterList= c("chrM"),select = FALSE,
                         sortBed = FALSE, uniqueBed = FALSE, minFregLen = 0,maxFregLen = 2e9){
    atacproc <- BedUtils$new(atacProc, bedInput = bedInput, bedOutput = bedOutput, report = report, reportPrefix = reportPrefix, mergePair = mergePair, downSample = downSample,
                             posOffset = posOffset, negOffset= negOffset, chrFilterList= chrFilterList,select = select,
                             sortBed = sortBed, uniqueBed = uniqueBed, minFregLen = minFregLen,maxFregLen = maxFregLen)
    atacproc$process()
    return(atacproc)
}
