LibComplexQC <-R6Class(
    classname = "LibComplexQC",
    inherit = ATACProc,
    public = list(
        initialize = function(atacProc,reportOutput=NULL,samInput=NULL,singleEnd = FALSE,subsample=TRUE,subsampleSize=4e6,editable=FALSE){
            super$initialize("LibComplexQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["samInput"]] <- atacProc$getParam("samOutput");
                regexProcName<-sprintf("(SAM|Sam|sam|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(SAM|Sam|sam)"
                private$singleEnd<-singleEnd
            }
            if(!is.null(samInput)){
                private$paramlist[["samInput"]] <- samInput;
            }
            if(is.null(reportOutput)){
                if(!is.null(private$paramlist[["samInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["samInput"]],regexProcName)
                    private$paramlist[["reportOutput"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
            }else{
                private$paramlist[["reportOutput"]] <- reportOutput;
            }

            private$paramlist[["subsample"]] <- subsample;
            private$paramlist[["subsampleSize"]] <- subsampleSize;



            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            if(!private$singleEnd){
                .sam2bed_merge_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportOutput"]],".tmp"),
                                    posOffset = 0, negOffset = 0,sortBed = !private$paramlist[["subsample"]],
                                    uniqueBed = FALSE, filterList = NULL,minFregLen = 0,maxFregLen = 1000000,saveExtLen = FALSE ,downSample=private$paramlist[["subsampleSize"]])
                print("test00")
            }else{
                .sam2bed_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportOutput"]],".tmp"),
                              posOffset = 0, negOffset = 0, sortBed = !private$paramlist[["subsample"]],uniqueBed = FALSE,  filterList = NULL,downSample=private$paramlist[["subsampleSize"]])
                print("test01")
            }
            print("test1")
            qcval<-.lib_complex_qc_call(bedfile=paste0(private$paramlist[["reportOutput"]],".tmp"), sortedBed=!private$paramlist[["subsample"]], max_reads=private$paramlist[["subsampleSize"]])
            print("test2")
            unlink(paste0(private$paramlist[["reportOutput"]],".tmp"))
            print(qcval)
            print(private$paramlist[["reportOutput"]])
            write.table(as.data.frame(qcval),file = private$paramlist[["reportOutput"]],quote=FALSE,sep="\t",row.names=FALSE)
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["samInput"]])){
                stop("samInput is required.")
            }


        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["samInput"]]);
            private$checkFileCreatable(private$paramlist[["reportOutput"]]);
        },
        getReportValImp = function(item){
            qcval <- as.list(read.table(file= private$paramlist[["reportOutput"]],header=TRUE))
            if(item == "report"){
                return(data.frame(Item=names(qcval),Value=as.character(qcval)))
            }else{
                return(qcval[[item]])
            }
        },
        getReportItemsImp = function(){
            return(c("report","NRF","PBC1","PBC2","one","two","total","reads"))
        }
    )


)
#' @name atacLibComplexQC
#' @aliases atacLibComplexQC
#' @aliases libComplexQC
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
#' @param atacProc \code{\link{ATACProc}} object scalar. 
#' It has to be the return value of upstream process:
#' \code{\link{atacBowtie2Mapping}}.
#' @param reportOutput \code{Character} scalar. 
#' The report file path 
#' @param samInput \code{Character} scalar. 
#' The SAM file input path.
#' @param singleEnd \code{Character} scalar. 
#' Single end data if TRUE. Paired end data if FALSE. 
#' @param subsample \code{Logical} scalar.
#' Down sample reads if TRUE, otherwise FALSE.
#' @param subsampleSize \code{Integer} scalar.
#' Down sample reads if the number is less than total number 
#' when \code{subsample} is TRUE
#' @details The parameter related to input and output file path
#' will be automatically 
#' obtained from \code{\link{ATACProc}} object(\code{atacProc}) or 
#' generated based on known parameters 
#' if their values are default(e.g. \code{NULL}).
#' Otherwise, the generated values will be overwrited.
#' If you want to use this function independently, 
#' \code{atacProc} should be set \code{NULL} 
#' or you can use \code{fregLenDistr} instead.
#' @return An invisible \code{\link{libComplexQC}} object scalar for downstream analysis.
#' @author Zheng Wei
#' @seealso 
#' \code{\link{atacSamToBed}} 
#' \code{\link{atacBedUtils}}

#' @rdname atacLibComplexQC
#' @export 
atacLibComplexQC<-function(atacProc,reportOutput=NULL,samInput=NULL,singleEnd = FALSE,subsample=TRUE,subsampleSize=4*10e6){
    libqc<-LibComplexQC$new(atacProc=atacProc, reportOutput=reportOutput,samInput=samInput,singleEnd = singleEnd,
                            subsample=subsample,subsampleSize=subsampleSize,editable=FALSE)
    libqc$process()
    invisible(libqc)
}
#' @rdname atacLibComplexQC
#' @export 
libComplexQC<-function(samInput, reportOutput=NULL,singleEnd = FALSE,subsample=TRUE,subsampleSize=4*10e6){
    libqc<-LibComplexQC$new(atacProc=NULL,reportOutput=reportOutput,samInput=samInput,singleEnd = singleEnd,
                            subsample=subsample,subsampleSize=subsampleSize,editable=FALSE)
    libqc$process()
    invisible(libqc)
}
