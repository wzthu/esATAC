LibComplexQC <-R6Class(
    classname = "LibComplexQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProc,reportPrefix=NULL,samInput=NULL,paired = FALSE,subsample=TRUE,subsampleSize=4e6,editable=FALSE){
            super$initialize("LibComplexQC",editable,list(arg1=atacProc))
            if(!is.null(atacProc)){
                private$paramlist[["samInput"]] <- atacProc$getParam("samOutput");
                regexProcName<-sprintf("(SAM|Sam|sam|%s)",atacProc$getProcName())
            }else{
                regexProcName<-"(SAM|Sam|sam)"
            }
            if(!is.null(samInput)){
                private$paramlist[["samInput"]] <- samInput;
            }
            if(is.null(reportPrefix)){
                if(!is.null(private$paramlist[["bedInput"]])){
                    prefix<-private$getBasenamePrefix(private$paramlist[["bedInput"]],regexProcName)
                    private$paramlist[["reportPrefix"]] <- file.path(.obtainConfigure("tmpdir"),paste0(prefix,".",self$getProcName(),".report"))
                }
                #private$paramlist[["reportPrefix"]] <- paste(private$paramlist[["samInput"]],".libcomplexreport",sep="");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }

            private$paramlist[["subsample"]] <- subsample;
            private$paramlist[["subsampleSize"]] <- subsampleSize;
            private$paramlist[["paired"]]<-paired


            private$paramValidation()
        }
    ),
    private = list(
        processing = function(){
            if(private$paramlist[["paired"]]){
                .sam2bed_merge_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportPrefix"]],".tmp"),
                                    posOffset = 0, negOffset = 0,sortBed = !private$paramlist[["subsample"]],
                                    uniqueBed = FALSE, filterList = NULL,minFregLen = 0,maxFregLen = 1000000,saveExtLen = FALSE ,downSample=private$paramlist[["subsampleSize"]])
            }else{
                .sam2bed_call(samfile = private$paramlist[["samInput"]], bedfile = paste0(private$paramlist[["reportPrefix"]],".tmp"),
                              posOffset = 0, negOffset = 0, sortBed = !private$paramlist[["subsample"]],uniqueBed = FALSE,  filterList = NULL,downSample=private$paramlist[["subsampleSize"]])
            }

            qcval<-.lib_complex_qc_call(bedfile=paste0(private$paramlist[["reportPrefix"]],".tmp"), sortedBed=!private$paramlist[["subsample"]], max_reads=private$paramlist[["subsampleSize"]])

            unlink(paste0(private$paramlist[["reportPrefix"]],".tmp"))
            private$paramlist[["qcval"]]<-qcval
            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)
        },
        checkRequireParam = function(){
            if(is.null(private$paramlist[["samInput"]])){
                stop("samInput is required.")
            }


        },
        checkAllPath = function(){
            private$checkFileExist(private$paramlist[["samInput"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
        }
    )


)

atacLibComplexQC<-function(atacProc,reportPrefix=NULL,samInput=NULL,paired = FALSE,subsample=TRUE,subsampleSize=4*10e6){
    libqc<-LibComplexQC$new(atacProc,reportPrefix=reportPrefix,samInput=samInput,paired = paired,
                            subsample=subsample,subsampleSize=subsampleSize,editable=FALSE)
    libqc$process()
    return(libqc)
}
