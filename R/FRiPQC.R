FRiPQC <-R6Class(
    classname = "FRiPQC",
    inherit = BaseProc,
    public = list(
        initialize = function(atacProcReads,atacProcPeak,reportPrefix=NULL,readsBedInput=NULL,peakBedInput,editable=FALSE){
            super$initialize("FRiPQC",editable,list(arg1=atacProcReads,arg2=atacProcPeak))
            if(!is.null(atacProcReads)){
                private$paramlist[["readsBedInput"]] <- atacProcReads$getParam("bedOutput");
                private$paramlist[["readsCount"]] <- atacProcReads$getParam("readsCount")
            }
            if(!is.null(atacProcPeak)){
                private$paramlist[["peakBedInput"]] <- atacProcReads$getParam("bedOutput");
            }
            if(!is.null(samInput)){
                private$paramlist[["samInput"]] <- samInput;
            }
            if(is.null(reportPrefix)){
                private$paramlist[["reportPrefix"]] <- paste0(private$paramlist[["samInput"]],".FripQCreport");
            }else{
                private$paramlist[["reportPrefix"]] <- reportPrefix;
            }
            private$checkFileExist(private$paramlist[["readsBedInput"]]);
            private$checkFileExist(private$paramlist[["peakBedInput"]]);
            private$checkPathExist(private$paramlist[["reportPrefix"]]);
            private$checkRequireParam();
        },
        processing = function(){
            super$processing()
            qcval=list();
            qcval[["peakReads"]]<-length(HelloRanges::R_bedtools_intersect(wa=TRUE,u=TRUE,a=private$paramlist[["readsBedInput"]],b=private$paramlist[["peakBedInput"]],bed = TRUE))
            if(is.null(private$paramlist[["readsCount"]])){
                qcval[["totalReads"]]<-R.utils::countLines(rivate$paramlist[["readsBedInput"]])
            }
            qcval[["FRiP"]]<-qcval[["peakReads"]]/qcval[["totalReads"]]
            #unlink(paste0(private$paramlist[["reportPrefix"]],".tmp"))
            private$paramlist[["qcval"]]<-qcval
            qcval<-as.matrix(qcval)
            write.table(qcval,file = paste0(private$paramlist[["reportPrefix"]],".txt"),sep="\t",quote = FALSE,col.names = FALSE)
            private$finish <- TRUE
        },
        setResultParam = function(fastqOutput1, fastqOutput2=NULL){
            super$setResultParam();
            private$paramlist[["fastqOutput1"]] <- fastqOutput1
            private$paramlist[["fastqOutput2"]] <- fastqOutput2
        }
    ),
    private = list(
        checkRequireParam = function(){
            if(private$editable){
                return();
            }
            if(is.null(private$paramlist[["readsBedInput"]])){
                stop("readsBedInput is required.")
            }
            if(is.null(private$paramlist[["peakBedInput"]])){
                stop("peakBedInput is required.")
            }
            if(is.null(private$paramlist[["readsCount"]])){
                stop("readsCount is required.")
            }


        }
    )


)
